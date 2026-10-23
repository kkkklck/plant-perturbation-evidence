from __future__ import annotations

from pathlib import Path
import hashlib
import json
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_lateral_root" / "GSE158761"
OUT = BASE / "plantvc_v1_blind"
META = BASE / "plantvc_v1_probe" / "metadata.tsv.gz"
RULES = ROOT / "PLANTVC_V1_GSE158761_BLIND_RULES_2026-10-20.md"
SOURCE = BASE / "GSE158761_all_cells_labeled.rds"
CORE = ROOT / "scripts" / "108_plantvc_v1.py"

ELIGIBLE = [
    ("Ambiguous Stele Cells", "ambiguous_stele"),
    ("Columella/ Root Cap", "columella_root_cap"),
    ("Epidermis", "epidermis"),
    ("Lateral Root Endodermis", "lateral_root_endodermis"),
    ("Mature Pericycle", "mature_pericycle"),
    ("Phloem", "phloem"),
    ("Phloem Pole Pericycle", "phloem_pole_pericycle"),
    ("Xylem", "xylem"),
    ("Xylem Pole Pericycle", "xylem_pole_pericycle"),
]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    if not OUT.exists():
        raise SystemExit("Freeze output directory missing")
    contexts_path = OUT / "contexts.tsv"
    used_path = OUT / "freeze_cells_used.tsv"
    withheld_path = OUT / "target_cells_withheld.tsv"
    for p in [contexts_path, used_path, withheld_path, RULES, SOURCE, CORE, META]:
        if not p.exists():
            raise SystemExit(f"Missing required file: {p}")

    meta = pd.read_csv(META, sep="\t")
    ct = pd.crosstab(meta["celltype"], meta["sample_number"])
    required_samples = [1, 2, 3, 5, 6, 7]
    derived = sorted(ct.index[(ct[required_samples] >= 20).all(axis=1)].astype(str).tolist())
    expected = sorted([x[0] for x in ELIGIBLE])
    if derived != expected:
        raise RuntimeError(f"Eligible cell-type mismatch: derived={derived} expected={expected}")

    used = pd.read_csv(used_path, sep="\t")
    withheld = pd.read_csv(withheld_path, sep="\t")
    if set(used["sample_number"].astype(int)) - {1, 2, 5, 6}:
        raise RuntimeError("Freeze used target sample numbers")
    if set(withheld["sample_number"].astype(int)) != {3, 7}:
        raise RuntimeError("Withheld target samples are not exactly 3 and 7")
    if set(used["cell"]).intersection(set(withheld["cell"])):
        raise RuntimeError("Used/withheld cell overlap detected")

    contexts = pd.read_csv(contexts_path, sep="\t")
    if sorted(contexts["celltype"].astype(str).tolist()) != expected:
        raise RuntimeError("contexts.tsv does not contain frozen eligible cell types")

    cfg_dir = OUT / "configs"
    cfg_dir.mkdir(exist_ok=True)
    ranking_root = OUT / "rankings"
    ranking_root.mkdir(exist_ok=True)

    records = []
    for celltype, slug in ELIGIBLE:
        evidence = OUT / f"evidence_{slug}.tsv.gz"
        if not evidence.exists():
            raise RuntimeError(f"Missing evidence file: {evidence}")
        ev = pd.read_csv(evidence, sep="\t")
        req = {"gene", "celltype", "technical_universe", "r8", "r8_repA", "r8_repB", "broad_shared_r8"}
        if not req.issubset(ev.columns):
            raise RuntimeError(f"Evidence columns missing for {celltype}: {sorted(req-set(ev.columns))}")
        if ev["celltype"].astype(str).nunique() != 1 or ev["celltype"].astype(str).iloc[0] != celltype:
            raise RuntimeError(f"Evidence celltype mismatch for {celltype}")

        outdir = ranking_root / slug
        cfg = {
            "task_name": f"GSE158761 8h -> 20h | {celltype}",
            "regime": "persistent",
            "input": str(evidence.relative_to(ROOT)).replace("\\", "/"),
            "output_dir": str(outdir.relative_to(ROOT)).replace("\\", "/"),
            "gene_column": "gene",
            "latest_column": "r8",
            "history_columns": ["r8_repA", "r8_repB"],
            "filter": {"column": "technical_universe", "equals": True},
            "top_k": [100, 250, 500, 1000],
        }
        cfg_path = cfg_dir / f"{slug}.json"
        cfg_path.write_text(json.dumps(cfg, indent=2, ensure_ascii=False), encoding="utf-8")
        subprocess.run([sys.executable, str(CORE), "--config", str(cfg_path)], cwd=ROOT, check=True)

        ranking = outdir / "ranking.tsv.gz"
        manifest = outdir / "manifest.json"
        if not ranking.exists() or not manifest.exists():
            raise RuntimeError(f"Unified PlantVC output missing for {celltype}")
        m = json.loads(manifest.read_text(encoding="utf-8"))
        if m.get("status") != "RANKED" or m.get("regime") != "persistent":
            raise RuntimeError(f"Unexpected PlantVC manifest for {celltype}: {m}")
        ntech = int(ev["technical_universe"].astype(bool).sum())
        if int(m["n_ranked_genes"]) != ntech:
            raise RuntimeError(f"Ranked gene count mismatch for {celltype}")

        records.append({
            "celltype": celltype,
            "slug": slug,
            "n_genes_all": len(ev),
            "n_technical_genes": ntech,
            "evidence_sha256": sha(evidence),
            "config_sha256": sha(cfg_path),
            "ranking_sha256": sha(ranking),
        })

    audit = pd.DataFrame(records)
    audit.to_csv(OUT / "freeze_context_audit.tsv", sep="\t", index=False)

    manifest = {
        "status": "FROZEN_BEFORE_20H_TARGET",
        "source_rds_sha256": sha(SOURCE),
        "rules_sha256": sha(RULES),
        "core_sha256": sha(CORE),
        "eligible_celltypes": [x[0] for x in ELIGIBLE],
        "pretarget_samples_used": [1, 2, 5, 6],
        "target_samples_withheld": [3, 7],
        "n_freeze_cells_used": int(len(used)),
        "n_target_cells_withheld": int(len(withheld)),
        "n_target_cells_used": 0,
        "primary_rule": "persistent: latest same-celltype 8h-vs-0h response predicts 20h-vs-0h ranking",
        "frozen_comparator": "leave-one-celltype-out mean 8h response across the other eight eligible cell types",
        "target_values_read": False,
        "contexts": records,
    }
    (OUT / "FREEZE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    print("GSE158761_FREEZE_AUDIT_SUCCESS")
    print(json.dumps({k: v for k, v in manifest.items() if k != "contexts"}, indent=2, ensure_ascii=False))
    print("\nCONTEXTS")
    print(audit[["celltype", "n_technical_genes", "ranking_sha256"]].to_string(index=False))
    print("\n20H TARGET VALUES HAVE NOT BEEN READ")


if __name__ == "__main__":
    main()
