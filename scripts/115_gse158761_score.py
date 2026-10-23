from __future__ import annotations

from pathlib import Path
import hashlib
import json
import math

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, hypergeom

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_lateral_root" / "GSE158761"
OUT = BASE / "plantvc_v1_blind"
RULES = ROOT / "PLANTVC_V1_GSE158761_BLIND_RULES_2026-10-20.md"
SOURCE = BASE / "GSE158761_all_cells_labeled.rds"
CORE = ROOT / "scripts" / "108_plantvc_v1.py"
KS = [100, 250, 500, 1000]

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


def safe_rho(a: np.ndarray, b: np.ndarray) -> float:
    x = spearmanr(a, b).statistic
    return float(x) if np.isfinite(x) else np.nan


def deterministic_order(score: np.ndarray, genes: np.ndarray) -> np.ndarray:
    # Primary descending score, deterministic ascending gene tie-break.
    return np.lexsort((genes.astype(str), -score))


def topk_metrics(score: np.ndarray, signed: np.ndarray, target: np.ndarray, genes: np.ndarray, k: int) -> dict:
    n = len(score)
    k = min(k, n)
    pred_order = deterministic_order(score, genes)
    targ_order = deterministic_order(np.abs(target), genes)
    pi = pred_order[:k]
    ti = targ_order[:k]
    inter = np.intersect1d(pi, ti)
    ov = len(inter)
    exp = k * k / n
    ranks = np.empty(n, int)
    ranks[pred_order] = np.arange(1, n + 1)
    pct = 1 - (ranks[ti] - 1) / max(n - 1, 1)
    return {
        "k": k,
        "overlap": ov,
        "random_expected_overlap": exp,
        "enrichment": ov / exp if exp else np.nan,
        "hypergeom_p": float(hypergeom.sf(ov - 1, n, k, k)),
        "sign_agreement_overlap": float(np.mean(np.sign(signed[inter]) == np.sign(target[inter]))) if ov else np.nan,
        "median_target_topk_predicted_percentile": float(np.median(pct)),
    }


def main() -> None:
    freeze_manifest_path = OUT / "FREEZE_MANIFEST.json"
    if not freeze_manifest_path.exists():
        raise SystemExit("Freeze manifest missing")
    fm = json.loads(freeze_manifest_path.read_text(encoding="utf-8"))
    if fm.get("status") != "FROZEN_BEFORE_20H_TARGET" or fm.get("target_values_read") is not False:
        raise RuntimeError("Freeze manifest is not a clean pretarget freeze")
    if sha(SOURCE) != fm["source_rds_sha256"]:
        raise RuntimeError("Source RDS hash changed after freeze")
    if sha(RULES) != fm["rules_sha256"]:
        raise RuntimeError("Frozen rules changed after freeze")
    if sha(CORE) != fm["core_sha256"]:
        raise RuntimeError("Unified PlantVC core changed after freeze")

    audit = pd.read_csv(OUT / "freeze_context_audit.tsv", sep="\t").set_index("slug")
    vector_rows = []
    top_rows = []
    rep_rows = []
    comparison_rows = []

    for celltype, slug in ELIGIBLE:
        evidence_path = OUT / f"evidence_{slug}.tsv.gz"
        target_path = OUT / f"target_{slug}.tsv.gz"
        ranking_path = OUT / "rankings" / slug / "ranking.tsv.gz"
        if not target_path.exists():
            raise RuntimeError(f"Target file missing for {celltype}")
        if sha(evidence_path) != audit.loc[slug, "evidence_sha256"]:
            raise RuntimeError(f"Evidence changed after freeze for {celltype}")
        if sha(ranking_path) != audit.loc[slug, "ranking_sha256"]:
            raise RuntimeError(f"Ranking changed after freeze for {celltype}")

        ev = pd.read_csv(evidence_path, sep="\t")
        ta = pd.read_csv(target_path, sep="\t")
        rk = pd.read_csv(ranking_path, sep="\t")
        if ev["gene"].astype(str).tolist() != ta["gene"].astype(str).tolist():
            raise RuntimeError(f"Gene order mismatch for {celltype}")

        tech = ev["technical_universe"].astype(bool).to_numpy()
        genes = ev.loc[tech, "gene"].astype(str).to_numpy()
        target = ta.loc[tech, "r20"].to_numpy(float)
        repA = ta.loc[tech, "r20_repA"].to_numpy(float)
        repB = ta.loc[tech, "r20_repB"].to_numpy(float)

        # PlantVC primary prediction comes from the frozen unified entrypoint output.
        rk2 = rk.set_index(rk["gene"].astype(str))
        if set(rk2.index) != set(genes):
            raise RuntimeError(f"Unified ranking technical universe mismatch for {celltype}")
        primary_signed = rk2.loc[genes, "predicted_response"].to_numpy(float)
        primary_score = rk2.loc[genes, "priority_score"].to_numpy(float)

        # Frozen comparator: leave-one-celltype-out shared 8h response.
        broad_signed = ev.loc[tech, "broad_shared_r8"].to_numpy(float)
        broad_score = np.abs(broad_signed)

        models = {
            "PlantVC_persistent_same_state": (primary_score, primary_signed),
            "frozen_broad_shared": (broad_score, broad_signed),
        }

        per_model_top100 = {}
        for model, (score, signed) in models.items():
            vector_rows.append({
                "celltype": celltype,
                "model": model,
                "n_genes": len(genes),
                "spearman_full_vector": safe_rho(signed, target),
                "sign_agreement_all": float(np.mean(np.sign(signed) == np.sign(target))),
            })
            for k in KS:
                row = topk_metrics(score, signed, target, genes, k)
                top_rows.append({"celltype": celltype, "model": model, **row})
                if k == 100:
                    per_model_top100[model] = row

        rep_rows.append({
            "celltype": celltype,
            "metric": "full_vector",
            "replicate_spearman": safe_rho(repA, repB),
            "replicate_sign_agreement": float(np.mean(np.sign(repA) == np.sign(repB))),
        })
        for k in KS:
            rr = topk_metrics(np.abs(repA), repA, repB, genes, k)
            rep_rows.append({"celltype": celltype, "metric": f"top{k}", **rr})

        p = next(r for r in vector_rows if r["celltype"] == celltype and r["model"] == "PlantVC_persistent_same_state")
        b = next(r for r in vector_rows if r["celltype"] == celltype and r["model"] == "frozen_broad_shared")
        comparison_rows.append({
            "celltype": celltype,
            "delta_rho_primary_minus_broad": p["spearman_full_vector"] - b["spearman_full_vector"],
            "delta_top100_overlap_primary_minus_broad": per_model_top100["PlantVC_persistent_same_state"]["overlap"] - per_model_top100["frozen_broad_shared"]["overlap"],
            "delta_top100_enrichment_primary_minus_broad": per_model_top100["PlantVC_persistent_same_state"]["enrichment"] - per_model_top100["frozen_broad_shared"]["enrichment"],
        })

    vectors = pd.DataFrame(vector_rows)
    tops = pd.DataFrame(top_rows)
    reps = pd.DataFrame(rep_rows)
    comps = pd.DataFrame(comparison_rows)
    vectors.to_csv(OUT / "blind_vector_results.tsv", sep="\t", index=False)
    tops.to_csv(OUT / "blind_topk_results.tsv", sep="\t", index=False)
    reps.to_csv(OUT / "target_replicate_reference.tsv", sep="\t", index=False)
    comps.to_csv(OUT / "primary_vs_broad.tsv", sep="\t", index=False)

    summary_rows = []
    for model in vectors["model"].unique():
        vv = vectors[vectors["model"] == model]
        row = {
            "model": model,
            "n_contexts": len(vv),
            "median_full_vector_rho": vv["spearman_full_vector"].median(),
            "median_sign_agreement_all": vv["sign_agreement_all"].median(),
        }
        for k in KS:
            tt = tops[(tops["model"] == model) & (tops["k"] == k)]
            row[f"median_top{k}_overlap"] = tt["overlap"].median()
            row[f"median_top{k}_enrichment"] = tt["enrichment"].median()
            row[f"median_top{k}_direction_agreement"] = tt["sign_agreement_overlap"].median()
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "blind_summary.tsv", sep="\t", index=False)

    win_summary = {
        "same_state_wins_rho": int((comps["delta_rho_primary_minus_broad"] > 0).sum()),
        "same_state_ties_rho": int((comps["delta_rho_primary_minus_broad"] == 0).sum()),
        "same_state_wins_top100_overlap": int((comps["delta_top100_overlap_primary_minus_broad"] > 0).sum()),
        "same_state_ties_top100_overlap": int((comps["delta_top100_overlap_primary_minus_broad"] == 0).sum()),
        "n_contexts": len(comps),
    }
    (OUT / "BLIND_RESULT.json").write_text(json.dumps({"summary": summary.to_dict("records"), "comparison": win_summary}, indent=2), encoding="utf-8")

    print("GSE158761_BLIND_SCORE_SUCCESS")
    print("\nSUMMARY")
    print(summary.round(4).to_string(index=False))
    print("\nPRIMARY VS BROAD")
    print(comps.round(4).to_string(index=False))
    print("\nWIN SUMMARY")
    print(json.dumps(win_summary, indent=2))


if __name__ == "__main__":
    main()
