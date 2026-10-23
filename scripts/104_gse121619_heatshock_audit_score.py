from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_heatshock" / "GSE121619"
OUT = BASE / "plantvc_v1_heatshock_blind"
RULES = ROOT / "PLANTVC_V1_HEATSHOCK_BLIND_RULES_2026-10-20.md"
SOURCE = BASE / "GSE121619_Control_Heatshock_cds.rds"
META = BASE / "GSE121619_Control_Heatshock_pData.tsv.gz"
TARGETS = {2, 3, 10}
DONORS = {5, 6, 8}
KS = [100, 250, 500, 1000]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def freeze_audit() -> None:
    pred_path = OUT / "frozen_prediction.tsv.gz"
    used_path = OUT / "freeze_cells_used.tsv"
    withheld_path = OUT / "target_heatshock_cells_withheld.tsv"
    tech_path = OUT / "technical_genes.txt"
    log_path = OUT / "freeze.log"
    for p in [pred_path, used_path, withheld_path, tech_path, log_path, RULES, SOURCE, META]:
        if not p.exists():
            raise RuntimeError(f"Missing freeze artifact: {p}")

    meta = pd.read_csv(META, sep="\t", compression="gzip")
    used = pd.read_csv(used_path, sep="\t")
    withheld = pd.read_csv(withheld_path, sep="\t")
    pred = pd.read_csv(pred_path, sep="\t", compression="gzip")

    leaked = used[(used["treatment"] == "HeatShock") & (used["cluster"].isin(TARGETS))]
    if len(leaked):
        raise RuntimeError(f"Target HeatShock leakage: {len(leaked)} cells")

    used_hs_clusters = set(used.loc[used.treatment == "HeatShock", "cluster"].astype(int))
    if used_hs_clusters != DONORS:
        raise RuntimeError(f"Unexpected HeatShock donor clusters in freeze: {used_hs_clusters}")

    expected_withheld = meta[(meta.treatment_id == "HeatShock") & (meta.louvain_component.isin(TARGETS))]
    if set(withheld["cell"]) != set(expected_withheld["cell"]):
        raise RuntimeError("Withheld target HeatShock cell set does not match metadata")

    if pred["gene"].duplicated().any():
        raise RuntimeError("Duplicate genes in frozen prediction")
    if pred["shared_prediction"].isna().any():
        raise RuntimeError("Missing values in shared prediction")

    log = log_path.read_text(encoding="utf-8", errors="replace")
    if "GSE121619_HEATSHOCK_FREEZE_SUCCESS" not in log:
        raise RuntimeError("Freeze success marker missing")
    if "TARGET_HEATSHOCK_VALUES_NOT_USED_FOR_ANY_PREDICTION" not in log:
        raise RuntimeError("Target firewall marker missing")

    manifest = {
        "status": "FROZEN_BEFORE_TARGET",
        "source_rds_sha256": sha256(SOURCE),
        "rules_sha256": sha256(RULES),
        "prediction_sha256": sha256(pred_path),
        "n_genes_all": int(len(pred)),
        "n_genes_technical": int(pred["technical"].sum()),
        "target_clusters_withheld": sorted(TARGETS),
        "donor_clusters_used": sorted(DONORS),
        "n_freeze_cells_used": int(len(used)),
        "n_heatshock_donor_cells_used": int((used.treatment == "HeatShock").sum()),
        "n_target_heatshock_cells_used": 0,
        "n_target_heatshock_cells_withheld": int(len(withheld)),
        "withheld_counts_by_cluster": {
            str(k): int(v)
            for k, v in withheld.groupby("cluster").size().sort_index().items()
        },
        "primary_prediction": "equal-weight mean HeatShock-Control log1p(CPM) response across donor clusters 5,6,8",
        "target_values_read": False,
    }
    (OUT / "FREEZE_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("GSE121619_HEATSHOCK_FREEZE_AUDIT_SUCCESS")
    print(json.dumps(manifest, indent=2))


def ordered_topk(genes: np.ndarray, values: np.ndarray, k: int) -> np.ndarray:
    # Deterministic: descending absolute value, then lexical gene ID.
    order = np.lexsort((genes.astype(str), -np.abs(values)))
    return order[: min(k, len(order))]


def score_one(genes: np.ndarray, pred: np.ndarray, target: np.ndarray, k: int) -> dict:
    n = len(genes)
    pi = ordered_topk(genes, pred, k)
    ti = ordered_topk(genes, target, k)
    pset = set(genes[pi])
    tset = set(genes[ti])
    overlap_genes = pset & tset
    overlap = len(overlap_genes)
    kk = min(k, n)
    expected = kk * kk / n
    enrich = overlap / expected if expected > 0 else np.nan
    pval = float(hypergeom.sf(overlap - 1, n, kk, kk))

    g_to_i = {g: i for i, g in enumerate(genes)}
    ov_idx = np.array([g_to_i[g] for g in overlap_genes], dtype=int)
    sign_agree = float(np.mean(np.sign(pred[ov_idx]) == np.sign(target[ov_idx]))) if len(ov_idx) else np.nan

    pred_order = np.lexsort((genes.astype(str), -np.abs(pred)))
    ranks = np.empty(n, dtype=int)
    ranks[pred_order] = np.arange(1, n + 1)
    target_top_ranks = ranks[ti]
    percentiles = 1.0 - (target_top_ranks - 1) / max(1, n - 1)

    return {
        "k": int(kk),
        "overlap": int(overlap),
        "random_expected_overlap": float(expected),
        "enrichment": float(enrich),
        "hypergeom_p": pval,
        "sign_agreement_overlap": sign_agree,
        "median_target_topk_predicted_percentile": float(np.median(percentiles)),
    }


def score() -> None:
    manifest_path = OUT / "FREEZE_MANIFEST.json"
    pred_path = OUT / "frozen_prediction.tsv.gz"
    target_path = OUT / "target_profiles.tsv.gz"
    if not manifest_path.exists():
        raise RuntimeError("Freeze manifest missing; target cannot be scored")
    if not target_path.exists():
        raise RuntimeError("Target export missing; run R script in target mode first")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    current_pred_hash = sha256(pred_path)
    if current_pred_hash != manifest["prediction_sha256"]:
        raise RuntimeError("Frozen prediction hash changed after target reveal")
    if sha256(RULES) != manifest["rules_sha256"]:
        raise RuntimeError("Frozen rule file changed after target reveal")

    pred = pd.read_csv(pred_path, sep="\t", compression="gzip")
    targ = pd.read_csv(target_path, sep="\t", compression="gzip")
    models = {
        "shared_v1": "shared_prediction",
        "pooled_donor": "pooled_donor_prediction",
        "donor_5": "donor_5",
        "donor_6": "donor_6",
        "donor_8": "donor_8",
    }

    metric_rows = []
    topk_rows = []
    for cluster in sorted(TARGETS):
        t = targ[targ.target_cluster == cluster][["gene", "target_response", "technical"]].copy()
        m = pred.merge(t, on="gene", how="inner", suffixes=("", "_target"))
        if len(m) != len(pred):
            raise RuntimeError(f"Gene mismatch for target cluster {cluster}")

        for universe in ["technical", "all"]:
            use = m["technical_target"].astype(bool).to_numpy() if universe == "technical" else np.ones(len(m), dtype=bool)
            mm = m.loc[use].copy()
            genes = mm["gene"].astype(str).to_numpy()
            y = mm["target_response"].to_numpy(float)
            for model, col in models.items():
                x = mm[col].to_numpy(float)
                rho = float(spearmanr(x, y).statistic)
                sign = float(np.mean(np.sign(x) == np.sign(y)))
                metric_rows.append({
                    "target_cluster": cluster,
                    "universe": universe,
                    "model": model,
                    "n_genes": len(mm),
                    "spearman_full_vector": rho,
                    "sign_agreement_all": sign,
                })
                for k in KS:
                    r = score_one(genes, x, y, k)
                    r.update({"target_cluster": cluster, "universe": universe, "model": model})
                    topk_rows.append(r)

    metrics = pd.DataFrame(metric_rows)
    topk = pd.DataFrame(topk_rows)
    metrics.to_csv(OUT / "blind_full_vector_metrics.tsv", sep="\t", index=False)
    topk.to_csv(OUT / "blind_topk_metrics.tsv", sep="\t", index=False)

    tech_shared = metrics[(metrics.universe == "technical") & (metrics.model == "shared_v1")]
    tech_top = topk[(topk.universe == "technical") & (topk.model == "shared_v1")]
    agg = tech_top.groupby("k").agg(
        median_overlap=("overlap", "median"),
        median_enrichment=("enrichment", "median"),
        median_sign_agreement_overlap=("sign_agreement_overlap", "median"),
        median_target_topk_predicted_percentile=("median_target_topk_predicted_percentile", "median"),
        max_hypergeom_p=("hypergeom_p", "max"),
    ).reset_index()
    agg.to_csv(OUT / "blind_shared_v1_aggregate.tsv", sep="\t", index=False)

    report = [
        "# PlantVC v1 acute heat-shock blind result",
        "",
        "Target states: clusters 2, 3, 10. Their HeatShock expression was excluded from the frozen prediction.",
        "Prediction: equal-weight mean response of donor states 5, 6, 8.",
        "",
        "## Full-vector descriptive result",
    ]
    for _, r in tech_shared.sort_values("target_cluster").iterrows():
        report.append(f"- cluster {int(r.target_cluster)}: rho={r.spearman_full_vector:.3f}, sign agreement={r.sign_agreement_all:.3f}")
    report += ["", "## Candidate-ranking utility"]
    for _, r in agg.iterrows():
        report.append(
            f"- top{int(r.k)}: median overlap={r.median_overlap:.1f}, median enrichment={r.median_enrichment:.2f}x, "
            f"median direction agreement={r.median_sign_agreement_overlap:.3f}, "
            f"median target percentile={r.median_target_topk_predicted_percentile:.3f}."
        )
    report += [
        "",
        "Interpretation is continuous. This benchmark tests acute cross-state candidate prioritization, not exact full-state simulation.",
    ]
    (OUT / "BLIND_RESULT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    print("GSE121619_HEATSHOCK_BLIND_SCORE_SUCCESS")
    print("\nTECHNICAL FULL VECTOR")
    print(metrics[metrics.universe.eq("technical")].to_string(index=False))
    print("\nSHARED V1 TECHNICAL TOP-K")
    print(tech_top.to_string(index=False))
    print("\nSHARED V1 AGGREGATE")
    print(agg.to_string(index=False))


def main() -> None:
    if len(sys.argv) != 2 or sys.argv[1] not in {"freeze", "score"}:
        raise SystemExit("Usage: python scripts/104_gse121619_heatshock_audit_score.py freeze|score")
    if sys.argv[1] == "freeze":
        freeze_audit()
    else:
        score()


if __name__ == "__main__":
    main()
