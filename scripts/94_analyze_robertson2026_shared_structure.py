from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026" / "pseudobulk_v0"
OUT = ROOT / "data" / "external_rice" / "Robertson2026" / "shared_structure_v0"
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(BASE / "group_metadata.tsv", sep="\t")
expr = pd.read_csv(BASE / "pseudobulk_log1p_cpm.tsv.gz", sep="\t", compression="gzip")
design = pd.read_csv(ROOT / "data" / "external_rice" / "Robertson2026" / "sample_design.tsv", sep="\t")
genes = expr.pop("gene").astype(str).to_numpy()
X = expr.to_numpy(dtype=float)  # genes x groups
cols = expr.columns.astype(str).tolist()
if cols != meta["group"].astype(str).tolist():
    raise RuntimeError("Expression columns do not match group_metadata order")

# Gene-universe filter is purely technical/descriptive: retain genes expressed in a reasonable
# fraction of pseudobulk groups. We report results for both all genes and the filtered set.
mean_expr = X.mean(axis=1)
nonzero_frac = (X > 0).mean(axis=1)
mask = (nonzero_frac >= 0.25) & (mean_expr >= 0.05)

conditions = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
control = "Well.Watered"
annotations = sorted(meta["cell_annotation"].unique())

# Use the GEO sample design to preserve the biological replicate / batch structure.
# Every stress replicate is contrasted against the well-watered replicate from the
# same experimental batch, then the matched response vectors are equal-weight averaged.
batch_by_sample = dict(zip(design["identity"].astype(str), design["batch"].astype(int)))
if set(meta["sample"].astype(str)) - set(batch_by_sample):
    missing = sorted(set(meta["sample"].astype(str)) - set(batch_by_sample))
    raise RuntimeError(f"Samples missing from sample_design.tsv: {missing}")

control_sample_by_batch = {}
for sample in meta.loc[meta.condition == control, "sample"].astype(str).unique():
    b = batch_by_sample[sample]
    if b in control_sample_by_batch:
        raise RuntimeError(f"Multiple well-watered samples for batch {b}")
    control_sample_by_batch[b] = sample

# Build batch-matched stress-minus-control response vectors per CellAnnotation.
def matched_responses(condition: str, annotation: str):
    t_idx = meta.index[(meta.condition == condition) & (meta.cell_annotation == annotation)].to_numpy()
    vecs = []
    pairs = []
    for ti in t_idx:
        stress_sample = str(meta.loc[ti, "sample"])
        batch = batch_by_sample[stress_sample]
        control_sample = control_sample_by_batch.get(batch)
        if control_sample is None:
            continue
        c_idx = meta.index[(meta["sample"].astype(str) == control_sample) &
                           (meta.cell_annotation == annotation)].to_numpy()
        if len(c_idx) != 1:
            continue
        vecs.append(X[:, ti] - X[:, c_idx[0]])
        pairs.append((stress_sample, control_sample, batch))
    return vecs, pairs

summary_rows = []
state_rows = []
pair_rows = []
response_matrices = {}

for cond in conditions:
    states = []
    responses = []
    for ann in annotations:
        vecs, _ = matched_responses(cond, ann)
        if len(vecs) == 0:
            continue
        states.append(ann)
        responses.append(np.vstack(vecs).mean(axis=0))
    R = np.vstack(responses)  # states x genes
    response_matrices[cond] = (states, R)

    for label, usemask in [("all", np.ones(len(genes), dtype=bool)), ("filtered", mask)]:
        M = R[:, usemask]
        u, s, vt = np.linalg.svd(M, full_matrices=False)
        energy = s**2
        frac = energy / energy.sum()
        shared = vt[0]
        # orient shared vector to align with mean response
        if np.dot(shared, M.mean(axis=0)) < 0:
            shared = -shared
        shared_unit = shared / np.linalg.norm(shared)

        # Pairwise correlations between cell-state responses
        prs = []
        for i in range(len(states)):
            for j in range(i + 1, len(states)):
                rho = spearmanr(M[i], M[j]).statistic
                prs.append(rho)
                if label == "filtered":
                    pair_rows.append({"condition": cond, "state_a": states[i], "state_b": states[j], "rho": rho})

        loo_rhos = []
        proj_fracs = []
        if label == "filtered":
            for i, st in enumerate(states):
                loo = np.delete(M, i, axis=0).mean(axis=0)
                rho = spearmanr(M[i], loo).statistic
                denom = float(np.dot(M[i], M[i]))
                proj = float(np.dot(M[i], shared_unit))
                proj_frac = (proj**2 / denom) if denom > 0 else np.nan
                loo_rhos.append(rho)
                proj_fracs.append(proj_frac)
                state_rows.append({
                    "condition": cond,
                    "cell_annotation": st,
                    "rho_to_LOO_shared": rho,
                    "projection_energy_fraction": proj_frac,
                    "response_rms": float(np.sqrt(np.mean(M[i] ** 2))),
                })

        summary_rows.append({
            "condition": cond,
            "gene_universe": label,
            "n_states": len(states),
            "n_genes": int(usemask.sum()),
            "rank1_energy_fraction": float(frac[0]),
            "rank2_cumulative_energy_fraction": float(frac[:2].sum()),
            "rank3_cumulative_energy_fraction": float(frac[:3].sum()),
            "median_pairwise_state_spearman": float(np.nanmedian(prs)),
            "mean_pairwise_state_spearman": float(np.nanmean(prs)),
            "median_state_rho_to_LOO_shared": float(np.nanmedian(loo_rhos)) if loo_rhos else np.nan,
            "mean_state_rho_to_LOO_shared": float(np.nanmean(loo_rhos)) if loo_rhos else np.nan,
            "median_projection_energy_fraction": float(np.nanmedian(proj_fracs)) if proj_fracs else np.nan,
            "mean_projection_energy_fraction": float(np.nanmean(proj_fracs)) if proj_fracs else np.nan,
        })

summary = pd.DataFrame(summary_rows)
states_df = pd.DataFrame(state_rows)
pairs_df = pd.DataFrame(pair_rows)
summary.to_csv(OUT / "condition_shared_structure_summary.tsv", sep="\t", index=False)
states_df.to_csv(OUT / "state_shared_structure_metrics.tsv", sep="\t", index=False)
pairs_df.to_csv(OUT / "pairwise_state_response_spearman.tsv", sep="\t", index=False)

# Stress-response similarity across conditions on the common filtered gene universe.
# Compare mean response vectors across cell states, not individual-cell annotations.
common_rows = []
for i, a in enumerate(conditions):
    sa, Ra = response_matrices[a]
    for b in conditions[i + 1:]:
        sb, Rb = response_matrices[b]
        common = sorted(set(sa).intersection(sb))
        ia = [sa.index(x) for x in common]
        ib = [sb.index(x) for x in common]
        ma = Ra[np.array(ia)][:, mask].mean(axis=0)
        mb = Rb[np.array(ib)][:, mask].mean(axis=0)
        common_rows.append({
            "condition_a": a,
            "condition_b": b,
            "n_common_states": len(common),
            "mean_response_spearman": float(spearmanr(ma, mb).statistic),
        })
pd.DataFrame(common_rows).to_csv(OUT / "condition_mean_response_similarity.tsv", sep="\t", index=False)

# Replicate reproducibility per stress x annotation using the same batch-matched
# stress-minus-control contrasts used in the primary descriptive decomposition.
rep_rows = []
for cond in conditions:
    for ann in annotations:
        vecs, labels = matched_responses(cond, ann)
        if len(vecs) < 2:
            continue
        vals = []
        for i in range(len(vecs)):
            for j in range(i + 1, len(vecs)):
                vals.append(spearmanr(vecs[i][mask], vecs[j][mask]).statistic)
        if vals:
            rep_rows.append({
                "condition": cond,
                "cell_annotation": ann,
                "n_matched_response_vectors": len(vecs),
                "matched_pairs": " | ".join(f"{s}-{c}:batch{b}" for s, c, b in labels),
                "n_pairwise_response_comparisons": len(vals),
                "median_replicate_response_spearman": float(np.nanmedian(vals)),
                "mean_replicate_response_spearman": float(np.nanmean(vals)),
            })
rep_df = pd.DataFrame(rep_rows)
rep_df.to_csv(OUT / "replicate_response_stability.tsv", sep="\t", index=False)

# Compact report, deliberately descriptive rather than thresholded PASS/FAIL.
filt = summary[summary.gene_universe == "filtered"].copy()
report = [
    "# Robertson 2026 rice stress shared-vs-state-specific structure",
    "",
    "Status: external descriptive replication / contrast. No biological hard threshold was preregistered.",
    "",
    "Response construction: each stress replicate is contrasted with the well-watered replicate from the same GEO batch; matched response vectors are then equal-weight averaged within condition x CellAnnotation.",
    "",
    f"Technical gene universe: {int(mask.sum())} / {len(genes)} genes (nonzero in >=25% groups and mean log1p(CPM) >=0.05).",
    "",
    "## Shared structure by condition",
]
for _, r in filt.iterrows():
    report += [
        f"- {r.condition}: rank1={r.rank1_energy_fraction:.3f}, rank1-3={r.rank3_cumulative_energy_fraction:.3f}, "
        f"median pairwise rho={r.median_pairwise_state_spearman:.3f}, median LOO-shared rho={r.median_state_rho_to_LOO_shared:.3f}, "
        f"median projection energy={r.median_projection_energy_fraction:.3f}."
    ]
if not rep_df.empty:
    report += ["", "## Replicate stability"]
    rr = rep_df.groupby("condition")["median_replicate_response_spearman"].median()
    for cond, val in rr.items():
        report.append(f"- {cond}: median across annotations={val:.3f}")
report += [
    "",
    "## Interpretation rule",
    "These values are read continuously. A weaker shared component is not a failure; it indicates more cell-state-specific response architecture. Stronger or weaker sharing across drought/salinity and severity is itself the result.",
]
(OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

payload = {
    "n_genes_total": int(len(genes)),
    "n_genes_filtered": int(mask.sum()),
    "n_annotations": int(len(annotations)),
    "conditions": conditions,
}
(OUT / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
print("ROBERTSON2026_SHARED_STRUCTURE_SUCCESS")
print(filt.to_string(index=False))
