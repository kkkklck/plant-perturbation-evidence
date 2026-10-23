from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026"
PB = BASE / "pseudobulk_v0"
OUT = BASE / "shared_structure_batch_matched_v0"
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(PB / "group_metadata.tsv", sep="\t")
design = pd.read_csv(BASE / "sample_design.tsv", sep="\t")
expr = pd.read_csv(PB / "pseudobulk_log1p_cpm.tsv.gz", sep="\t", compression="gzip")
genes = expr.pop("gene").astype(str).to_numpy()
X = expr.to_numpy(float)  # genes x groups
if expr.columns.astype(str).tolist() != meta["group"].astype(str).tolist():
    raise RuntimeError("Expression columns do not match group metadata")

# technical gene universe; same as unpaired descriptive analysis
mean_expr = X.mean(axis=1)
nonzero_frac = (X > 0).mean(axis=1)
mask = (nonzero_frac >= 0.25) & (mean_expr >= 0.05)

# normalize sample names / construct exact batch matched control map
control_rows = design[design.stress_type == "control"].copy()
control_by_batch = dict(zip(control_rows.batch.astype(int), control_rows.identity.astype(str)))

stress_conditions = [
    ("Mild.Drought", "MildDrought"),
    ("Moderate.Drought", "ModerateDrought"),
    ("Mild.Salinity", "MildSalinity"),
    ("Moderate.Salinity", "ModerateSalinity"),
]
annotations = sorted(meta.cell_annotation.unique())

# quick exact lookup group -> column index
col_index = {g: i for i, g in enumerate(meta["group"].astype(str))}

def vec(sample: str, annotation: str) -> np.ndarray:
    key = f"{sample}||{annotation}"
    if key not in col_index:
        raise KeyError(key)
    return X[:, col_index[key]]

summary_rows = []
state_rows = []
rep_rows = []
pair_rows = []
mean_response_by_condition = {}

for cond_label, design_condition in stress_conditions:
    stress_samples = design[design.condition == design_condition].copy().sort_values("batch")
    if len(stress_samples) != 2:
        raise RuntimeError(f"Expected 2 stress samples for {design_condition}, got {len(stress_samples)}")

    states = []
    avg_responses = []
    rep_responses_by_state = {}

    for ann in annotations:
        paired = []
        pair_labels = []
        for _, row in stress_samples.iterrows():
            stress_sample = str(row.identity)
            batch = int(row.batch)
            control_sample = control_by_batch[batch]
            delta = vec(stress_sample, ann) - vec(control_sample, ann)
            paired.append(delta)
            pair_labels.append((stress_sample, control_sample, batch))
        paired = np.vstack(paired)  # 2 x genes
        avg = paired.mean(axis=0)
        states.append(ann)
        avg_responses.append(avg)
        rep_responses_by_state[ann] = paired
        rho_rep = spearmanr(paired[0, mask], paired[1, mask]).statistic
        rep_rows.append({
            "condition": cond_label,
            "cell_annotation": ann,
            "stress_sample_1": pair_labels[0][0],
            "control_sample_1": pair_labels[0][1],
            "batch_1": pair_labels[0][2],
            "stress_sample_2": pair_labels[1][0],
            "control_sample_2": pair_labels[1][1],
            "batch_2": pair_labels[1][2],
            "matched_replicate_response_spearman": rho_rep,
        })

    R = np.vstack(avg_responses)  # states x genes
    mean_response_by_condition[cond_label] = R[:, mask].mean(axis=0)

    for universe, usemask in [("all", np.ones(len(genes), dtype=bool)), ("filtered", mask)]:
        M = R[:, usemask]
        u, s, vt = np.linalg.svd(M, full_matrices=False)
        frac = (s**2) / np.sum(s**2)
        shared = vt[0].copy()
        if np.dot(shared, M.mean(axis=0)) < 0:
            shared *= -1
        shared /= np.linalg.norm(shared)

        pair_rhos = []
        loo_rhos = []
        proj_fracs = []
        for i in range(len(states)):
            for j in range(i + 1, len(states)):
                rho = spearmanr(M[i], M[j]).statistic
                pair_rhos.append(rho)
                if universe == "filtered":
                    pair_rows.append({"condition": cond_label, "state_a": states[i], "state_b": states[j], "rho": rho})
        if universe == "filtered":
            for i, st in enumerate(states):
                loo = np.delete(M, i, axis=0).mean(axis=0)
                rho = spearmanr(M[i], loo).statistic
                denom = float(np.dot(M[i], M[i]))
                p = float(np.dot(M[i], shared))
                pf = (p*p/denom) if denom > 0 else np.nan
                loo_rhos.append(rho)
                proj_fracs.append(pf)
                state_rows.append({
                    "condition": cond_label,
                    "cell_annotation": st,
                    "rho_to_LOO_shared": rho,
                    "projection_energy_fraction": pf,
                    "response_rms": float(np.sqrt(np.mean(M[i]**2))),
                })

        summary_rows.append({
            "condition": cond_label,
            "gene_universe": universe,
            "n_states": len(states),
            "n_genes": int(usemask.sum()),
            "rank1_energy_fraction": float(frac[0]),
            "rank2_cumulative_energy_fraction": float(frac[:2].sum()),
            "rank3_cumulative_energy_fraction": float(frac[:3].sum()),
            "median_pairwise_state_spearman": float(np.nanmedian(pair_rhos)),
            "mean_pairwise_state_spearman": float(np.nanmean(pair_rhos)),
            "median_state_rho_to_LOO_shared": float(np.nanmedian(loo_rhos)) if loo_rhos else np.nan,
            "mean_state_rho_to_LOO_shared": float(np.nanmean(loo_rhos)) if loo_rhos else np.nan,
            "median_projection_energy_fraction": float(np.nanmedian(proj_fracs)) if proj_fracs else np.nan,
            "mean_projection_energy_fraction": float(np.nanmean(proj_fracs)) if proj_fracs else np.nan,
        })

summary = pd.DataFrame(summary_rows)
states_df = pd.DataFrame(state_rows)
rep_df = pd.DataFrame(rep_rows)
pair_df = pd.DataFrame(pair_rows)
summary.to_csv(OUT / "condition_shared_structure_summary.tsv", sep="\t", index=False)
states_df.to_csv(OUT / "state_shared_structure_metrics.tsv", sep="\t", index=False)
rep_df.to_csv(OUT / "matched_replicate_response_stability.tsv", sep="\t", index=False)
pair_df.to_csv(OUT / "pairwise_state_response_spearman.tsv", sep="\t", index=False)

# Mean-response similarity between stress conditions
sim_rows = []
labels = [x[0] for x in stress_conditions]
for i, a in enumerate(labels):
    for b in labels[i+1:]:
        sim_rows.append({
            "condition_a": a,
            "condition_b": b,
            "mean_response_spearman": float(spearmanr(mean_response_by_condition[a], mean_response_by_condition[b]).statistic),
        })
pd.DataFrame(sim_rows).to_csv(OUT / "condition_mean_response_similarity.tsv", sep="\t", index=False)

# Compare paired result with previous all-control-average sensitivity result
old_path = BASE / "shared_structure_v0" / "condition_shared_structure_summary.tsv"
comparison = None
if old_path.exists():
    old = pd.read_csv(old_path, sep="\t")
    old = old[old.gene_universe == "filtered"].copy()
    new = summary[summary.gene_universe == "filtered"].copy()
    cols = [
        "condition", "rank1_energy_fraction", "rank3_cumulative_energy_fraction",
        "median_pairwise_state_spearman", "median_state_rho_to_LOO_shared",
        "median_projection_energy_fraction"
    ]
    comparison = old[cols].merge(new[cols], on="condition", suffixes=("_all_control_avg", "_batch_matched"))
    comparison.to_csv(OUT / "sensitivity_vs_all_control_average.tsv", sep="\t", index=False)

# Report continuous values only, no biological hard threshold.
f = summary[summary.gene_universe == "filtered"]
report = [
    "# Robertson 2026 batch-matched shared-vs-state-specific stress structure",
    "",
    "Design-aware analysis: each stress replicate is differenced against the well-watered replicate from the same GEO batch, then the two matched response vectors are equal-weight averaged per cell annotation.",
    "",
    "No biological hard threshold is used; all structure metrics are interpreted continuously.",
    "",
    f"Technical gene universe: {int(mask.sum())} / {len(genes)} genes.",
    "",
    "## Shared structure",
]
for _, r in f.iterrows():
    report.append(
        f"- {r.condition}: rank1={r.rank1_energy_fraction:.3f}, rank1-3={r.rank3_cumulative_energy_fraction:.3f}, "
        f"median pairwise rho={r.median_pairwise_state_spearman:.3f}, median LOO-shared rho={r.median_state_rho_to_LOO_shared:.3f}, "
        f"median projection energy={r.median_projection_energy_fraction:.3f}."
    )
report += ["", "## Matched replicate stability"]
for cond, sub in rep_df.groupby("condition"):
    report.append(
        f"- {cond}: median={sub.matched_replicate_response_spearman.median():.3f}, "
        f"mean={sub.matched_replicate_response_spearman.mean():.3f}, range=[{sub.matched_replicate_response_spearman.min():.3f}, {sub.matched_replicate_response_spearman.max():.3f}]."
    )
report += [
    "",
    "## Scope",
    "This rice experiment is a multi-day soil stress atlas, not a minutes-scale acute perturbation. It is used as an external contrast for the fraction of shared versus cell-state-specific response, not as a direct replication of the Arabidopsis 15-min rewatering temporal regime.",
]
(OUT / "REPORT.md").write_text("\n".join(report)+"\n", encoding="utf-8")

print("ROBERTSON2026_BATCH_MATCHED_SUCCESS")
print(f.to_string(index=False))
print("\nMATCHED_REPLICATE_STABILITY")
print(rep_df.groupby("condition").matched_replicate_response_spearman.agg(["median","mean","min","max"]).to_string())
