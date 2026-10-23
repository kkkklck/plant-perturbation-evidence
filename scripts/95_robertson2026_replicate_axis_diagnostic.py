from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026" / "pseudobulk_v0"
OUT = ROOT / "data" / "external_rice" / "Robertson2026" / "shared_structure_v0"
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(BASE / "group_metadata.tsv", sep="\t")
design = pd.read_csv(ROOT / "data" / "external_rice" / "Robertson2026" / "sample_design.tsv", sep="\t")
expr = pd.read_csv(BASE / "pseudobulk_log1p_cpm.tsv.gz", sep="\t", compression="gzip")
genes = expr.pop("gene").astype(str).to_numpy()
X = expr.to_numpy(dtype=float)
if expr.columns.astype(str).tolist() != meta["group"].astype(str).tolist():
    raise RuntimeError("Expression columns do not match group_metadata order")

# Same technical universe used in script 94.
mean_expr = X.mean(axis=1)
nonzero_frac = (X > 0).mean(axis=1)
mask = (nonzero_frac >= 0.25) & (mean_expr >= 0.05)

conditions = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
control = "Well.Watered"
annotations = sorted(meta["cell_annotation"].astype(str).unique())

batch_by_sample = dict(zip(design["identity"].astype(str), design["batch"].astype(int)))
control_sample_by_batch: dict[int, str] = {}
for sample in meta.loc[meta.condition == control, "sample"].astype(str).unique():
    b = batch_by_sample[sample]
    if b in control_sample_by_batch:
        raise RuntimeError(f"Multiple control samples for batch {b}")
    control_sample_by_batch[b] = sample

# Fast exact lookup: each sample x annotation group should be unique.
idx_lookup: dict[tuple[str, str], int] = {}
for i, row in meta.iterrows():
    key = (str(row["sample"]), str(row["cell_annotation"]))
    if key in idx_lookup:
        raise RuntimeError(f"Duplicate group: {key}")
    idx_lookup[key] = int(i)


def rho(a: np.ndarray, b: np.ndarray) -> float:
    return float(spearmanr(a, b).statistic)


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    den = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / den) if den > 0 else np.nan


def rank1_metrics(M: np.ndarray) -> dict[str, float | np.ndarray]:
    u, s, vt = np.linalg.svd(M, full_matrices=False)
    v = vt[0].copy()
    consensus = M.mean(axis=0)
    if np.dot(v, consensus) < 0:
        v *= -1
        u[:, 0] *= -1
    energy = s**2
    loading_share = u[:, 0] ** 2  # sums to 1
    return {
        "axis": v,
        "rank1_fraction": float(energy[0] / energy.sum()),
        "max_loading_share": float(loading_share.max()),
        "effective_states": float(1.0 / np.sum(loading_share**2)),
        "rank1_vs_consensus_rho": rho(v, consensus),
        "rank1_vs_consensus_cosine": cosine(v, consensus),
    }

rows = []
state_rows = []
for cond in conditions:
    stress_samples = design.loc[design["condition"].eq(cond.replace(".", "")), "identity"].astype(str).tolist()
    # sample_design condition labels omit punctuation; fall back by identity prefix if needed.
    if len(stress_samples) != 2:
        stress_samples = sorted(meta.loc[meta.condition == cond, "sample"].astype(str).unique())
    if len(stress_samples) != 2:
        raise RuntimeError(f"Expected 2 stress samples for {cond}, got {stress_samples}")

    replicate_mats = []
    pair_labels = []
    for ss in stress_samples:
        b = batch_by_sample[ss]
        cs = control_sample_by_batch[b]
        mat = []
        for ann in annotations:
            ti = idx_lookup[(ss, ann)]
            ci = idx_lookup[(cs, ann)]
            mat.append((X[:, ti] - X[:, ci])[mask])
        replicate_mats.append(np.vstack(mat))
        pair_labels.append(f"{ss}-{cs}:batch{b}")

    A, B = replicate_mats
    A1 = rank1_metrics(A)
    B1 = rank1_metrics(B)
    consA = A.mean(axis=0)
    consB = B.mean(axis=0)

    state_rhos = []
    for i, ann in enumerate(annotations):
        r = rho(A[i], B[i])
        state_rhos.append(r)
        state_rows.append({
            "condition": cond,
            "cell_annotation": ann,
            "replicate_response_rho": r,
            "replicate1_response_rms": float(np.sqrt(np.mean(A[i] ** 2))),
            "replicate2_response_rms": float(np.sqrt(np.mean(B[i] ** 2))),
        })

    avg = (A + B) / 2.0
    avg1 = rank1_metrics(avg)
    avg_cons = avg.mean(axis=0)
    loo = []
    consensus_proj = []
    cunit = avg_cons / np.linalg.norm(avg_cons)
    for i in range(len(annotations)):
        loo_mean = np.delete(avg, i, axis=0).mean(axis=0)
        loo.append(rho(avg[i], loo_mean))
        den = float(np.dot(avg[i], avg[i]))
        p = float(np.dot(avg[i], cunit))
        consensus_proj.append((p * p / den) if den > 0 else np.nan)

    rows.append({
        "condition": cond,
        "matched_pair_1": pair_labels[0],
        "matched_pair_2": pair_labels[1],
        "n_states": len(annotations),
        "n_genes": int(mask.sum()),
        "median_state_replicate_rho": float(np.nanmedian(state_rhos)),
        "mean_state_replicate_rho": float(np.nanmean(state_rhos)),
        "flattened_response_replicate_rho": rho(A.ravel(), B.ravel()),
        "consensus_response_replicate_rho": rho(consA, consB),
        "consensus_response_replicate_cosine": cosine(consA, consB),
        "rank1_axis_replicate_rho": rho(A1["axis"], B1["axis"]),
        "rank1_axis_replicate_cosine": cosine(A1["axis"], B1["axis"]),
        "rep1_rank1_fraction": A1["rank1_fraction"],
        "rep2_rank1_fraction": B1["rank1_fraction"],
        "rep1_max_rank1_loading_share": A1["max_loading_share"],
        "rep2_max_rank1_loading_share": B1["max_loading_share"],
        "rep1_rank1_effective_states": A1["effective_states"],
        "rep2_rank1_effective_states": B1["effective_states"],
        "rep1_rank1_vs_consensus_rho": A1["rank1_vs_consensus_rho"],
        "rep2_rank1_vs_consensus_rho": B1["rank1_vs_consensus_rho"],
        "avg_rank1_fraction": avg1["rank1_fraction"],
        "avg_max_rank1_loading_share": avg1["max_loading_share"],
        "avg_rank1_effective_states": avg1["effective_states"],
        "avg_rank1_vs_consensus_rho": avg1["rank1_vs_consensus_rho"],
        "avg_rank1_vs_consensus_cosine": avg1["rank1_vs_consensus_cosine"],
        "avg_median_LOO_shared_rho": float(np.nanmedian(loo)),
        "avg_median_consensus_projection_energy": float(np.nanmedian(consensus_proj)),
    })

out = pd.DataFrame(rows)
states = pd.DataFrame(state_rows)
out.to_csv(OUT / "replicate_axis_diagnostic.tsv", sep="\t", index=False)
states.to_csv(OUT / "state_replicate_diagnostic.tsv", sep="\t", index=False)

print("ROBERTSON2026_REPLICATE_AXIS_DIAGNOSTIC_SUCCESS")
print(out.to_string(index=False))
