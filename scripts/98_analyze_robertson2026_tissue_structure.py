from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026" / "tissue_pseudobulk_v0"
DESIGN = pd.read_csv(ROOT / "data" / "external_rice" / "Robertson2026" / "sample_design.tsv", sep="\t")
OUT = ROOT / "data" / "external_rice" / "Robertson2026" / "tissue_shared_structure_v0"
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(BASE / "group_metadata.tsv", sep="\t")
expr = pd.read_csv(BASE / "pseudobulk_log1p_cpm.tsv.gz", sep="\t", compression="gzip")
genes = expr.pop("gene").astype(str).to_numpy()
X = expr.to_numpy(dtype=float)
if expr.columns.astype(str).tolist() != meta["group"].astype(str).tolist():
    raise RuntimeError("Expression columns do not match tissue group metadata")

mean_expr = X.mean(axis=1)
nonzero_frac = (X > 0).mean(axis=1)
mask = (nonzero_frac >= 0.25) & (mean_expr >= 0.05)

conditions = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
control = "Well.Watered"
tissues = sorted(meta["tissue_system"].astype(str).unique())
batch_by_sample = dict(zip(DESIGN["identity"].astype(str), DESIGN["batch"].astype(int)))
control_by_batch = {}
for sample in meta.loc[meta.condition == control, "sample"].astype(str).unique():
    control_by_batch[batch_by_sample[sample]] = sample

idx = {(str(r["sample"]), str(r["tissue_system"])): int(i) for i, r in meta.iterrows()}

def rho(a, b):
    return float(spearmanr(a, b).statistic)

def cosine(a, b):
    den = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / den) if den > 0 else np.nan

def rank1(M):
    u, s, vt = np.linalg.svd(M, full_matrices=False)
    v = vt[0].copy()
    cons = M.mean(axis=0)
    if np.dot(v, cons) < 0:
        v *= -1
        u[:, 0] *= -1
    e = s**2
    load = u[:, 0] ** 2
    return {
        "axis": v,
        "frac": float(e[0] / e.sum()),
        "max_load": float(load.max()),
        "effective_states": float(1.0 / np.sum(load**2)),
        "rho_consensus": rho(v, cons),
        "cos_consensus": cosine(v, cons),
    }

rows = []
state_rows = []
for cond in conditions:
    stress_samples = sorted(meta.loc[meta.condition == cond, "sample"].astype(str).unique())
    if len(stress_samples) != 2:
        raise RuntimeError(f"Expected 2 stress samples for {cond}")
    reps = []
    pair_labels = []
    for ss in stress_samples:
        b = batch_by_sample[ss]
        cs = control_by_batch[b]
        mat = []
        for ts in tissues:
            mat.append((X[:, idx[(ss, ts)]] - X[:, idx[(cs, ts)]])[mask])
        reps.append(np.vstack(mat))
        pair_labels.append(f"{ss}-{cs}:batch{b}")
    A, B = reps
    avg = (A + B) / 2.0
    avg_r1 = rank1(avg)
    A1, B1 = rank1(A), rank1(B)
    consA, consB = A.mean(axis=0), B.mean(axis=0)
    pairwise = []
    loo = []
    cunit = avg.mean(axis=0)
    cunit = cunit / np.linalg.norm(cunit)
    proj = []
    for i, ts in enumerate(tissues):
        rr = rho(A[i], B[i])
        state_rows.append({
            "condition": cond,
            "tissue_system": ts,
            "replicate_response_rho": rr,
            "avg_response_rms": float(np.sqrt(np.mean(avg[i]**2))),
        })
        for j in range(i+1, len(tissues)):
            pairwise.append(rho(avg[i], avg[j]))
        loo_mean = np.delete(avg, i, axis=0).mean(axis=0)
        loo.append(rho(avg[i], loo_mean))
        den = float(np.dot(avg[i], avg[i]))
        p = float(np.dot(avg[i], cunit))
        proj.append((p*p/den) if den > 0 else np.nan)
    rows.append({
        "condition": cond,
        "n_tissues": len(tissues),
        "n_genes": int(mask.sum()),
        "median_tissue_replicate_rho": float(np.median([rho(A[i], B[i]) for i in range(len(tissues))])),
        "consensus_response_replicate_rho": rho(consA, consB),
        "consensus_response_replicate_cosine": cosine(consA, consB),
        "rank1_axis_replicate_rho": rho(A1["axis"], B1["axis"]),
        "rank1_axis_replicate_cosine": cosine(A1["axis"], B1["axis"]),
        "avg_rank1_fraction": avg_r1["frac"],
        "avg_rank1_effective_states": avg_r1["effective_states"],
        "avg_max_rank1_loading_share": avg_r1["max_load"],
        "avg_rank1_vs_consensus_rho": avg_r1["rho_consensus"],
        "median_pairwise_tissue_rho": float(np.median(pairwise)),
        "median_tissue_LOO_shared_rho": float(np.median(loo)),
        "median_consensus_projection_energy": float(np.median(proj)),
        "matched_pair_1": pair_labels[0],
        "matched_pair_2": pair_labels[1],
    })

res = pd.DataFrame(rows)
st = pd.DataFrame(state_rows)
res.to_csv(OUT / "tissue_structure_summary.tsv", sep="\t", index=False)
st.to_csv(OUT / "tissue_replicate_metrics.tsv", sep="\t", index=False)
print("ROBERTSON2026_TISSUE_STRUCTURE_SUCCESS")
print(res.to_string(index=False))
print("\nTISSUE_REPLICATE_METRICS")
print(st.to_string(index=False))
