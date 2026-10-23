from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026" / "tissue_pseudobulk_v0"
OUT = ROOT / "data" / "external_rice" / "Robertson2026" / "tissue_batch_adjusted_v0"
OUT.mkdir(parents=True, exist_ok=True)

design = pd.read_csv(ROOT / "data" / "external_rice" / "Robertson2026" / "sample_design.tsv", sep="\t")
meta = pd.read_csv(BASE / "group_metadata.tsv", sep="\t")
expr = pd.read_csv(BASE / "pseudobulk_log1p_cpm.tsv.gz", sep="\t", compression="gzip")
genes = expr.pop("gene").astype(str).to_numpy()
X = expr.to_numpy(dtype=float)
if expr.columns.astype(str).tolist() != meta["group"].astype(str).tolist():
    raise RuntimeError("Expression columns do not match group metadata")

# Technical expression universe only; no condition-effect thresholding.
mean_expr = X.mean(axis=1)
nonzero_frac = (X > 0).mean(axis=1)
mask = (nonzero_frac >= 0.25) & (mean_expr >= 0.05)
genes_f = genes[mask]

conditions = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
control = "Well.Watered"
tissues = sorted(meta["tissue_system"].astype(str).unique())
samples = design["identity"].astype(str).tolist()

# Canonical sample metadata from the Seurat object, plus batch from GEO design.
sc = meta[["sample", "condition"]].drop_duplicates().set_index("sample")
batch_by_sample = dict(zip(design["identity"].astype(str), design["batch"].astype(int)))
for s in samples:
    if s not in sc.index:
        raise RuntimeError(f"Sample missing from tissue pseudobulk metadata: {s}")

# group lookup
idx = {(str(r["sample"]), str(r["tissue_system"])): int(i) for i, r in meta.iterrows()}

# Design matrix: intercept + batch2 + batch3 + four condition indicators.
def make_design(keep_samples: list[str]):
    rows = []
    for s in keep_samples:
        cond = str(sc.loc[s, "condition"])
        batch = batch_by_sample[s]
        rows.append([
            1.0,
            1.0 if batch == 2 else 0.0,
            1.0 if batch == 3 else 0.0,
            1.0 if cond == "Mild.Drought" else 0.0,
            1.0 if cond == "Moderate.Drought" else 0.0,
            1.0 if cond == "Mild.Salinity" else 0.0,
            1.0 if cond == "Moderate.Salinity" else 0.0,
        ])
    D = np.asarray(rows, dtype=float)
    if np.linalg.matrix_rank(D) < D.shape[1]:
        raise RuntimeError(f"Design rank deficient for samples: {keep_samples}")
    return D

coef_index = {
    "Mild.Drought": 3,
    "Moderate.Drought": 4,
    "Mild.Salinity": 5,
    "Moderate.Salinity": 6,
}

def fit_effects(keep_samples: list[str]):
    D = make_design(keep_samples)
    pinv = np.linalg.pinv(D)
    out = {c: [] for c in conditions}
    for ts in tissues:
        Y = np.vstack([X[:, idx[(s, ts)]][mask] for s in keep_samples])  # samples x genes
        B = pinv @ Y  # coefficients x genes
        for c in conditions:
            out[c].append(B[coef_index[c]].copy())
    return {c: np.vstack(v) for c, v in out.items()}  # tissues x genes


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

full = fit_effects(samples)
summary_rows = []
state_rows = []
jack_rows = []

for cond in conditions:
    M = full[cond]
    r1 = rank1(M)
    pair = []
    loo = []
    cons = M.mean(axis=0)
    cunit = cons / np.linalg.norm(cons)
    proj = []
    for i, ts in enumerate(tissues):
        for j in range(i + 1, len(tissues)):
            pair.append(rho(M[i], M[j]))
        loo_mean = np.delete(M, i, axis=0).mean(axis=0)
        loo_rho = rho(M[i], loo_mean)
        den = float(np.dot(M[i], M[i]))
        p = float(np.dot(M[i], cunit))
        proj_frac = p*p/den if den > 0 else np.nan
        loo.append(loo_rho)
        proj.append(proj_frac)
        state_rows.append({
            "condition": cond,
            "tissue_system": ts,
            "effect_rms": float(np.sqrt(np.mean(M[i]**2))),
            "rho_to_LOO_shared": loo_rho,
            "projection_energy_fraction": proj_frac,
        })

    # Leave-one-stress-sample jackknife: each estimate keeps the other replicate,
    # all controls, and the other stress conditions, allowing batch terms to be re-estimated.
    stress_samples = [s for s in samples if str(sc.loc[s, "condition"]) == cond]
    if len(stress_samples) != 2:
        raise RuntimeError(f"Expected two stress samples for {cond}")
    jk = []
    labels = []
    for omit in stress_samples:
        kept = [s for s in samples if s != omit]
        eff = fit_effects(kept)[cond]
        jk.append(eff)
        labels.append(f"omit:{omit}")
    A, B = jk
    state_jk = [rho(A[i], B[i]) for i in range(len(tissues))]
    consA, consB = A.mean(axis=0), B.mean(axis=0)
    A1, B1 = rank1(A), rank1(B)
    for i, ts in enumerate(tissues):
        jack_rows.append({
            "condition": cond,
            "tissue_system": ts,
            "jackknife_effect_rho": state_jk[i],
            "omit_a": labels[0],
            "omit_b": labels[1],
        })

    summary_rows.append({
        "condition": cond,
        "n_tissues": len(tissues),
        "n_genes": int(mask.sum()),
        "rank1_energy_fraction": r1["frac"],
        "rank1_effective_tissues": r1["effective_states"],
        "max_rank1_loading_share": r1["max_load"],
        "rank1_vs_consensus_rho": r1["rho_consensus"],
        "median_pairwise_tissue_rho": float(np.median(pair)),
        "median_tissue_LOO_shared_rho": float(np.median(loo)),
        "median_consensus_projection_energy": float(np.median(proj)),
        "median_tissue_jackknife_rho": float(np.median(state_jk)),
        "consensus_jackknife_rho": rho(consA, consB),
        "consensus_jackknife_cosine": cosine(consA, consB),
        "rank1_axis_jackknife_rho": rho(A1["axis"], B1["axis"]),
        "rank1_axis_jackknife_cosine": cosine(A1["axis"], B1["axis"]),
    })

summary = pd.DataFrame(summary_rows)
states = pd.DataFrame(state_rows)
jack = pd.DataFrame(jack_rows)
summary.to_csv(OUT / "batch_adjusted_tissue_summary.tsv", sep="\t", index=False)
states.to_csv(OUT / "batch_adjusted_tissue_state_metrics.tsv", sep="\t", index=False)
jack.to_csv(OUT / "batch_adjusted_tissue_jackknife.tsv", sep="\t", index=False)

# Save the full batch-adjusted coefficient matrices for later gene-level interpretation.
for cond in conditions:
    M = full[cond]
    df = pd.DataFrame(M.T, columns=tissues)
    df.insert(0, "gene", genes_f)
    safe = cond.lower().replace(".", "_")
    df.to_csv(OUT / f"effect_{safe}.tsv.gz", sep="\t", index=False, compression="gzip")

print("ROBERTSON2026_BATCH_ADJUSTED_TISSUE_SUCCESS")
print(summary.to_string(index=False))
print("\nSTATE_METRICS")
print(states.to_string(index=False))
