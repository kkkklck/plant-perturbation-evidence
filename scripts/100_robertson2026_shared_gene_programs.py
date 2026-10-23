from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026" / "tissue_batch_adjusted_v0"
OUT = BASE
conditions = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
files = {c: BASE / f"effect_{c.lower().replace('.', '_')}.tsv.gz" for c in conditions}

mats = {}
means = {}
rank1_axes = {}
for c, p in files.items():
    df = pd.read_csv(p, sep="\t", compression="gzip")
    genes = df.pop("gene").astype(str)
    M = df.to_numpy(dtype=float).T  # tissues x genes
    mats[c] = (genes.to_numpy(), M, df.columns.astype(str).tolist())
    means[c] = M.mean(axis=0)
    u, s, vt = np.linalg.svd(M, full_matrices=False)
    v = vt[0].copy()
    if np.dot(v, means[c]) < 0:
        v *= -1
    rank1_axes[c] = v

# Verify identical gene order.
g0 = mats[conditions[0]][0]
for c in conditions[1:]:
    if not np.array_equal(g0, mats[c][0]):
        raise RuntimeError(f"Gene order differs for {c}")

pair_rows = []
for i, a in enumerate(conditions):
    for b in conditions[i+1:]:
        pair_rows.append({
            "condition_a": a,
            "condition_b": b,
            "mean_effect_spearman": float(spearmanr(means[a], means[b]).statistic),
            "rank1_axis_spearman": float(spearmanr(rank1_axes[a], rank1_axes[b]).statistic),
            "rank1_axis_cosine": float(np.dot(rank1_axes[a], rank1_axes[b]) / (np.linalg.norm(rank1_axes[a]) * np.linalg.norm(rank1_axes[b]))),
        })
pair_df = pd.DataFrame(pair_rows)
pair_df.to_csv(OUT / "condition_effect_similarity.tsv", sep="\t", index=False)

# Gene-level sharedness within each condition: mean effect, tissue SD, sign consistency,
# and residual RMS around the across-tissue mean. No hard biological threshold is used.
gene_rows = []
top_sets = {}
for c in conditions:
    genes, M, tissues = mats[c]
    mu = M.mean(axis=0)
    sd = M.std(axis=0, ddof=0)
    resid_rms = np.sqrt(np.mean((M - mu[None, :])**2, axis=0))
    pos = (M > 0).sum(axis=0)
    neg = (M < 0).sum(axis=0)
    sign_cons = np.maximum(pos, neg) / M.shape[0]
    same_dir_all = np.maximum(pos, neg) == M.shape[0]
    # Rank by absolute mean effect among 4/4 same-direction genes for an interpretable shared list.
    order = np.argsort(-np.abs(mu))
    rank = np.empty_like(order)
    rank[order] = np.arange(1, len(order)+1)
    for j, g in enumerate(genes):
        gene_rows.append({
            "condition": c,
            "gene": g,
            "mean_effect": float(mu[j]),
            "tissue_sd": float(sd[j]),
            "residual_rms": float(resid_rms[j]),
            "sign_consistency": float(sign_cons[j]),
            "same_direction_4of4": bool(same_dir_all[j]),
            "abs_mean_rank_all_genes": int(rank[j]),
        })
    same_idx = np.where(same_dir_all)[0]
    same_order = same_idx[np.argsort(-np.abs(mu[same_idx]))]
    top_sets[c] = set(genes[same_order[:100]])
    top = pd.DataFrame({
        "gene": genes[same_order[:200]],
        "mean_effect": mu[same_order[:200]],
        "residual_rms": resid_rms[same_order[:200]],
        "sign_consistency": sign_cons[same_order[:200]],
    })
    top.to_csv(OUT / f"shared_genes_{c.lower().replace('.', '_')}_top200.tsv", sep="\t", index=False)

pd.DataFrame(gene_rows).to_csv(OUT / "gene_sharedness_all_conditions.tsv.gz", sep="\t", index=False, compression="gzip")

overlap_rows = []
for i, a in enumerate(conditions):
    for b in conditions[i+1:]:
        inter = top_sets[a] & top_sets[b]
        overlap_rows.append({
            "condition_a": a,
            "condition_b": b,
            "top100_shared_gene_overlap": len(inter),
            "jaccard_top100": len(inter) / len(top_sets[a] | top_sets[b]),
        })
overlap_df = pd.DataFrame(overlap_rows)
overlap_df.to_csv(OUT / "shared_gene_top100_overlap.tsv", sep="\t", index=False)

print("ROBERTSON2026_SHARED_GENE_PROGRAMS_SUCCESS")
print("\nCONDITION EFFECT SIMILARITY")
print(pair_df.to_string(index=False))
print("\nTOP100 SHARED-GENE OVERLAP")
print(overlap_df.to_string(index=False))
for c in conditions:
    p = OUT / f"shared_genes_{c.lower().replace('.', '_')}_top200.tsv"
    df = pd.read_csv(p, sep="\t").head(12)
    print(f"\n{c} TOP SHARED GENES")
    print(df.to_string(index=False))
