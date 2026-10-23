#!/usr/bin/env python3
"""Low-rank structure of matched-batch DESeq2 pseudobulk effects in Robertson 2026.

Compares modeled stress-condition coefficients with the nuisance batch coefficient
from the same state/condition models. No biological pass/fail threshold is used.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026"
INP = BASE / "deseq2_effects_v0" / "CellAnnotation_DESeq2_effects.tsv.gz"
OUT = BASE / "deseq2_effects_v0"
CONDITIONS = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
EFFECTS = {
    "condition_raw": "condition_lfc_raw",
    "condition_shrunk": "condition_lfc_shrunk",
    "batch_raw": "batch_lfc_raw",
    "batch_shrunk": "batch_lfc_shrunk",
}


def rho(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    r = spearmanr(a, b).statistic
    return float(r) if np.isfinite(r) else np.nan


def cosine(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    return float(np.dot(a,b)/(na*nb)) if na > 0 and nb > 0 else np.nan


def structure(R):
    gram = R @ R.T
    vals = np.linalg.eigvalsh(gram)
    vals = np.sort(np.maximum(vals, 0))[::-1]
    frac = vals / vals.sum()
    cum = np.cumsum(frac)
    shared = R.mean(axis=0)
    sh2 = float(np.dot(shared, shared))
    pair = []
    loo = []
    proj = []
    for i in range(R.shape[0]):
        y = R[i]
        lo = np.mean(np.delete(R, i, axis=0), axis=0)
        loo.append(rho(y, lo))
        alpha = float(np.dot(y, shared)/sh2) if sh2 > 0 else np.nan
        p = alpha*shared if np.isfinite(alpha) else np.full_like(y, np.nan)
        e = float(np.dot(y,y))
        re = float(np.dot(y-p,y-p)/e) if e > 0 else np.nan
        proj.append(1-re if np.isfinite(re) else np.nan)
        for j in range(i+1,R.shape[0]):
            pair.append(rho(R[i],R[j]))
    return {
        "rank1": float(frac[0]),
        "rank2_cum": float(cum[min(1,len(cum)-1)]),
        "rank3_cum": float(cum[min(2,len(cum)-1)]),
        "median_pairwise": float(np.nanmedian(pair)),
        "median_loo": float(np.nanmedian(loo)),
        "median_projection": float(np.nanmedian(proj)),
        "shared": shared,
    }


def make_matrix(x: pd.DataFrame, value_col: str):
    states = sorted(x.state.astype(str).unique())
    genes = sorted(x.gene.astype(str).unique())
    p = x.pivot(index="state", columns="gene", values=value_col).reindex(index=states, columns=genes)
    arr = p.to_numpy(float)
    finite = np.isfinite(arr).all(axis=0)
    return states, np.asarray(genes)[finite].tolist(), arr[:,finite]


def main():
    if not INP.exists():
        raise FileNotFoundError(INP)
    x = pd.read_csv(INP, sep="\t", compression="gzip")
    states = sorted(x.state.astype(str).unique())
    if len(states) != 14:
        raise RuntimeError(f"expected 14 states, got {len(states)}")

    rows = []
    shared_vectors = {}
    gene_rows = []
    for cond in CONDITIONS:
        z = x[x.condition == cond].copy()
        for effect_name, col in EFFECTS.items():
            sts, genes, R = make_matrix(z, col)
            if sts != states:
                raise RuntimeError(f"{cond}/{effect_name}: state mismatch")
            S = structure(R)
            rows.append({
                "condition": cond,
                "effect": effect_name,
                "n_states": len(sts),
                "n_finite_genes": len(genes),
                "rank1_energy_fraction": S["rank1"],
                "rank2_cumulative_energy_fraction": S["rank2_cum"],
                "rank3_cumulative_energy_fraction": S["rank3_cum"],
                "median_pairwise_state_spearman": S["median_pairwise"],
                "median_rho_shared_LOO": S["median_loo"],
                "median_projection_energy_fraction": S["median_projection"],
            })
            shared_vectors[(cond,effect_name)] = (genes,S["shared"])

            # Descriptive shared-gene metrics only for condition effects.
            if effect_name.startswith("condition_"):
                mean = R.mean(axis=0)
                sign_pos = (R>0).mean(axis=0)
                sign_neg = (R<0).mean(axis=0)
                sign_cons = np.maximum(sign_pos,sign_neg)
                sd = R.std(axis=0,ddof=1)
                for g,m,sig,sc in zip(genes,mean,sd,sign_cons):
                    gene_rows.append({
                        "condition":cond,"effect":effect_name,"gene":g,
                        "mean_lfc":float(m),"sd_across_states":float(sig),
                        "sign_consistency":float(sc),
                    })

    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "DESeq2_structure_summary.tsv", sep="\t", index=False)
    pd.DataFrame(gene_rows).to_csv(OUT / "DESeq2_condition_gene_sharedness.tsv.gz", sep="\t", index=False, compression="gzip")

    # Same-model stress-vs-batch comparison.
    comp = []
    for cond in CONDITIONS:
        for suffix in ["raw","shrunk"]:
            c = summary[(summary.condition==cond)&(summary.effect==f"condition_{suffix}")].iloc[0]
            b = summary[(summary.condition==cond)&(summary.effect==f"batch_{suffix}")].iloc[0]
            comp.append({
                "condition":cond,"estimator":suffix,
                "condition_rank1":c.rank1_energy_fraction,
                "batch_rank1":b.rank1_energy_fraction,
                "rank1_difference_condition_minus_batch":c.rank1_energy_fraction-b.rank1_energy_fraction,
                "condition_median_LOO":c.median_rho_shared_LOO,
                "batch_median_LOO":b.median_rho_shared_LOO,
                "LOO_difference_condition_minus_batch":c.median_rho_shared_LOO-b.median_rho_shared_LOO,
                "condition_median_pairwise":c.median_pairwise_state_spearman,
                "batch_median_pairwise":b.median_pairwise_state_spearman,
            })
    comp = pd.DataFrame(comp)
    comp.to_csv(OUT / "DESeq2_condition_vs_batch.tsv", sep="\t", index=False)

    # Cross-condition similarity requires a common finite-gene intersection per estimator.
    cross = []
    for suffix in ["raw","shrunk"]:
        for i,a in enumerate(CONDITIONS):
            for b in CONDITIONS[i+1:]:
                ga,va = shared_vectors[(a,f"condition_{suffix}")]
                gb,vb = shared_vectors[(b,f"condition_{suffix}")]
                da = dict(zip(ga,va)); db = dict(zip(gb,vb))
                common = sorted(set(da)&set(db))
                xa = np.array([da[g] for g in common],float)
                xb = np.array([db[g] for g in common],float)
                cross.append({
                    "estimator":suffix,"condition_a":a,"condition_b":b,
                    "n_common_genes":len(common),
                    "spearman_shared_vectors":rho(xa,xb),
                    "cosine_shared_vectors":cosine(xa,xb),
                })
    cross = pd.DataFrame(cross)
    cross.to_csv(OUT / "DESeq2_cross_condition_shared_similarity.tsv",sep="\t",index=False)

    # Top descriptive shared genes from shrinkage (not a significance claim).
    gg = pd.DataFrame(gene_rows)
    top = []
    for cond in CONDITIONS:
        g = gg[(gg.condition==cond)&(gg.effect=="condition_shrunk")].copy()
        g["abs_mean_lfc"] = g.mean_lfc.abs()
        g = g.sort_values(["sign_consistency","abs_mean_lfc"],ascending=[False,False]).head(50)
        top.append(g)
    pd.concat(top,ignore_index=True).to_csv(OUT / "DESeq2_descriptive_top_shared_genes.tsv",sep="\t",index=False)

    compact={"status":"MODEL_BASED_CONTINUOUS_STRUCTURE_NO_BIOLOGICAL_PASS_FAIL_GATE","condition_vs_batch":comp.to_dict(orient="records")}
    (OUT/"DESeq2_structure_summary.json").write_text(json.dumps(compact,indent=2),encoding="utf-8")

    lines=[
        "# Robertson 2026 DESeq2 modeled stress-effect structure","",
        "Matched batches are modeled with design `~ batch + condition`. Both raw and normal-shrunken LFCs are reported. No biological pass/fail threshold is used.","",
        "## Stress condition coefficient versus nuisance batch coefficient","",
        "| condition | estimator | stress rank1 | batch rank1 | delta rank1 | stress median LOO | batch median LOO | delta LOO |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in comp.itertuples(index=False):
        lines.append(f"| {r.condition} | {r.estimator} | {r.condition_rank1:.3f} | {r.batch_rank1:.3f} | {r.rank1_difference_condition_minus_batch:+.3f} | {r.condition_median_LOO:.3f} | {r.batch_median_LOO:.3f} | {r.LOO_difference_condition_minus_batch:+.3f} |")
    lines += ["","## Cross-condition shared-vector similarity (condition effects)",""]
    for r in cross.itertuples(index=False):
        lines.append(f"- {r.estimator}: {r.condition_a} vs {r.condition_b}: rho={r.spearman_shared_vectors:.3f}, cosine={r.cosine_shared_vectors:.3f}, genes={r.n_common_genes}")
    (OUT/"DESEQ2_REPORT.md").write_text("\n".join(lines)+"\n",encoding="utf-8")

    print("ROBERTSON2026_DESEQ2_STRUCTURE_SUCCESS")
    print(summary.to_string(index=False))
    print("\nCONDITION VS BATCH")
    print(comp.to_string(index=False))
    print("\nCROSS")
    print(cross.to_string(index=False))


if __name__ == "__main__":
    main()
