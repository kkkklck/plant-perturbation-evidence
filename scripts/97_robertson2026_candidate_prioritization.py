from __future__ import annotations

from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026"
PB = BASE / "pseudobulk_v0"
DE = BASE / "deseq2_v0"
OUT = BASE / "candidate_prioritization_v0"
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(PB / "group_metadata.tsv", sep="\t")
design = pd.read_csv(BASE / "sample_design.tsv", sep="\t")
expr = pd.read_csv(PB / "pseudobulk_log1p_cpm.tsv.gz", sep="\t", compression="gzip")
genes = expr.pop("gene").astype(str).to_numpy()
X = expr.to_numpy(float)
if expr.columns.astype(str).tolist() != meta.group.astype(str).tolist():
    raise RuntimeError("group order mismatch")
col_index = {g:i for i,g in enumerate(meta.group.astype(str))}
annotations = sorted(meta.cell_annotation.unique())
control_by_batch = dict(zip(design.loc[design.stress_type=="control","batch"].astype(int), design.loc[design.stress_type=="control","identity"].astype(str)))

sig = pd.read_csv(DE / "deseq2_fdr05_genes.tsv.gz", sep="\t", compression="gzip")
# R wrote TRUE/FALSE; normalize robustly.
if sig["strong"].dtype != bool:
    sig["strong"] = sig["strong"].astype(str).str.upper().eq("TRUE")

stress_map = [
    ("Mild.Drought","MildDrought"),
    ("Moderate.Drought","ModerateDrought"),
    ("Mild.Salinity","MildSalinity"),
    ("Moderate.Salinity","ModerateSalinity"),
]

def vec(sample, ann):
    return X[:, col_index[f"{sample}||{ann}"]]

def matrices(mode, design_cond):
    stress_samples = design.loc[design.condition==design_cond,"identity"].astype(str).tolist()
    R=[]; C=[]
    for ann in annotations:
        if mode=="author_like":
            ctrls=design.loc[design.stress_type=="control","identity"].astype(str).tolist()
            ctrl=np.mean([vec(c,ann) for c in ctrls],axis=0)
            stress=np.mean([vec(s,ann) for s in stress_samples],axis=0)
            R.append(stress-ctrl); C.append(ctrl)
        elif mode=="batch_matched":
            deltas=[]; controls=[]
            for s in stress_samples:
                row=design.loc[design.identity==s].iloc[0]
                c=control_by_batch[int(row.batch)]
                deltas.append(vec(s,ann)-vec(c,ann)); controls.append(vec(c,ann))
            R.append(np.mean(deltas,axis=0)); C.append(np.mean(controls,axis=0))
        else:
            raise ValueError(mode)
    return np.vstack(R),np.vstack(C)

def auroc(y,score):
    y=np.asarray(y,dtype=bool); score=np.asarray(score,float)
    n1=int(y.sum()); n0=int((~y).sum())
    if n1==0 or n0==0: return np.nan
    r=rankdata(score,method="average")
    return float((r[y].sum()-n1*(n1+1)/2)/(n1*n0))

def average_precision(y,score):
    y=np.asarray(y,dtype=bool); score=np.asarray(score,float)
    npos=int(y.sum())
    if npos==0: return np.nan
    order=np.argsort(-score,kind="mergesort")
    yy=y[order].astype(int)
    cum=np.cumsum(yy)
    ranks=np.arange(1,len(yy)+1)
    return float(np.sum((cum/ranks)*yy)/npos)

def top_metrics(y,score,frac):
    y=np.asarray(y,dtype=bool); n=len(y); npos=int(y.sum())
    k=max(1,int(math.ceil(frac*n)))
    order=np.argsort(-score,kind="mergesort")[:k]
    hits=int(y[order].sum())
    recall=hits/npos if npos else np.nan
    precision=hits/k
    prevalence=npos/n
    enrichment=(precision/prevalence) if prevalence>0 else np.nan
    return hits,k,recall,precision,enrichment

rows=[]
for mode in ["author_like","batch_matched"]:
    for cond_label,design_cond in stress_map:
        R,C=matrices(mode,design_cond)
        for i,ann in enumerate(annotations):
            truth=sig[(sig["mode"]==mode)&(sig["condition"]==cond_label)&(sig["annotation"]==ann)&sig["strong"]].copy()
            pos_genes=set(truth.gene.astype(str))
            y=np.array([g in pos_genes for g in genes],dtype=bool)
            npos=int(y.sum())
            if npos==0:
                continue

            other=np.delete(R,i,axis=0)
            signed_mean=other.mean(axis=0)
            pred={
                "loo_abs_signed_mean":np.abs(signed_mean),
                "loo_mean_abs":np.mean(np.abs(other),axis=0),
                "target_control_expression":C[i],
            }
            lfc_map=dict(zip(truth.gene.astype(str),truth.log2FoldChange.astype(float)))
            pos_idx=np.where(y)[0]
            actual_sign=np.sign([lfc_map[genes[j]] for j in pos_idx])
            shared_sign=np.sign(signed_mean[pos_idx])
            nonzero=actual_sign!=0
            sign_agreement=float(np.mean(shared_sign[nonzero]==actual_sign[nonzero])) if nonzero.any() else np.nan

            for name,score in pred.items():
                rec={
                    "mode":mode,"condition":cond_label,"cell_annotation":ann,"predictor":name,
                    "n_genes":len(genes),"n_true_strong_deg":npos,"prevalence":npos/len(genes),
                    "auroc":auroc(y,score),"average_precision":average_precision(y,score),
                    "shared_direction_agreement_on_true_deg":sign_agreement if name=="loo_abs_signed_mean" else np.nan,
                }
                for frac in [0.05,0.10,0.20]:
                    hits,k,recall,precision,enrich=top_metrics(y,score,frac)
                    tag=str(int(frac*100))
                    rec[f"top{tag}_k"]=k
                    rec[f"top{tag}_hits"]=hits
                    rec[f"top{tag}_recall"]=recall
                    rec[f"top{tag}_precision"]=precision
                    rec[f"top{tag}_enrichment"]=enrich
                rows.append(rec)

res=pd.DataFrame(rows)
res.to_csv(OUT/"per_state_candidate_prioritization.tsv",sep="\t",index=False)

# Aggregate continuously; no pass/fail threshold.
agg=res.groupby(["mode","condition","predictor"]).agg(
    n_scoreable_states=("cell_annotation","nunique"),
    median_n_true_deg=("n_true_strong_deg","median"),
    median_auroc=("auroc","median"),
    median_ap=("average_precision","median"),
    median_top5_recall=("top5_recall","median"),
    median_top5_enrichment=("top5_enrichment","median"),
    median_top10_recall=("top10_recall","median"),
    median_top10_enrichment=("top10_enrichment","median"),
    median_top20_recall=("top20_recall","median"),
    median_top20_enrichment=("top20_enrichment","median"),
    median_direction_agreement=("shared_direction_agreement_on_true_deg","median"),
).reset_index()
agg.to_csv(OUT/"candidate_prioritization_summary.tsv",sep="\t",index=False)

# Compact cross-condition summary for the main LOO signed-shared predictor.
main=agg[agg.predictor=="loo_abs_signed_mean"].copy()
print("ROBERTSON2026_CANDIDATE_PRIORITIZATION_SUCCESS")
print(main.to_string(index=False))
