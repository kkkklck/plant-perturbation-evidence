from __future__ import annotations

from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data"/"external_rice"/"Robertson2026"
PB=BASE/"pseudobulk_v0"
DE=BASE/"deseq2_v0"
OUT=BASE/"cross_condition_prioritization_v0"
OUT.mkdir(parents=True,exist_ok=True)

meta=pd.read_csv(PB/"group_metadata.tsv",sep="\t")
design=pd.read_csv(BASE/"sample_design.tsv",sep="\t")
expr=pd.read_csv(PB/"pseudobulk_log1p_cpm.tsv.gz",sep="\t",compression="gzip")
genes=expr.pop("gene").astype(str).to_numpy(); X=expr.to_numpy(float)
if expr.columns.astype(str).tolist()!=meta.group.astype(str).tolist(): raise RuntimeError("group order mismatch")
col_index={g:i for i,g in enumerate(meta.group.astype(str))}
annotations=sorted(meta.cell_annotation.unique())
control_by_batch=dict(zip(design.loc[design.stress_type=="control","batch"].astype(int),design.loc[design.stress_type=="control","identity"].astype(str)))

sig=pd.read_csv(DE/"deseq2_fdr05_genes.tsv.gz",sep="\t",compression="gzip")
if sig["strong"].dtype!=bool: sig["strong"]=sig["strong"].astype(str).str.upper().eq("TRUE")

cond_map={
    "Mild.Drought":"MildDrought",
    "Moderate.Drought":"ModerateDrought",
    "Mild.Salinity":"MildSalinity",
    "Moderate.Salinity":"ModerateSalinity",
}
cond_meta={
    "Mild.Drought":("drought","mild"),
    "Moderate.Drought":("drought","moderate"),
    "Mild.Salinity":("salinity","mild"),
    "Moderate.Salinity":("salinity","moderate"),
}

def vec(sample,ann): return X[:,col_index[f"{sample}||{ann}"]]

def batch_matched_matrix(label):
    dcond=cond_map[label]
    ss=design.loc[design.condition==dcond,"identity"].astype(str).tolist()
    R=[]; C=[]
    for ann in annotations:
        deltas=[]; ctrls=[]
        for s in ss:
            row=design.loc[design.identity==s].iloc[0]
            c=control_by_batch[int(row.batch)]
            deltas.append(vec(s,ann)-vec(c,ann)); ctrls.append(vec(c,ann))
        R.append(np.mean(deltas,axis=0)); C.append(np.mean(ctrls,axis=0))
    return np.vstack(R),np.vstack(C)

R={}; C={}
for c in cond_map: R[c],C[c]=batch_matched_matrix(c)

def auroc(y,s):
    y=np.asarray(y,bool); n1=int(y.sum()); n0=int((~y).sum())
    if not n1 or not n0:return np.nan
    r=rankdata(s,method="average")
    return float((r[y].sum()-n1*(n1+1)/2)/(n1*n0))

def ap(y,s):
    y=np.asarray(y,bool); n=int(y.sum())
    if not n:return np.nan
    o=np.argsort(-s,kind="mergesort"); z=y[o].astype(int); c=np.cumsum(z); ranks=np.arange(1,len(z)+1)
    return float(np.sum((c/ranks)*z)/n)

def top(y,s,f=.10):
    y=np.asarray(y,bool); n=len(y); np1=int(y.sum()); k=max(1,math.ceil(f*n)); o=np.argsort(-s,kind="mergesort")[:k]
    hit=int(y[o].sum()); rec=hit/np1 if np1 else np.nan; prec=hit/k; prev=np1/n; enr=prec/prev if prev else np.nan
    return hit,k,rec,prec,enr

rows=[]
conditions=list(cond_map)
for target in conditions:
    ttype,tsev=cond_meta[target]
    for source in conditions:
        if source==target: continue
        stype,ssev=cond_meta[source]
        if stype==ttype:
            relation="same_stress_different_severity"
        elif ssev==tsev:
            relation="cross_stress_same_severity"
        else:
            relation="cross_stress_cross_severity"
        # global similarity is descriptive context, independent of target cell annotation.
        cond_rho=float(spearmanr(R[source].mean(axis=0),R[target].mean(axis=0)).statistic)
        for i,ann in enumerate(annotations):
            truth=sig[(sig["mode"]=="batch_matched")&(sig["condition"]==target)&(sig["annotation"]==ann)&sig["strong"]].copy()
            pos=set(truth.gene.astype(str)); y=np.array([g in pos for g in genes],bool); npos=int(y.sum())
            if npos==0: continue
            source_same=R[source][i]
            source_loo=np.delete(R[source],i,axis=0).mean(axis=0)
            source_all=R[source].mean(axis=0)
            pred={
                "source_same_state":np.abs(source_same),
                "source_LOO_shared":np.abs(source_loo),
                "source_all_state_shared":np.abs(source_all),
                "target_control_expression":C[target][i],
            }
            lfc=dict(zip(truth.gene.astype(str),truth.log2FoldChange.astype(float)))
            idx=np.where(y)[0]; actual=np.sign([lfc[genes[j]] for j in idx])
            signs={
                "source_same_state":np.sign(source_same[idx]),
                "source_LOO_shared":np.sign(source_loo[idx]),
                "source_all_state_shared":np.sign(source_all[idx]),
            }
            for name,score in pred.items():
                hit,k,rec,prec,enr=top(y,score,.10)
                signagr=np.nan
                if name in signs:
                    nz=actual!=0; signagr=float(np.mean(signs[name][nz]==actual[nz])) if nz.any() else np.nan
                rows.append({
                    "source_condition":source,"target_condition":target,"relation":relation,
                    "condition_mean_response_spearman":cond_rho,"cell_annotation":ann,"predictor":name,
                    "n_true_strong_deg":npos,"auroc":auroc(y,score),"average_precision":ap(y,score),
                    "top10_recall":rec,"top10_enrichment":enr,"direction_agreement_on_true_deg":signagr,
                })

res=pd.DataFrame(rows)
res.to_csv(OUT/"per_state_cross_condition_prioritization.tsv",sep="\t",index=False)
agg=res.groupby(["source_condition","target_condition","relation","predictor"]).agg(
    n_scoreable_states=("cell_annotation","nunique"),
    condition_mean_response_spearman=("condition_mean_response_spearman","first"),
    median_true_deg=("n_true_strong_deg","median"),
    median_auroc=("auroc","median"),
    median_ap=("average_precision","median"),
    median_top10_recall=("top10_recall","median"),
    median_top10_enrichment=("top10_enrichment","median"),
    median_direction_agreement=("direction_agreement_on_true_deg","median"),
).reset_index()
agg.to_csv(OUT/"cross_condition_prioritization_summary.tsv",sep="\t",index=False)

main=agg[agg.predictor.isin(["source_same_state","source_LOO_shared","target_control_expression"])].copy()
print("ROBERTSON2026_CROSS_CONDITION_PRIORITIZATION_SUCCESS")
print(main.sort_values(["target_condition","source_condition","predictor"]).to_string(index=False))
