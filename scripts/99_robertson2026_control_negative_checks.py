from __future__ import annotations

from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data"/"external_rice"/"Robertson2026"
PB=BASE/"pseudobulk_v0"; DE=BASE/"deseq2_v0"; OUT=BASE/"cross_condition_prioritization_v0"
OUT.mkdir(parents=True,exist_ok=True)
meta=pd.read_csv(PB/"group_metadata.tsv",sep="\t")
expr=pd.read_csv(PB/"pseudobulk_log1p_cpm.tsv.gz",sep="\t",compression="gzip")
genes=expr.pop("gene").astype(str).to_numpy(); X=expr.to_numpy(float)
col_index={g:i for i,g in enumerate(meta.group.astype(str))}
annotations=sorted(meta.cell_annotation.unique())
sig=pd.read_csv(DE/"deseq2_fdr05_genes.tsv.gz",sep="\t",compression="gzip")
if sig["strong"].dtype!=bool:sig["strong"]=sig["strong"].astype(str).str.upper().eq("TRUE")
controls=["Well.Watered.R1","Well.Watered.R2","Well.Watered.R3"]
targets=["Mild.Drought","Moderate.Drought","Mild.Salinity","Moderate.Salinity"]

def vec(s,a):return X[:,col_index[f"{s}||{a}"]]
def auroc(y,s):
 y=np.asarray(y,bool);n1=int(y.sum());n0=int((~y).sum())
 if not n1 or not n0:return np.nan
 r=rankdata(s,method="average");return float((r[y].sum()-n1*(n1+1)/2)/(n1*n0))
def top10(y,s):
 y=np.asarray(y,bool);n=len(y);npos=int(y.sum());k=math.ceil(.1*n);o=np.argsort(-s,kind="mergesort")[:k];h=int(y[o].sum());rec=h/npos if npos else np.nan;prev=npos/n;prec=h/k;return rec,(prec/prev if prev else np.nan)

rows=[]
for target in targets:
 for ann in annotations:
  truth=sig[(sig["mode"]=="batch_matched")&(sig["condition"]==target)&(sig["annotation"]==ann)&sig["strong"]]
  pos=set(truth.gene.astype(str));y=np.array([g in pos for g in genes],bool)
  if y.sum()==0:continue
  v=[vec(c,ann) for c in controls]
  scores={
   "WW2_minus_WW1":np.abs(v[1]-v[0]),
   "WW3_minus_WW1":np.abs(v[2]-v[0]),
   "WW3_minus_WW2":np.abs(v[2]-v[1]),
   "mean_abs_WW_pair_difference":np.mean(np.vstack([np.abs(v[1]-v[0]),np.abs(v[2]-v[0]),np.abs(v[2]-v[1])]),axis=0),
  }
  for name,score in scores.items():
   rec,enr=top10(y,score)
   rows.append({"target_condition":target,"cell_annotation":ann,"negative_control":name,"n_true_deg":int(y.sum()),"auroc":auroc(y,score),"top10_recall":rec,"top10_enrichment":enr})
res=pd.DataFrame(rows)
res.to_csv(OUT/"control_only_negative_checks.tsv",sep="\t",index=False)
agg=res.groupby(["target_condition","negative_control"]).agg(n_states=("cell_annotation","nunique"),median_auroc=("auroc","median"),median_top10_recall=("top10_recall","median"),median_top10_enrichment=("top10_enrichment","median")).reset_index()
agg.to_csv(OUT/"control_only_negative_check_summary.tsv",sep="\t",index=False)
print("ROBERTSON2026_CONTROL_NEGATIVE_CHECK_SUCCESS")
print(agg.to_string(index=False))
