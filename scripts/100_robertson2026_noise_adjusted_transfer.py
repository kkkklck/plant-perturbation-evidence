from __future__ import annotations

from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data"/"external_rice"/"Robertson2026"
PB=BASE/"pseudobulk_v0"; DE=BASE/"deseq2_v0"; OUT=BASE/"cross_condition_prioritization_v0"
meta=pd.read_csv(PB/"group_metadata.tsv",sep="\t")
design=pd.read_csv(BASE/"sample_design.tsv",sep="\t")
expr=pd.read_csv(PB/"pseudobulk_log1p_cpm.tsv.gz",sep="\t",compression="gzip")
genes=expr.pop("gene").astype(str).to_numpy(); X=expr.to_numpy(float)
col_index={g:i for i,g in enumerate(meta.group.astype(str))}; annotations=sorted(meta.cell_annotation.unique())
sig=pd.read_csv(DE/"deseq2_fdr05_genes.tsv.gz",sep="\t",compression="gzip")
if sig["strong"].dtype!=bool:sig["strong"]=sig["strong"].astype(str).str.upper().eq("TRUE")
controls=["Well.Watered.R1","Well.Watered.R2","Well.Watered.R3"]
control_by_batch=dict(zip(design.loc[design.stress_type=="control","batch"].astype(int),design.loc[design.stress_type=="control","identity"].astype(str)))
cond_map={"Mild.Drought":"MildDrought","Moderate.Drought":"ModerateDrought","Mild.Salinity":"MildSalinity","Moderate.Salinity":"ModerateSalinity"}

def vec(s,a):return X[:,col_index[f"{s}||{a}"]]
def response(label):
 ss=design.loc[design.condition==cond_map[label],"identity"].astype(str).tolist(); R=[]
 for ann in annotations:
  ds=[]
  for s in ss:
   row=design.loc[design.identity==s].iloc[0];c=control_by_batch[int(row.batch)];ds.append(vec(s,ann)-vec(c,ann))
  R.append(np.mean(ds,axis=0))
 return np.vstack(R)
R={c:response(c) for c in cond_map}

def auroc(y,s):
 y=np.asarray(y,bool);n1=int(y.sum());n0=int((~y).sum())
 if not n1 or not n0:return np.nan
 r=rankdata(s,method="average");return float((r[y].sum()-n1*(n1+1)/2)/(n1*n0))
def top10(y,s):
 y=np.asarray(y,bool);n=len(y);npos=int(y.sum());k=math.ceil(.1*n);o=np.argsort(-s,kind="mergesort")[:k];h=int(y[o].sum());rec=h/npos;prec=h/k;prev=npos/n;return rec,prec/prev

rows=[]; negdir=[]
for target in cond_map:
 for i,ann in enumerate(annotations):
  truth=sig[(sig["mode"]=="batch_matched")&(sig["condition"]==target)&(sig["annotation"]==ann)&sig["strong"]].copy()
  if len(truth)==0:continue
  pos=set(truth.gene.astype(str));y=np.array([g in pos for g in genes],bool)
  lfc=dict(zip(truth.gene.astype(str),truth.log2FoldChange.astype(float)));idx=np.where(y)[0];actual=np.sign([lfc[genes[j]] for j in idx])

  # Control variability estimated only from OTHER cell annotations.
  noise_vectors=[]
  signed_control_vectors=[]
  for j,a2 in enumerate(annotations):
   if j==i:continue
   w=[vec(c,a2) for c in controls]
   noise_vectors += [np.abs(w[1]-w[0]),np.abs(w[2]-w[0]),np.abs(w[2]-w[1])]
  noise=np.median(np.vstack(noise_vectors),axis=0)
  positive=noise[noise>0]
  floor=float(np.quantile(positive,0.10)) if len(positive) else 1e-6
  denom=np.maximum(noise,floor)

  # Signed LOO control-only contrasts for direction negative controls.
  for cname,(c1,c2) in {"WW2-WW1":(controls[1],controls[0]),"WW3-WW1":(controls[2],controls[0]),"WW3-WW2":(controls[2],controls[1])}.items():
   cv=np.mean([vec(c1,a2)-vec(c2,a2) for j,a2 in enumerate(annotations) if j!=i],axis=0)
   predsign=np.sign(cv[idx]);nz=actual!=0
   negdir.append({"target_condition":target,"cell_annotation":ann,"control_contrast":cname,"direction_agreement":float(np.mean(predsign[nz]==actual[nz]))})

  for source in cond_map:
   if source==target:continue
   shared=np.delete(R[source],i,axis=0).mean(axis=0)
   scores={
    "raw_source_LOO_abs":np.abs(shared),
    "noise_adjusted_ratio":np.abs(shared)/denom,
    "noise_adjusted_excess":np.maximum(np.abs(shared)-noise,0),
    "control_variability":noise,
   }
   predsign=np.sign(shared[idx]);nz=actual!=0;diragr=float(np.mean(predsign[nz]==actual[nz]))
   for name,score in scores.items():
    rec,enr=top10(y,score)
    rows.append({"source_condition":source,"target_condition":target,"cell_annotation":ann,"predictor":name,"n_true_deg":int(y.sum()),"auroc":auroc(y,score),"top10_recall":rec,"top10_enrichment":enr,"direction_agreement":diragr if name.startswith("raw_") or name.startswith("noise_adjusted") else np.nan})

res=pd.DataFrame(rows);res.to_csv(OUT/"noise_adjusted_cross_condition_per_state.tsv",sep="\t",index=False)
agg=res.groupby(["source_condition","target_condition","predictor"]).agg(n_states=("cell_annotation","nunique"),median_auroc=("auroc","median"),median_top10_recall=("top10_recall","median"),median_top10_enrichment=("top10_enrichment","median"),median_direction_agreement=("direction_agreement","median")).reset_index()
agg.to_csv(OUT/"noise_adjusted_cross_condition_summary.tsv",sep="\t",index=False)
nd=pd.DataFrame(negdir);nd.to_csv(OUT/"control_only_direction_negative.tsv",sep="\t",index=False)
nda=nd.groupby(["target_condition","control_contrast"]).direction_agreement.agg(["median","mean","min","max"]).reset_index();nda.to_csv(OUT/"control_only_direction_negative_summary.tsv",sep="\t",index=False)
print("ROBERTSON2026_NOISE_ADJUSTED_TRANSFER_SUCCESS")
print(agg[agg.predictor.isin(["raw_source_LOO_abs","noise_adjusted_ratio","control_variability"])].sort_values(["target_condition","source_condition","predictor"]).to_string(index=False))
print("\nCONTROL_ONLY_DIRECTION_NEGATIVE")
print(nda.to_string(index=False))
