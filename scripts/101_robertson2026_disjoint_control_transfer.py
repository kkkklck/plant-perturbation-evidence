from __future__ import annotations
from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"data"/"external_rice"/"Robertson2026"; PB=BASE/"pseudobulk_v0"; DE=BASE/"deseq2_v0"; OUT=BASE/"cross_condition_prioritization_v0"
meta=pd.read_csv(PB/"group_metadata.tsv",sep="\t"); design=pd.read_csv(BASE/"sample_design.tsv",sep="\t")
expr=pd.read_csv(PB/"pseudobulk_log1p_cpm.tsv.gz",sep="\t",compression="gzip"); genes=expr.pop("gene").astype(str).to_numpy(); X=expr.to_numpy(float)
idx={g:i for i,g in enumerate(meta.group.astype(str))}; anns=sorted(meta.cell_annotation.unique())
sig=pd.read_csv(DE/"deseq2_fdr05_genes.tsv.gz",sep="\t",compression="gzip")
if sig["strong"].dtype!=bool:sig["strong"]=sig["strong"].astype(str).str.upper().eq("TRUE")
controls={1:"Well.Watered.R1",2:"Well.Watered.R2",3:"Well.Watered.R3"}
cond_map={"Mild.Drought":"MildDrought","Moderate.Drought":"ModerateDrought","Mild.Salinity":"MildSalinity","Moderate.Salinity":"ModerateSalinity"}

def vec(s,a):return X[:,idx[f"{s}||{a}"]]
def auroc(y,s):
 y=np.asarray(y,bool);n1=int(y.sum());n0=int((~y).sum());r=rankdata(s,method="average");return float((r[y].sum()-n1*(n1+1)/2)/(n1*n0)) if n1 and n0 else np.nan
def top10(y,s):
 n=len(y);npos=int(y.sum());k=math.ceil(.1*n);o=np.argsort(-s,kind="mergesort")[:k];h=int(y[o].sum());return h/npos,(h/k)/(npos/n)

def target_batches(label):
 ss=design.loc[design.condition==cond_map[label]]
 return set(ss.batch.astype(int))

def disjoint_control_for_target(label):
 unused=sorted(set(controls)-target_batches(label))
 if len(unused)!=1:raise RuntimeError((label,unused))
 return controls[unused[0]],unused[0]

rows=[]
for target in cond_map:
 dc,db=disjoint_control_for_target(target)
 for source in cond_map:
  if source==target:continue
  ss=design.loc[design.condition==cond_map[source],"identity"].astype(str).tolist()
  Rsrc=[]
  for a in anns:
   stress=np.mean([vec(s,a) for s in ss],axis=0)
   Rsrc.append(stress-vec(dc,a))
  Rsrc=np.vstack(Rsrc)
  for i,a in enumerate(anns):
   truth=sig[(sig["mode"]=="batch_matched")&(sig["condition"]==target)&(sig["annotation"]==a)&sig["strong"]]
   if len(truth)==0:continue
   pos=set(truth.gene.astype(str));y=np.array([g in pos for g in genes],bool);lfc=dict(zip(truth.gene.astype(str),truth.log2FoldChange.astype(float)));gi=np.where(y)[0];actual=np.sign([lfc[genes[j]] for j in gi])
   preds={"disjoint_control_source_same_state":Rsrc[i],"disjoint_control_source_LOO_shared":np.delete(Rsrc,i,axis=0).mean(axis=0)}
   for name,signed in preds.items():
    score=np.abs(signed);rec,enr=top10(y,score);ps=np.sign(signed[gi]);nz=actual!=0
    rows.append({"source_condition":source,"target_condition":target,"disjoint_control":dc,"disjoint_control_batch":db,"cell_annotation":a,"predictor":name,"n_true_deg":int(y.sum()),"auroc":auroc(y,score),"top10_recall":rec,"top10_enrichment":enr,"direction_agreement":float(np.mean(ps[nz]==actual[nz]))})
res=pd.DataFrame(rows);res.to_csv(OUT/"disjoint_control_cross_condition_per_state.tsv",sep="\t",index=False)
agg=res.groupby(["source_condition","target_condition","disjoint_control","predictor"]).agg(n_states=("cell_annotation","nunique"),median_auroc=("auroc","median"),median_top10_recall=("top10_recall","median"),median_top10_enrichment=("top10_enrichment","median"),median_direction_agreement=("direction_agreement","median")).reset_index()
agg.to_csv(OUT/"disjoint_control_cross_condition_summary.tsv",sep="\t",index=False)
print("ROBERTSON2026_DISJOINT_CONTROL_TRANSFER_SUCCESS")
print(agg.to_string(index=False))
