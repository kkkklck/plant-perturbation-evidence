from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026"
PB = BASE / "pseudobulk_v0"
DE = BASE / "deseq2_v0"
OUT = BASE / "deg_focused_structure_v0"
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(PB / "group_metadata.tsv", sep="\t")
design = pd.read_csv(BASE / "sample_design.tsv", sep="\t")
expr = pd.read_csv(PB / "pseudobulk_log1p_cpm.tsv.gz", sep="\t", compression="gzip")
genes = expr.pop("gene").astype(str).to_numpy()
X = expr.to_numpy(float)
if expr.columns.astype(str).tolist() != meta.group.astype(str).tolist():
    raise RuntimeError("group order mismatch")

membership = pd.read_csv(DE / "condition_deg_union_membership.tsv", sep="\t")
annotations = sorted(meta.cell_annotation.unique())
col_index = {g:i for i,g in enumerate(meta.group.astype(str))}
control_by_batch = dict(zip(design.loc[design.stress_type=="control","batch"].astype(int), design.loc[design.stress_type=="control","identity"].astype(str)))

stress_map = [
    ("Mild.Drought","MildDrought"),
    ("Moderate.Drought","ModerateDrought"),
    ("Mild.Salinity","MildSalinity"),
    ("Moderate.Salinity","ModerateSalinity"),
]

def vec(sample, ann):
    return X[:, col_index[f"{sample}||{ann}"]]

def response_matrix(mode, cond_label, design_cond):
    stress_samples = design.loc[design.condition==design_cond,"identity"].astype(str).tolist()
    R=[]; paired_by_state={}
    for ann in annotations:
        if mode == "author_like":
            ctrl_samples = design.loc[design.stress_type=="control","identity"].astype(str).tolist()
            d = np.mean([vec(s,ann) for s in stress_samples],axis=0) - np.mean([vec(s,ann) for s in ctrl_samples],axis=0)
            R.append(d)
        elif mode == "batch_matched":
            pairs=[]
            for s in stress_samples:
                row=design.loc[design.identity==s].iloc[0]
                c=control_by_batch[int(row.batch)]
                pairs.append(vec(s,ann)-vec(c,ann))
            pairs=np.vstack(pairs)
            R.append(pairs.mean(axis=0)); paired_by_state[ann]=pairs
        else:
            raise ValueError(mode)
    return np.vstack(R), paired_by_state

def metrics(M):
    if M.shape[1] < 2:
        return dict(rank1=np.nan,rank3=np.nan,median_pairwise=np.nan,median_loo=np.nan,median_proj=np.nan)
    u,s,vt=np.linalg.svd(M,full_matrices=False)
    frac=s*s/np.sum(s*s)
    shared=vt[0].copy()
    if np.dot(shared,M.mean(axis=0))<0: shared*=-1
    shared/=np.linalg.norm(shared)
    prs=[]; loos=[]; projs=[]
    for i in range(M.shape[0]):
        for j in range(i+1,M.shape[0]):
            prs.append(spearmanr(M[i],M[j]).statistic)
        loo=np.delete(M,i,axis=0).mean(axis=0)
        loos.append(spearmanr(M[i],loo).statistic)
        den=float(np.dot(M[i],M[i])); p=float(np.dot(M[i],shared))
        projs.append(p*p/den if den>0 else np.nan)
    return dict(
        rank1=float(frac[0]),
        rank3=float(frac[:min(3,len(frac))].sum()),
        median_pairwise=float(np.nanmedian(prs)),
        median_loo=float(np.nanmedian(loos)),
        median_proj=float(np.nanmedian(projs)),
    )

rows=[]; rep_rows=[]
for mode in ["author_like","batch_matched"]:
    for cond_label,design_cond in stress_map:
        R,pairs=response_matrix(mode,cond_label,design_cond)
        m=membership[(membership["mode"]==mode)&(membership["condition"]==cond_label)].copy()
        for panel_name,sel in [
            ("union_DEG",m.n_annotations_strong>=1),
            ("shared_ge2",m.n_annotations_strong>=2),
            ("shared_ge3",m.n_annotations_strong>=3),
        ]:
            pg=set(m.loc[sel,"gene"].astype(str))
            mask=np.array([g in pg for g in genes])
            if mask.sum()==0: continue
            met=metrics(R[:,mask])
            rows.append({"mode":mode,"condition":cond_label,"panel":panel_name,"n_genes":int(mask.sum()),**met})
            if mode=="batch_matched":
                vals=[]
                for ann in annotations:
                    P=pairs[ann][:,mask]
                    vals.append(spearmanr(P[0],P[1]).statistic)
                rep_rows.append({
                    "condition":cond_label,"panel":panel_name,"n_genes":int(mask.sum()),
                    "median_matched_replicate_rho":float(np.nanmedian(vals)),
                    "mean_matched_replicate_rho":float(np.nanmean(vals)),
                    "min_matched_replicate_rho":float(np.nanmin(vals)),
                    "max_matched_replicate_rho":float(np.nanmax(vals)),
                })

out=pd.DataFrame(rows)
rep=pd.DataFrame(rep_rows)
out.to_csv(OUT/"deg_focused_structure.tsv",sep="\t",index=False)
rep.to_csv(OUT/"deg_focused_matched_replicate_stability.tsv",sep="\t",index=False)

print("ROBERTSON2026_DEG_FOCUSED_SUCCESS")
print(out.to_string(index=False))
print("\nBATCH_MATCHED_REPLICATE_STABILITY_ON_DEG_PANELS")
print(rep.to_string(index=False))
