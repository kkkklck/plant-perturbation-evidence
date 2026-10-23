from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026" / "tissuesystem_pseudobulk_v0"
DESIGN = ROOT / "data" / "external_rice" / "Robertson2026" / "sample_design.tsv"
OUT = ROOT / "data" / "external_rice" / "Robertson2026" / "tissuesystem_structure_v0"
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(BASE / "group_metadata.tsv", sep="\t")
expr = pd.read_csv(BASE / "tissuesystem_log1p_cpm.tsv.gz", sep="\t")
design = pd.read_csv(DESIGN, sep="\t")
genes = expr.pop("gene").astype(str).to_numpy()
X = expr.to_numpy(float)
if list(expr.columns.astype(str)) != list(meta.group.astype(str)):
    raise RuntimeError("Expression columns do not match metadata")

mean_expr = X.mean(axis=1)
nonzero_frac = (X > 0).mean(axis=1)
mask = (nonzero_frac >= 0.25) & (mean_expr >= 0.05)

sample_batch = dict(zip(design.identity.astype(str), design.batch.astype(int)))
conditions = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
control = "Well.Watered"
states = sorted(meta.tissue_system.unique())
ctrl_by_batch = {}
for s in meta.loc[meta.condition.eq(control), "sample"].astype(str).unique():
    ctrl_by_batch[sample_batch[s]] = s

def idx(sample, state):
    z = meta.index[(meta["sample"].astype(str) == sample) & (meta.tissue_system == state)].to_numpy()
    return z[0] if len(z) == 1 else None

def matched_vecs(cond, state):
    out=[]; labels=[]
    for ti in meta.index[(meta.condition.eq(cond)) & (meta.tissue_system.eq(state))].to_numpy():
        ss = str(meta.loc[ti,"sample"]); b=sample_batch[ss]; cs=ctrl_by_batch[b]; ci=idx(cs,state)
        if ci is None: continue
        out.append((X[:,ti]-X[:,ci])[mask]); labels.append((ss,cs,b))
    return out, labels

# baseline noise from 3 controls per tissue system
noise = {}
ctrl_samples = sorted(meta.loc[meta.condition.eq(control),"sample"].astype(str).unique())
for st in states:
    vals=[]; rhos=[]
    for i in range(len(ctrl_samples)):
        for j in range(i+1,len(ctrl_samples)):
            a=idx(ctrl_samples[i],st); b=idx(ctrl_samples[j],st)
            va=X[:,a][mask]; vb=X[:,b][mask]
            vals.append(float(np.sqrt(np.mean((va-vb)**2))))
            rhos.append(float(spearmanr(va,vb).statistic))
    noise[st]=(float(np.median(vals)), float(np.median(rhos)))

summary=[]; state_rows=[]
for cond in conditions:
    reps_by_state={}
    avg=[]; used=[]
    for st in states:
        vs, labs = matched_vecs(cond,st)
        if not vs: continue
        reps_by_state[st]=(vs,labs)
        avg.append(np.vstack(vs).mean(axis=0)); used.append(st)
    R=np.vstack(avg)
    u,s,vt=np.linalg.svd(R,full_matrices=False)
    e=s*s; frac=e/e.sum(); axis=vt[0]
    if np.dot(axis,R.mean(axis=0))<0: axis=-axis
    prs=[spearmanr(R[i],R[j]).statistic for i in range(len(used)) for j in range(i+1,len(used))]
    loo=[]; proj=[]
    for i,st in enumerate(used):
        l=np.delete(R,i,axis=0).mean(axis=0)
        rr=float(spearmanr(R[i],l).statistic); loo.append(rr)
        denom=float(np.dot(R[i],R[i])); p=float(np.dot(R[i],axis)); pf=(p*p/denom) if denom else np.nan; proj.append(pf)
        vs,labs=reps_by_state[st]
        rep_rho=float(spearmanr(vs[0],vs[1]).statistic) if len(vs)>=2 else np.nan
        signal=float(np.mean([np.sqrt(np.mean(v*v)) for v in vs]))
        nr,cr=noise[st]
        state_rows.append(dict(condition=cond,tissue_system=st,rho_to_LOO_shared=rr,projection_energy_fraction=pf,replicate_response_rho=rep_rho,response_rms=signal,baseline_control_control_rms=nr,response_to_noise_ratio=signal/nr if nr>0 else np.nan,baseline_control_expression_rho=cr))
    # replicate-specific consensus + rank1 stability across the four tissue systems
    # all conditions have two matched replicate responses per state
    R1=np.vstack([reps_by_state[st][0][0] for st in used]); R2=np.vstack([reps_by_state[st][0][1] for st in used])
    c1=R1.mean(axis=0); c2=R2.mean(axis=0)
    _,ss1,v1=np.linalg.svd(R1,full_matrices=False); _,ss2,v2=np.linalg.svd(R2,full_matrices=False)
    a1=v1[0]; a2=v2[0]
    if np.dot(a1,c1)<0:a1=-a1
    if np.dot(a2,c2)<0:a2=-a2
    summary.append(dict(condition=cond,n_states=len(used),n_genes=int(mask.sum()),rank1_energy_fraction=float(frac[0]),rank3_cumulative_energy_fraction=float(frac[:3].sum()),median_pairwise_state_spearman=float(np.nanmedian(prs)),median_state_rho_to_LOO_shared=float(np.nanmedian(loo)),median_projection_energy_fraction=float(np.nanmedian(proj)),median_state_replicate_rho=float(np.nanmedian([spearmanr(reps_by_state[st][0][0],reps_by_state[st][0][1]).statistic for st in used])),consensus_response_replicate_rho=float(spearmanr(c1,c2).statistic),rank1_axis_replicate_rho=float(spearmanr(a1,a2).statistic),median_response_to_noise=float(np.nanmedian([r["response_to_noise_ratio"] for r in state_rows if r["condition"]==cond]))))

sdf=pd.DataFrame(summary); stf=pd.DataFrame(state_rows)
sdf.to_csv(OUT/"summary.tsv",sep="\t",index=False)
stf.to_csv(OUT/"state_metrics.tsv",sep="\t",index=False)
print("ROBERTSON2026_TISSUESYSTEM_STRUCTURE_SUCCESS")
print(sdf.to_string(index=False))
print("\nSTATE_METRICS")
print(stf.to_string(index=False))
