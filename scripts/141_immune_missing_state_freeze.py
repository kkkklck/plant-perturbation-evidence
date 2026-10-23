from __future__ import annotations
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data' / 'prospective_task' / 'immune_missing_state'
SRC = BASE / 'source_only'
TGT = BASE / 'target_only'
OUT = BASE / 'blind_v0_1'
OUT.mkdir(parents=True, exist_ok=True)
RULE = ROOT / '免疫刺激_未见状态_任务结构盲测规则_2026-10-21.md'
TASK_RULE = ROOT / '任务结构优先证据选择规则_v0.1_冻结_2026-10-21.md'

STIMS = ['flg22','Pep1','Chitin','SCOOP10','SCREW1']
SOURCE_STATES = [
    'Rootcap_Columella',
    'Meristem_Epidermis','Meristem_Endodermis','Meristem_Stele',
    'Elongation_Epidermis','Elongation_Cortex','Elongation_Stele',
    'Maturation_Epidermis','Maturation_Cortex','Maturation_Endodermis',
]
TARGET_STATES = ['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']
THRESH = {'rank1_energy_min':0.50,'median_pairwise_rho_min':0.35,'effective_state_fraction_min':0.50}


def sha256(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''):
            h.update(b)
    return h.hexdigest()


def read_state(p: Path) -> pd.Series:
    d=pd.read_csv(p,sep='\t',usecols=['gene','avg_logFC'])
    d['gene']=d['gene'].astype(str)
    d['avg_logFC']=pd.to_numeric(d['avg_logFC'],errors='raise')
    d=d.drop_duplicates('gene',keep='first').set_index('gene')['avg_logFC']
    if not np.isfinite(d.to_numpy(float)).all():
        raise RuntimeError(f'nonfinite values in {p}')
    return d


def geometry(M: np.ndarray) -> dict:
    # M = genes x source states; identical to earlier sharedness diagnostic convention.
    u,s,vt=np.linalg.svd(M,full_matrices=False)
    e=(s*s)/(s*s).sum()
    p=vt[0]**2
    p=p/p.sum()
    prs=[]; loo=[]
    for i in range(M.shape[1]):
        for j in range(i+1,M.shape[1]):
            prs.append(float(spearmanr(M[:,i],M[:,j]).statistic))
        oth=np.mean(np.delete(M,i,axis=1),axis=1)
        loo.append(float(spearmanr(M[:,i],oth).statistic))
    eff=float(1/np.sum(p*p))
    return {
        'n_genes':int(M.shape[0]),
        'n_source_states':int(M.shape[1]),
        'rank1_energy':float(e[0]),
        'rank1_3_energy':float(e[:min(3,len(e))].sum()),
        'median_pairwise_rho':float(np.nanmedian(prs)),
        'median_loo_shared_rho':float(np.nanmedian(loo)),
        'max_rank1_loading_share':float(np.max(p)),
        'effective_states':eff,
        'effective_state_fraction':eff/M.shape[1],
    }

if TGT.exists() and any(TGT.rglob('*')):
    raise RuntimeError('target_only exists before freeze; abort')
if not SRC.exists():
    raise FileNotFoundError(SRC)

script_paths=[ROOT/'scripts'/'140_immune_missing_state_source_firewall.R',
              ROOT/'scripts'/'141_immune_missing_state_freeze.py',
              ROOT/'scripts'/'142_immune_missing_state_target_firewall.R',
              ROOT/'scripts'/'143_immune_missing_state_score.py']
for p in [RULE,TASK_RULE,*script_paths]:
    if not p.exists(): raise FileNotFoundError(p)

routes={}; source_hashes={}; prediction_hashes={}
for stim in STIMS:
    vals={}
    source_hashes[stim]={}
    for state in SOURCE_STATES:
        p=SRC/stim/f'{state}.tsv.gz'
        if not p.exists(): raise FileNotFoundError(p)
        source_hashes[stim][state]=sha256(p)
        vals[state]=read_state(p)
    common=set(vals[SOURCE_STATES[0]].index)
    for state in SOURCE_STATES[1:]:
        common &= set(vals[state].index)
    genes=sorted(common)
    if len(genes)<1000:
        raise RuntimeError(f'too few common source genes for {stim}: {len(genes)}')
    M=np.column_stack([vals[s].reindex(genes).to_numpy(float) for s in SOURCE_STATES])
    g=geometry(M)
    allowed=(g['rank1_energy']>=THRESH['rank1_energy_min'] and
             g['median_pairwise_rho']>=THRESH['median_pairwise_rho_min'] and
             g['effective_state_fraction']>=THRESH['effective_state_fraction_min'])
    pred=np.mean(M,axis=1)
    df=pd.DataFrame({'gene':genes,'predicted_response':pred})
    df['priority_score']=df['predicted_response'].abs()
    df=df.sort_values(['priority_score','gene'],ascending=[False,True],kind='mergesort').reset_index(drop=True)
    df['priority_rank']=np.arange(1,len(df)+1)
    df['predicted_direction']=np.where(df['predicted_response']>0,'up',np.where(df['predicted_response']<0,'down','flat'))
    df['gene_level_output_allowed']=bool(allowed)
    op=OUT/f'FROZEN_{stim}_PREDICTION.tsv.gz'
    df.to_csv(op,sep='\t',index=False,compression='gzip')
    prediction_hashes[stim]=sha256(op)
    routes[stim]={
        'task_type':'missing_state',
        'route':'ALLOW_SHARED_GENE_LEVEL' if allowed else 'REJECT_GENE_LEVEL',
        'gene_level_output_allowed':bool(allowed),
        'source_geometry':g,
        'thresholds':THRESH,
        'n_frozen_genes':len(genes),
        'prediction_file':str(op.relative_to(ROOT)),
        'prediction_sha256':prediction_hashes[stim],
    }

manifest={
    'status':'FROZEN_BEFORE_TARGET_REVEAL',
    'rule_version':'task_structure_v0.1',
    'task_type':'missing_state',
    'stimuli':STIMS,
    'source_states':SOURCE_STATES,
    'target_states_withheld':TARGET_STATES,
    'shared_thresholds':THRESH,
    'routes':routes,
    'source_file_sha256':source_hashes,
    'rule_sha256':sha256(RULE),
    'task_rule_sha256':sha256(TASK_RULE),
    'script_sha256':{p.name:sha256(p) for p in script_paths},
    'target_files_present':False,
    'target_values_read_by_prediction':False,
    'upstream_monolithic_rds_firewall_limitation':True,
    'upstream_target_values_used_for_any_decision':False,
}
(OUT/'FREEZE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('IMMUNE_MISSING_STATE_FREEZE_SUCCESS')
for stim in STIMS:
    r=routes[stim]
    g=r['source_geometry']
    print(f"{stim}\t{r['route']}\tgenes={r['n_frozen_genes']}\trank1={g['rank1_energy']:.4f}\tpair={g['median_pairwise_rho']:.4f}\teff_frac={g['effective_state_fraction']:.4f}\tloo={g['median_loo_shared_rho']:.4f}")
print('TARGET_FILES_PRESENT=0')
