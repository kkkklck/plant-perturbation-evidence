from __future__ import annotations
from pathlib import Path
import hashlib, io, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'data'/'prospective_router'/'GSE46205_salt'
SRC = BASE/'source_1h_samples'
OUT = BASE/'router_blind'
OUT.mkdir(parents=True, exist_ok=True)
RULES = ROOT/'GSE46205_自动路由独立前瞻盲测规则_2026-10-21.md'
ROUTER = ROOT/'来源结构路由规则_v0.2_冻结_2026-10-21.md'
TARGET_DIR = BASE/'target_3h_samples'

GROUPS = {
    'epidermis': {
        'control':['GSM1126230','GSM1126231','GSM1126232'],
        'salt1h':['GSM1126236','GSM1126237','GSM1126238'],
        'target3h':['GSM1126239','GSM1126240','GSM1126241'],
    },
    'cortex': {
        'control':['GSM1126254','GSM1126255','GSM1126256'],
        'salt1h':['GSM1126260','GSM1126261','GSM1126262'],
        'target3h':['GSM1126263','GSM1126264','GSM1126265'],
    },
    'endodermis': {
        'control':['GSM1126278','GSM1126279','GSM1126280'],
        'salt1h':['GSM1126284','GSM1126285','GSM1126286'],
        'target3h':['GSM1126287','GSM1126288','GSM1126289'],
    },
    'stele': {
        'control':['GSM1126302','GSM1126303','GSM1126304'],
        'salt1h':['GSM1126308','GSM1126309','GSM1126310'],
        'target3h':['GSM1126311','GSM1126312','GSM1126313'],
    },
}

THR = {
    'shared': {'rank1_energy_min':0.50,'median_pairwise_rho_min':0.35,'effective_state_fraction_min':0.50},
    'state_specific': {'rank1_energy_max':0.40,'median_pairwise_rho_max':0.20},
}

def sha(p:Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for x in iter(lambda:f.read(1<<20),b''): h.update(x)
    return h.hexdigest()

def read_soft(p:Path):
    lines=[]; inside=False
    with p.open('r',encoding='utf-8',errors='strict') as f:
        for line in f:
            if line.startswith('!sample_table_begin'): inside=True; continue
            if line.startswith('!sample_table_end'): break
            if inside: lines.append(line)
    d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE'])
    gsm=p.stem.split('.')[0]
    d=d.rename(columns={'ID_REF':'probe','VALUE':gsm})
    d['probe']=d.probe.astype(str); d[gsm]=pd.to_numeric(d[gsm],errors='raise')
    return d

# Hard firewall: target directory/files must not exist before freeze.
if TARGET_DIR.exists() and any(TARGET_DIR.iterdir()):
    raise RuntimeError('3h target files present before freeze')
for ct,g in GROUPS.items():
    for gsm in g['target3h']:
        if (BASE/f'{gsm}.soft.txt').exists(): raise RuntimeError(f'target present before freeze: {gsm}')

# Read all 24 allowed source arrays and require exact common probe set.
all_ids=[]
for g in GROUPS.values(): all_ids += g['control'] + g['salt1h']
z=None; source_hash={}
for gsm in all_ids:
    p=SRC/f'{gsm}.soft.txt'
    if not p.exists(): raise FileNotFoundError(p)
    source_hash[gsm]=sha(p)
    q=read_soft(p)
    z=q if z is None else z.merge(q,on='probe',how='inner',validate='one_to_one')
if z is None: raise RuntimeError('no source data')
if len(z)<20000: raise RuntimeError(f'unexpectedly small common probe universe: {len(z)}')
for gsm in all_ids:
    v=z[gsm].to_numpy(float)
    if np.any(~np.isfinite(v)) or np.any(v<=0): raise RuntimeError(f'invalid values {gsm}')
    z[gsm]=np.log2(v)

cts=list(GROUPS)
R=[]
for ct in cts:
    g=GROUPS[ct]
    ctrl=z[g['control']].mean(axis=1).to_numpy(float)
    salt=z[g['salt1h']].mean(axis=1).to_numpy(float)
    R.append(salt-ctrl)
M=np.column_stack(R)  # probes x states

# Geometry uses states x probes.
A=M.T
u,s,vt=np.linalg.svd(A,full_matrices=False)
e=s*s/(s*s).sum()
p=u[:,0]**2; p=p/p.sum()
pairs=[]; loos=[]
for i in range(A.shape[0]):
    for j in range(i+1,A.shape[0]): pairs.append(float(spearmanr(A[i],A[j]).statistic))
    oth=np.mean(np.delete(A,i,axis=0),axis=0); loos.append(float(spearmanr(A[i],oth).statistic))
geom={
    'rank1_energy':float(e[0]),
    'rank1_3_energy':float(e[:3].sum()),
    'median_pairwise_rho':float(np.nanmedian(pairs)),
    'median_loo_shared_rho':float(np.nanmedian(loos)),
    'max_rank1_loading_share':float(p.max()),
    'effective_states':float(1/np.sum(p*p)),
    'effective_state_fraction':float((1/np.sum(p*p))/len(cts)),
}
shared=(geom['rank1_energy']>=0.50 and geom['median_pairwise_rho']>=0.35 and geom['effective_state_fraction']>=0.50)
state=(geom['rank1_energy']<=0.40 and geom['median_pairwise_rho']<=0.20)
if shared: route='强共享型'; primary='shared_other3_1h'
elif state: route='强状态特异型'; primary='same_state_1h'
else: route='路线歧义'; primary=None

rows=[]
for j,ct in enumerate(cts):
    same=M[:,j]
    shared_other=np.mean(np.delete(M,j,axis=1),axis=1)
    q=pd.DataFrame({'probe':z.probe.astype(str),'target_celltype':ct,'same_state_1h':same,'shared_other3_1h':shared_other})
    q['same_state_1h_abs_rank']=pd.Series(np.abs(same)).rank(method='first',ascending=False).to_numpy(int)
    q['shared_other3_1h_abs_rank']=pd.Series(np.abs(shared_other)).rank(method='first',ascending=False).to_numpy(int)
    if primary is None:
        q['primary_prediction']=np.nan; q['primary_abs_rank']=np.nan
    else:
        q['primary_prediction']=q[primary]
        q['primary_abs_rank']=q[f'{primary}_abs_rank']
    rows.append(q)
out=pd.concat(rows,ignore_index=True)
pred=OUT/'FROZEN_PREDICTIONS.tsv.gz'; out.to_csv(pred,sep='\t',index=False,compression='gzip')
manifest={
    'status':'FROZEN_BEFORE_3H_TARGET',
    'router_version':'v0.2',
    'route':route,
    'primary_predictor':primary,
    'source_geometry':geom,
    'thresholds':THR,
    'n_celltypes':len(cts),
    'celltypes':cts,
    'n_probes_per_celltype':int(len(z)),
    'target_samples':{ct:g['target3h'] for ct,g in GROUPS.items()},
    'target_files_present':False,
    'target_values_read':False,
    'rules_sha256':sha(RULES),
    'router_rules_sha256':sha(ROUTER),
    'source_sample_sha256':source_hash,
    'prediction_sha256':sha(pred),
    'freeze_script_sha256':sha(Path(__file__)),
    'candidate_unit':'ATH1 probe set'
}
(OUT/'FREEZE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'ROUTE_DECISION.json').write_text(json.dumps({'route':route,'primary_predictor':primary,'source_geometry':geom},ensure_ascii=False,indent=2),encoding='utf-8')
print('GSE46205_ROUTER_FREEZE_SUCCESS'); print(json.dumps(manifest,ensure_ascii=False,indent=2)); print('3H TARGET FILES NOT PRESENT OR READ')
