from __future__ import annotations

from pathlib import Path
import hashlib
import io
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, hypergeom

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_router'/'Dinneny2008_salt'
PRE=BASE/'pretarget_samples'
TGT=BASE/'target_4h_samples'
OUT=BASE/'router_blind'
PRED=OUT/'FROZEN_PREDICTIONS.tsv.gz'
MAN=OUT/'FREEZE_MANIFEST.json'
TARGETS=['GSM184931','GSM184932']
BASELINE=['GSM184925','GSM184926']
MODELS=['shared_1h','persistent_1h','whole_0p5h']
KS=[100,250,500,1000]


def sha256(path:Path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for ch in iter(lambda:f.read(1<<20),b''):h.update(ch)
    return h.hexdigest()


def read_soft(path:Path):
    lines=[];inside=False
    with path.open('r',encoding='utf-8',errors='strict') as f:
        for line in f:
            if line.startswith('!sample_table_begin'):inside=True;continue
            if line.startswith('!sample_table_end'):break
            if inside:lines.append(line)
    d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE'])
    d=d.rename(columns={'ID_REF':'probe','VALUE':path.stem.split('.')[0]})
    d['probe']=d.probe.astype(str);d.iloc[:,1]=pd.to_numeric(d.iloc[:,1],errors='raise')
    return d

with MAN.open(encoding='utf-8') as f: man=json.load(f)
if man.get('status')!='FROZEN_BEFORE_4H_TARGET':raise RuntimeError('invalid freeze status')
if sha256(PRED)!=man['prediction_sha256']:raise RuntimeError('frozen prediction hash changed')
# The prospective router must remain ambiguous; scoring cannot silently change the route.
if man.get('route')!='路线歧义' or man.get('primary_predictor') is not None:
    raise RuntimeError('frozen route is not the expected ambiguous decision')

pred=pd.read_csv(PRED,sep='\t')
# Rebuild baseline only from pretarget files already hashed in freeze manifest.
allx=None
for gsm in BASELINE:
    p=PRE/f'{gsm}.soft.txt'
    if sha256(p)!=man['pretarget_sample_sha256'][gsm]:raise RuntimeError(f'pretarget file changed: {gsm}')
    q=read_soft(p)
    allx=q if allx is None else allx.merge(q,on='probe',how='inner',validate='one_to_one')
for gsm in TARGETS:
    p=TGT/f'{gsm}.soft.txt'
    if not p.exists():raise FileNotFoundError(p)
    q=read_soft(p)
    allx=allx.merge(q,on='probe',how='inner',validate='one_to_one')

z=pred.merge(allx,on='probe',how='inner',validate='one_to_one')
if len(z)!=len(pred):raise RuntimeError(f'frozen probe loss {len(pred)} -> {len(z)}')
for c in BASELINE+TARGETS:
    v=z[c].to_numpy(float)
    if np.any(~np.isfinite(v)) or np.any(v<=0):raise RuntimeError(f'invalid target/base values in {c}')
    z[c]=np.log2(v)
base=z[BASELINE].mean(axis=1).to_numpy(float)
rep={gsm:z[gsm].to_numpy(float)-base for gsm in TARGETS}
target=np.mean(np.column_stack([rep[x] for x in TARGETS]),axis=1)
z['target_4h']=target
N=len(z)
truth_rank=pd.Series(np.abs(target)).rank(method='first',ascending=False).to_numpy(int)

vrows=[];trows=[]
for model in MODELS:
    pv=z[model].to_numpy(float)
    rho=float(spearmanr(pv,target).statistic)
    mask=(pv!=0)&(target!=0)&np.isfinite(pv)&np.isfinite(target)
    sign=float(np.mean(np.sign(pv[mask])==np.sign(target[mask]))) if mask.any() else np.nan
    vrows.append({'model':model,'n_probes':N,'spearman_full_vector':rho,'sign_agreement_all_nonzero':sign})
    pr=pd.Series(np.abs(pv)).rank(method='first',ascending=False).to_numpy(int)
    for k in KS:
        kk=min(k,N);ps=set(np.where(pr<=kk)[0]);ts=set(np.where(truth_rank<=kk)[0]);ov=sorted(ps&ts);n=len(ov)
        exp=kk*kk/N
        trows.append({'model':model,'k':kk,'overlap':n,'random_expected_overlap':exp,'enrichment':n/exp,
                      'hypergeom_p':float(hypergeom.sf(n-1,N,kk,kk)),
                      'sign_agreement_overlap':float(np.mean(np.sign(pv[ov])==np.sign(target[ov]))) if ov else np.nan})

vdf=pd.DataFrame(vrows);tdf=pd.DataFrame(trows)
vdf.to_csv(OUT/'blind_4h_vector_results.tsv',sep='\t',index=False)
tdf.to_csv(OUT/'blind_4h_topk_results.tsv',sep='\t',index=False)
# Target replicate reference
va,vb=rep[TARGETS[0]],rep[TARGETS[1]]
rho_rep=float(spearmanr(va,vb).statistic)
rep_rows=[]
for k in [100,500]:
    ra=pd.Series(np.abs(va)).rank(method='first',ascending=False).to_numpy(int)
    rb=pd.Series(np.abs(vb)).rank(method='first',ascending=False).to_numpy(int)
    ov=len(set(np.where(ra<=k)[0])&set(np.where(rb<=k)[0]))
    rep_rows.append({'metric':f'top{k}','overlap':ov,'enrichment':ov/(k*k/N)})
pd.DataFrame(rep_rows).to_csv(OUT/'target_4h_replicate_topk.tsv',sep='\t',index=False)

best=vdf.sort_values('spearman_full_vector',ascending=False).iloc[0]
top100=tdf[tdf.k==100].sort_values('overlap',ascending=False)
summary={
    'freeze_route':man['route'],'freeze_primary_predictor':man['primary_predictor'],
    'n_probes':N,'target_replicate_rho':rho_rep,
    'model_full_rho':dict(zip(vdf.model,vdf.spearman_full_vector)),
    'model_top100_overlap':dict(zip(top100.model,top100.overlap.astype(int))),
    'model_top100_enrichment':dict(zip(top100.model,top100.enrichment)),
    'posthoc_best_full_vector':str(best['model']),
    'interpretation':'Calibration of ambiguous region only; no pre-reveal winner was selected.'
}
(OUT/'BLIND_4H_RESULT.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print('DINNENY2008_4H_SCORE_SUCCESS')
print(json.dumps(summary,ensure_ascii=False,indent=2))
print('\nVECTOR')
print(vdf.to_string(index=False))
print('\nTOPK')
print(tdf.to_string(index=False))
