from __future__ import annotations
from pathlib import Path
import hashlib,io,json
import numpy as np,pandas as pd

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_router'/'Dinneny2008_iron'
PRE=BASE/'pretarget_samples'; H48=BASE/'target_48h_samples'; OUT=BASE/'router_72h_blind'
OUT.mkdir(parents=True,exist_ok=True)
OLD=BASE/'router_blind'/'FROZEN_PREDICTIONS.tsv.gz'
RULE=ROOT/'来源结构路由规则_v0.2_冻结_2026-10-21.md'
BASELINE=['GSM265461','GSM265462']; H48S=['GSM265471','GSM265472']; H72S=['GSM265473','GSM265474']

def sha(p:Path):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(1<<20),b''): h.update(x)
 return h.hexdigest()

def read_soft(p:Path):
 lines=[];inside=False
 with p.open('r',encoding='utf-8',errors='strict') as f:
  for line in f:
   if line.startswith('!sample_table_begin'): inside=True; continue
   if line.startswith('!sample_table_end'): break
   if inside: lines.append(line)
 d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE'])
 gsm=p.stem.split('.')[0]; d=d.rename(columns={'ID_REF':'probe','VALUE':gsm}); d['probe']=d.probe.astype(str); d[gsm]=pd.to_numeric(d[gsm],errors='raise'); return d

# Hard target firewall: 72h files must not exist at freeze time.
for gsm in H72S:
 for p in [BASE/'target_72h_samples'/f'{gsm}.soft.txt', BASE/f'{gsm}.soft.txt']:
  if p.exists(): raise RuntimeError(f'72h target present before freeze: {p}')

old=pd.read_csv(OLD,sep='\t')
z=None
for gsm in BASELINE:
 p=PRE/f'{gsm}.soft.txt'; q=read_soft(p); z=q if z is None else z.merge(q,on='probe',validate='one_to_one')
for gsm in H48S:
 p=H48/f'{gsm}.soft.txt'
 if not p.exists(): raise FileNotFoundError(p)
 z=z.merge(read_soft(p),on='probe',validate='one_to_one')
z=old.merge(z,on='probe',how='inner',validate='one_to_one')
if len(z)!=len(old): raise RuntimeError(f'probe loss {len(old)}->{len(z)}')
for c in BASELINE+H48S:
 v=z[c].to_numpy(float)
 if np.any(~np.isfinite(v)) or np.any(v<=0): raise RuntimeError(f'invalid values {c}')
 z[c]=np.log2(v)
base=z[BASELINE].mean(axis=1).to_numpy(float)
r48=z[H48S].mean(axis=1).to_numpy(float)-base

out=old[['probe','shared_24h','persistent_24h','whole_12h','whole_6h','whole_3h']].copy()
out['persistent_48h']=r48
for c in ['shared_24h','persistent_24h','whole_12h','whole_6h','whole_3h','persistent_48h']:
 out[f'{c}_abs_rank']=pd.Series(np.abs(out[c].to_numpy(float))).rank(method='first',ascending=False).to_numpy(int)
out['primary_prediction']=out['persistent_48h']
out['primary_abs_rank']=out['persistent_48h_abs_rank']
pred=OUT/'FROZEN_72H_PREDICTIONS.tsv.gz'; out.to_csv(pred,sep='\t',index=False,compression='gzip')
manifest={
 'status':'FROZEN_BEFORE_72H_TARGET',
 'router_version':'v0.2',
 'route':'强状态特异型',
 'primary_predictor':'persistent_48h',
 'source_geometry_24h':{'rank1_energy':0.3716401632264418,'median_pairwise_rho':0.18651657864508348,'effective_state_fraction':0.5994434638991472},
 'n_frozen_probes':int(len(out)),
 'target_72h_samples':H72S,
 'target_72h_downloaded':False,
 'target_values_read':False,
 'router_rules_sha256':sha(RULE),
 'old_prediction_sha256':sha(OLD),
 'prediction_sha256':sha(pred),
 'baseline_sha256':{g:sha(PRE/f'{g}.soft.txt') for g in BASELINE},
 'history_48h_sha256':{g:sha(H48/f'{g}.soft.txt') for g in H48S},
 'candidate_unit':'ATH1 probe set'
}
(OUT/'FREEZE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('DINNENY2008_IRON_72H_FREEZE_SUCCESS'); print(json.dumps(manifest,ensure_ascii=False,indent=2)); print('72H TARGET FILES NOT PRESENT OR READ')
