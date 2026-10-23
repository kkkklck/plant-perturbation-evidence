from __future__ import annotations

from pathlib import Path
import gzip, io, hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_router'/'Dinneny2008_iron'
FE_MATRIX=BASE/'GSE10501_series_matrix.txt.gz'
STD_MATRIX=ROOT/'data'/'prospective_router'/'Dinneny2008_salt'/'GSE7641_series_matrix.txt.gz'
PRE=BASE/'pretarget_samples'
OUT=BASE/'router_blind';OUT.mkdir(parents=True,exist_ok=True)
RULES=ROOT/'Dinneny2008_缺铁自动路由前瞻盲测规则_2026-10-20.md'
ROUTER_RULES=ROOT/'来源结构路由规则_v0.1_冻结_2026-10-20.md'

STD={
'epidermis_lrc':['GSM184889','GSM184890','GSM184891'],
'columella':['GSM184892','GSM184893','GSM184894'],
'cortex':['GSM184895','GSM184896','GSM184897'],
'endodermis_qc':['GSM184898','GSM184899','GSM184900'],
'stele':['GSM184901','GSM184902','GSM184903'],
}
FE={
'epidermis_lrc':['GSM266662','GSM266663','GSM266664','GSM266665'],
'columella':['GSM266666','GSM266667','GSM266668'],
'cortex':['GSM266669','GSM266670','GSM266671'],
'endodermis_qc':['GSM266672','GSM266673','GSM266674'],
'stele':['GSM266675','GSM266676','GSM266677'],
}
WHOLE={
'baseline':['GSM265461','GSM265462'],
'3h':['GSM265463','GSM265464'],
'6h':['GSM265465','GSM265466'],
'12h':['GSM265467','GSM265468'],
'24h':['GSM265469','GSM265470'],
}
TARGETS=['GSM265471','GSM265472','GSM265473','GSM265474']

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for ch in iter(lambda:f.read(1<<20),b''):h.update(ch)
 return h.hexdigest()

def read_matrix(path):
 lines=[];inside=False
 with gzip.open(path,'rt',encoding='utf-8',errors='strict') as f:
  for line in f:
   if line.startswith('!series_matrix_table_begin'):inside=True;continue
   if line.startswith('!series_matrix_table_end'):break
   if inside:lines.append(line)
 d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t')
 d=d.rename(columns={d.columns[0]:'probe'});d['probe']=d.probe.astype(str).str.strip().str.strip('"')
 return d

def read_soft(path):
 lines=[];inside=False
 with path.open('r',encoding='utf-8',errors='strict') as f:
  for line in f:
   if line.startswith('!sample_table_begin'):inside=True;continue
   if line.startswith('!sample_table_end'):break
   if inside:lines.append(line)
 d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE'])
 gsm=path.stem.split('.')[0];d=d.rename(columns={'ID_REF':'probe','VALUE':gsm});d['probe']=d.probe.astype(str);d[gsm]=pd.to_numeric(d[gsm],errors='raise')
 return d

def geometry(M):
 X=np.asarray(M,float).T
 u,s,vt=np.linalg.svd(X,full_matrices=False);e=s*s/(s*s).sum();p=u[:,0]**2;p=p/p.sum()
 prs=[];loo=[]
 for i in range(X.shape[0]):
  for j in range(i+1,X.shape[0]):prs.append(float(spearmanr(X[i],X[j]).statistic))
  other=X[np.arange(X.shape[0])!=i].mean(axis=0);loo.append(float(spearmanr(X[i],other).statistic))
 eff=float(1/np.sum(p*p))
 return {'rank1_energy':float(e[0]),'rank1_3_energy':float(e[:min(3,len(e))].sum()),
 'median_pairwise_rho':float(np.median(prs)),'median_loo_shared_rho':float(np.median(loo)),
 'max_rank1_loading_share':float(p.max()),'effective_states':eff,'effective_state_fraction':eff/X.shape[0]}

def main():
 # No target files or full GSE10502 matrix/family allowed before freeze.
 for gsm in TARGETS:
  if any(BASE.rglob(f'*{gsm}*')): raise RuntimeError(f'TARGET PRESENT BEFORE FREEZE {gsm}')
 for bad in ['GSE10502_series_matrix.txt.gz','GSE10502_family.soft.gz']:
  if (BASE/bad).exists():raise RuntimeError(f'TARGET-CONTAINING FILE PRESENT {bad}')
 std=read_matrix(STD_MATRIX);fe=read_matrix(FE_MATRIX)
 scols=[x for v in STD.values() for x in v];fcols=[x for v in FE.values() for x in v]
 for c in scols:
  if c not in std.columns:raise RuntimeError(f'missing std {c}')
  std[c]=pd.to_numeric(std[c],errors='raise')
 for c in fcols:
  if c not in fe.columns:raise RuntimeError(f'missing fe {c}')
  fe[c]=pd.to_numeric(fe[c],errors='raise')
 z=std[['probe']+scols].merge(fe[['probe']+fcols],on='probe',how='inner',validate='one_to_one')
 whole=None;hashes={}
 for group in WHOLE.values():
  for gsm in group:
   p=PRE/f'{gsm}.soft.txt';hashes[gsm]=sha(p);q=read_soft(p)
   whole=q if whole is None else whole.merge(q,on='probe',how='inner',validate='one_to_one')
 z=z.merge(whole,on='probe',how='inner',validate='one_to_one')
 z=z[~z.probe.str.upper().str.startswith('AFFX-')].copy()
 allcols=scols+fcols+[x for v in WHOLE.values() for x in v]
 valid=np.all(np.isfinite(z[allcols].to_numpy(float))&(z[allcols].to_numpy(float)>0),axis=1)
 z=z.loc[valid].reset_index(drop=True)
 for c in allcols:z[c]=np.log2(z[c].to_numpy(float))
 states=list(STD);R=[]
 for st in states:
  r=z[FE[st]].mean(axis=1).to_numpy(float)-z[STD[st]].mean(axis=1).to_numpy(float)
  z[f'resp_{st}']=r;R.append(r)
 R=np.column_stack(R);geom=geometry(R)
 base=z[WHOLE['baseline']].mean(axis=1).to_numpy(float)
 for t in ['3h','6h','12h','24h']:
  z[f'whole_{t}']=z[WHOLE[t]].mean(axis=1).to_numpy(float)-base
 z['shared_24h']=R.mean(axis=1);z['persistent_24h']=z['whole_24h']
 is_shared=geom['rank1_energy']>=.50 and geom['median_pairwise_rho']>=.35 and geom['effective_state_fraction']>=.50
 is_state=geom['rank1_energy']<=.40 and geom['median_pairwise_rho']<=.15 and geom['effective_state_fraction']<=.30
 if is_shared:route='强共享型';primary='shared_24h'
 elif is_state:route='强状态特异型';primary='persistent_24h'
 else:route='路线歧义';primary=None
 for name in ['shared_24h','persistent_24h','whole_12h','whole_6h','whole_3h']:
  z[f'{name}_abs_rank']=pd.Series(np.abs(z[name])).rank(method='first',ascending=False).astype(int)
 z['primary_prediction']=z[primary] if primary else np.nan;z['primary_abs_rank']=z[f'{primary}_abs_rank'] if primary else np.nan
 keep=['probe']+[f'resp_{s}' for s in states]+['shared_24h','persistent_24h','whole_12h','whole_6h','whole_3h']+[f'{x}_abs_rank' for x in ['shared_24h','persistent_24h','whole_12h','whole_6h','whole_3h']]+['primary_prediction','primary_abs_rank']
 pp=OUT/'FROZEN_PREDICTIONS.tsv.gz';z[keep].to_csv(pp,sep='\t',index=False,compression='gzip')
 info={'status':'FROZEN_BEFORE_48H_TARGET','route':route,'primary_predictor':primary,'source_geometry':geom,
 'thresholds':{'shared':{'rank1_energy_min':.50,'median_pairwise_rho_min':.35,'effective_state_fraction_min':.50},'state_specific':{'rank1_energy_max':.40,'median_pairwise_rho_max':.15,'effective_state_fraction_max':.30}},
 'n_frozen_probes':len(z),'target_48h_samples':['GSM265471','GSM265472'],'target_72h_samples':['GSM265473','GSM265474'],
 'target_values_read':False,'target_files_present':False,'candidate_unit':'ATH1 probe set'}
 (OUT/'ROUTE_DECISION.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf-8')
 man={**info,'rules_sha256':sha(RULES),'router_rules_sha256':sha(ROUTER_RULES),'fe_source_sha256':sha(FE_MATRIX),'standard_source_sha256':sha(STD_MATRIX),'pretarget_sample_sha256':hashes,'prediction_sha256':sha(pp),'freeze_script_sha256':sha(Path(__file__))}
 (OUT/'FREEZE_MANIFEST.json').write_text(json.dumps(man,ensure_ascii=False,indent=2),encoding='utf-8')
 print('DINNENY2008_IRON_ROUTER_FREEZE_SUCCESS');print(json.dumps(info,ensure_ascii=False,indent=2));print('48H/72H TARGET FILES NOT PRESENT OR READ')
if __name__=='__main__':main()
