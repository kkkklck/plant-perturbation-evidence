from __future__ import annotations
from pathlib import Path
import hashlib,io,json
import numpy as np,pandas as pd
from scipy.stats import spearmanr,hypergeom
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_router'/'Dinneny2008_iron';PRE=BASE/'pretarget_samples';TGT=BASE/'target_48h_samples';OUT=BASE/'router_blind'
PRED=OUT/'FROZEN_PREDICTIONS.tsv.gz';MAN=OUT/'FREEZE_MANIFEST.json';TARGETS=['GSM265471','GSM265472'];BASELINE=['GSM265461','GSM265462']
MODELS=['shared_24h','persistent_24h','whole_12h','whole_6h','whole_3h'];KS=[100,250,500,1000]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(1<<20),b''):h.update(x)
 return h.hexdigest()
def read_soft(p):
 lines=[];inside=False
 with p.open('r',encoding='utf-8',errors='strict') as f:
  for line in f:
   if line.startswith('!sample_table_begin'):inside=True;continue
   if line.startswith('!sample_table_end'):break
   if inside:lines.append(line)
 d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE']);gsm=p.stem.split('.')[0]
 d=d.rename(columns={'ID_REF':'probe','VALUE':gsm});d['probe']=d.probe.astype(str);d[gsm]=pd.to_numeric(d[gsm],errors='raise');return d
man=json.loads(MAN.read_text(encoding='utf-8'))
if man['status']!='FROZEN_BEFORE_48H_TARGET' or man['route']!='路线歧义' or man['primary_predictor'] is not None:raise RuntimeError('freeze decision changed')
if sha(PRED)!=man['prediction_sha256']:raise RuntimeError('prediction hash changed')
pred=pd.read_csv(PRED,sep='\t')
z=None
for gsm in BASELINE:
 p=PRE/f'{gsm}.soft.txt'
 if sha(p)!=man['pretarget_sample_sha256'][gsm]:raise RuntimeError(f'pretarget changed {gsm}')
 q=read_soft(p);z=q if z is None else z.merge(q,on='probe',validate='one_to_one')
for gsm in TARGETS:
 p=TGT/f'{gsm}.soft.txt'
 if not p.exists():raise FileNotFoundError(p)
 z=z.merge(read_soft(p),on='probe',validate='one_to_one')
z=pred.merge(z,on='probe',how='inner',validate='one_to_one')
if len(z)!=len(pred):raise RuntimeError(f'probe loss {len(pred)}->{len(z)}')
for c in BASELINE+TARGETS:
 v=z[c].to_numpy(float)
 if np.any(~np.isfinite(v)) or np.any(v<=0):raise RuntimeError(f'invalid {c}')
 z[c]=np.log2(v)
base=z[BASELINE].mean(axis=1).to_numpy(float);rep={x:z[x].to_numpy(float)-base for x in TARGETS};target=np.mean(np.column_stack(list(rep.values())),axis=1);N=len(z)
truth_rank=pd.Series(np.abs(target)).rank(method='first',ascending=False).to_numpy(int)
vr=[];tr=[]
for model in MODELS:
 pv=z[model].to_numpy(float);mask=(pv!=0)&(target!=0)&np.isfinite(pv)&np.isfinite(target)
 vr.append({'model':model,'n_probes':N,'spearman_full_vector':float(spearmanr(pv,target).statistic),'sign_agreement_all_nonzero':float(np.mean(np.sign(pv[mask])==np.sign(target[mask])))})
 pr=pd.Series(np.abs(pv)).rank(method='first',ascending=False).to_numpy(int)
 for k in KS:
  kk=min(k,N);ov=sorted(set(np.where(pr<=kk)[0])&set(np.where(truth_rank<=kk)[0]));n=len(ov);exp=kk*kk/N
  tr.append({'model':model,'k':kk,'overlap':n,'random_expected_overlap':exp,'enrichment':n/exp,'hypergeom_p':float(hypergeom.sf(n-1,N,kk,kk)),'sign_agreement_overlap':float(np.mean(np.sign(pv[ov])==np.sign(target[ov]))) if ov else np.nan})
vdf=pd.DataFrame(vr);tdf=pd.DataFrame(tr);vdf.to_csv(OUT/'blind_48h_vector_results.tsv',sep='\t',index=False);tdf.to_csv(OUT/'blind_48h_topk_results.tsv',sep='\t',index=False)
ra=pd.Series(np.abs(rep[TARGETS[0]])).rank(method='first',ascending=False).to_numpy(int);rb=pd.Series(np.abs(rep[TARGETS[1]])).rank(method='first',ascending=False).to_numpy(int)
rep_rho=float(spearmanr(rep[TARGETS[0]],rep[TARGETS[1]]).statistic);rep100=len(set(np.where(ra<=100)[0])&set(np.where(rb<=100)[0]))
summary={'freeze_route':man['route'],'freeze_primary_predictor':None,'n_probes':N,'target_replicate_rho':rep_rho,'target_replicate_top100_overlap':rep100,'model_full_rho':dict(zip(vdf.model,vdf.spearman_full_vector)),'model_top100_overlap':dict(zip(tdf[tdf.k==100].model,tdf[tdf.k==100].overlap.astype(int))),'posthoc_best_full_vector':str(vdf.sort_values('spearman_full_vector',ascending=False).model.iloc[0]),'interpretation':'Second prospective ambiguous-region calibration; no pre-reveal winner selected.'}
(OUT/'BLIND_48H_RESULT.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print('DINNENY2008_IRON_48H_SCORE_SUCCESS');print(json.dumps(summary,ensure_ascii=False,indent=2));print('\nVECTOR');print(vdf.to_string(index=False));print('\nTOPK');print(tdf.to_string(index=False))
