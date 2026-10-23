from pathlib import Path
import hashlib,io,json
import numpy as np,pandas as pd
from scipy.stats import spearmanr,hypergeom
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_task'/'sulfur';SRC=BASE/'source_0_3h';TGT=BASE/'target_12h_samples';OUT=BASE/'task_12h_blind';PRED=OUT/'FROZEN_12H_PREDICTION.tsv.gz';MAN=OUT/'FREEZE_MANIFEST.json'
ZERO=['GSM744800','GSM744801'];TARGET=['GSM744804','GSM744805'];KS=[100,250,500,1000]
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
 d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE']);gsm=p.stem.split('.')[0];d=d.rename(columns={'ID_REF':'probe','VALUE':gsm});d['probe']=d.probe.astype(str);d[gsm]=pd.to_numeric(d[gsm],errors='raise');return d
man=json.loads(MAN.read_text(encoding='utf-8'))
if man['status']!='FROZEN_BEFORE_12H_TARGET' or man['rule_version']!='task_structure_v0.1' or man['primary_predictor']!='same_state_3h':raise RuntimeError('freeze changed')
if sha(PRED)!=man['prediction_sha256']:raise RuntimeError('prediction hash changed')
pred=pd.read_csv(PRED,sep='\t');z=None
# Baseline 0h is source and was frozen; verify unchanged.
for g in ZERO:
 p=SRC/f'{g}.soft.txt'
 if sha(p)!=man['source_sha256'][g]:raise RuntimeError(f'baseline changed {g}')
 q=read_soft(p);z=q if z is None else z.merge(q,on='probe',validate='one_to_one')
for g in TARGET:
 p=TGT/f'{g}.soft.txt'
 if not p.exists():raise FileNotFoundError(p)
 z=z.merge(read_soft(p),on='probe',validate='one_to_one')
z=pred.merge(z,on='probe',how='inner',validate='one_to_one')
if len(z)!=len(pred):raise RuntimeError(f'probe loss {len(pred)}->{len(z)}')
for c in ZERO+TARGET:
 v=z[c].to_numpy(float)
 if np.any(~np.isfinite(v)) or np.any(v<=0):raise RuntimeError(f'invalid {c}')
 z[c]=np.log2(v)
base=z[ZERO].mean(axis=1).to_numpy(float);reps={g:z[g].to_numpy(float)-base for g in TARGET};truth=np.mean(np.column_stack(list(reps.values())),axis=1);pv=z.same_state_3h.to_numpy(float);N=len(z)
rho=float(spearmanr(pv,truth).statistic);mask=(pv!=0)&(truth!=0);sign=float(np.mean(np.sign(pv[mask])==np.sign(truth[mask])));pr=pd.Series(np.abs(pv)).rank(method='first',ascending=False).to_numpy(int);trank=pd.Series(np.abs(truth)).rank(method='first',ascending=False).to_numpy(int)
rows=[]
for k in KS:
 kk=min(k,N);ov=sorted(set(np.where(pr<=kk)[0])&set(np.where(trank<=kk)[0]));n=len(ov);exp=kk*kk/N
 rows.append({'k':kk,'overlap':n,'random_expected_overlap':exp,'enrichment':n/exp,'hypergeom_p':float(hypergeom.sf(n-1,N,kk,kk)),'sign_agreement_overlap':float(np.mean(np.sign(pv[ov])==np.sign(truth[ov]))) if ov else np.nan})
tdf=pd.DataFrame(rows);tdf.to_csv(OUT/'blind_12h_topk_results.tsv',sep='\t',index=False)
va,vb=reps[TARGET[0]],reps[TARGET[1]];ra=pd.Series(np.abs(va)).rank(method='first',ascending=False).to_numpy(int);rb=pd.Series(np.abs(vb)).rank(method='first',ascending=False).to_numpy(int);rep_rho=float(spearmanr(va,vb).statistic);rep100=len(set(np.where(ra<=100)[0])&set(np.where(rb<=100)[0]))
pd.DataFrame([{'model':'same_state_3h','n_probes':N,'spearman_full_vector':rho,'sign_agreement_all_nonzero':sign}]).to_csv(OUT/'blind_12h_vector_results.tsv',sep='\t',index=False)
pd.DataFrame([{'rep_a':TARGET[0],'rep_b':TARGET[1],'spearman_full_vector':rep_rho,'top100_overlap':rep100}]).to_csv(OUT/'target_12h_replicate_reference.tsv',sep='\t',index=False)
top=tdf.set_index('k');summary={'freeze_task_type':man['task_type'],'freeze_primary_predictor':man['primary_predictor'],'rule_version':man['rule_version'],'n_probes':N,'target_replicate_rho':rep_rho,'target_replicate_top100_overlap':rep100,'primary_full_rho':rho,'primary_sign_agreement_all_nonzero':sign,'primary_top100_overlap':int(top.loc[100,'overlap']),'primary_top100_enrichment':float(top.loc[100,'enrichment']),'primary_top100_direction':float(top.loc[100,'sign_agreement_overlap']) if pd.notna(top.loc[100,'sign_agreement_overlap']) else None}
(OUT/'BLIND_12H_RESULT.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print('SULFUR_12H_TASK_SCORE_SUCCESS');print(json.dumps(summary,ensure_ascii=False,indent=2));print('\nTOPK');print(tdf.to_string(index=False))
