from __future__ import annotations
from pathlib import Path
import hashlib, io, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, hypergeom

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_router'/'GSE46205_salt'; SRC=BASE/'source_1h_samples'; TGT=BASE/'target_8h_samples'; OUT=BASE/'task_8h_blind'
PRED=OUT/'FROZEN_8H_PREDICTIONS.tsv.gz'; MAN=OUT/'FREEZE_MANIFEST.json'; KS=[100,250,500,1000]
GROUPS={
'epidermis':{'control':['GSM1126230','GSM1126231','GSM1126232'],'target':['GSM1126242','GSM1126243','GSM1126244']},
'cortex':{'control':['GSM1126254','GSM1126255','GSM1126256'],'target':['GSM1126266','GSM1126267','GSM1126268']},
'endodermis':{'control':['GSM1126278','GSM1126279','GSM1126280'],'target':['GSM1126290','GSM1126291','GSM1126292']},
'stele':{'control':['GSM1126302','GSM1126303','GSM1126304'],'target':['GSM1126314','GSM1126315','GSM1126316']},
}
MODELS=['same_state_3h','shared_other3_3h']

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

man=json.loads(MAN.read_text(encoding='utf-8'))
if man['status']!='FROZEN_BEFORE_8H_TARGET' or man['rule_version']!='task_structure_v0.1' or man['task_type']!='future_time' or man['primary_predictor']!='same_state_3h': raise RuntimeError('freeze decision changed')
if sha(PRED)!=man['prediction_sha256']: raise RuntimeError('prediction hash changed')
pred=pd.read_csv(PRED,sep='\t')
vector=[]; topk=[]; reps=[]
for ct,g in GROUPS.items():
 qpred=pred[pred.target_celltype==ct].copy(); z=None
 for gsm in g['control']:
  p=SRC/f'{gsm}.soft.txt'
  if sha(p)!=man['source_control_sha256'][gsm]: raise RuntimeError(f'control changed {gsm}')
  q=read_soft(p); z=q if z is None else z.merge(q,on='probe',validate='one_to_one')
 for gsm in g['target']:
  p=TGT/f'{gsm}.soft.txt'
  if not p.exists(): raise FileNotFoundError(p)
  z=z.merge(read_soft(p),on='probe',validate='one_to_one')
 z=qpred.merge(z,on='probe',how='inner',validate='one_to_one')
 if len(z)!=len(qpred): raise RuntimeError(f'probe loss {ct}: {len(qpred)}->{len(z)}')
 for c in g['control']+g['target']:
  v=z[c].to_numpy(float)
  if np.any(~np.isfinite(v)) or np.any(v<=0): raise RuntimeError(f'invalid {c}')
  z[c]=np.log2(v)
 base=z[g['control']].mean(axis=1).to_numpy(float)
 r={gsm:z[gsm].to_numpy(float)-base for gsm in g['target']}
 truth=np.mean(np.column_stack(list(r.values())),axis=1); N=len(z)
 truth_rank=pd.Series(np.abs(truth)).rank(method='first',ascending=False).to_numpy(int)
 ids=g['target']
 for i in range(3):
  for j in range(i+1,3):
   a,b=ids[i],ids[j]; va,vb=r[a],r[b]
   ra=pd.Series(np.abs(va)).rank(method='first',ascending=False).to_numpy(int); rb=pd.Series(np.abs(vb)).rank(method='first',ascending=False).to_numpy(int)
   reps.append({'celltype':ct,'rep_a':a,'rep_b':b,'spearman_full_vector':float(spearmanr(va,vb).statistic),'top100_overlap':len(set(np.where(ra<=100)[0])&set(np.where(rb<=100)[0]))})
 for model in MODELS:
  pv=z[model].to_numpy(float); mask=(pv!=0)&(truth!=0)&np.isfinite(pv)&np.isfinite(truth)
  vector.append({'celltype':ct,'model':model,'n_probes':N,'spearman_full_vector':float(spearmanr(pv,truth).statistic),'sign_agreement_all_nonzero':float(np.mean(np.sign(pv[mask])==np.sign(truth[mask])))})
  pr=pd.Series(np.abs(pv)).rank(method='first',ascending=False).to_numpy(int)
  for k in KS:
   kk=min(k,N); ov=sorted(set(np.where(pr<=kk)[0])&set(np.where(truth_rank<=kk)[0])); n=len(ov); exp=kk*kk/N
   topk.append({'celltype':ct,'model':model,'k':kk,'overlap':n,'random_expected_overlap':exp,'enrichment':n/exp,'hypergeom_p':float(hypergeom.sf(n-1,N,kk,kk)),'sign_agreement_overlap':float(np.mean(np.sign(pv[ov])==np.sign(truth[ov]))) if ov else np.nan})
vdf=pd.DataFrame(vector); tdf=pd.DataFrame(topk); rdf=pd.DataFrame(reps)
vdf.to_csv(OUT/'blind_8h_vector_results.tsv',sep='\t',index=False); tdf.to_csv(OUT/'blind_8h_topk_results.tsv',sep='\t',index=False); rdf.to_csv(OUT/'target_8h_replicate_reference.tsv',sep='\t',index=False)
p=vdf[vdf.model=='same_state_3h'].set_index('celltype'); a=vdf[vdf.model=='shared_other3_3h'].set_index('celltype')
tp=tdf[(tdf.model=='same_state_3h')&(tdf.k==100)].set_index('celltype'); ta=tdf[(tdf.model=='shared_other3_3h')&(tdf.k==100)].set_index('celltype')
summary={'freeze_task_type':man['task_type'],'freeze_primary_predictor':man['primary_predictor'],'rule_version':man['rule_version'],'n_celltypes':4,'n_probes_per_celltype':int(man['n_probes_per_celltype']),'shared_geometry_3h_descriptive_only':man['shared_geometry_3h_descriptive_only'],'target_replicate_rho_median':float(rdf.spearman_full_vector.median()),'target_replicate_top100_overlap_median':float(rdf.top100_overlap.median()),'primary_full_rho_by_celltype':{ct:float(p.loc[ct,'spearman_full_vector']) for ct in GROUPS},'alternative_full_rho_by_celltype':{ct:float(a.loc[ct,'spearman_full_vector']) for ct in GROUPS},'primary_full_rho_median':float(p.spearman_full_vector.median()),'alternative_full_rho_median':float(a.spearman_full_vector.median()),'primary_full_vector_wins':int(sum(p.loc[ct,'spearman_full_vector']>a.loc[ct,'spearman_full_vector'] for ct in GROUPS)),'primary_top100_overlap_by_celltype':{ct:int(tp.loc[ct,'overlap']) for ct in GROUPS},'alternative_top100_overlap_by_celltype':{ct:int(ta.loc[ct,'overlap']) for ct in GROUPS},'primary_top100_wins':int(sum(tp.loc[ct,'overlap']>ta.loc[ct,'overlap'] for ct in GROUPS)),'primary_top100_ties':int(sum(tp.loc[ct,'overlap']==ta.loc[ct,'overlap'] for ct in GROUPS)),'primary_top100_enrichment_median':float(tp.enrichment.median()),'primary_top100_direction_median':float(tp.sign_agreement_overlap.median())}
(OUT/'BLIND_8H_RESULT.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print('GSE46205_8H_TASK_SCORE_SUCCESS'); print(json.dumps(summary,ensure_ascii=False,indent=2)); print('\nVECTOR'); print(vdf.to_string(index=False)); print('\nTOPK'); print(tdf.to_string(index=False)); print('\nREPLICATES'); print(rdf.to_string(index=False))
