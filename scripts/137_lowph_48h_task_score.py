from __future__ import annotations
from pathlib import Path
import hashlib,io,json
import numpy as np,pandas as pd
from scipy.stats import spearmanr,hypergeom

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_task'/'low_pH'; SRC=BASE/'whole_root_24h_source'; TGT=BASE/'target_48h_samples'; OUT=BASE/'task_48h_blind'
PRED=OUT/'FROZEN_48H_PREDICTIONS.tsv.gz'; MAN=OUT/'FREEZE_MANIFEST.json'
STD24=['GSM744782','GSM744783']; STD48=['GSM744784','GSM744785']; LOW48=['GSM744798','GSM744799']
MODELS=['persistent_whole_24h','shared_celltypes_24h']; KS=[100,250,500,1000]

def sha(p:Path):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(1<<20),b''):h.update(x)
 return h.hexdigest()

def read_soft(p:Path):
 lines=[];inside=False
 with p.open('r',encoding='utf-8',errors='strict') as f:
  for line in f:
   if line.startswith('!sample_table_begin'):inside=True;continue
   if line.startswith('!sample_table_end'):break
   if inside:lines.append(line)
 d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE'])
 gsm=p.stem.split('.')[0];d=d.rename(columns={'ID_REF':'probe','VALUE':gsm});d['probe']=d.probe.astype(str);d[gsm]=pd.to_numeric(d[gsm],errors='raise');return d

man=json.loads(MAN.read_text(encoding='utf-8'))
if man['status']!='FROZEN_BEFORE_48H_TARGET' or man['rule_version']!='task_structure_v0.1' or man['task_type']!='future_time' or man['primary_predictor']!='persistent_whole_24h':raise RuntimeError('freeze decision changed')
if sha(PRED)!=man['prediction_sha256']:raise RuntimeError('prediction hash changed')
pred=pd.read_csv(PRED,sep='\t')
z=None
# 48h response must use time-matched 48h standard controls, all withheld together with low-pH targets.
for gsm in STD48+LOW48:
 p=TGT/f'{gsm}.soft.txt'
 if not p.exists():raise FileNotFoundError(p)
 q=read_soft(p);z=q if z is None else z.merge(q,on='probe',validate='one_to_one')
z=pred.merge(z,on='probe',how='inner',validate='one_to_one')
if len(z)!=len(pred):raise RuntimeError(f'probe loss {len(pred)}->{len(z)}')
for c in STD48+LOW48:
 v=z[c].to_numpy(float)
 if np.any(~np.isfinite(v)) or np.any(v<=0):raise RuntimeError(f'invalid {c}')
 z[c]=np.log2(v)
base=z[STD48].mean(axis=1).to_numpy(float)
reps={g:z[g].to_numpy(float)-base for g in LOW48}
target=np.mean(np.column_stack(list(reps.values())),axis=1);N=len(z);truth_rank=pd.Series(np.abs(target)).rank(method='first',ascending=False).to_numpy(int)
vr=[];tr=[]
for model in MODELS:
 pv=z[model].to_numpy(float);mask=(pv!=0)&(target!=0)&np.isfinite(pv)&np.isfinite(target)
 vr.append({'model':model,'n_probes':N,'spearman_full_vector':float(spearmanr(pv,target).statistic),'sign_agreement_all_nonzero':float(np.mean(np.sign(pv[mask])==np.sign(target[mask])))})
 pr=pd.Series(np.abs(pv)).rank(method='first',ascending=False).to_numpy(int)
 for k in KS:
  kk=min(k,N);ov=sorted(set(np.where(pr<=kk)[0])&set(np.where(truth_rank<=kk)[0]));n=len(ov);exp=kk*kk/N
  tr.append({'model':model,'k':kk,'overlap':n,'random_expected_overlap':exp,'enrichment':n/exp,'hypergeom_p':float(hypergeom.sf(n-1,N,kk,kk)),'sign_agreement_overlap':float(np.mean(np.sign(pv[ov])==np.sign(target[ov]))) if ov else np.nan})
vdf=pd.DataFrame(vr);tdf=pd.DataFrame(tr);vdf.to_csv(OUT/'blind_48h_vector_results.tsv',sep='\t',index=False);tdf.to_csv(OUT/'blind_48h_topk_results.tsv',sep='\t',index=False)
# replicate reference: each low-pH replicate relative to the common mean 48h standard control.
a,b=LOW48;va,vb=reps[a],reps[b];ra=pd.Series(np.abs(va)).rank(method='first',ascending=False).to_numpy(int);rb=pd.Series(np.abs(vb)).rank(method='first',ascending=False).to_numpy(int)
rdf=pd.DataFrame([{'rep_a':a,'rep_b':b,'spearman_full_vector':float(spearmanr(va,vb).statistic),'top100_overlap':len(set(np.where(ra<=100)[0])&set(np.where(rb<=100)[0])),'top500_overlap':len(set(np.where(ra<=500)[0])&set(np.where(rb<=500)[0]))}]);rdf.to_csv(OUT/'target_48h_replicate_reference.tsv',sep='\t',index=False)
p=vdf.set_index('model');tk=tdf[tdf.k==100].set_index('model')
summary={'freeze_task_type':man['task_type'],'freeze_primary_predictor':man['primary_predictor'],'rule_version':man['rule_version'],'n_probes':N,'source_shared_geometry_descriptive_only':man['celltype_shared_geometry_descriptive_only'],'target_replicate_rho':float(rdf.spearman_full_vector.iloc[0]),'target_replicate_top100_overlap':int(rdf.top100_overlap.iloc[0]),'primary_full_rho':float(p.loc['persistent_whole_24h','spearman_full_vector']),'alternative_full_rho':float(p.loc['shared_celltypes_24h','spearman_full_vector']),'primary_top100_overlap':int(tk.loc['persistent_whole_24h','overlap']),'alternative_top100_overlap':int(tk.loc['shared_celltypes_24h','overlap']),'primary_top100_enrichment':float(tk.loc['persistent_whole_24h','enrichment']),'alternative_top100_enrichment':float(tk.loc['shared_celltypes_24h','enrichment']),'primary_top100_direction':float(tk.loc['persistent_whole_24h','sign_agreement_overlap']) if pd.notna(tk.loc['persistent_whole_24h','sign_agreement_overlap']) else None,'primary_wins_full_vector':bool(p.loc['persistent_whole_24h','spearman_full_vector']>p.loc['shared_celltypes_24h','spearman_full_vector']),'primary_wins_top100':bool(tk.loc['persistent_whole_24h','overlap']>tk.loc['shared_celltypes_24h','overlap'])}
(OUT/'BLIND_48H_RESULT.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print('LOWPH_48H_TASK_SCORE_SUCCESS');print(json.dumps(summary,ensure_ascii=False,indent=2));print('\nVECTOR');print(vdf.to_string(index=False));print('\nTOPK');print(tdf.to_string(index=False));print('\nREPLICATE');print(rdf.to_string(index=False))
