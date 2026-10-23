from __future__ import annotations
from pathlib import Path
import hashlib, io, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_router'/'GSE46205_salt'
SRC=BASE/'source_1h_samples'; H3=BASE/'target_3h_samples'; OUT=BASE/'task_8h_blind'
OUT.mkdir(parents=True,exist_ok=True)
RULE=ROOT/'任务结构优先证据选择规则_v0.1_冻结_2026-10-21.md'
TARGET_DIR=BASE/'target_8h_samples'
GROUPS={
'epidermis':{'control':['GSM1126230','GSM1126231','GSM1126232'],'h3':['GSM1126239','GSM1126240','GSM1126241'],'h8':['GSM1126242','GSM1126243','GSM1126244']},
'cortex':{'control':['GSM1126254','GSM1126255','GSM1126256'],'h3':['GSM1126263','GSM1126264','GSM1126265'],'h8':['GSM1126266','GSM1126267','GSM1126268']},
'endodermis':{'control':['GSM1126278','GSM1126279','GSM1126280'],'h3':['GSM1126287','GSM1126288','GSM1126289'],'h8':['GSM1126290','GSM1126291','GSM1126292']},
'stele':{'control':['GSM1126302','GSM1126303','GSM1126304'],'h3':['GSM1126311','GSM1126312','GSM1126313'],'h8':['GSM1126314','GSM1126315','GSM1126316']},
}

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

# 8h target must not be present before freeze.
if TARGET_DIR.exists() and any(TARGET_DIR.iterdir()): raise RuntimeError('8h target already present')
for g in GROUPS.values():
 for gsm in g['h8']:
  if (BASE/f'{gsm}.soft.txt').exists(): raise RuntimeError(f'8h target present: {gsm}')

# Build 3h same-state response for all four states using previously revealed 3h values.
z=None; src_hash={}; h3_hash={}
for ct,g in GROUPS.items():
 for gsm in g['control']:
  p=SRC/f'{gsm}.soft.txt'; src_hash[gsm]=sha(p); q=read_soft(p); z=q if z is None else z.merge(q,on='probe',validate='one_to_one')
 for gsm in g['h3']:
  p=H3/f'{gsm}.soft.txt'; h3_hash[gsm]=sha(p); z=z.merge(read_soft(p),on='probe',validate='one_to_one')
for c in [x for g in GROUPS.values() for x in (g['control']+g['h3'])]:
 v=z[c].to_numpy(float)
 if np.any(~np.isfinite(v)) or np.any(v<=0): raise RuntimeError(f'invalid {c}')
 z[c]=np.log2(v)
cts=list(GROUPS); R=[]
for ct in cts:
 g=GROUPS[ct]; base=z[g['control']].mean(axis=1).to_numpy(float); h3=z[g['h3']].mean(axis=1).to_numpy(float); R.append(h3-base)
M=np.column_stack(R)  # probes x states

# Descriptive 3h shared geometry; task rule says this may not override same-state time history.
A=M.T; u,s,vt=np.linalg.svd(A,full_matrices=False); e=s*s/(s*s).sum(); p=u[:,0]**2; p=p/p.sum()
pairs=[];loo=[]
for i in range(A.shape[0]):
 for j in range(i+1,A.shape[0]): pairs.append(float(spearmanr(A[i],A[j]).statistic))
 oth=np.mean(np.delete(A,i,axis=0),axis=0);loo.append(float(spearmanr(A[i],oth).statistic))
eff=float(1/np.sum(p*p))
geom={'rank1_energy':float(e[0]),'rank1_3_energy':float(e[:3].sum()),'median_pairwise_rho':float(np.median(pairs)),'median_loo_shared_rho':float(np.median(loo)),'max_rank1_loading_share':float(p.max()),'effective_states':eff,'effective_state_fraction':eff/len(cts)}

rows=[]
for j,ct in enumerate(cts):
 same=M[:,j]; shared=np.mean(np.delete(M,j,axis=1),axis=1)
 q=pd.DataFrame({'probe':z.probe.astype(str),'target_celltype':ct,'same_state_3h':same,'shared_other3_3h':shared})
 q['same_state_3h_abs_rank']=pd.Series(np.abs(same)).rank(method='first',ascending=False).to_numpy(int)
 q['shared_other3_3h_abs_rank']=pd.Series(np.abs(shared)).rank(method='first',ascending=False).to_numpy(int)
 q['primary_prediction']=same; q['primary_abs_rank']=q['same_state_3h_abs_rank']
 rows.append(q)
out=pd.concat(rows,ignore_index=True)
pred=OUT/'FROZEN_8H_PREDICTIONS.tsv.gz'; out.to_csv(pred,sep='\t',index=False,compression='gzip')
manifest={'status':'FROZEN_BEFORE_8H_TARGET','rule_version':'task_structure_v0.1','task_type':'future_time','decision':'最近同状态历史','primary_predictor':'same_state_3h','alternative_predictor':'shared_other3_3h','shared_geometry_3h_descriptive_only':geom,'n_celltypes':4,'celltypes':cts,'n_probes_per_celltype':int(len(z)),'target_8h_samples':{ct:g['h8'] for ct,g in GROUPS.items()},'target_values_read':False,'target_files_present':False,'rule_sha256':sha(RULE),'source_control_sha256':src_hash,'history_3h_sha256':h3_hash,'prediction_sha256':sha(pred),'freeze_script_sha256':sha(Path(__file__)),'candidate_unit':'ATH1 probe set'}
(OUT/'FREEZE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('GSE46205_8H_TASK_FREEZE_SUCCESS');print(json.dumps(manifest,ensure_ascii=False,indent=2));print('8H TARGET FILES NOT PRESENT OR READ')
