from pathlib import Path
import hashlib,io,json
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_task'/'sulfur';SRC=BASE/'source_0_3h';OUT=BASE/'task_12h_blind';OUT.mkdir(parents=True,exist_ok=True)
RULE=ROOT/'任务结构优先证据选择规则_v0.1_冻结_2026-10-21.md';PROT=ROOT/'缺硫_3h到12h_任务结构独立盲测规则_2026-10-21.md'
ZERO=['GSM744800','GSM744801'];H3=['GSM744802','GSM744803'];TARGET=['GSM744804','GSM744805'];TDIR=BASE/'target_12h_samples'
def sha(p):
 h=hashlib.sha256();
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
if TDIR.exists() and any(TDIR.iterdir()):raise RuntimeError('12h target present before freeze')
for g in TARGET:
 if (BASE/f'{g}.soft.txt').exists():raise RuntimeError(f'target present {g}')
z=None;hs={}
for g in ZERO+H3:
 p=SRC/f'{g}.soft.txt';hs[g]=sha(p);q=read_soft(p);z=q if z is None else z.merge(q,on='probe',validate='one_to_one')
for c in ZERO+H3:
 v=z[c].to_numpy(float)
 if np.any(~np.isfinite(v)) or np.any(v<=0):raise RuntimeError(f'invalid {c}')
 z[c]=np.log2(v)
pred=z[H3].mean(axis=1).to_numpy(float)-z[ZERO].mean(axis=1).to_numpy(float)
out=pd.DataFrame({'probe':z.probe.astype(str),'same_state_3h':pred});out['same_state_3h_abs_rank']=pd.Series(np.abs(pred)).rank(method='first',ascending=False).to_numpy(int);out['primary_prediction']=pred;out['primary_abs_rank']=out.same_state_3h_abs_rank
pp=OUT/'FROZEN_12H_PREDICTION.tsv.gz';out.to_csv(pp,sep='\t',index=False,compression='gzip')
man={'status':'FROZEN_BEFORE_12H_TARGET','rule_version':'task_structure_v0.1','task_type':'future_time','decision':'最近同状态历史','primary_predictor':'same_state_3h','n_probes':int(len(out)),'target_samples':TARGET,'target_files_present':False,'target_values_read':False,'rule_sha256':sha(RULE),'protocol_sha256':sha(PROT),'source_sha256':hs,'prediction_sha256':sha(pp),'freeze_script_sha256':sha(Path(__file__)),'candidate_unit':'ATH1 probe set'}
(OUT/'FREEZE_MANIFEST.json').write_text(json.dumps(man,ensure_ascii=False,indent=2),encoding='utf-8');print('SULFUR_12H_TASK_FREEZE_SUCCESS');print(json.dumps(man,ensure_ascii=False,indent=2));print('12H TARGET FILES NOT PRESENT OR READ')
