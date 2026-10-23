from __future__ import annotations
from pathlib import Path
import hashlib, io, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_task'/'low_pH'
CELL=BASE/'GSE30095_series_matrix.txt.gz'
WHOLE=BASE/'whole_root_24h_source'
OUT=BASE/'task_48h_blind'; OUT.mkdir(parents=True,exist_ok=True)
RULE=ROOT/'任务结构优先证据选择规则_v0.1_冻结_2026-10-21.md'
PROTOCOL=ROOT/'低pH_48h_任务结构独立盲测规则_2026-10-21.md'
TARGET_DIR=BASE/'target_48h_samples'

CELL_GROUPS={
 'columella':{'low':['GSM744726','GSM744727','GSM744728'],'std':['GSM744729','GSM744730','GSM744731']},
 'epidermis':{'low':['GSM744732','GSM744733','GSM744734'],'std':['GSM744735','GSM744736','GSM744737']},
 'cortex':{'low':['GSM744738','GSM744739','GSM744740'],'std':['GSM744741','GSM744742','GSM744743']},
 'endodermis':{'low':['GSM744744','GSM744745','GSM744746'],'std':['GSM744747','GSM744748','GSM744749']},
 'stele':{'low':['GSM744750','GSM744751','GSM744752'],'std':['GSM744753','GSM744754','GSM744755']},
}
WHOLE_STD=['GSM744782','GSM744783']; WHOLE_LOW=['GSM744796','GSM744797']
TARGETS=['GSM744784','GSM744785','GSM744798','GSM744799']

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
 gsm=p.stem.split('.')[0]; d=d.rename(columns={'ID_REF':'probe','VALUE':gsm}); d['probe']=d.probe.astype(str); d[gsm]=pd.to_numeric(d[gsm],errors='raise'); return d

if TARGET_DIR.exists() and any(TARGET_DIR.iterdir()): raise RuntimeError('48h target present before freeze')
for gsm in TARGETS:
 if (BASE/f'{gsm}.soft.txt').exists(): raise RuntimeError(f'target present before freeze {gsm}')

# The entire GSE30095 matrix is allowed source data: all samples are 24h cell-layer source profiles.
cell=pd.read_csv(CELL,sep='\t',compression='gzip',comment='!',quotechar='"')
cell=cell.rename(columns={'ID_REF':'probe'}); cell['probe']=cell.probe.astype(str)
needed=[x for g in CELL_GROUPS.values() for x in (g['low']+g['std'])]
miss=[c for c in needed if c not in cell.columns]
if miss: raise RuntimeError(f'missing cell source columns {miss[:5]}')
for c in needed:
 cell[c]=pd.to_numeric(cell[c],errors='raise')
 v=cell[c].to_numpy(float)
 if np.any(~np.isfinite(v)) or np.any(v<=0): raise RuntimeError(f'invalid cell source {c}')
 cell[c]=np.log2(v)

R=[]; ctnames=list(CELL_GROUPS)
for ct in ctnames:
 g=CELL_GROUPS[ct]
 R.append(cell[g['low']].mean(axis=1).to_numpy(float)-cell[g['std']].mean(axis=1).to_numpy(float))
M=np.column_stack(R)
shared=M.mean(axis=1)

# Descriptive source geometry only; future-time task rule prevents it overriding same-state history.
A=M.T; u,s,vt=np.linalg.svd(A,full_matrices=False); e=s*s/(s*s).sum(); p=u[:,0]**2;p=p/p.sum();pairs=[];loo=[]
for i in range(A.shape[0]):
 for j in range(i+1,A.shape[0]): pairs.append(float(spearmanr(A[i],A[j]).statistic))
 oth=np.mean(np.delete(A,i,axis=0),axis=0); loo.append(float(spearmanr(A[i],oth).statistic))
eff=float(1/np.sum(p*p))
geom={'rank1_energy':float(e[0]),'rank1_3_energy':float(e[:3].sum()),'median_pairwise_rho':float(np.median(pairs)),'median_loo_shared_rho':float(np.median(loo)),'max_rank1_loading_share':float(p.max()),'effective_states':eff,'effective_state_fraction':eff/len(ctnames)}

z=None; whole_hash={}
for gsm in WHOLE_STD+WHOLE_LOW:
 pth=WHOLE/f'{gsm}.soft.txt'
 if not pth.exists(): raise FileNotFoundError(pth)
 whole_hash[gsm]=sha(pth);q=read_soft(pth);z=q if z is None else z.merge(q,on='probe',validate='one_to_one')
for c in WHOLE_STD+WHOLE_LOW:
 v=z[c].to_numpy(float)
 if np.any(~np.isfinite(v)) or np.any(v<=0): raise RuntimeError(f'invalid whole source {c}')
 z[c]=np.log2(v)
persistent=z[WHOLE_LOW].mean(axis=1).to_numpy(float)-z[WHOLE_STD].mean(axis=1).to_numpy(float)
whole=pd.DataFrame({'probe':z.probe.astype(str),'persistent_whole_24h':persistent})
shared_df=pd.DataFrame({'probe':cell.probe.astype(str),'shared_celltypes_24h':shared})
out=whole.merge(shared_df,on='probe',how='inner',validate='one_to_one')
if len(out)<20000: raise RuntimeError(f'unexpected probe universe {len(out)}')
for c in ['persistent_whole_24h','shared_celltypes_24h']:
 out[f'{c}_abs_rank']=pd.Series(np.abs(out[c].to_numpy(float))).rank(method='first',ascending=False).to_numpy(int)
out['primary_prediction']=out['persistent_whole_24h'];out['primary_abs_rank']=out['persistent_whole_24h_abs_rank']
pred=OUT/'FROZEN_48H_PREDICTIONS.tsv.gz';out.to_csv(pred,sep='\t',index=False,compression='gzip')
manifest={'status':'FROZEN_BEFORE_48H_TARGET','rule_version':'task_structure_v0.1','task_type':'future_time','decision':'最近同状态历史','primary_predictor':'persistent_whole_24h','alternative_predictor':'shared_celltypes_24h','celltype_shared_geometry_descriptive_only':geom,'n_celltypes':len(ctnames),'celltypes':ctnames,'n_probes':int(len(out)),'target_samples':TARGETS,'target_files_present':False,'target_values_read':False,'rule_sha256':sha(RULE),'protocol_sha256':sha(PROTOCOL),'cell_source_matrix_sha256':sha(CELL),'whole_24h_source_sha256':whole_hash,'prediction_sha256':sha(pred),'freeze_script_sha256':sha(Path(__file__)),'candidate_unit':'ATH1 probe set'}
(OUT/'FREEZE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('LOWPH_48H_TASK_FREEZE_SUCCESS');print(json.dumps(manifest,ensure_ascii=False,indent=2));print('48H TARGET FILES NOT PRESENT OR READ')
