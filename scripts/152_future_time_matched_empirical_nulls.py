from __future__ import annotations
from pathlib import Path
import io, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data'/'publication_controls_v1'/'future_time_matched_nulls'
OUT.mkdir(parents=True,exist_ok=True)
SEED=20261021
N_PERM=10000
K=100
rng=np.random.default_rng(SEED)


def ranks_abs(v):
    order=np.argsort(-np.abs(v),kind='mergesort'); r=np.empty(len(v),int); r[order]=np.arange(1,len(v)+1); return r


def quantile_bins(x, q=5):
    # robust quantile bins; duplicates='drop' may reduce bin count
    return pd.qcut(pd.Series(x).rank(method='average'),q=q,labels=False,duplicates='drop').to_numpy(int)


def matched_null(candidate_idx, truth_top_idx, feat1, feat2, n_perm=N_PERM):
    b1=quantile_bins(feat1,5); b2=quantile_bins(feat2,5)
    strata=list(zip(b1,b2))
    cand_counts={}
    for i in candidate_idx:
        cand_counts[strata[i]]=cand_counts.get(strata[i],0)+1
    pools={s:np.array([i for i,z in enumerate(strata) if z==s],dtype=int) for s in cand_counts}
    truth=set(map(int,truth_top_idx))
    null=np.empty(n_perm,int)
    for p in range(n_perm):
        sel=[]
        for s,n in cand_counts.items():
            pool=pools[s]
            replace=len(pool)<n
            sel.extend(rng.choice(pool,size=n,replace=replace).tolist())
        null[p]=len(set(sel)&truth)
    obs=len(set(map(int,candidate_idx))&truth)
    return obs,null,cand_counts


def summarize(dataset,context,n,obs,null):
    return {'dataset':dataset,'context':context,'n_genes':n,'k':K,'observed_overlap':int(obs),
            'null_mean':float(np.mean(null)),'null_sd':float(np.std(null,ddof=1)),
            'null_q95':float(np.quantile(null,.95)),'null_q99':float(np.quantile(null,.99)),
            'matched_enrichment':float(obs/np.mean(null)) if np.mean(null)>0 else np.inf,
            'empirical_p_one_sided':float((1+np.sum(null>=obs))/(len(null)+1))}

rows=[]; null_frames=[]

# ---------- GSE153103 salt 6h and 12h ----------
WB=ROOT/'data/external_salt/GSE153103/GSE153103_1_raw_reads_counts.xlsx'
PRE_ALL=['RAS048','RAS049','RAS050','RAS051','RAS056','RAS057','RAS058','RAS059']
h=pd.read_excel(WB,sheet_name='raw_counts',nrows=0); gene_col=h.columns[0]
d=pd.read_excel(WB,sheet_name='raw_counts',usecols=[gene_col]+PRE_ALL).rename(columns={gene_col:'gene'})
d['gene']=d['gene'].astype(str)
for c in PRE_ALL: d[c]=pd.to_numeric(d[c],errors='raise')
lib=d[PRE_ALL].sum(axis=0)
L=pd.DataFrame({'gene':d.gene})
for c in PRE_ALL:
    L[c]=np.log2(d[c].to_numpy(float)/float(lib[c])*1e6+1)
base_expr=((L['RAS048']+L['RAS056'])/2).to_numpy(float)
r1=((L['RAS049']-L['RAS048'])+(L['RAS057']-L['RAS056']))/2
r3=((L['RAS050']-L['RAS048'])+(L['RAS058']-L['RAS056']))/2
r6=((L['RAS051']-L['RAS048'])+(L['RAS059']-L['RAS056']))/2

for target_time,blind_dir,pred_col,target_col,hist_cols in [
    (6,'plantvc_v1_6h_blind','pred_baseline_3h','r6',[r1.to_numpy(float),r3.to_numpy(float)]),
    (12,'plantvc_v1_1_12h_blind','pred_v1_1','r12',[r1.to_numpy(float),r3.to_numpy(float),r6.to_numpy(float)])]:
    B=ROOT/'data/external_salt/GSE153103'/blind_dir
    fr=pd.read_csv(B/'FROZEN_PRETARGET_RANKING.tsv.gz',sep='\t')
    tg=pd.read_csv(B/f'target_{target_time}h_response.tsv.gz',sep='\t')
    z=fr.merge(tg[['gene',target_col]],on='gene',validate='one_to_one')
    if z.gene.astype(str).tolist()!=L.gene.astype(str).tolist(): raise RuntimeError('salt gene order mismatch')
    mask=z.technical_universe.to_numpy(bool)
    pred=z.loc[mask,pred_col].to_numpy(float); truth=z.loc[mask,target_col].to_numpy(float)
    f1=base_expr[mask]
    H=np.column_stack([x[mask] for x in hist_cols])
    f2=np.std(np.column_stack([np.zeros(mask.sum()),H]),axis=1)
    pr=ranks_abs(pred); tr=ranks_abs(truth); ci=np.where(pr<=K)[0]; ti=np.where(tr<=K)[0]
    obs,null,_=matched_null(ci,ti,f1,f2)
    rows.append(summarize(f'GSE153103_salt_{target_time}h','bulk',mask.sum(),obs,null))
    null_frames.append(pd.DataFrame({'dataset':f'GSE153103_salt_{target_time}h','context':'bulk','null_overlap':null}))

# ---------- GSE46205 8h, four cell types ----------
BASE=ROOT/'data/prospective_router/GSE46205_salt'; SRC=BASE/'source_1h_samples'; H3=BASE/'target_3h_samples'; H8=BASE/'target_8h_samples'
GROUPS={
'epidermis':{'control':['GSM1126230','GSM1126231','GSM1126232'],'h1':['GSM1126236','GSM1126237','GSM1126238'],'h3':['GSM1126239','GSM1126240','GSM1126241'],'h8':['GSM1126242','GSM1126243','GSM1126244']},
'cortex':{'control':['GSM1126254','GSM1126255','GSM1126256'],'h1':['GSM1126260','GSM1126261','GSM1126262'],'h3':['GSM1126263','GSM1126264','GSM1126265'],'h8':['GSM1126266','GSM1126267','GSM1126268']},
'endodermis':{'control':['GSM1126278','GSM1126279','GSM1126280'],'h1':['GSM1126284','GSM1126285','GSM1126286'],'h3':['GSM1126287','GSM1126288','GSM1126289'],'h8':['GSM1126290','GSM1126291','GSM1126292']},
'stele':{'control':['GSM1126302','GSM1126303','GSM1126304'],'h1':['GSM1126308','GSM1126309','GSM1126310'],'h3':['GSM1126311','GSM1126312','GSM1126313'],'h8':['GSM1126314','GSM1126315','GSM1126316']}}

def read_soft(p):
    lines=[]; inside=False
    with p.open('r',encoding='utf-8') as f:
        for line in f:
            if line.startswith('!sample_table_begin'): inside=True; continue
            if line.startswith('!sample_table_end'): break
            if inside: lines.append(line)
    d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE'])
    gsm=p.stem.split('.')[0]; d=d.rename(columns={'ID_REF':'probe','VALUE':gsm}); d[gsm]=pd.to_numeric(d[gsm],errors='raise'); return d

for ct,g in GROUPS.items():
    z=None
    for phase,folder in [('control',SRC),('h1',SRC),('h3',H3),('h8',H8)]:
        for gsm in g[phase]:
            q=read_soft(folder/f'{gsm}.soft.txt'); z=q if z is None else z.merge(q,on='probe',validate='one_to_one')
    for c in g['control']+g['h1']+g['h3']+g['h8']:
        z[c]=np.log2(z[c].to_numpy(float))
    base=z[g['control']].mean(axis=1).to_numpy(float)
    r1=z[g['h1']].mean(axis=1).to_numpy(float)-base
    r3=z[g['h3']].mean(axis=1).to_numpy(float)-base
    r8=z[g['h8']].mean(axis=1).to_numpy(float)-base
    feat2=np.std(np.column_stack([np.zeros(len(z)),r1,r3]),axis=1)
    pr=ranks_abs(r3); tr=ranks_abs(r8); ci=np.where(pr<=K)[0]; ti=np.where(tr<=K)[0]
    obs,null,_=matched_null(ci,ti,base,feat2)
    rows.append(summarize('GSE46205_salt_8h',ct,len(z),obs,null))
    null_frames.append(pd.DataFrame({'dataset':'GSE46205_salt_8h','context':ct,'null_overlap':null}))

res=pd.DataFrame(rows)
# BH correction over these six primary matched-null tests.
p=res.empirical_p_one_sided.to_numpy(float); order=np.argsort(p); q=np.empty(len(p)); prev=1.0
for rank,idx in reversed(list(enumerate(order,start=1))):
    val=min(prev,p[idx]*len(p)/rank); q[idx]=val; prev=val
res['q_bh_6_tests']=q
res.to_csv(OUT/'future_time_matched_null_summary.tsv',sep='\t',index=False)
pd.concat(null_frames,ignore_index=True).to_csv(OUT/'future_time_matched_null_distributions.tsv.gz',sep='\t',index=False,compression='gzip')
summary={'seed':SEED,'n_permutations_each':N_PERM,'n_tests':len(res),'all_q_lt_0_05':bool((res.q_bh_6_tests<.05).all()),'results':res.to_dict(orient='records')}
(OUT/'SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print('FUTURE_TIME_MATCHED_NULLS_SUCCESS')
print(res.round(5).to_string(index=False))
