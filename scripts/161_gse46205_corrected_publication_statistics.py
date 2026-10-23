from __future__ import annotations

from pathlib import Path
from itertools import product
import io
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data' / 'prospective_router' / 'GSE46205_salt'
SRC = BASE / 'source_1h_samples'
H3 = BASE / 'target_3h_samples'
H8 = BASE / 'target_8h_samples'
CORR = ROOT / 'data' / 'publication_controls_v1' / 'gse46205_preprocessing_correction'
OUT = ROOT / 'data' / 'publication_controls_v1' / 'gse46205_corrected_publication_stats'
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261021
B = 10000
K = 100
rng = np.random.default_rng(SEED)

GROUPS = {
'epidermis':{'control':['GSM1126230','GSM1126231','GSM1126232'],'h1':['GSM1126236','GSM1126237','GSM1126238'],'h3':['GSM1126239','GSM1126240','GSM1126241'],'h8':['GSM1126242','GSM1126243','GSM1126244']},
'cortex':{'control':['GSM1126254','GSM1126255','GSM1126256'],'h1':['GSM1126260','GSM1126261','GSM1126262'],'h3':['GSM1126263','GSM1126264','GSM1126265'],'h8':['GSM1126266','GSM1126267','GSM1126268']},
'endodermis':{'control':['GSM1126278','GSM1126279','GSM1126280'],'h1':['GSM1126284','GSM1126285','GSM1126286'],'h3':['GSM1126287','GSM1126288','GSM1126289'],'h8':['GSM1126290','GSM1126291','GSM1126292']},
'stele':{'control':['GSM1126302','GSM1126303','GSM1126304'],'h1':['GSM1126308','GSM1126309','GSM1126310'],'h3':['GSM1126311','GSM1126312','GSM1126313'],'h8':['GSM1126314','GSM1126315','GSM1126316']},
}


def read_soft(p: Path) -> pd.DataFrame:
    lines=[]; inside=False
    with p.open('r',encoding='utf-8',errors='strict') as f:
        for line in f:
            if line.startswith('!sample_table_begin'): inside=True; continue
            if line.startswith('!sample_table_end'): break
            if inside: lines.append(line)
    d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE'])
    gsm=p.stem.split('.')[0]
    d=d.rename(columns={'ID_REF':'probe','VALUE':gsm})
    d['probe']=d['probe'].astype(str); d[gsm]=pd.to_numeric(d[gsm],errors='raise')
    return d


def qbin(x, q=5):
    return pd.qcut(pd.Series(np.asarray(x,float)).rank(method='average'), q=q, labels=False, duplicates='drop').to_numpy(int)


def matched_null(candidate_idx, truth_idx, feat1, feat2):
    b1=qbin(feat1,5); b2=qbin(feat2,5); strata=list(zip(b1,b2))
    counts={}
    for i in candidate_idx: counts[strata[i]]=counts.get(strata[i],0)+1
    pools={s:np.array([i for i,z in enumerate(strata) if z==s],int) for s in counts}
    truth=set(map(int,truth_idx)); null=np.empty(B,int)
    for b in range(B):
        chosen=[]
        for s,n in counts.items():
            chosen.extend(rng.choice(pools[s],size=n,replace=False).tolist())
        null[b]=len(set(chosen)&truth)
    obs=len(set(map(int,candidate_idx))&truth)
    return obs,null


def bh(pvals):
    p=np.asarray(pvals,float); n=len(p); order=np.argsort(p); ranked=p[order]
    q=ranked*n/np.arange(1,n+1); q=np.minimum.accumulate(q[::-1])[::-1]; q=np.clip(q,0,1)
    out=np.empty(n,float); out[order]=q; return out


def exact_signflip(diffs):
    d=np.asarray(diffs,float); obs=float(np.mean(d)); vals=[]
    for signs in product([-1.0,1.0],repeat=len(d)):
        vals.append(float(np.mean(d*np.asarray(signs))))
    vals=np.asarray(vals)
    return obs,float(np.mean(vals>=obs-1e-12)),len(vals)

# Native deposited GCRMA values are already log2 normalized; no extra transform.
z=None
for ct,g in GROUPS.items():
    for phase,folder in [('control',SRC),('h1',SRC),('h3',H3),('h8',H8)]:
        for gsm in g[phase]:
            q=read_soft(folder/f'{gsm}.soft.txt'); z=q if z is None else z.merge(q,on='probe',validate='one_to_one')

resp={1:{},3:{},8:{}}; baseline={}
for ct,g in GROUPS.items():
    base=z[g['control']].mean(axis=1).to_numpy(float); baseline[ct]=base
    for h,key in [(1,'h1'),(3,'h3'),(8,'h8')]:
        resp[h][ct]=z[g[key]].mean(axis=1).to_numpy(float)-base

# Corrected matched nulls for four GSE46205 3h->8h contexts.
null_rows=[]; null_frames=[]
for ct in GROUPS:
    pred=resp[3][ct]; truth=resp[8][ct]; r1=resp[1][ct]
    pr=np.argsort(-np.abs(pred),kind='mergesort')[:K]; tr=np.argsort(-np.abs(truth),kind='mergesort')[:K]
    history_disp=np.std(np.column_stack([np.zeros(len(pred)),r1,pred]),axis=1)
    obs,null=matched_null(pr,tr,baseline[ct],history_disp)
    null_rows.append({
        'dataset':'GSE46205_salt_8h','context':ct,'n_genes':len(pred),'k':K,'observed_overlap':obs,
        'null_mean':float(null.mean()),'null_sd':float(null.std(ddof=1)),'null_q95':float(np.quantile(null,.95)),
        'null_q99':float(np.quantile(null,.99)),'matched_enrichment':float(obs/null.mean()),
        'empirical_p_one_sided':float((1+np.sum(null>=obs))/(B+1)),
        'matching_variables':'native GCRMA control expression quintile + 0h/1h/3h history dispersion quintile',
        'preprocessing_status':'corrected_post_reveal_native_GCRMA_log2',
    })
    null_frames.append(pd.DataFrame({'dataset':'GSE46205_salt_8h','context':ct,'null_overlap':null}))
new_null=pd.DataFrame(null_rows)
new_null.to_csv(OUT/'gse46205_corrected_matched_null.tsv',sep='\t',index=False)
pd.concat(null_frames,ignore_index=True).to_csv(OUT/'gse46205_corrected_matched_null_distributions.tsv.gz',sep='\t',index=False,compression='gzip')

# Replace corrupted GSE46205 rows in the prior 16-test table and recalculate BH across all future-time tests.
old=pd.read_csv(ROOT/'data/publication_controls_v1/future_time_matched_nulls_complete/future_time_matched_null_complete_summary.tsv',sep='\t')
keep=old[old.dataset!='GSE46205_salt_8h'].copy()
# normalize columns
allcols=sorted(set(keep.columns)|set(new_null.columns))
for c in allcols:
    if c not in keep: keep[c]=np.nan
    if c not in new_null: new_null[c]=np.nan
combined=pd.concat([keep[allcols],new_null[allcols]],ignore_index=True)
combined['q_bh_all_future_time_tests_corrected']=bh(combined.empirical_p_one_sided.to_numpy(float))
combined.to_csv(OUT/'future_time_matched_null_complete_corrected.tsv',sep='\t',index=False)

# Corrected paired same-state vs shared statistics for both 1->3 and 3->8.
route=pd.read_csv(CORR/'corrected_route_comparison.tsv',sep='\t')
pair_rows=[]; detail=[]
for target_h in [3,8]:
    q=route[route.target_h==target_h]
    for metric in ['spearman_full_vector','top100_overlap']:
        w=q.pivot(index='celltype',columns='model',values=metric)
        d=w['same_state_history']-w['shared_other3']
        obs,p,nperm=exact_signflip(d)
        pair_rows.append({
            'dataset':f'GSE46205_salt_{target_h}h_corrected','metric':metric,'n_contexts':len(d),
            'mean_paired_difference':obs,'median_paired_difference':float(np.median(d)),
            'wins':int((d>0).sum()),'ties':int((d==0).sum()),'losses':int((d<0).sum()),
            'exact_signflip_p_one_sided':p,'n_permutations':nperm,
            'preprocessing_status':'corrected_post_reveal_native_GCRMA_log2',
        })
        for ct,x in d.items(): detail.append({'target_h':target_h,'celltype':ct,'metric':metric,'paired_difference':x})
pd.DataFrame(pair_rows).to_csv(OUT/'corrected_paired_context_inference.tsv',sep='\t',index=False)
pd.DataFrame(detail).to_csv(OUT/'corrected_paired_context_differences.tsv',sep='\t',index=False)

summary={
    'status':'POST_REVEAL_CORRECTED_PUBLICATION_STATISTICS',
    'seed':SEED,'n_permutations_matched_null':B,
    'gse46205_matched_null_results':new_null.to_dict(orient='records'),
    'paired_inference':pd.DataFrame(pair_rows).to_dict(orient='records'),
    'n_future_time_tests_after_replacement':int(len(combined)),
    'n_future_time_q_lt_0_05':int((combined.q_bh_all_future_time_tests_corrected<.05).sum()),
    'interpretation':'Corrected GSE46205 numerical analyses preserve the full-vector route conclusion; 1h->3h Top100 becomes 3 wins + 1 tie, while 3h->8h remains 4/4 wins. These corrected results are post-reveal implementation audits, not untouched prospective scores.'
}
(OUT/'SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print('GSE46205_CORRECTED_PUBLICATION_STATS_SUCCESS')
print('\nMATCHED NULL')
print(new_null.round(5).to_string(index=False))
print('\nPAIRED')
print(pd.DataFrame(pair_rows).round(5).to_string(index=False))
print('\nFUTURE TIME BH:')
print('n=',len(combined),'q<.05=',int((combined.q_bh_all_future_time_tests_corrected<.05).sum()))
