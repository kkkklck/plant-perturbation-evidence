from __future__ import annotations

from pathlib import Path
import io
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'publication_controls_v1' / 'future_time_strong_baselines'
OUT.mkdir(parents=True, exist_ok=True)
K = 100


def ranks_abs(v: np.ndarray) -> np.ndarray:
    order = np.argsort(-np.abs(v), kind='mergesort')
    r = np.empty(len(v), dtype=int)
    r[order] = np.arange(1, len(v) + 1)
    return r


def score(name: str, dataset: str, context: str, pred: np.ndarray, truth: np.ndarray) -> dict:
    n = len(truth)
    pr = ranks_abs(pred)
    tr = ranks_abs(truth)
    pi = np.where(pr <= min(K,n))[0]
    ti = np.where(tr <= min(K,n))[0]
    ov = np.intersect1d(pi, ti)
    expected = min(K,n)**2 / n
    return {
        'dataset': dataset,
        'context': context,
        'model': name,
        'n_genes': n,
        'spearman_full_vector': float(spearmanr(pred, truth).statistic),
        'top100_overlap': int(len(ov)),
        'top100_enrichment': float(len(ov) / expected),
        'top100_direction_agreement': float(np.mean(np.sign(pred[ov]) == np.sign(truth[ov]))) if len(ov) else np.nan,
    }


def linear_predict(times: list[float], matrix: np.ndarray, target_time: float) -> np.ndarray:
    # matrix genes x timepoints; ordinary least squares independently per gene
    x = np.asarray(times, dtype=float)
    X = np.column_stack([np.ones(len(x)), x])
    beta = np.linalg.pinv(X) @ matrix.T  # 2 x genes
    return (np.array([1.0, target_time]) @ beta).astype(float)


rows = []

# -----------------------------------------------------------------------------
# GSE153103 salt: 0h baseline response = 0, observed response histories at 1/3/6h.
# Evaluate 6h and 12h targets on exactly the original technical universe.
# -----------------------------------------------------------------------------
base6 = ROOT / 'data' / 'external_salt' / 'GSE153103' / 'plantvc_v1_6h_blind'
f6 = pd.read_csv(base6 / 'FROZEN_PRETARGET_RANKING.tsv.gz', sep='\t')
t6 = pd.read_csv(base6 / 'target_6h_response.tsv.gz', sep='\t')
z6 = f6[f6['technical_universe'].astype(bool)].merge(t6[['gene','r6']], on='gene', validate='one_to_one')
r1 = z6['r1'].to_numpy(float); r3 = z6['r3'].to_numpy(float); truth6 = z6['r6'].to_numpy(float)
models6 = {
    'latest_3h': r3,
    'older_1h': r1,
    'history_mean_1h_3h': (r1+r3)/2,
    'linear_last2_to_6h': r3 + ((6-3)/(3-1))*(r3-r1),
    'linear_with_baseline_0_1_3_to_6h': linear_predict([0,1,3], np.column_stack([np.zeros(len(r1)), r1, r3]), 6),
}
for nm,p in models6.items(): rows.append(score(nm,'GSE153103_salt_6h','bulk',p,truth6))

base12 = ROOT / 'data' / 'external_salt' / 'GSE153103' / 'plantvc_v1_1_12h_blind'
f12 = pd.read_csv(base12 / 'FROZEN_PRETARGET_RANKING.tsv.gz', sep='\t')
t12 = pd.read_csv(base12 / 'target_12h_response.tsv.gz', sep='\t')
z12 = f12[f12['technical_universe'].astype(bool)].merge(t12[['gene','r12']], on='gene', validate='one_to_one')
r1 = z12['r1'].to_numpy(float); r3 = z12['r3'].to_numpy(float); r6 = z12['r6'].to_numpy(float); truth12 = z12['r12'].to_numpy(float)
models12 = {
    'latest_6h': r6,
    'older_3h': r3,
    'oldest_1h': r1,
    'history_mean_1h_3h_6h': (r1+r3+r6)/3,
    'linear_last2_3h_6h_to_12h': r6 + ((12-6)/(6-3))*(r6-r3),
    'linear_all_0_1_3_6_to_12h': linear_predict([0,1,3,6], np.column_stack([np.zeros(len(r1)),r1,r3,r6]), 12),
}
for nm,p in models12.items(): rows.append(score(nm,'GSE153103_salt_12h','bulk',p,truth12))

# -----------------------------------------------------------------------------
# GSE46205 salt 8h: reconstruct 1h and 3h same-state responses from frozen/revealed samples,
# compare latest history to earlier history, mean, linear extrapolation and cross-state mean.
# -----------------------------------------------------------------------------
BASE = ROOT / 'data' / 'prospective_router' / 'GSE46205_salt'
SRC = BASE / 'source_1h_samples'
H3 = BASE / 'target_3h_samples'
H8 = BASE / 'target_8h_samples'
GROUPS = {
'epidermis': {'control':['GSM1126230','GSM1126231','GSM1126232'], 'h1':['GSM1126236','GSM1126237','GSM1126238'], 'h3':['GSM1126239','GSM1126240','GSM1126241'], 'h8':['GSM1126242','GSM1126243','GSM1126244']},
'cortex': {'control':['GSM1126254','GSM1126255','GSM1126256'], 'h1':['GSM1126260','GSM1126261','GSM1126262'], 'h3':['GSM1126263','GSM1126264','GSM1126265'], 'h8':['GSM1126266','GSM1126267','GSM1126268']},
'endodermis': {'control':['GSM1126278','GSM1126279','GSM1126280'], 'h1':['GSM1126284','GSM1126285','GSM1126286'], 'h3':['GSM1126287','GSM1126288','GSM1126289'], 'h8':['GSM1126290','GSM1126291','GSM1126292']},
'stele': {'control':['GSM1126302','GSM1126303','GSM1126304'], 'h1':['GSM1126308','GSM1126309','GSM1126310'], 'h3':['GSM1126311','GSM1126312','GSM1126313'], 'h8':['GSM1126314','GSM1126315','GSM1126316']},
}

def read_soft(p: Path) -> pd.DataFrame:
    lines=[]; inside=False
    with p.open('r',encoding='utf-8',errors='strict') as f:
        for line in f:
            if line.startswith('!sample_table_begin'):
                inside=True; continue
            if line.startswith('!sample_table_end'):
                break
            if inside: lines.append(line)
    d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE'])
    gsm=p.stem.split('.')[0]
    d=d.rename(columns={'ID_REF':'probe','VALUE':gsm})
    d['probe']=d['probe'].astype(str)
    d[gsm]=pd.to_numeric(d[gsm],errors='raise')
    return d

# Merge all required samples once.
z=None
for ct,g in GROUPS.items():
    for phase, folder in [('control',SRC),('h1',SRC),('h3',H3),('h8',H8)]:
        for gsm in g[phase]:
            q=read_soft(folder/f'{gsm}.soft.txt')
            z=q if z is None else z.merge(q,on='probe',validate='one_to_one')
all_samples=[x for g in GROUPS.values() for phase in ['control','h1','h3','h8'] for x in g[phase]]
for c in all_samples:
    v=z[c].to_numpy(float)
    if np.any(v<=0) or np.any(~np.isfinite(v)): raise RuntimeError(f'invalid microarray value {c}')
    z[c]=np.log2(v)

r1_by={}; r3_by={}; r8_by={}
for ct,g in GROUPS.items():
    base=z[g['control']].mean(axis=1).to_numpy(float)
    r1_by[ct]=z[g['h1']].mean(axis=1).to_numpy(float)-base
    r3_by[ct]=z[g['h3']].mean(axis=1).to_numpy(float)-base
    r8_by[ct]=z[g['h8']].mean(axis=1).to_numpy(float)-base
for ct in GROUPS:
    r1=r1_by[ct]; r3=r3_by[ct]; truth=r8_by[ct]
    shared=np.mean(np.column_stack([r3_by[o] for o in GROUPS if o!=ct]),axis=1)
    models={
        'latest_3h':r3,
        'older_1h':r1,
        'history_mean_1h_3h':(r1+r3)/2,
        'linear_last2_to_8h':r3 + ((8-3)/(3-1))*(r3-r1),
        'linear_with_baseline_0_1_3_to_8h':linear_predict([0,1,3],np.column_stack([np.zeros(len(r1)),r1,r3]),8),
        'shared_other_states_3h':shared,
    }
    for nm,p in models.items(): rows.append(score(nm,'GSE46205_salt_8h',ct,p,truth))

res=pd.DataFrame(rows)
res.to_csv(OUT/'future_time_baseline_results.tsv',sep='\t',index=False)

# Aggregate by dataset/model with equal weight to cell-type contexts where applicable.
agg=(res.groupby(['dataset','model'],as_index=False)
     .agg(n_contexts=('context','size'),median_rho=('spearman_full_vector','median'),mean_rho=('spearman_full_vector','mean'),
          median_top100_overlap=('top100_overlap','median'),mean_top100_overlap=('top100_overlap','mean'),
          median_top100_enrichment=('top100_enrichment','median'),median_direction=('top100_direction_agreement','median')))
agg.to_csv(OUT/'future_time_baseline_aggregate.tsv',sep='\t',index=False)

# Per-dataset rank of models by median/mean rho and top100 overlap.
summary={}
for ds,q in agg.groupby('dataset'):
    qq=q.sort_values(['median_rho','median_top100_overlap'],ascending=False)
    summary[ds]={
        'best_by_median_rho':str(qq.iloc[0]['model']),
        'latest_history_row':q[q['model'].str.startswith('latest_')].to_dict(orient='records'),
        'all_models':q.to_dict(orient='records'),
    }
(OUT/'SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')

print('FUTURE_TIME_STRONG_BASELINES_SUCCESS')
print('\nAGGREGATE')
print(agg.round(4).to_string(index=False))
print('\nGSE46205_CONTEXT_DETAIL')
print(res[res.dataset=='GSE46205_salt_8h'].round(4).to_string(index=False))
