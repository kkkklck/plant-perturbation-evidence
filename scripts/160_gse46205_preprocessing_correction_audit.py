from __future__ import annotations

from pathlib import Path
import io
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data' / 'prospective_router' / 'GSE46205_salt'
SRC = BASE / 'source_1h_samples'
H3 = BASE / 'target_3h_samples'
H8 = BASE / 'target_8h_samples'
OUT = ROOT / 'data' / 'publication_controls_v1' / 'gse46205_preprocessing_correction'
OUT.mkdir(parents=True, exist_ok=True)
K = 100

GROUPS = {
    'epidermis': {
        'control':['GSM1126230','GSM1126231','GSM1126232'],
        'h1':['GSM1126236','GSM1126237','GSM1126238'],
        'h3':['GSM1126239','GSM1126240','GSM1126241'],
        'h8':['GSM1126242','GSM1126243','GSM1126244'],
    },
    'cortex': {
        'control':['GSM1126254','GSM1126255','GSM1126256'],
        'h1':['GSM1126260','GSM1126261','GSM1126262'],
        'h3':['GSM1126263','GSM1126264','GSM1126265'],
        'h8':['GSM1126266','GSM1126267','GSM1126268'],
    },
    'endodermis': {
        'control':['GSM1126278','GSM1126279','GSM1126280'],
        'h1':['GSM1126284','GSM1126285','GSM1126286'],
        'h3':['GSM1126287','GSM1126288','GSM1126289'],
        'h8':['GSM1126290','GSM1126291','GSM1126292'],
    },
    'stele': {
        'control':['GSM1126302','GSM1126303','GSM1126304'],
        'h1':['GSM1126308','GSM1126309','GSM1126310'],
        'h3':['GSM1126311','GSM1126312','GSM1126313'],
        'h8':['GSM1126314','GSM1126315','GSM1126316'],
    },
}


def read_soft(p: Path) -> pd.DataFrame:
    lines = []
    inside = False
    processing = []
    with p.open('r', encoding='utf-8', errors='strict') as f:
        for line in f:
            if line.startswith('!Sample_data_processing'):
                processing.append(line.strip())
            if line.startswith('!sample_table_begin'):
                inside = True
                continue
            if line.startswith('!sample_table_end'):
                break
            if inside:
                lines.append(line)
    d = pd.read_csv(io.StringIO(''.join(lines)), sep='\t', usecols=['ID_REF','VALUE'])
    gsm = p.stem.split('.')[0]
    d = d.rename(columns={'ID_REF':'probe','VALUE':gsm})
    d['probe'] = d['probe'].astype(str)
    d[gsm] = pd.to_numeric(d[gsm], errors='raise')
    return d


def topk(pred: np.ndarray, truth: np.ndarray, k: int = K) -> dict:
    n = len(pred)
    kk = min(k, n)
    po = np.argsort(-np.abs(pred), kind='mergesort')[:kk]
    to = np.argsort(-np.abs(truth), kind='mergesort')[:kk]
    ov = np.intersect1d(po, to)
    exp = kk * kk / n
    return {
        'top100_overlap': int(len(ov)),
        'top100_random_expected': float(exp),
        'top100_enrichment': float(len(ov) / exp),
        'top100_direction': float(np.mean(np.sign(pred[ov]) == np.sign(truth[ov]))) if len(ov) else np.nan,
    }


def geometry(M: np.ndarray) -> dict:
    # M: probes x states; responses are differences on the deposited GCRMA log2 scale.
    A = M.T
    u, s, vt = np.linalg.svd(A, full_matrices=False)
    energy = s*s/(s*s).sum()
    load = u[:,0]**2
    load = load/load.sum()
    pairs = []
    loos = []
    for i in range(A.shape[0]):
        for j in range(i+1, A.shape[0]):
            pairs.append(float(spearmanr(A[i], A[j]).statistic))
        oth = np.mean(np.delete(A, i, axis=0), axis=0)
        loos.append(float(spearmanr(A[i], oth).statistic))
    eff = float(1/np.sum(load*load))
    return {
        'rank1_energy': float(energy[0]),
        'rank1_3_energy': float(energy[:3].sum()),
        'median_pairwise_rho': float(np.nanmedian(pairs)),
        'median_loo_shared_rho': float(np.nanmedian(loos)),
        'max_rank1_loading_share': float(load.max()),
        'effective_states': eff,
        'effective_state_fraction': eff / A.shape[0],
    }


# Read every required sample exactly once. Deposited VALUE is already normalized GCRMA log2 signal.
z = None
for ct, g in GROUPS.items():
    for phase, folder in [('control', SRC), ('h1', SRC), ('h3', H3), ('h8', H8)]:
        for gsm in g[phase]:
            q = read_soft(folder / f'{gsm}.soft.txt')
            z = q if z is None else z.merge(q, on='probe', validate='one_to_one')

all_samples = [x for g in GROUPS.values() for ph in ['control','h1','h3','h8'] for x in g[ph]]
for c in all_samples:
    v = z[c].to_numpy(float)
    if np.any(~np.isfinite(v)):
        raise RuntimeError(f'non-finite {c}')

responses = {1:{}, 3:{}, 8:{}}
for ct, g in GROUPS.items():
    base = z[g['control']].mean(axis=1).to_numpy(float)
    for h, key in [(1,'h1'), (3,'h3'), (8,'h8')]:
        responses[h][ct] = z[g[key]].mean(axis=1).to_numpy(float) - base

cts = list(GROUPS)
M1 = np.column_stack([responses[1][ct] for ct in cts])
M3 = np.column_stack([responses[3][ct] for ct in cts])
geom1 = geometry(M1)
geom3 = geometry(M3)

rows = []
for target_h, source_h in [(3,1),(8,3)]:
    Msrc = M1 if source_h == 1 else M3
    for j, ct in enumerate(cts):
        same = Msrc[:,j]
        shared = np.mean(np.delete(Msrc, j, axis=1), axis=1)
        truth = responses[target_h][ct]
        for model, pred in [('same_state_history', same), ('shared_other3', shared)]:
            row = {
                'target_h': target_h,
                'source_h': source_h,
                'celltype': ct,
                'model': model,
                'n_probes': len(truth),
                'spearman_full_vector': float(spearmanr(pred, truth).statistic),
                **topk(pred, truth),
            }
            rows.append(row)

# Corrected 8h strong baselines using 1h/3h response on the native deposited log2 scale.
def linear_predict(times: list[float], matrix: np.ndarray, target_time: float) -> np.ndarray:
    x = np.asarray(times, float)
    X = np.column_stack([np.ones(len(x)), x])
    beta = np.linalg.pinv(X) @ matrix.T
    return np.array([1.0, target_time]) @ beta

base_rows = []
for ct in cts:
    r1 = responses[1][ct]
    r3 = responses[3][ct]
    y = responses[8][ct]
    shared3 = np.mean(np.column_stack([responses[3][x] for x in cts if x != ct]), axis=1)
    models = {
        'latest_3h': r3,
        'older_1h': r1,
        'history_mean_1h_3h': (r1+r3)/2,
        'linear_last2_to_8h': r3 + ((8-3)/(3-1))*(r3-r1),
        'linear_with_baseline_0_1_3_to_8h': linear_predict([0,1,3], np.column_stack([np.zeros(len(r1)),r1,r3]), 8),
        'shared_other_states_3h': shared3,
    }
    for model, pred in models.items():
        base_rows.append({
            'dataset':'GSE46205_salt_8h_corrected',
            'celltype':ct,
            'model':model,
            'spearman_full_vector':float(spearmanr(pred,y).statistic),
            **topk(pred,y),
        })

res = pd.DataFrame(rows)
res.to_csv(OUT/'corrected_route_comparison.tsv', sep='\t', index=False)
baseline_df = pd.DataFrame(base_rows)
baseline_df.to_csv(OUT/'corrected_8h_strong_baselines.tsv', sep='\t', index=False)
agg = (baseline_df.groupby('model', as_index=False)
       .agg(n_contexts=('celltype','nunique'), median_rho=('spearman_full_vector','median'), mean_rho=('spearman_full_vector','mean'),
            median_top100_overlap=('top100_overlap','median'), mean_top100_overlap=('top100_overlap','mean'),
            median_top100_enrichment=('top100_enrichment','median'), median_top100_direction=('top100_direction','median')))
agg.to_csv(OUT/'corrected_8h_strong_baseline_aggregate.tsv', sep='\t', index=False)

# Compare corrected values against historical double-log outputs where available.
historical = []
for target_h, folder, vf_name, tf_name in [
    (3, BASE/'router_blind', 'blind_3h_vector_results.tsv', 'blind_3h_topk_results.tsv'),
    (8, BASE/'task_8h_blind', 'blind_8h_vector_results.tsv', 'blind_8h_topk_results.tsv'),
]:
    vf = folder / vf_name
    tf = folder / tf_name
    if vf.exists() and tf.exists():
        v = pd.read_csv(vf, sep='\t')
        t = pd.read_csv(tf, sep='\t')
        historical.append({'target_h':target_h,'vector_columns':v.columns.tolist(),'topk_columns':t.columns.tolist()})

summary = {
    'status': 'POST_REVEAL_PREPROCESSING_CORRECTION_AUDIT',
    'reason': 'GEO GSE46205 deposited VALUE is already GCRMA normalized log2 signal; historical scripts 131-135 applied an additional log2 transform.',
    'correction': 'Use deposited VALUE directly; responses are treatment mean minus matched control mean on that native log2 scale.',
    'historical_freeze_files_modified': False,
    'geometry_1h_corrected': geom1,
    'geometry_3h_corrected': geom3,
    'route_comparison': {},
    'prospective_status_note': 'Because the implementation error was discovered after target reveal, corrected analyses must not be described as untouched prospective freeze-reveal evidence. The original task split/rule was frozen, but corrected numerical results are post-reveal implementation audits.',
}
for target_h in [3,8]:
    q = res[res.target_h == target_h]
    w = q.pivot(index='celltype', columns='model', values='spearman_full_vector')
    wt = q.pivot(index='celltype', columns='model', values='top100_overlap')
    summary['route_comparison'][f'{target_h}h'] = {
        'same_state_wins_rho': int((w['same_state_history'] > w['shared_other3']).sum()),
        'same_state_wins_top100': int((wt['same_state_history'] > wt['shared_other3']).sum()),
        'same_state_ties_top100': int((wt['same_state_history'] == wt['shared_other3']).sum()),
        'median_same_state_rho': float(w['same_state_history'].median()),
        'median_shared_rho': float(w['shared_other3'].median()),
        'median_same_state_top100': float(wt['same_state_history'].median()),
        'median_shared_top100': float(wt['shared_other3'].median()),
    }

(OUT/'SUMMARY.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

print('GSE46205_PREPROCESSING_CORRECTION_AUDIT_SUCCESS')
print(json.dumps(summary, ensure_ascii=False, indent=2))
print('\nCORRECTED ROUTE COMPARISON')
print(res.round(5).to_string(index=False))
print('\nCORRECTED 8H BASELINE AGGREGATE')
print(agg.round(5).to_string(index=False))
