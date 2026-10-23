from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data' / 'external_heatshock_crossstudy' / 'GSE290176'
HEAT_OUT = BASE / 'crossstudy_blind'
TARGET = BASE / 'GSE290176_counts_HS_Ath.txt.gz'
TRUE_PRED = HEAT_OUT / 'frozen_crossstudy_prediction.tsv.gz'
OUT = ROOT / 'data' / 'publication_controls_v1' / 'heatshock_crossstudy_negative_controls'
OUT.mkdir(parents=True, exist_ok=True)

CK = [f'Col0-CK{i}' for i in range(1, 7)]
HS = [f'Col0-HS{i}' for i in range(1, 4)]
SEED = 20261021
N_PERM = 10000
K = 100


def norm_gene(x: str) -> str:
    return str(x).strip().strip('"').upper()


def read_target() -> pd.Series:
    counts = pd.read_csv(TARGET, sep='\t', compression='gzip')
    req = ['Geneid'] + CK + HS
    miss = [c for c in req if c not in counts.columns]
    if miss:
        raise RuntimeError(f'missing target columns: {miss}')
    counts = counts[req].copy()
    counts['gene'] = counts['Geneid'].map(norm_gene)
    for c in CK + HS:
        counts[c] = pd.to_numeric(counts[c], errors='coerce').fillna(0.0)
    counts = counts.groupby('gene', as_index=False)[CK + HS].sum()
    libs = counts[CK + HS].sum(axis=0)
    lc = pd.DataFrame({'gene': counts['gene']})
    for c in CK + HS:
        cpm = counts[c].to_numpy(float) / float(libs[c]) * 1e6
        lc[c] = np.log2(cpm + 1.0)
    lc['target_response'] = lc[HS].mean(axis=1) - lc[CK].mean(axis=1)
    return lc.set_index('gene')['target_response']


def read_program(path: Path, value_col: str, *, filter_col: str | None = None) -> pd.Series:
    d = pd.read_csv(path, sep='\t')
    d['gene'] = d['gene'].map(norm_gene)
    if filter_col is not None:
        if filter_col not in d.columns:
            raise RuntimeError(f'{filter_col} missing from {path}')
        d = d[d[filter_col].astype(bool)].copy()
    d[value_col] = pd.to_numeric(d[value_col], errors='coerce')
    d = d[np.isfinite(d[value_col])].drop_duplicates('gene', keep='first')
    return d.set_index('gene')[value_col]


def ranks_abs(v: np.ndarray) -> np.ndarray:
    # stable descending rank by absolute value
    order = np.argsort(-np.abs(v), kind='mergesort')
    r = np.empty(len(v), dtype=int)
    r[order] = np.arange(1, len(v) + 1)
    return r


def score(program: np.ndarray, truth: np.ndarray, k: int = K) -> dict:
    if len(program) != len(truth):
        raise ValueError('length mismatch')
    n = len(program)
    rho = float(spearmanr(program, truth).statistic)
    pr = ranks_abs(program)
    tr = ranks_abs(truth)
    kk = min(k, n)
    pi = np.where(pr <= kk)[0]
    ti = np.where(tr <= kk)[0]
    ov = np.intersect1d(pi, ti, assume_unique=False)
    overlap = int(len(ov))
    expected = kk * kk / n
    sign = float(np.mean(np.sign(program[ov]) == np.sign(truth[ov]))) if overlap else np.nan
    return {
        'n_genes': n,
        'spearman_full_vector': rho,
        'top100_overlap': overlap,
        'top100_random_expected': expected,
        'top100_enrichment': overlap / expected,
        'top100_direction_agreement': sign,
    }


target = read_target()
true_heat = read_program(TRUE_PRED, 'shared_prediction')

programs: dict[str, pd.Series] = {}
# Wrong stress history from an independent Arabidopsis salt time course.
salt = pd.read_csv(ROOT / 'data' / 'external_salt' / 'GSE153103' / 'plantvc_v1_1_12h_blind' / 'FROZEN_PRETARGET_RANKING.tsv.gz', sep='\t')
salt['gene'] = salt['gene'].map(norm_gene)
salt = salt[salt['technical_universe'].astype(bool)].drop_duplicates('gene')
programs['salt_6h_response'] = salt.set_index('gene')['r6'].astype(float)
programs['salt_3h_response'] = salt.set_index('gene')['r3'].astype(float)

# Wrong immune programs: use exactly the already frozen source-state shared responses.
for stim in ['Pep1', 'SCOOP10', 'SCREW1', 'flg22', 'Chitin']:
    p = ROOT / 'data' / 'prospective_task' / 'immune_missing_state' / 'blind_v0_1' / f'FROZEN_{stim}_PREDICTION.tsv.gz'
    programs[f'immume_{stim}'] = read_program(p, 'predicted_response')

rows = []
for name, wrong in programs.items():
    genes = sorted(set(target.index) & set(true_heat.index) & set(wrong.index))
    if len(genes) < 1000:
        raise RuntimeError(f'too few common genes for {name}: {len(genes)}')
    tv = target.reindex(genes).to_numpy(float)
    hv = true_heat.reindex(genes).to_numpy(float)
    wv = wrong.reindex(genes).to_numpy(float)
    hs = score(hv, tv)
    ws = score(wv, tv)
    row = {
        'wrong_program': name,
        'common_universe_genes': len(genes),
        'true_heat_rho': hs['spearman_full_vector'],
        'wrong_program_rho': ws['spearman_full_vector'],
        'delta_rho_true_minus_wrong': hs['spearman_full_vector'] - ws['spearman_full_vector'],
        'true_heat_top100_overlap': hs['top100_overlap'],
        'wrong_program_top100_overlap': ws['top100_overlap'],
        'delta_top100_overlap_true_minus_wrong': hs['top100_overlap'] - ws['top100_overlap'],
        'true_heat_top100_enrichment': hs['top100_enrichment'],
        'wrong_program_top100_enrichment': ws['top100_enrichment'],
        'true_heat_top100_direction': hs['top100_direction_agreement'],
        'wrong_program_top100_direction': ws['top100_direction_agreement'],
    }
    rows.append(row)

pairwise = pd.DataFrame(rows)
pairwise.to_csv(OUT / 'wrong_perturbation_pairwise.tsv', sep='\t', index=False)

# Exact-distribution label-shuffle null on the full frozen heat universe.
genes = sorted(set(target.index) & set(true_heat.index))
tv = target.reindex(genes).to_numpy(float)
hv = true_heat.reindex(genes).to_numpy(float)
obs = score(hv, tv)
N = len(genes)
truth_rank = ranks_abs(tv)
truth_top = set(np.where(truth_rank <= K)[0].tolist())
# Spearman can be computed as Pearson correlation of rank vectors. Precompute target signed ranks.
truth_signed_rank = pd.Series(tv).rank(method='average').to_numpy(float)
true_signed_rank = pd.Series(hv).rank(method='average').to_numpy(float)
obs_rho = float(np.corrcoef(true_signed_rank, truth_signed_rank)[0,1])
obs_top = obs['top100_overlap']
rng = np.random.default_rng(SEED)
null_rho = np.empty(N_PERM, dtype=float)
null_top = np.empty(N_PERM, dtype=int)
base_idx = np.arange(N)
for i in range(N_PERM):
    perm = rng.permutation(base_idx)
    rp = true_signed_rank[perm]
    null_rho[i] = float(np.corrcoef(rp, truth_signed_rank)[0,1])
    # absolute top100 identities are permuted with gene labels
    perm_top = set(perm[np.argsort(-np.abs(hv[perm]), kind='mergesort')[:K]].tolist())
    # The previous line would select original top genes after permutation incorrectly.
    # Correct mapping: program values hv[perm] are assigned to positions 0..N-1.
    assigned_top_positions = set(np.argsort(-np.abs(hv[perm]), kind='mergesort')[:K].tolist())
    null_top[i] = len(assigned_top_positions & truth_top)

emp = {
    'seed': SEED,
    'n_permutations': N_PERM,
    'n_genes': N,
    'observed_full_rho': obs_rho,
    'null_rho_mean': float(np.mean(null_rho)),
    'null_rho_sd': float(np.std(null_rho, ddof=1)),
    'empirical_p_rho_one_sided': float((1 + np.sum(null_rho >= obs_rho)) / (N_PERM + 1)),
    'observed_top100_overlap': int(obs_top),
    'null_top100_mean': float(np.mean(null_top)),
    'null_top100_q95': float(np.quantile(null_top, 0.95)),
    'null_top100_q99': float(np.quantile(null_top, 0.99)),
    'empirical_p_top100_one_sided': float((1 + np.sum(null_top >= obs_top)) / (N_PERM + 1)),
}
(OUT / 'gene_label_shuffle_summary.json').write_text(json.dumps(emp, ensure_ascii=False, indent=2), encoding='utf-8')
pd.DataFrame({'rho': null_rho, 'top100_overlap': null_top}).to_csv(OUT / 'gene_label_shuffle_null.tsv.gz', sep='\t', index=False, compression='gzip')

summary = {
    'primary_full_universe': obs,
    'shuffle_null': emp,
    'wrong_programs_n': int(len(pairwise)),
    'true_heat_wins_rho_pairwise': int((pairwise['delta_rho_true_minus_wrong'] > 0).sum()),
    'true_heat_wins_top100_pairwise': int((pairwise['delta_top100_overlap_true_minus_wrong'] > 0).sum()),
    'true_heat_ties_top100_pairwise': int((pairwise['delta_top100_overlap_true_minus_wrong'] == 0).sum()),
}
(OUT / 'SUMMARY.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

print('HEATSHOCK_CROSSSTUDY_NEGATIVE_CONTROLS_SUCCESS')
print('\nPAIRWISE_WRONG_PROGRAMS')
print(pairwise.round(4).to_string(index=False))
print('\nSHUFFLE_NULL')
print(json.dumps(emp, ensure_ascii=False, indent=2))
