from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'publication_controls_v1' / 'heatshock_matched_nulls'
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261021
N_PERM = 10000
K = 100
rng = np.random.default_rng(SEED)


def quantile_bins(x: np.ndarray, q: int) -> np.ndarray:
    s = pd.Series(np.asarray(x, float))
    return pd.qcut(s.rank(method='average'), q=q, labels=False, duplicates='drop').to_numpy(int)


def matched_null(candidate_idx: np.ndarray, truth_idx: np.ndarray, features: list[tuple[np.ndarray, int]]) -> np.ndarray:
    bins = [quantile_bins(x, q) for x, q in features]
    strata = list(zip(*bins))
    counts: dict[tuple[int, ...], int] = {}
    for i in candidate_idx:
        s = strata[int(i)]
        counts[s] = counts.get(s, 0) + 1
    pools = {s: np.array([i for i, z in enumerate(strata) if z == s], dtype=int) for s in counts}
    truth = set(map(int, truth_idx))
    null = np.empty(N_PERM, int)
    for b in range(N_PERM):
        selected: list[int] = []
        for s, n in counts.items():
            pool = pools[s]
            if len(pool) < n:
                raise RuntimeError(f'stratum {s}: need {n}, have {len(pool)}')
            selected.extend(rng.choice(pool, size=n, replace=False).tolist())
        null[b] = len(set(selected) & truth)
    return null


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q_sorted = p[order] * n / np.arange(1, n + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty(n, float)
    q[order] = np.clip(q_sorted, 0, 1)
    return q


def score_case(dataset: str, context: str, pred: np.ndarray, truth: np.ndarray, donors: np.ndarray) -> tuple[dict, pd.DataFrame]:
    pred = np.asarray(pred, float)
    truth = np.asarray(truth, float)
    donors = np.asarray(donors, float)
    if donors.ndim != 2 or donors.shape[0] != len(pred):
        raise RuntimeError('invalid donor matrix')
    cand = np.argsort(-np.abs(pred), kind='mergesort')[:K]
    truth_idx = np.argsort(-np.abs(truth), kind='mergesort')[:K]
    donor_sd = np.std(donors, axis=1)
    null = matched_null(cand, truth_idx, [(np.abs(pred), 20), (donor_sd, 5)])
    obs = int(len(np.intersect1d(cand, truth_idx)))
    mean = float(np.mean(null))
    row = {
        'dataset': dataset,
        'context': context,
        'n_genes': int(len(pred)),
        'k': K,
        'observed_overlap': obs,
        'null_mean': mean,
        'null_sd': float(np.std(null, ddof=1)),
        'null_q95': float(np.quantile(null, 0.95)),
        'null_q99': float(np.quantile(null, 0.99)),
        'matched_enrichment': float(obs / mean) if mean > 0 else np.inf,
        'empirical_p_one_sided': float((1 + np.sum(null >= obs)) / (N_PERM + 1)),
        'matching_variables': '|shared heat response| ventile + donor-response SD quintile',
    }
    nf = pd.DataFrame({'dataset': dataset, 'context': context, 'null_overlap': null})
    return row, nf


def missing_state_cases() -> tuple[list[dict], list[pd.DataFrame]]:
    base = ROOT / 'data' / 'external_heatshock' / 'GSE121619' / 'plantvc_v1_heatshock_blind'
    pred = pd.read_csv(base / 'frozen_prediction.tsv.gz', sep='\t')
    pred = pred[pred['technical'].astype(bool)].copy()
    pred['gene'] = pred['gene'].astype(str)
    genes = pred['gene'].tolist()
    pv = pred['shared_prediction'].to_numpy(float)
    donors = pred[['donor_5', 'donor_6', 'donor_8']].to_numpy(float)

    target = pd.read_csv(base / 'target_profiles.tsv.gz', sep='\t')
    target = target[target['technical'].astype(bool)].copy()
    target['gene'] = target['gene'].astype(str)
    rows: list[dict] = []
    nulls: list[pd.DataFrame] = []
    for cl in [2, 3, 10]:
        q = target[target['target_cluster'] == cl].set_index('gene').reindex(genes)
        if q['target_response'].isna().any():
            raise RuntimeError(f'missing target genes cluster {cl}')
        row, nf = score_case('GSE121619_heatshock_missing_state', f'cluster_{cl}', pv,
                             q['target_response'].to_numpy(float), donors)
        rows.append(row)
        nulls.append(nf)
    return rows, nulls


def cross_experiment_case() -> tuple[dict, pd.DataFrame]:
    base = ROOT / 'data' / 'external_heatshock_crossstudy' / 'GSE290176'
    pred = pd.read_csv(base / 'crossstudy_blind' / 'frozen_crossstudy_prediction.tsv.gz', sep='\t')
    pred['gene'] = pred['gene'].astype(str).str.upper()
    donors = pred[['donor_5', 'donor_6', 'donor_8']].to_numpy(float)
    pv = pred['shared_prediction'].to_numpy(float)

    counts = pd.read_csv(base / 'GSE290176_counts_HS_Ath.txt.gz', sep='\t')
    ck = [f'Col0-CK{i}' for i in range(1, 7)]
    hs = [f'Col0-HS{i}' for i in range(1, 4)]
    counts = counts[['Geneid'] + ck + hs].copy()
    counts['gene'] = counts['Geneid'].astype(str).str.strip().str.strip('"').str.upper()
    for c in ck + hs:
        counts[c] = pd.to_numeric(counts[c], errors='coerce').fillna(0.0)
    counts = counts.groupby('gene', as_index=False)[ck + hs].sum()
    libs = counts[ck + hs].sum(axis=0)
    for c in ck + hs:
        counts[c] = np.log2(counts[c].to_numpy(float) / float(libs[c]) * 1e6 + 1)
    counts['target_response'] = counts[hs].mean(axis=1) - counts[ck].mean(axis=1)
    z = pred.merge(counts[['gene', 'target_response']], on='gene', validate='one_to_one')
    if len(z) != len(pred):
        raise RuntimeError('cross-study target merge lost genes')
    # Preserve prediction order when merging.
    truth = z['target_response'].to_numpy(float)
    return score_case('GSE290176_heatshock_crossstudy', 'whole_seedling', pv, truth, donors)


def main() -> None:
    rows, nulls = missing_state_cases()
    row, nf = cross_experiment_case()
    rows.append(row)
    nulls.append(nf)
    res = pd.DataFrame(rows)
    res['q_bh_4_heatshock_tests'] = bh(res['empirical_p_one_sided'].to_numpy(float))
    res.to_csv(OUT / 'heatshock_matched_null_summary.tsv', sep='\t', index=False)
    pd.concat(nulls, ignore_index=True).to_csv(OUT / 'heatshock_matched_null_distributions.tsv.gz', sep='\t', index=False, compression='gzip')
    summary = {
        'seed': SEED,
        'n_permutations_each': N_PERM,
        'n_tests': int(len(res)),
        'all_q_lt_0_05': bool((res['q_bh_4_heatshock_tests'] < 0.05).all()),
        'matching': '|shared heat response| ventile + donor-response SD quintile; both source-only',
        'results': res.to_dict(orient='records'),
    }
    (OUT / 'SUMMARY.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print('HEATSHOCK_MATCHED_NULLS_SUCCESS')
    print(res.round(5).to_string(index=False))


if __name__ == '__main__':
    main()
