from __future__ import annotations

from pathlib import Path
import io
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'publication_controls_v1' / 'future_time_matched_nulls_complete'
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261021
N_PERM = 10000
K = 100
rng = np.random.default_rng(SEED)


def quantile_bins(x: np.ndarray, q: int) -> np.ndarray:
    s = pd.Series(np.asarray(x, float))
    return pd.qcut(s.rank(method='average'), q=q, labels=False, duplicates='drop').to_numpy(int)


def matched_null(candidate_idx: np.ndarray, truth_idx: np.ndarray, features: list[tuple[np.ndarray, int]], n_perm: int = N_PERM) -> np.ndarray:
    bins = [quantile_bins(x, q) for x, q in features]
    strata = list(zip(*bins))
    counts: dict[tuple[int, ...], int] = {}
    for i in candidate_idx:
        s = strata[int(i)]
        counts[s] = counts.get(s, 0) + 1
    pools = {s: np.array([i for i, z in enumerate(strata) if z == s], dtype=int) for s in counts}
    truth = set(map(int, truth_idx))
    null = np.empty(n_perm, int)
    for b in range(n_perm):
        selected: list[int] = []
        for s, n in counts.items():
            pool = pools[s]
            if len(pool) < n:
                # Should be rare because the candidate set itself is contained in the pool;
                # keep behavior explicit rather than silently changing strata.
                raise RuntimeError(f'stratum {s}: need {n}, have {len(pool)}')
            selected.extend(rng.choice(pool, size=n, replace=False).tolist())
        null[b] = len(set(selected) & truth)
    return null


def summarize(dataset: str, context: str, n_genes: int, observed: int, null: np.ndarray,
              matching: str, source: str) -> dict:
    mean = float(np.mean(null))
    return {
        'dataset': dataset,
        'context': context,
        'n_genes': int(n_genes),
        'k': K,
        'observed_overlap': int(observed),
        'null_mean': mean,
        'null_sd': float(np.std(null, ddof=1)),
        'null_q95': float(np.quantile(null, 0.95)),
        'null_q99': float(np.quantile(null, 0.99)),
        'matched_enrichment': float(observed / mean) if mean > 0 else np.inf,
        'empirical_p_one_sided': float((1 + np.sum(null >= observed)) / (len(null) + 1)),
        'matching_variables': matching,
        'source': source,
    }


def read_soft(path: Path) -> pd.DataFrame:
    lines: list[str] = []
    inside = False
    with path.open('r', encoding='utf-8', errors='strict') as f:
        for line in f:
            if line.startswith('!sample_table_begin'):
                inside = True
                continue
            if line.startswith('!sample_table_end'):
                break
            if inside:
                lines.append(line)
    d = pd.read_csv(io.StringIO(''.join(lines)), sep='\t', usecols=['ID_REF', 'VALUE'])
    gsm = path.stem.split('.')[0]
    d = d.rename(columns={'ID_REF': 'probe', 'VALUE': gsm})
    d['probe'] = d['probe'].astype(str)
    d[gsm] = pd.to_numeric(d[gsm], errors='raise')
    return d


def sulfur_row() -> tuple[dict, pd.DataFrame]:
    base = ROOT / 'data' / 'prospective_task' / 'sulfur'
    src = base / 'source_0_3h'
    blind = base / 'task_12h_blind'
    zero = ['GSM744800', 'GSM744801']
    h3 = ['GSM744802', 'GSM744803']

    z = None
    for gsm in zero + h3:
        q = read_soft(src / f'{gsm}.soft.txt')
        z = q if z is None else z.merge(q, on='probe', validate='one_to_one')
    for c in zero + h3:
        v = z[c].to_numpy(float)
        if np.any(v <= 0) or np.any(~np.isfinite(v)):
            raise RuntimeError(f'invalid sulfur source values: {c}')
        z[c] = np.log2(v)

    fr = pd.read_csv(blind / 'FROZEN_12H_PREDICTION.tsv.gz', sep='\t')
    # Target table was not persisted as a single file by the original scorer, so reconstruct it from the two revealed 12h samples.
    tgt_dir = base / 'target_12h_samples'
    target_ids = ['GSM744804', 'GSM744805']
    for gsm in target_ids:
        z = z.merge(read_soft(tgt_dir / f'{gsm}.soft.txt'), on='probe', validate='one_to_one')
    for c in target_ids:
        v = z[c].to_numpy(float)
        if np.any(v <= 0) or np.any(~np.isfinite(v)):
            raise RuntimeError(f'invalid sulfur target values: {c}')
        z[c] = np.log2(v)

    z = fr.merge(z, on='probe', validate='one_to_one')
    baseline = z[zero].mean(axis=1).to_numpy(float)
    h3_resp = np.column_stack([z[c].to_numpy(float) - baseline for c in h3])
    source_dispersion = np.std(h3_resp, axis=1)
    target = z[target_ids].mean(axis=1).to_numpy(float) - baseline
    pred = z['same_state_3h'].to_numpy(float)

    cand = np.argsort(-np.abs(pred), kind='mergesort')[:K]
    truth = np.argsort(-np.abs(target), kind='mergesort')[:K]
    observed = len(np.intersect1d(cand, truth))
    null = matched_null(cand, truth, [(baseline, 5), (source_dispersion, 5)])
    row = summarize(
        'sulfur_3h_to_12h', 'bulk', len(z), observed, null,
        '0h mean expression quintile + 3h source-replicate response-SD quintile',
        'scripts/156_future_time_remaining_matched_nulls.py'
    )
    return row, pd.DataFrame({'dataset': row['dataset'], 'context': row['context'], 'null_overlap': null})


def lateral_root_rows() -> tuple[list[dict], list[pd.DataFrame]]:
    base = ROOT / 'data' / 'external_lateral_root' / 'GSE158761' / 'plantvc_v1_blind'
    audit = pd.read_csv(base / 'freeze_context_audit.tsv', sep='\t')
    rows: list[dict] = []
    null_frames: list[pd.DataFrame] = []
    for r in audit.itertuples(index=False):
        ev = pd.read_csv(base / f'evidence_{r.slug}.tsv.gz', sep='\t')
        tg = pd.read_csv(base / f'target_{r.slug}.tsv.gz', sep='\t')
        rk = pd.read_csv(base / 'rankings' / r.slug / 'ranking.tsv.gz', sep='\t')
        tech = ev['technical_universe'].astype(bool).to_numpy()
        genes = ev.loc[tech, 'gene'].astype(str).to_numpy()
        if ev['gene'].astype(str).tolist() != tg['gene'].astype(str).tolist():
            raise RuntimeError(f'lateral-root gene mismatch: {r.slug}')
        rk2 = rk.set_index(rk['gene'].astype(str)).loc[genes]
        if not np.array_equal(rk2['gene'].astype(str).to_numpy(), genes):
            raise RuntimeError(f'lateral-root ranking mismatch: {r.slug}')

        pred = rk2['predicted_response'].to_numpy(float)
        truth = tg.loc[tech, 'r20'].to_numpy(float)
        r8 = ev.loc[tech, 'r8'].to_numpy(float)
        rep_a = ev.loc[tech, 'r8_repA'].to_numpy(float)
        rep_b = ev.loc[tech, 'r8_repB'].to_numpy(float)
        rep_disagreement = np.abs(rep_a - rep_b)

        cand = np.where(rk2['priority_rank'].to_numpy(int) <= K)[0]
        truth_idx = np.argsort(-np.abs(truth), kind='mergesort')[:K]
        observed = len(np.intersect1d(cand, truth_idx))
        # Conservative null: preserve source-response magnitude very finely, plus replicate disagreement.
        # This asks whether the exact frozen top genes persist to 20h beyond merely being strong 8h responders.
        null = matched_null(cand, truth_idx, [(np.abs(r8), 20), (rep_disagreement, 5)])
        row = summarize(
            'GSE158761_lateral_root_8h_to_20h', str(r.celltype), len(truth), observed, null,
            '|8h same-state response| ventile + 8h replicate-disagreement quintile',
            'scripts/156_future_time_remaining_matched_nulls.py'
        )
        rows.append(row)
        null_frames.append(pd.DataFrame({'dataset': row['dataset'], 'context': row['context'], 'null_overlap': null}))
    return rows, null_frames


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n, float)
    out[order] = q
    return out


def main() -> None:
    existing_path = ROOT / 'data' / 'publication_controls_v1' / 'future_time_matched_nulls' / 'future_time_matched_null_summary.tsv'
    existing = pd.read_csv(existing_path, sep='\t')
    existing = existing.drop(columns=[c for c in existing.columns if c.startswith('q_bh_')], errors='ignore')
    existing['matching_variables'] = 'see original script 152; baseline expression + source temporal variability'
    existing['source'] = 'scripts/152_future_time_matched_empirical_nulls.py'

    srow, snull = sulfur_row()
    lr_rows, lr_nulls = lateral_root_rows()
    added = pd.DataFrame([srow] + lr_rows)
    complete = pd.concat([existing, added], ignore_index=True, sort=False)
    complete['q_bh_all_future_time_tests'] = bh(complete['empirical_p_one_sided'].to_numpy(float))
    complete.to_csv(OUT / 'future_time_matched_null_complete_summary.tsv', sep='\t', index=False)

    all_nulls = [snull] + lr_nulls
    pd.concat(all_nulls, ignore_index=True).to_csv(OUT / 'new_null_distributions.tsv.gz', sep='\t', index=False, compression='gzip')

    summary = {
        'seed': SEED,
        'n_permutations_each_new_test': N_PERM,
        'n_existing_tests_reused': int(len(existing)),
        'n_new_tests': int(len(added)),
        'n_total_future_time_tests': int(len(complete)),
        'n_q_lt_0_05': int((complete['q_bh_all_future_time_tests'] < 0.05).sum()),
        'new_results': added.to_dict(orient='records'),
    }
    (OUT / 'SUMMARY.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print('FUTURE_TIME_REMAINING_MATCHED_NULLS_SUCCESS')
    print(json.dumps({k: v for k, v in summary.items() if k != 'new_results'}, ensure_ascii=False, indent=2))
    print('\nNEW TESTS')
    print(complete[complete['source'].str.contains('156_')][['dataset','context','observed_overlap','null_mean','null_q99','matched_enrichment','empirical_p_one_sided','q_bh_all_future_time_tests']].round(5).to_string(index=False))


if __name__ == '__main__':
    main()
