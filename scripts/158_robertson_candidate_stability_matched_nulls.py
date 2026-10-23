from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ROB = ROOT / 'data' / 'external_rice' / 'Robertson2026'
OUT = ROB / 'candidate_recovery_v0' / 'matched_null_v1'
OUT.mkdir(parents=True, exist_ok=True)

DESIGN = pd.read_csv(ROB / 'sample_design.tsv', sep='\t')
META = pd.read_csv(ROB / 'tissuesystem_pseudobulk_v0' / 'group_metadata.tsv', sep='\t')
EXPR = pd.read_csv(ROB / 'tissuesystem_pseudobulk_v0' / 'tissuesystem_log1p_cpm.tsv.gz', sep='\t')
GENES_ALL = EXPR.pop('gene').astype(str).to_numpy()
X_ALL = EXPR.to_numpy(float)

SAMPLE_BATCH = dict(zip(DESIGN.identity.astype(str), DESIGN.batch.astype(int)))
CONDITIONS = ['Mild.Drought', 'Mild.Salinity', 'Moderate.Drought', 'Moderate.Salinity']
CONTROL = 'Well.Watered'
K = 100
B = 10000
SEED = 20261021
rng = np.random.default_rng(SEED)

# Match the original technical universe used for candidate stability.
MASK = ((X_ALL > 0).mean(axis=1) >= 0.25) & (X_ALL.mean(axis=1) >= 0.05)
GENES = GENES_ALL[MASK]
X = X_ALL[MASK]


def qbin(x: np.ndarray, q: int) -> np.ndarray:
    s = pd.Series(np.asarray(x, float))
    return pd.qcut(s.rank(method='average'), q=q, labels=False, duplicates='drop').to_numpy(int)


def bh(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    qq = ranked * n / np.arange(1, n + 1)
    qq = np.minimum.accumulate(qq[::-1])[::-1]
    qq = np.clip(qq, 0, 1)
    out = np.empty(n, float)
    out[order] = qq
    return out


def control_by_batch_for_state(state: str) -> dict[int, np.ndarray]:
    out = {}
    rows = META.index[(META.condition.eq(CONTROL)) & (META.tissue_system.astype(str).eq(state))].to_numpy()
    for i in rows:
        sample = str(META.loc[i, 'sample'])
        batch = SAMPLE_BATCH[sample]
        out[batch] = X[:, i]
    return out


def treatment_vectors(cond: str, state: str) -> tuple[list[np.ndarray], list[int]]:
    ctrls = control_by_batch_for_state(state)
    vecs, batches = [], []
    rows = META.index[(META.condition.eq(cond)) & (META.tissue_system.astype(str).eq(state))].to_numpy()
    for i in rows:
        sample = str(META.loc[i, 'sample'])
        batch = SAMPLE_BATCH[sample]
        if batch not in ctrls:
            continue
        vecs.append(X[:, i] - ctrls[batch])
        batches.append(batch)
    return vecs, batches


def matched_random_set(counts: dict[int, int], pools: dict[int, np.ndarray], rng: np.random.Generator) -> set[int]:
    chosen = []
    for s, c in counts.items():
        pool = pools[s]
        if c > len(pool):
            raise RuntimeError(f'stratum {s}: need {c}, have {len(pool)}')
        chosen.extend(rng.choice(pool, size=c, replace=False).tolist())
    return set(chosen)


rows = []
null_frames = []
states = sorted(META.tissue_system.astype(str).unique())
for cond in CONDITIONS:
    for state in states:
        vecs, batches = treatment_vectors(cond, state)
        if len(vecs) != 2:
            continue
        v1, v2 = vecs
        b1, b2 = batches

        # Source-independent matching variables: baseline expression and control-batch instability.
        ctrls = control_by_batch_for_state(state)
        if b1 not in ctrls or b2 not in ctrls:
            raise RuntimeError(f'missing controls for {cond} / {state}')
        c1, c2 = ctrls[b1], ctrls[b2]
        baseline_mean = (c1 + c2) / 2.0
        baseline_instability = np.abs(c1 - c2)
        expr_bin = qbin(baseline_mean, 5)
        instab_bin = qbin(baseline_instability, 5)
        strata = expr_bin * 5 + instab_bin

        top1 = np.argsort(-np.abs(v1), kind='mergesort')[:K]
        top2 = np.argsort(-np.abs(v2), kind='mergesort')[:K]
        observed = len(set(top1.tolist()) & set(top2.tolist()))

        def counts_for(top: np.ndarray) -> dict[int, int]:
            u, c = np.unique(strata[top], return_counts=True)
            return {int(a): int(b) for a, b in zip(u, c)}

        cnt1 = counts_for(top1)
        cnt2 = counts_for(top2)
        needed = sorted(set(cnt1) | set(cnt2))
        pools = {int(s): np.where(strata == s)[0] for s in needed}

        null = np.empty(B, dtype=int)
        for i in range(B):
            a = matched_random_set(cnt1, pools, rng)
            b = matched_random_set(cnt2, pools, rng)
            null[i] = len(a & b)

        p = (1 + int(np.sum(null >= observed))) / (B + 1)
        rows.append({
            'condition': cond,
            'tissue_system': state,
            'n_genes': int(len(GENES)),
            'k': K,
            'observed_overlap': int(observed),
            'null_mean': float(null.mean()),
            'null_sd': float(null.std(ddof=1)),
            'null_q95': float(np.quantile(null, 0.95)),
            'null_q99': float(np.quantile(null, 0.99)),
            'matched_enrichment': float(observed / null.mean()) if null.mean() > 0 else np.inf,
            'empirical_p_one_sided': float(p),
            'matching_variables': 'control mean expression quintile + control batch-instability quintile',
        })
        null_frames.append(pd.DataFrame({
            'condition': cond,
            'tissue_system': state,
            'null_overlap': null,
        }))

res = pd.DataFrame(rows)
res['q_bh_16_tests'] = bh(res['empirical_p_one_sided'].to_numpy(float))
res.to_csv(OUT / 'tissuesystem_top100_matched_null.tsv', sep='\t', index=False)
pd.concat(null_frames, ignore_index=True).to_csv(OUT / 'tissuesystem_top100_matched_null_distributions.tsv.gz', sep='\t', index=False, compression='gzip')

agg = (res.groupby('condition', as_index=False)
       .agg(n_tissue_systems=('tissue_system', 'nunique'),
            median_observed_overlap=('observed_overlap', 'median'),
            median_null_mean=('null_mean', 'median'),
            median_matched_enrichment=('matched_enrichment', 'median'),
            min_matched_enrichment=('matched_enrichment', 'min'),
            n_q_lt_0_05=('q_bh_16_tests', lambda x: int((x < 0.05).sum())),
            max_q=('q_bh_16_tests', 'max')))
agg.to_csv(OUT / 'condition_summary.tsv', sep='\t', index=False)

summary = {
    'seed': SEED,
    'n_permutations_each': B,
    'n_tests': int(len(res)),
    'all_16_q_lt_0_05': bool((res['q_bh_16_tests'] < 0.05).all()),
    'matching_variables': ['control mean expression', 'absolute control difference between matched batches'],
    'purpose': 'Test whether cross-replicate Top100 stability exceeds what is expected from genes with similar baseline expression and baseline stability.',
    'condition_summary': agg.to_dict(orient='records'),
}
(OUT / 'SUMMARY.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

print('ROBERTSON_MATCHED_NULL_SUCCESS')
print(res.round(5).to_string(index=False))
print('\nCONDITION SUMMARY')
print(agg.round(5).to_string(index=False))
