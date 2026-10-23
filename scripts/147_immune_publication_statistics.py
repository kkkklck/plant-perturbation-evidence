from __future__ import annotations

from pathlib import Path
from itertools import combinations, permutations
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data' / 'prospective_task' / 'immune_missing_state'
SRC = BASE / 'source_only'
TGT = BASE / 'target_only'
BLIND = BASE / 'blind_v0_1'
OUT = BASE / 'publication_stats_v1'
OUT.mkdir(parents=True, exist_ok=True)

STIMS = ['flg22', 'Pep1', 'Chitin', 'SCOOP10', 'SCREW1']
SOURCE_STATES = [
    'Rootcap_Columella',
    'Meristem_Epidermis','Meristem_Endodermis','Meristem_Stele',
    'Elongation_Epidermis','Elongation_Cortex','Elongation_Stele',
    'Maturation_Epidermis','Maturation_Cortex','Maturation_Endodermis',
]
TARGET_STATES = ['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']
SEED = 20261021
B_MATCHED = 10000
B_BOOT = 50000
PRIMARY_K = 100


def bh_qvalues(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n, float)
    out[order] = q
    return out


def quantile_bin(x: pd.Series, q: int = 5) -> pd.Series:
    # Rank first to make qcut deterministic under ties.
    r = x.rank(method='first', pct=True)
    b = np.floor(np.minimum(r.to_numpy(float), 0.999999) * q).astype(int)
    return pd.Series(b, index=x.index)


def read_source_state(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, sep='\t', usecols=['gene','avg_logFC','pct.1','pct.2','p_val_adj'])
    d['gene'] = d['gene'].astype(str)
    for c in ['avg_logFC','pct.1','pct.2','p_val_adj']:
        d[c] = pd.to_numeric(d[c], errors='coerce')
    d = d.drop_duplicates('gene', keep='first').set_index('gene')
    return d


def read_target(path: Path) -> pd.Series:
    d = pd.read_csv(path, sep='\t', usecols=['gene','avg_logFC'])
    d['gene'] = d['gene'].astype(str)
    d['avg_logFC'] = pd.to_numeric(d['avg_logFC'], errors='coerce')
    return d.drop_duplicates('gene', keep='first').set_index('gene')['avg_logFC']


def empirical_matched_null(top_idx: np.ndarray, truth_idx: set[int], strata: np.ndarray,
                           n_genes: int, rng: np.random.Generator, B: int) -> tuple[np.ndarray, dict]:
    # Preserve the observed candidate set's composition across target-independent source-derived strata.
    unique, counts = np.unique(strata[top_idx], return_counts=True)
    pools = {s: np.where(strata == s)[0] for s in unique}
    null = np.zeros(B, dtype=int)
    for b in range(B):
        chosen = []
        for s, c in zip(unique, counts):
            pool = pools[s]
            if c > len(pool):
                raise RuntimeError(f'stratum {s}: need {c}, have {len(pool)}')
            chosen.extend(rng.choice(pool, size=int(c), replace=False).tolist())
        null[b] = len(set(chosen) & truth_idx)
    obs = len(set(top_idx.tolist()) & truth_idx)
    p = (1 + int(np.sum(null >= obs))) / (B + 1)
    info = {
        'observed_overlap': int(obs),
        'null_mean': float(np.mean(null)),
        'null_median': float(np.median(null)),
        'null_q025': float(np.quantile(null, 0.025)),
        'null_q975': float(np.quantile(null, 0.975)),
        'empirical_enrichment_vs_matched_mean': float(obs / np.mean(null)) if np.mean(null) > 0 else np.inf,
        'empirical_p_one_sided': float(p),
    }
    return null, info


def exact_spearman_permutation(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan, np.nan
    obs = float(spearmanr(x, y).statistic)
    vals = []
    for perm in permutations(y.tolist()):
        vals.append(float(spearmanr(x, np.asarray(perm, float)).statistic))
    vals = np.asarray(vals)
    p = float(np.mean(np.abs(vals) >= abs(obs) - 1e-12))
    return obs, p


def hierarchical_bootstrap(target_df: pd.DataFrame, allowed_map: dict[str, bool], outcome: str,
                           rng: np.random.Generator, B: int) -> dict:
    allowed_stims = [s for s in STIMS if allowed_map[s]]
    rejected_stims = [s for s in STIMS if not allowed_map[s]]
    by_stim = {s: target_df[target_df['stimulus'] == s][outcome].dropna().to_numpy(float) for s in STIMS}
    diffs = np.empty(B, float)
    for b in range(B):
        group_vals = {}
        for label, stims in [('allow', allowed_stims), ('reject', rejected_stims)]:
            sampled_stims = rng.choice(stims, size=len(stims), replace=True)
            stim_means = []
            for s in sampled_stims:
                arr = by_stim[s]
                sample = rng.choice(arr, size=len(arr), replace=True)
                stim_means.append(float(np.mean(sample)))
            group_vals[label] = float(np.mean(stim_means))
        diffs[b] = group_vals['allow'] - group_vals['reject']

    # Stimulus is the gate unit. Exact label permutation therefore operates on 5 stimulus-level summaries.
    stim_summary = target_df.groupby('stimulus', as_index=False)[outcome].mean().set_index('stimulus')[outcome]
    obs_allow = float(stim_summary[allowed_stims].mean())
    obs_reject = float(stim_summary[rejected_stims].mean())
    obs_diff = obs_allow - obs_reject
    perm_diffs = []
    for allow_combo in combinations(STIMS, len(allowed_stims)):
        allow_combo = list(allow_combo)
        reject_combo = [s for s in STIMS if s not in allow_combo]
        perm_diffs.append(float(stim_summary[allow_combo].mean() - stim_summary[reject_combo].mean()))
    perm_diffs = np.asarray(perm_diffs)
    p_one = float(np.mean(perm_diffs >= obs_diff - 1e-12))

    return {
        'outcome': outcome,
        'n_stimuli': len(STIMS),
        'n_allowed_stimuli': len(allowed_stims),
        'n_rejected_stimuli': len(rejected_stims),
        'allowed_stimulus_equalweight_mean': obs_allow,
        'rejected_stimulus_equalweight_mean': obs_reject,
        'difference_allow_minus_reject': obs_diff,
        'bootstrap_ci025': float(np.quantile(diffs, 0.025)),
        'bootstrap_ci975': float(np.quantile(diffs, 0.975)),
        'bootstrap_prob_difference_gt0': float(np.mean(diffs > 0)),
        'exact_stimulus_label_permutation_p_one_sided': p_one,
        'exact_permutation_n': int(len(perm_diffs)),
    }


man = json.loads((BLIND / 'FREEZE_MANIFEST.json').read_text(encoding='utf-8'))
vector = pd.read_csv(BLIND / 'blind_vector_results.tsv', sep='\t')
topk = pd.read_csv(BLIND / 'blind_topk_results.tsv', sep='\t')
top100 = topk[topk['k'] == PRIMARY_K].copy()

allowed_map = {s: bool(man['routes'][s]['gene_level_output_allowed']) for s in STIMS}
rng = np.random.default_rng(SEED)

# ---------- Source-only features and confounder candidates ----------
source_rows = []
source_feature_frames: dict[str, pd.DataFrame] = {}
for stim in STIMS:
    pred = pd.read_csv(BLIND / f'FROZEN_{stim}_PREDICTION.tsv.gz', sep='\t')
    genes = pred['gene'].astype(str).tolist()
    state_tables = {s: read_source_state(SRC / stim / f'{s}.tsv.gz') for s in SOURCE_STATES}
    M = np.column_stack([state_tables[s].reindex(genes)['avg_logFC'].to_numpy(float) for s in SOURCE_STATES])
    pct2 = np.column_stack([state_tables[s].reindex(genes)['pct.2'].to_numpy(float) for s in SOURCE_STATES])
    padj = np.column_stack([state_tables[s].reindex(genes)['p_val_adj'].to_numpy(float) for s in SOURCE_STATES])

    feat = pd.DataFrame({
        'gene': genes,
        'median_control_detection': np.nanmedian(pct2, axis=1),
        'response_dispersion_sd': np.nanstd(M, axis=1),
        'mean_abs_source_response': np.nanmean(np.abs(M), axis=1),
        'shared_response': np.nanmean(M, axis=1),
        'shared_abs_response': np.abs(np.nanmean(M, axis=1)),
        'source_significant_fraction': np.nanmean((padj < 0.05) & (np.abs(M) >= 0.25), axis=1),
    })
    feat['detect_bin'] = quantile_bin(feat['median_control_detection'], 5).to_numpy(int)
    feat['disp_bin'] = quantile_bin(feat['response_dispersion_sd'], 5).to_numpy(int)
    feat['matched_stratum'] = feat['detect_bin'] * 5 + feat['disp_bin']
    source_feature_frames[stim] = feat

    geom = man['routes'][stim]['source_geometry']
    gate_margin = min(
        geom['rank1_energy'] / man['shared_thresholds']['rank1_energy_min'],
        geom['median_pairwise_rho'] / man['shared_thresholds']['median_pairwise_rho_min'],
        geom['effective_state_fraction'] / man['shared_thresholds']['effective_state_fraction_min'],
    )
    source_rows.append({
        'stimulus': stim,
        'allowed': allowed_map[stim],
        'rank1_energy': geom['rank1_energy'],
        'median_pairwise_rho': geom['median_pairwise_rho'],
        'effective_state_fraction': geom['effective_state_fraction'],
        'gate_margin_min_normalized': gate_margin,
        'shared_rms': float(np.sqrt(np.mean(feat['shared_response'].to_numpy(float) ** 2))),
        'median_abs_shared_response': float(np.median(feat['shared_abs_response'])),
        'top100_mean_abs_shared_response': float(np.mean(np.sort(feat['shared_abs_response'].to_numpy(float))[-100:])),
        'mean_source_significant_fraction': float(np.mean(feat['source_significant_fraction'])),
        'median_control_detection': float(np.median(feat['median_control_detection'])),
        'n_frozen_genes': int(len(feat)),
    })
source_metrics = pd.DataFrame(source_rows)

# ---------- Matched empirical null for top100 ----------
emp_rows = []
null_dir = OUT / 'matched_null_top100'
null_dir.mkdir(exist_ok=True)
for stim in STIMS:
    pred = pd.read_csv(BLIND / f'FROZEN_{stim}_PREDICTION.tsv.gz', sep='\t')
    genes = pred['gene'].astype(str).tolist()
    N = len(genes)
    top_idx = np.where(pred['priority_rank'].to_numpy(int) <= PRIMARY_K)[0]
    feat = source_feature_frames[stim].set_index('gene').reindex(genes)
    strata = feat['matched_stratum'].to_numpy(int)
    for target in TARGET_STATES:
        ts = read_target(TGT / stim / f'{target}.tsv.gz')
        tv = np.array([float(ts[g]) if g in ts.index else 0.0 for g in genes], dtype=float)
        truth_order = np.lexsort((np.array(genes), -np.abs(tv)))
        truth_idx = set(truth_order[:PRIMARY_K].tolist())
        null, info = empirical_matched_null(top_idx, truth_idx, strata, N, rng, B_MATCHED)
        np.savez_compressed(null_dir / f'{stim}__{target}__k100.npz', overlap=null)
        emp_rows.append({
            'stimulus': stim,
            'target_state': target,
            'gene_level_output_allowed': allowed_map[stim],
            'k': PRIMARY_K,
            'match_variables': 'source_control_detection_quintile + source_response_dispersion_quintile',
            'n_genes': N,
            **info,
        })
emp = pd.DataFrame(emp_rows)
emp['empirical_q_bh_20_targets'] = bh_qvalues(emp['empirical_p_one_sided'].to_numpy(float))
emp.to_csv(OUT / 'top100_matched_empirical_null.tsv', sep='\t', index=False)

# ---------- Target-level joined table ----------
joined = vector.merge(top100[['stimulus','target_state','overlap','enrichment','direction_agreement_overlap']],
                      on=['stimulus','target_state'], how='left')
joined = joined.merge(emp[['stimulus','target_state','null_mean','null_q025','null_q975',
                           'empirical_enrichment_vs_matched_mean','empirical_p_one_sided','empirical_q_bh_20_targets']],
                      on=['stimulus','target_state'], how='left')
joined.to_csv(OUT / 'target_level_publication_stats.tsv', sep='\t', index=False)

# ---------- Hierarchical bootstrap and exact stimulus-level label permutation ----------
boot_input = joined.rename(columns={
    'spearman_frozen_universe': 'full_rho',
    'overlap': 'top100_overlap',
    'enrichment': 'top100_hypergeom_enrichment',
    'empirical_enrichment_vs_matched_mean': 'top100_matched_enrichment',
})
boot_results = []
for outcome in ['full_rho','top100_overlap','top100_matched_enrichment']:
    boot_results.append(hierarchical_bootstrap(boot_input, allowed_map, outcome, rng, B_BOOT))
pd.DataFrame(boot_results).to_csv(OUT / 'hierarchical_bootstrap_gate_vs_reject.tsv', sep='\t', index=False)

# ---------- Stimulus-level confounder analysis ----------
stim_perf = joined.groupby('stimulus', as_index=False).agg(
    median_target_full_rho=('spearman_frozen_universe','median'),
    median_target_top100_overlap=('overlap','median'),
    median_target_matched_enrichment=('empirical_enrichment_vs_matched_mean','median'),
    median_target_empirical_p=('empirical_p_one_sided','median'),
)
stim = source_metrics.merge(stim_perf, on='stimulus', how='left')
stim.to_csv(OUT / 'stimulus_level_source_vs_target.tsv', sep='\t', index=False)

predictors = [
    'rank1_energy','median_pairwise_rho','effective_state_fraction','gate_margin_min_normalized',
    'shared_rms','median_abs_shared_response','top100_mean_abs_shared_response',
    'mean_source_significant_fraction','median_control_detection','n_frozen_genes'
]
outcomes = ['median_target_full_rho','median_target_top100_overlap','median_target_matched_enrichment']
corr_rows = []
for xcol in predictors:
    x = stim[xcol].to_numpy(float)
    for ycol in outcomes:
        y = stim[ycol].to_numpy(float)
        rho, p = exact_spearman_permutation(x, y)
        corr_rows.append({
            'source_predictor': xcol,
            'target_outcome': ycol,
            'spearman_rho': rho,
            'exact_two_sided_permutation_p_n5': p,
        })
corr = pd.DataFrame(corr_rows)
corr['q_bh_within_outcome'] = np.nan
for ycol in outcomes:
    m = corr['target_outcome'] == ycol
    vals = corr.loc[m, 'exact_two_sided_permutation_p_n5'].to_numpy(float)
    finite = np.isfinite(vals)
    q = np.full(len(vals), np.nan, float)
    if finite.any():
        q[finite] = bh_qvalues(vals[finite])
    corr.loc[m, 'q_bh_within_outcome'] = q
corr.to_csv(OUT / 'stimulus_level_confounder_correlations.tsv', sep='\t', index=False)

# ---------- Summary ----------
summary = {
    'analysis': 'immune_missing_state_publication_statistics_v1',
    'seed': SEED,
    'matched_null_draws_per_target': B_MATCHED,
    'hierarchical_bootstrap_draws': B_BOOT,
    'independent_gate_units': 5,
    'allowed_stimuli': [s for s in STIMS if allowed_map[s]],
    'rejected_stimuli': [s for s in STIMS if not allowed_map[s]],
    'matched_null_definition': 'source-only control-detection quintile x source-response-dispersion quintile',
    'multiple_testing': 'Benjamini-Hochberg across 20 target-state top100 empirical tests',
    'stimulus_level_inference_note': 'gate is constant within stimulus; effective independent unit is stimulus, not target state',
    'exact_gate_label_permutation_resolution': 0.1,
}
(OUT / 'ANALYSIS_MANIFEST.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

print('IMMUNE_PUBLICATION_STATS_SUCCESS')
print('\nHIERARCHICAL_BOOTSTRAP')
print(pd.DataFrame(boot_results).round(4).to_string(index=False))
print('\nSTIMULUS_LEVEL')
print(stim.round(4).to_string(index=False))
print('\nTOP100_MATCHED_EMPIRICAL')
print(emp[['stimulus','target_state','observed_overlap','null_mean','empirical_enrichment_vs_matched_mean','empirical_p_one_sided','empirical_q_bh_20_targets']].round(4).to_string(index=False))
print('\nCONFOUNDER_CORRELATIONS_TOP')
print(corr.sort_values(['target_outcome','exact_two_sided_permutation_p_n5']).groupby('target_outcome').head(5).round(4).to_string(index=False))
