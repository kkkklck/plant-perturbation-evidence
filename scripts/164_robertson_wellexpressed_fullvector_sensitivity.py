from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
ROB = ROOT / 'data' / 'external_rice' / 'Robertson2026'
OUT = ROOT / 'data' / 'publication_controls_v2' / 'robertson_wellexpressed_fullvector'
OUT.mkdir(parents=True, exist_ok=True)

DESIGN = pd.read_csv(ROB / 'sample_design.tsv', sep='\t')
META = pd.read_csv(ROB / 'tissuesystem_pseudobulk_v0' / 'group_metadata.tsv', sep='\t')
EXPR = pd.read_csv(ROB / 'tissuesystem_pseudobulk_v0' / 'tissuesystem_log1p_cpm.tsv.gz', sep='\t')
GENES_ALL = EXPR.pop('gene').astype(str).to_numpy()
X_ALL = EXPR.to_numpy(float)

SAMPLE_BATCH = dict(zip(DESIGN.identity.astype(str), DESIGN.batch.astype(int)))
CONDITIONS = ['Mild.Drought', 'Moderate.Drought', 'Mild.Salinity', 'Moderate.Salinity']
CONTROL = 'Well.Watered'

# Original technical universe used in the paper/candidate analysis.
TECH = ((X_ALL > 0).mean(axis=1) >= 0.25) & (X_ALL.mean(axis=1) >= 0.05)


def safe_rho(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or np.nanstd(a) == 0 or np.nanstd(b) == 0:
        return np.nan
    return float(spearmanr(a, b).statistic)


def controls_for_state(state: str) -> dict[int, np.ndarray]:
    out: dict[int, np.ndarray] = {}
    idx = META.index[(META.condition.eq(CONTROL)) & (META.tissue_system.astype(str).eq(state))].to_numpy()
    for i in idx:
        sample = str(META.loc[i, 'sample'])
        out[SAMPLE_BATCH[sample]] = X_ALL[:, i]
    return out


def treatment_vectors(cond: str, state: str):
    ctrls = controls_for_state(state)
    vecs, batches = [], []
    idx = META.index[(META.condition.eq(cond)) & (META.tissue_system.astype(str).eq(state))].to_numpy()
    for i in idx:
        sample = str(META.loc[i, 'sample'])
        batch = SAMPLE_BATCH[sample]
        if batch not in ctrls:
            continue
        vecs.append(X_ALL[:, i] - ctrls[batch])
        batches.append(batch)
    return vecs, batches, ctrls


rows = []
states = sorted(META.tissue_system.astype(str).unique())
for cond in CONDITIONS:
    for state in states:
        vecs, batches, ctrls = treatment_vectors(cond, state)
        if len(vecs) != 2:
            continue
        v1, v2 = vecs
        b1, b2 = batches
        c1, c2 = ctrls[b1], ctrls[b2]
        baseline = (c1 + c2) / 2.0

        tech_idx = np.where(TECH)[0]
        baseline_tech = baseline[tech_idx]
        positive = baseline_tech[baseline_tech > 0]
        med_positive = float(np.median(positive)) if len(positive) else np.nan
        q50 = float(np.quantile(baseline_tech, 0.50))
        q75 = float(np.quantile(baseline_tech, 0.75))

        masks = {
            'technical_universe': TECH,
            'control_mean_gt_0.5': TECH & (baseline >= 0.5),
            'control_mean_gt_1.0': TECH & (baseline >= 1.0),
            'control_top50pct_within_technical': TECH & (baseline >= q50),
            'control_top25pct_within_technical': TECH & (baseline >= q75),
        }

        for universe, mask in masks.items():
            n = int(mask.sum())
            rho = safe_rho(v1[mask], v2[mask])
            rows.append({
                'condition': cond,
                'tissue_system': state,
                'universe': universe,
                'n_genes': n,
                'replicate_spearman_rho': rho,
                'control_baseline_q50_technical': q50,
                'control_baseline_q75_technical': q75,
                'control_baseline_median_positive_technical': med_positive,
            })

res = pd.DataFrame(rows)
res.to_csv(OUT / 'context_level_fullvector_sensitivity.tsv', sep='\t', index=False)

summary = (res.groupby(['universe', 'condition'], as_index=False)
           .agg(n_contexts=('tissue_system', 'nunique'),
                median_n_genes=('n_genes', 'median'),
                median_rho=('replicate_spearman_rho', 'median'),
                min_rho=('replicate_spearman_rho', 'min'),
                max_rho=('replicate_spearman_rho', 'max')))
summary.to_csv(OUT / 'condition_summary.tsv', sep='\t', index=False)

overall = (res.groupby('universe', as_index=False)
           .agg(n_contexts=('replicate_spearman_rho', 'size'),
                median_n_genes=('n_genes', 'median'),
                median_rho=('replicate_spearman_rho', 'median'),
                mean_rho=('replicate_spearman_rho', 'mean'),
                q25_rho=('replicate_spearman_rho', lambda x: float(np.quantile(x, 0.25))),
                q75_rho=('replicate_spearman_rho', lambda x: float(np.quantile(x, 0.75))),
                n_rho_lt_0_3=('replicate_spearman_rho', lambda x: int((x < 0.3).sum())),
                n_rho_lt_0_5=('replicate_spearman_rho', lambda x: int((x < 0.5).sum()))))
overall.to_csv(OUT / 'overall_summary.tsv', sep='\t', index=False)

# Direct paired changes relative to the original technical universe.
wide = res.pivot_table(index=['condition', 'tissue_system'], columns='universe', values='replicate_spearman_rho').reset_index()
base = wide['technical_universe']
paired_rows = []
for u in [c for c in wide.columns if c not in ['condition', 'tissue_system', 'technical_universe']]:
    delta = wide[u] - base
    paired_rows.append({
        'universe': u,
        'median_delta_rho_vs_technical': float(np.nanmedian(delta)),
        'mean_delta_rho_vs_technical': float(np.nanmean(delta)),
        'n_improved': int(np.sum(delta > 0)),
        'n_worsened': int(np.sum(delta < 0)),
        'n_total': int(np.sum(np.isfinite(delta))),
    })
pd.DataFrame(paired_rows).to_csv(OUT / 'paired_delta_vs_technical.tsv', sep='\t', index=False)

payload = {
    'purpose': 'Sensitivity test for the reviewer concern that low-expression genes alone drive limited full-vector replicate concordance in Robertson rice stress data.',
    'technical_universe_definition': 'nonzero in >=25% pseudobulk groups and global mean log1p(CPM) >=0.05',
    'additional_filters_are_control_only': True,
    'n_contexts': int((res.universe == 'technical_universe').sum()),
    'overall': overall.to_dict(orient='records'),
}
(OUT / 'SUMMARY.json').write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')

print('ROBERTSON_WELLEXPRESSED_FULLVECTOR_SENSITIVITY_SUCCESS')
print('\nOVERALL')
print(overall.round(4).to_string(index=False))
print('\nPAIRED DELTAS')
print(pd.DataFrame(paired_rows).round(4).to_string(index=False))
print('\nBY CONDITION')
print(summary.round(4).to_string(index=False))
