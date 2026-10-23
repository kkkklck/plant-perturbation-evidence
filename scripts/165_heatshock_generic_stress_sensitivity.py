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
OUT = ROOT / 'data' / 'publication_controls_v2' / 'heatshock_generic_stress_sensitivity'
OUT.mkdir(parents=True, exist_ok=True)

CK = [f'Col0-CK{i}' for i in range(1, 7)]
HS = [f'Col0-HS{i}' for i in range(1, 4)]
K = 100


def norm_gene(x: str) -> str:
    return str(x).strip().strip('"').upper()


def read_target() -> pd.Series:
    counts = pd.read_csv(TARGET, sep='\t', compression='gzip')
    req = ['Geneid'] + CK + HS
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


def read_program(path: Path, value_col: str) -> pd.Series:
    d = pd.read_csv(path, sep='\t')
    d['gene'] = d['gene'].map(norm_gene)
    d[value_col] = pd.to_numeric(d[value_col], errors='coerce')
    d = d[np.isfinite(d[value_col])].drop_duplicates('gene', keep='first')
    return d.set_index('gene')[value_col]


def ranks_abs(v: np.ndarray) -> np.ndarray:
    order = np.argsort(-np.abs(v), kind='mergesort')
    r = np.empty(len(v), dtype=int)
    r[order] = np.arange(1, len(v) + 1)
    return r


def score(program: np.ndarray, truth: np.ndarray, k: int = K) -> dict:
    rho = float(spearmanr(program, truth).statistic)
    pr = ranks_abs(program)
    tr = ranks_abs(truth)
    pi = np.where(pr <= min(k, len(program)))[0]
    ti = np.where(tr <= min(k, len(truth)))[0]
    ov = np.intersect1d(pi, ti)
    return {
        'n_genes': int(len(program)),
        'spearman_full_vector': rho,
        'top100_overlap': int(len(ov)),
        'top100_direction_agreement': float(np.mean(np.sign(program[ov]) == np.sign(truth[ov]))) if len(ov) else np.nan,
    }


target = read_target()
heat = read_program(TRUE_PRED, 'shared_prediction')

programs: dict[str, pd.Series] = {}
salt = pd.read_csv(ROOT / 'data' / 'external_salt' / 'GSE153103' / 'plantvc_v1_1_12h_blind' / 'FROZEN_PRETARGET_RANKING.tsv.gz', sep='\t')
salt['gene'] = salt['gene'].map(norm_gene)
salt = salt[salt['technical_universe'].astype(bool)].drop_duplicates('gene')
programs['salt_6h'] = salt.set_index('gene')['r6'].astype(float)
programs['salt_3h'] = salt.set_index('gene')['r3'].astype(float)
for stim in ['Pep1', 'SCOOP10', 'SCREW1', 'flg22', 'Chitin']:
    p = ROOT / 'data' / 'prospective_task' / 'immune_missing_state' / 'blind_v0_1' / f'FROZEN_{stim}_PREDICTION.tsv.gz'
    programs[stim] = read_program(p, 'predicted_response')

# Use a single common universe so every generic-stress definition is source-only and directly comparable.
common = set(target.index) & set(heat.index)
for s in programs.values():
    common &= set(s.index)
genes = np.array(sorted(common), dtype=object)
if len(genes) < 1000:
    raise RuntimeError(f'too few genes in all-program common universe: {len(genes)}')

tv = target.reindex(genes).to_numpy(float)
hv = heat.reindex(genes).to_numpy(float)
wrong_mat = np.column_stack([programs[k].reindex(genes).to_numpy(float) for k in programs])
wrong_names = list(programs)

# Per-wrong-program top candidate membership and percentile scores.
wrong_ranks = np.column_stack([ranks_abs(wrong_mat[:, j]) for j in range(wrong_mat.shape[1])])
wrong_top100 = wrong_ranks <= 100
wrong_top250 = wrong_ranks <= 250
freq100 = wrong_top100.sum(axis=1)
freq250 = wrong_top250.sum(axis=1)
mean_abs_percentile = np.mean(1.0 - (wrong_ranks - 1) / max(len(genes) - 1, 1), axis=1)

# Data-driven generic-stress sets, all defined without heat target outcomes.
union_top100 = freq100 >= 1
recurrent_top100 = freq100 >= 2
recurrent_top250 = freq250 >= 2
panstress_top100_idx = np.argsort(-mean_abs_percentile, kind='mergesort')[:100]
panstress_top250_idx = np.argsort(-mean_abs_percentile, kind='mergesort')[:250]
panstress_top100 = np.zeros(len(genes), dtype=bool); panstress_top100[panstress_top100_idx] = True
panstress_top250 = np.zeros(len(genes), dtype=bool); panstress_top250[panstress_top250_idx] = True

def original_overlap_indices(program: np.ndarray, truth: np.ndarray, k: int = 100) -> np.ndarray:
    pi = np.argsort(-np.abs(program), kind='mergesort')[:k]
    ti = np.argsort(-np.abs(truth), kind='mergesort')[:k]
    return np.intersect1d(pi, ti)

orig_overlap = original_overlap_indices(hv, tv, K)

sets = {
    'none': np.zeros(len(genes), dtype=bool),
    'exclude_union_wrong_top100': union_top100,
    'exclude_recurrent_wrong_top100_ge2': recurrent_top100,
    'exclude_recurrent_wrong_top250_ge2': recurrent_top250,
    'exclude_panstress_top100_mean_rank': panstress_top100,
    'exclude_panstress_top250_mean_rank': panstress_top250,
}

rows = []
for name, excl in sets.items():
    keep = ~excl
    sc = score(hv[keep], tv[keep], K)
    rows.append({
        'sensitivity': name,
        'excluded_genes': int(excl.sum()),
        'remaining_genes': int(keep.sum()),
        **sc,
    })

sens = pd.DataFrame(rows)
sens.to_csv(OUT / 'heat_transfer_after_generic_stress_exclusion.tsv', sep='\t', index=False)

# How much of the original heat Top100 overlap is generic by each source-only definition?
overlap_rows = []
for name, excl in sets.items():
    if name == 'none':
        continue
    n = int(excl[orig_overlap].sum())
    overlap_rows.append({
        'generic_definition': name,
        'original_heat_overlap_n': int(len(orig_overlap)),
        'overlap_genes_flagged_generic': n,
        'fraction_original_overlap_flagged_generic': float(n / len(orig_overlap)) if len(orig_overlap) else np.nan,
    })
overlap_df = pd.DataFrame(overlap_rows)
overlap_df.to_csv(OUT / 'original_heat_overlap_generic_fraction.tsv', sep='\t', index=False)

# A generic-stress candidate baseline: rank genes only by how consistently they rank strongly across non-heat programs.
# This tests candidate identity, not direction, because source assays differ in scale and biology.
truth_top = set(np.argsort(-np.abs(tv), kind='mergesort')[:K].tolist())
generic_top100 = set(panstress_top100_idx.tolist())
generic_overlap = len(truth_top & generic_top100)
heat_top = set(np.argsort(-np.abs(hv), kind='mergesort')[:K].tolist())
heat_overlap = len(truth_top & heat_top)

baseline = {
    'common_universe_genes': int(len(genes)),
    'wrong_programs': wrong_names,
    'generic_score': 'mean absolute-rank percentile across 7 non-heat programs',
    'generic_top100_target_overlap': int(generic_overlap),
    'true_heat_top100_target_overlap_same_common_universe': int(heat_overlap),
    'true_minus_generic_overlap': int(heat_overlap - generic_overlap),
}
(OUT / 'generic_stress_candidate_baseline.json').write_text(json.dumps(baseline, ensure_ascii=False, indent=2), encoding='utf-8')

# Save gene-level flags for auditability.
gene_df = pd.DataFrame({
    'gene': genes,
    'nonheat_top100_frequency': freq100,
    'nonheat_top250_frequency': freq250,
    'nonheat_mean_abs_rank_percentile': mean_abs_percentile,
    'generic_union_top100': union_top100,
    'generic_recurrent_top100_ge2': recurrent_top100,
    'generic_recurrent_top250_ge2': recurrent_top250,
    'generic_panstress_top100': panstress_top100,
    'generic_panstress_top250': panstress_top250,
    'heat_abs_rank': ranks_abs(hv),
    'target_abs_rank': ranks_abs(tv),
})
gene_df.to_csv(OUT / 'gene_level_generic_stress_flags.tsv.gz', sep='\t', index=False, compression='gzip')

summary = {
    'purpose': 'Test whether cross-study heat Top100 recovery is explained mainly by genes that are broadly strong responders across unrelated non-heat perturbations.',
    'all_generic_definitions_use_source_only_nonheat_programs': True,
    'baseline': baseline,
    'sensitivities': sens.to_dict(orient='records'),
    'original_overlap_generic_fraction': overlap_df.to_dict(orient='records'),
}
(OUT / 'SUMMARY.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

print('HEATSHOCK_GENERIC_STRESS_SENSITIVITY_SUCCESS')
print('\nGENERIC BASELINE')
print(json.dumps(baseline, ensure_ascii=False, indent=2))
print('\nEXCLUSION SENSITIVITIES')
print(sens.round(4).to_string(index=False))
print('\nORIGINAL OVERLAP FLAGGED GENERIC')
print(overlap_df.round(4).to_string(index=False))
