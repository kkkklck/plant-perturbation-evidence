from __future__ import annotations

from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'paper_figures_v2'
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


sources: list[Path] = []
outputs: list[dict] = []


def source(rel: str) -> Path:
    p = ROOT / rel
    if not p.exists():
        raise FileNotFoundError(p)
    sources.append(p)
    return p


def emit(name: str, rows: list[dict]) -> None:
    cols = ['section', 'dataset', 'context', 'model', 'metric', 'value', 'aux', 'note']
    df = pd.DataFrame(rows)
    for c in cols:
        if c not in df.columns:
            df[c] = ''
    df = df[cols]
    path = OUT / name
    df.to_csv(path, sep='\t', index=False)
    outputs.append({'path': str(path.relative_to(ROOT)).replace('\\', '/'), 'rows': int(len(df)), 'sha256': sha256(path)})


def add(rows, section, dataset, context, model, metric, value, aux='', note=''):
    rows.append({
        'section': section,
        'dataset': dataset,
        'context': str(context),
        'model': str(model),
        'metric': str(metric),
        'value': value,
        'aux': aux,
        'note': note,
    })


# =============================================================================
# Figure 1: candidate stability and target redefinition
# =============================================================================
fig1 = []

state = pd.read_csv(source('data/external_rice/Robertson2026/tissuesystem_structure_v0/state_metrics.tsv'), sep='\t')
for _, r in state.iterrows():
    add(fig1, 'rice_full_vector', 'Robertson2026', f"{r['condition']}|{r['tissue_system']}", 'biological_replicates',
        'replicate_response_spearman', float(r['replicate_response_rho']))

cand = pd.read_csv(source('data/external_rice/Robertson2026/candidate_recovery_v0/candidate_overlap_significance.tsv'), sep='\t')
cand = cand[(cand['level'] == 'TissueSystem') & (cand['k'] == 100)]
for _, r in cand.iterrows():
    ctx = f"{r['condition']}|{r['state']}"
    add(fig1, 'rice_top100_random_expectation', 'Robertson2026', ctx, 'Top100_response_candidates', 'observed_overlap', int(r['overlap']))
    add(fig1, 'rice_top100_random_expectation', 'Robertson2026', ctx, 'Top100_response_candidates', 'random_expected_overlap', float(r['expected_overlap']))
    add(fig1, 'rice_top100_random_expectation', 'Robertson2026', ctx, 'Top100_response_candidates', 'random_enrichment', float(r['enrichment']))

matched = pd.read_csv(source('data/external_rice/Robertson2026/candidate_recovery_v0/matched_null_v1/tissuesystem_top100_matched_null.tsv'), sep='\t')
for _, r in matched.iterrows():
    ctx = f"{r['condition']}|{r['tissue_system']}"
    for metric in ['observed_overlap', 'null_mean', 'null_q99', 'matched_enrichment', 'empirical_p_one_sided', 'q_bh_16_tests']:
        add(fig1, 'rice_top100_matched_null', 'Robertson2026', ctx, 'Top100_response_candidates', metric, float(r[metric]),
            note=str(r['matching_variables']))

matched_sum = pd.read_csv(source('data/external_rice/Robertson2026/candidate_recovery_v0/matched_null_v1/condition_summary.tsv'), sep='\t')
for _, r in matched_sum.iterrows():
    for metric in ['median_observed_overlap', 'median_null_mean', 'median_matched_enrichment', 'min_matched_enrichment', 'n_q_lt_0_05', 'max_q']:
        add(fig1, 'rice_condition_summary', 'Robertson2026', r['condition'], 'Top100_response_candidates', metric, float(r[metric]))

emit('fig1_candidate_stability_and_problem_definition.tsv', fig1)


# =============================================================================
# Figure 2: future time
# =============================================================================
fig2 = []

old_future = pd.read_csv(source('data/paper_figures_v1/fig3_future_time.tsv'), sep='\t')
for _, r in old_future.iterrows():
    # Historical GSE46205 rows used a second log2 on already-log2 GCRMA values and are superseded below.
    if str(r['dataset']).startswith('GSE46205'):
        continue
    add(fig2, 'frozen_reveal_primary_results', r['dataset'], r['context'], r['model'], r['metric'], float(r['value']))

corr_route = pd.read_csv(source('data/publication_controls_v1/gse46205_preprocessing_correction/corrected_route_comparison.tsv'), sep='\t')
for _, r in corr_route.iterrows():
    ds = f"GSE46205_salt_{int(r['target_h'])}h_corrected"
    for metric in ['spearman_full_vector','top100_overlap','top100_enrichment','top100_direction']:
        add(fig2, 'post_reveal_preprocessing_corrected_results', ds, r['celltype'], r['model'], metric, float(r[metric]),
            note='Deposited GCRMA VALUE used directly; corrected after target reveal, so these scores are not labeled untouched prospective results.')

base = pd.read_csv(source('data/publication_baseline_suite_v2_corrected/baseline_suite_aggregate.tsv'), sep='\t')
base = base[base['task_family'] == 'future_time']
for _, r in base.iterrows():
    for metric in ['median_rho', 'mean_rho', 'median_top100_overlap', 'mean_top100_overlap', 'median_top100_enrichment', 'median_top100_direction']:
        v = r[metric]
        if pd.notna(v):
            add(fig2, 'strong_baseline_aggregate', r['dataset'], 'aggregate', r['model'], metric, float(v),
                aux=f"n_contexts={int(r['n_contexts'])}; class={r['model_class']}; pretarget={bool(r['pretarget_eligible'])}")

nulls = pd.read_csv(source('data/publication_controls_v1/gse46205_corrected_publication_stats/future_time_matched_null_complete_corrected.tsv'), sep='\t')
for _, r in nulls.iterrows():
    for metric in ['observed_overlap', 'null_mean', 'null_q99', 'matched_enrichment', 'empirical_p_one_sided', 'q_bh_all_future_time_tests_corrected']:
        v = r[metric]
        if pd.notna(v):
            add(fig2, 'matched_empirical_null', r['dataset'], r['context'], 'primary_future_time_prediction', metric, float(v),
                note=str(r['matching_variables']))

paired = pd.read_csv(source('data/publication_stats_v1/context_level/paired_context_inference.tsv'), sep='\t')
paired = paired[~paired['dataset'].astype(str).str.startswith('GSE46205')]
paired_corr = pd.read_csv(source('data/publication_controls_v1/gse46205_corrected_publication_stats/corrected_paired_context_inference.tsv'), sep='\t')
paired = pd.concat([paired, paired_corr], ignore_index=True, sort=False)
for _, r in paired.iterrows():
    for metric in ['mean_paired_difference', 'median_paired_difference', 'wins', 'ties', 'losses', 'exact_signflip_p_one_sided']:
        add(fig2, 'paired_context_inference', r['dataset'], r['metric'], 'same_state_minus_shared', metric, float(r[metric]),
            aux=f"n_contexts={int(r['n_contexts'])}; permutations={int(r['n_permutations'])}")

# Corrected source-sharing values. Historical frozen files are retained in the audit archive but not plotted.
add(fig2, 'source_structure_context', 'GSE46205_salt_3h_test_corrected', '1h_source', 'cross_state_source_structure', 'rank1_energy_fraction', 0.7289499091483366,
    note='Post-reveal preprocessing correction using deposited GCRMA log2 values directly; original task split unchanged.')
add(fig2, 'source_structure_context', 'GSE46205_salt_8h_test_corrected', '3h_source', 'cross_state_source_structure', 'rank1_energy_fraction', 0.744648398399414,
    note='Post-reveal preprocessing correction using deposited GCRMA log2 values directly; rule choice had been fixed before 8h target reveal, but corrected numerical score is post-reveal.')

emit('fig2_future_time.tsv', fig2)


# =============================================================================
# Figure 3: missing state and cross-experiment transfer
# =============================================================================
fig3 = []

old_missing = pd.read_csv(source('data/paper_figures_v1/fig4_missing_state_transfer.tsv'), sep='\t')
# Keep heatshock strict holdout and independent cross-experiment rows here; immune reliability goes to Fig.4.
for _, r in old_missing.iterrows():
    if r['section'] in {'heatshock_missing_state', 'cross_experiment_heatshock'}:
        add(fig3, 'frozen_reveal_primary_results', r['section'], f"{r['stimulus']}|{r['target']}", r['route'], r['metric'], float(r['value']))

heat_base = pd.read_csv(source('data/publication_baselines_v1/heatshock/heatshock_baseline_summary.tsv'), sep='\t')
for _, r in heat_base.iterrows():
    for metric in ['median_rho', 'median_sign', 'median_top100_overlap', 'median_top100_enrichment', 'median_top100_direction']:
        v = r[metric]
        if pd.notna(v):
            add(fig3, 'heatshock_strong_baselines', r['dataset'], 'aggregate', r['model'], metric, float(v),
                aux=f"n_contexts={int(r['n_contexts'])}")

heat_null = pd.read_csv(source('data/publication_controls_v1/heatshock_matched_nulls/heatshock_matched_null_summary.tsv'), sep='\t')
for _, r in heat_null.iterrows():
    for metric in ['observed_overlap', 'null_mean', 'null_q99', 'matched_enrichment', 'empirical_p_one_sided', 'q_bh_4_heatshock_tests']:
        add(fig3, 'heatshock_matched_null', r['dataset'], r['context'], 'shared_heat_program', metric, float(r[metric]),
            note=str(r['matching_variables']))

wrong = pd.read_csv(source('data/publication_controls_v1/heatshock_crossstudy_negative_controls/wrong_perturbation_pairwise.tsv'), sep='\t')
for _, r in wrong.iterrows():
    ctx = r['wrong_program']
    for metric in ['true_heat_rho', 'wrong_program_rho', 'delta_rho_true_minus_wrong', 'true_heat_top100_overlap', 'wrong_program_top100_overlap', 'delta_top100_overlap_true_minus_wrong']:
        add(fig3, 'cross_experiment_wrong_perturbation_controls', 'GSE290176_crossstudy', ctx, 'true_heat_vs_wrong_program', metric, float(r[metric]),
            aux=f"common_genes={int(r['common_universe_genes'])}")

emit('fig3_missing_state_and_cross_experiment.tsv', fig3)


# =============================================================================
# Figure 4: reliability stratification and boundary
# =============================================================================
fig4 = []

# Frozen immune results and source geometry.
for _, r in old_missing.iterrows():
    if r['section'] in {'immune_missing_state', 'immune_source_geometry'}:
        add(fig4, 'immune_frozen_reveal', r['stimulus'], r['target'], r['route'], r['metric'], float(r['value']))

stim = pd.read_csv(source('data/prospective_task/immune_missing_state/publication_stats_v1/stimulus_level_source_vs_target.tsv'), sep='\t')
for _, r in stim.iterrows():
    route = 'higher_reliability_frozen_group' if bool(r['allowed']) else 'lower_reliability_frozen_group'
    for metric in ['rank1_energy', 'median_pairwise_rho', 'effective_state_fraction', 'shared_rms', 'top100_mean_abs_shared_response',
                   'median_target_full_rho', 'median_target_top100_overlap', 'median_target_matched_enrichment']:
        add(fig4, 'stimulus_level_source_vs_target', 'immune_missing_state', r['stimulus'], route, metric, float(r[metric]),
            aux=f"n_frozen_genes={int(r['n_frozen_genes'])}")

boot = pd.read_csv(source('data/prospective_task/immune_missing_state/publication_stats_v1/hierarchical_bootstrap_gate_vs_reject.tsv'), sep='\t')
for _, r in boot.iterrows():
    for metric in ['allowed_stimulus_equalweight_mean', 'rejected_stimulus_equalweight_mean', 'difference_allow_minus_reject',
                   'bootstrap_ci025', 'bootstrap_ci975', 'bootstrap_prob_difference_gt0',
                   'exact_stimulus_label_permutation_p_one_sided']:
        add(fig4, 'hierarchical_stimulus_inference', 'immune_missing_state', r['outcome'], 'higher_minus_lower_reliability', metric, float(r[metric]),
            aux=f"n_stimuli={int(r['n_stimuli'])}; exact_permutations={int(r['exact_permutation_n'])}")

emp = pd.read_csv(source('data/prospective_task/immune_missing_state/publication_stats_v1/top100_matched_empirical_null.tsv'), sep='\t')
for _, r in emp.iterrows():
    route = 'higher_reliability_frozen_group' if bool(r['gene_level_output_allowed']) else 'lower_reliability_frozen_group'
    ctx = f"{r['stimulus']}|{r['target_state']}"
    for metric in ['observed_overlap', 'null_mean', 'null_q975', 'empirical_enrichment_vs_matched_mean', 'empirical_p_one_sided', 'empirical_q_bh_20_targets']:
        add(fig4, 'immune_matched_null', 'immune_missing_state', ctx, route, metric, float(r[metric]), note=str(r['match_variables']))

conf = pd.read_csv(source('data/prospective_task/immune_missing_state/publication_stats_v1/stimulus_level_confounder_correlations.tsv'), sep='\t')
for _, r in conf.iterrows():
    add(fig4, 'source_metric_confounding', 'immune_missing_state', r['source_predictor'], r['target_outcome'], 'spearman_rho', float(r['spearman_rho']))
    add(fig4, 'source_metric_confounding', 'immune_missing_state', r['source_predictor'], r['target_outcome'], 'exact_permutation_p_n5', float(r['exact_two_sided_permutation_p_n5']))
    add(fig4, 'source_metric_confounding', 'immune_missing_state', r['source_predictor'], r['target_outcome'], 'q_bh_within_outcome', float(r['q_bh_within_outcome']))

emit('fig4_reliability_stratification_and_boundary.tsv', fig4)


manifest = {
    'status': 'PAPER_FIGURE_TABLES_V2_TOP_TIER_CLOSURE',
    'output_dir': 'data/paper_figures_v2',
    'main_figure_count': 4,
    'story': [
        'candidate stability and prediction-target redefinition',
        'future time: latest same-state history as robust default',
        'missing state / experiment: cross-state shared evidence transfer',
        'task structure first; source capability controls output granularity and uncertainty boundary',
    ],
    'source_files': [
        {'path': str(p.relative_to(ROOT)).replace('\\', '/'), 'sha256': sha256(p)}
        for p in sorted(set(sources), key=lambda x: str(x))
    ],
    'outputs': outputs,
    'important_boundary': 'Immune source geometry is retained as the historically frozen reliability rule, but the v2 story does not claim geometry is independent of response amplitude.',
}
(OUT / 'FIGURE_TABLE_MANIFEST.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')

print('PAPER_FIGURE_TABLES_V2_SUCCESS')
print(json.dumps(manifest, ensure_ascii=False, indent=2))
