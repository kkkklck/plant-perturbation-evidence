from __future__ import annotations

from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'publication_baseline_suite_v1'
OUT.mkdir(parents=True, exist_ok=True)
K = 100


def _top100(pred: np.ndarray, truth: np.ndarray) -> tuple[float, float, float]:
    n = len(truth)
    k = min(K, n)
    p = np.argsort(-np.abs(pred), kind='mergesort')[:k]
    t = np.argsort(-np.abs(truth), kind='mergesort')[:k]
    ov = np.intersect1d(p, t)
    expected = k * k / n
    enrich = len(ov) / expected if expected else np.nan
    direction = float(np.mean(np.sign(pred[ov]) == np.sign(truth[ov]))) if len(ov) else np.nan
    return float(len(ov)), float(enrich), direction


def _row(task_family: str, dataset: str, context: str, model: str, model_class: str,
         n_genes: int, rho: float | None, top100_overlap: float | None,
         top100_enrichment: float | None, top100_direction: float | None,
         pretarget_eligible: bool, source_file: str, note: str = '') -> dict:
    return {
        'task_family': task_family,
        'dataset': dataset,
        'context': context,
        'model': model,
        'model_class': model_class,
        'pretarget_eligible': bool(pretarget_eligible),
        'n_genes': int(n_genes),
        'spearman_full_vector': rho,
        'top100_overlap': top100_overlap,
        'top100_enrichment': top100_enrichment,
        'top100_direction_agreement': top100_direction,
        'source_file': source_file,
        'note': note,
    }


def add_random_expectation(rows: list[dict], task_family: str, dataset: str, context: str,
                           n_genes: int, source_file: str) -> None:
    k = min(K, n_genes)
    expected = k * k / n_genes
    rows.append(_row(
        task_family, dataset, context, 'random_rank_expectation', 'random_expectation',
        n_genes, 0.0, float(expected), 1.0, np.nan, True, source_file,
        'Analytic random-ranking expectation; rho shown as null expectation 0, not one arbitrary random draw.'
    ))


def load_future_time(rows: list[dict]) -> None:
    p = ROOT / 'data' / 'publication_controls_v1' / 'future_time_strong_baselines' / 'future_time_baseline_results.tsv'
    d = pd.read_csv(p, sep='\t')
    for r in d.itertuples(index=False):
        rows.append(_row(
            'future_time', str(r.dataset), str(r.context), str(r.model), 'strong_baseline',
            int(r.n_genes), float(r.spearman_full_vector), float(r.top100_overlap),
            float(r.top100_enrichment), float(r.top100_direction_agreement), True,
            str(p.relative_to(ROOT)),
            'Uses only history available before the target time.'
        ))
    for (ds, ctx), q in d.groupby(['dataset', 'context'], sort=False):
        add_random_expectation(rows, 'future_time', str(ds), str(ctx), int(q.iloc[0].n_genes), str(p.relative_to(ROOT)))


def load_missing_state(rows: list[dict]) -> None:
    p = ROOT / 'data' / 'publication_controls_v1' / 'missing_state_strong_baselines' / 'missing_state_baseline_results.tsv'
    d = pd.read_csv(p, sep='\t')
    for r in d.itertuples(index=False):
        rows.append(_row(
            'missing_state', str(r.dataset), str(r.context), str(r.model), 'strong_baseline',
            int(r.n_genes), float(r.spearman_full_vector), float(r.top100_overlap),
            float(r.top100_enrichment), float(r.top100_direction_agreement), True,
            str(p.relative_to(ROOT)),
            'Target state excluded from source construction.'
        ))
    for (ds, ctx), q in d.groupby(['dataset', 'context'], sort=False):
        add_random_expectation(rows, 'missing_state', str(ds), str(ctx), int(q.iloc[0].n_genes), str(p.relative_to(ROOT)))


def load_cross_experiment(rows: list[dict], negative_rows: list[dict]) -> None:
    base = ROOT / 'data' / 'external_heatshock_crossstudy' / 'GSE290176' / 'crossstudy_blind'
    vp = base / 'blind_vector_results.tsv'
    tp = base / 'blind_topk_results.tsv'
    v = pd.read_csv(vp, sep='\t')
    t = pd.read_csv(tp, sep='\t')
    t100 = t[t['k'] == 100].set_index('model')
    for r in v.itertuples(index=False):
        z = t100.loc[r.model]
        rows.append(_row(
            'cross_experiment', 'GSE121619_to_GSE290176_heatshock', 'whole_seedling_target', str(r.model),
            'matched_perturbation_baseline', int(r.n_genes), float(r.spearman_full_vector),
            float(z['overlap']), float(z['enrichment']), float(z['sign_agreement_overlap']), True,
            str(vp.relative_to(ROOT)),
            'Independent-study heat-shock transfer; target control expression was not used during freezing.'
        ))
    add_random_expectation(rows, 'cross_experiment', 'GSE121619_to_GSE290176_heatshock',
                           'whole_seedling_target', int(v.iloc[0].n_genes), str(vp.relative_to(ROOT)))

    npair = ROOT / 'data' / 'publication_controls_v1' / 'heatshock_crossstudy_negative_controls' / 'wrong_perturbation_pairwise.tsv'
    d = pd.read_csv(npair, sep='\t')
    for r in d.itertuples(index=False):
        negative_rows.append({
            'task_family': 'cross_experiment',
            'dataset': 'GSE121619_to_GSE290176_heatshock',
            'control_type': 'wrong_perturbation',
            'control_name': str(r.wrong_program),
            'common_universe_genes': int(r.common_universe_genes),
            'true_heat_rho': float(r.true_heat_rho),
            'control_rho': float(r.wrong_program_rho),
            'delta_rho_true_minus_control': float(r.delta_rho_true_minus_wrong),
            'true_heat_top100_overlap': float(r.true_heat_top100_overlap),
            'control_top100_overlap': float(r.wrong_program_top100_overlap),
            'delta_top100_true_minus_control': float(r.delta_top100_overlap_true_minus_wrong),
            'source_file': str(npair.relative_to(ROOT)),
        })


def load_combination(rows: list[dict]) -> None:
    base = ROOT / 'data' / 'external_combination' / 'GSE253757' / 'combination_blind'
    fp = base / 'FROZEN_COMBINATION_PREDICTIONS.tsv.gz'
    fr = pd.read_csv(fp, sep='\t')
    tech = fr['technical_universe'].astype(bool).to_numpy()
    parts = {
        'pat1_path1': ['resp_pat1', 'resp_path1'],
        'pat1_path2': ['resp_pat1', 'resp_path2'],
        'path1_path2': ['resp_path1', 'resp_path2'],
        'pat1_path1_path2': ['resp_pat1', 'resp_path1', 'resp_path2'],
    }
    for target, cols in parts.items():
        tp = base / f'TARGET_{target}.tsv.gz'
        td = pd.read_csv(tp, sep='\t')
        if not fr['gene'].astype(str).equals(td['gene'].astype(str)):
            raise RuntimeError(f'gene order mismatch for {target}')
        mask = tech & td['technical_universe'].astype(bool).to_numpy()
        truth = td.loc[mask, 'target_response'].to_numpy(float)
        components = {c: fr.loc[mask, c].to_numpy(float) for c in cols}
        sum_pred = fr.loc[mask, f'pred_{target}'].to_numpy(float)
        mean_pred = np.mean(np.column_stack(list(components.values())), axis=1)
        per_gene_strongest = fr.loc[mask, f'baseline_strongest_{target}'].to_numpy(float)
        rms = {c: float(np.sqrt(np.mean(v * v))) for c, v in components.items()}
        selected = max(rms, key=rms.get)

        models: list[tuple[str, str, np.ndarray, str]] = [
            ('simple_sum', 'composition', sum_pred, 'Frozen signed sum of constituent single-mutant responses.'),
            ('component_mean', 'composition', mean_pred, 'Arithmetic mean of constituent responses; rank metrics are scale-invariant relative to the sum.'),
            ('strongest_single_pretarget_rms', 'single_baseline', components[selected], f'Global single member selected before target reveal by source-only RMS: {selected}.'),
            ('strongest_constituent_per_gene', 'single_baseline', per_gene_strongest, 'Per-gene strongest constituent using source-only constituent responses.'),
        ]
        for c, v in components.items():
            models.append((f'single_{c.replace("resp_", "")}', 'single_baseline', v, 'Individual constituent single-mutant response.'))

        for model_name, model_class, pred, note in models:
            rho = float(spearmanr(pred, truth).statistic)
            ov, enrich, direction = _top100(pred, truth)
            rows.append(_row(
                'combination', 'GSE253757_PAT_family', target, model_name, model_class,
                len(truth), rho, ov, enrich, direction, True, str(fp.relative_to(ROOT)), note
            ))
        add_random_expectation(rows, 'combination', 'GSE253757_PAT_family', target, len(truth), str(fp.relative_to(ROOT)))


def aggregate(long: pd.DataFrame) -> pd.DataFrame:
    q = long.copy()
    out = (q.groupby(['task_family', 'dataset', 'model', 'model_class', 'pretarget_eligible'], dropna=False, as_index=False)
           .agg(n_contexts=('context', 'size'),
                median_rho=('spearman_full_vector', 'median'),
                mean_rho=('spearman_full_vector', 'mean'),
                median_top100_overlap=('top100_overlap', 'median'),
                mean_top100_overlap=('top100_overlap', 'mean'),
                median_top100_enrichment=('top100_enrichment', 'median'),
                median_top100_direction=('top100_direction_agreement', 'median')))
    return out


def main() -> None:
    rows: list[dict] = []
    negative_rows: list[dict] = []
    load_future_time(rows)
    load_missing_state(rows)
    load_cross_experiment(rows, negative_rows)
    load_combination(rows)

    long = pd.DataFrame(rows)
    long.to_csv(OUT / 'baseline_suite_long.tsv', sep='\t', index=False)
    agg = aggregate(long)
    agg.to_csv(OUT / 'baseline_suite_aggregate.tsv', sep='\t', index=False)
    pd.DataFrame(negative_rows).to_csv(OUT / 'negative_controls_cross_experiment.tsv', sep='\t', index=False)

    # Compact manuscript-facing summary: compare only models that cover every target context in that dataset.
    # This prevents a constituent-specific model that exists for only a subset of combination targets from
    # being mislabeled as the dataset-wide winner.
    manuscript = []
    elig = agg[(agg['pretarget_eligible']) & (agg['model'] != 'random_rank_expectation')].copy()
    for (task, ds), g in elig.groupby(['task_family', 'dataset'], sort=False):
        required_contexts = int(long[(long['task_family'] == task) & (long['dataset'] == ds)]['context'].nunique())
        full = g[g['n_contexts'] == required_contexts].copy()
        if full.empty:
            raise RuntimeError(f'no full-coverage baseline for {task} / {ds}')
        best = full.sort_values(['median_rho', 'median_top100_overlap'], ascending=False).iloc[0]
        manuscript.append({
            'task_family': task,
            'dataset': ds,
            'best_model_by_median_rho': best['model'],
            'median_rho': float(best['median_rho']),
            'median_top100_overlap': float(best['median_top100_overlap']),
            'n_contexts': int(best['n_contexts']),
            'required_contexts': required_contexts,
        })
    manuscript = pd.DataFrame(manuscript)
    manuscript.to_csv(OUT / 'manuscript_baseline_winners.tsv', sep='\t', index=False)

    manifest = {
        'suite_version': 'publication_baseline_suite_v1',
        'k_primary': K,
        'principles': [
            'Only target-compatible simple baselines are compared.',
            'Main baseline rows use information available before target reveal.',
            'Random ranking is represented by its analytic expectation rather than one arbitrary draw.',
            'Cross-experiment wrong-perturbation controls are stored separately because each uses a matched common gene universe.',
            'Combination component mean is intentionally retained even though positive rescaling makes its rank metrics identical to the simple sum.'
        ],
        'n_rows_long': int(len(long)),
        'n_rows_aggregate': int(len(agg)),
        'n_negative_controls': int(len(negative_rows)),
        'task_families': sorted(long['task_family'].unique().tolist()),
        'outputs': [
            'baseline_suite_long.tsv',
            'baseline_suite_aggregate.tsv',
            'negative_controls_cross_experiment.tsv',
            'manuscript_baseline_winners.tsv',
        ],
    }
    (OUT / 'MANIFEST.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')

    print('PUBLICATION_BASELINE_SUITE_SUCCESS')
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    print('\nMANUSCRIPT BASELINE WINNERS')
    print(manuscript.round(4).to_string(index=False))


if __name__ == '__main__':
    main()
