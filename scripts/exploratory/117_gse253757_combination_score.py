from __future__ import annotations

from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data' / 'external_combination' / 'GSE253757'
RAW = BASE / 'GSE253757_RawData_withAnnotations.xls.gz'
RULES = ROOT / '组合扰动_GSE253757_盲测规则_2026-10-20.md'
OUT = BASE / 'combination_blind'
FREEZE_SCRIPT = ROOT / 'scripts' / '116_gse253757_combination_freeze.py'
PRED = OUT / 'FROZEN_COMBINATION_PREDICTIONS.tsv.gz'
MANIFEST = OUT / 'FREEZE_MANIFEST.json'

ANNOT = ['TAIR_locus','Description','EntrezID','EntrezIDs','Gene_symbol','Ensembl_gene_id']
CONTROL = ['1','2','3']
TARGETS = {
    'pat1_path1': ['13','14','15'],
    'pat1_path2': ['16','17','18'],
    'path1_path2': ['19','20','21'],
    'pat1_path1_path2': ['22','23','24'],
}
PARTS = {
    'pat1_path1': ['resp_pat1','resp_path1'],
    'pat1_path2': ['resp_pat1','resp_path2'],
    'path1_path2': ['resp_path1','resp_path2'],
    'pat1_path1_path2': ['resp_pat1','resp_path1','resp_path2'],
}
KS = [100, 250, 500, 1000]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def topk(score: np.ndarray, signed: np.ndarray, target: np.ndarray, k: int) -> dict:
    n = len(score)
    k = min(k, n)
    pi = np.argpartition(-score, k - 1)[:k]
    ti = np.argpartition(-np.abs(target), k - 1)[:k]
    inter = np.intersect1d(pi, ti)
    ov = len(inter)
    exp = k * k / n
    order = np.argsort(-score, kind='mergesort')
    ranks = np.empty(n, int)
    ranks[order] = np.arange(1, n + 1)
    pct = 1 - (ranks[ti] - 1) / max(n - 1, 1)
    return {
        'k': k,
        'overlap': ov,
        'random_expected_overlap': exp,
        'enrichment': ov / exp if exp else np.nan,
        'hypergeom_p': float(hypergeom.sf(ov - 1, n, k, k)),
        'sign_agreement_overlap': float(np.mean(np.sign(signed[inter]) == np.sign(target[inter]))) if ov else np.nan,
        'median_target_topk_predicted_percentile': float(np.median(pct)),
    }


def main() -> None:
    m = json.loads(MANIFEST.read_text(encoding='utf-8'))
    if m.get('status') != 'FROZEN_BEFORE_COMBINATION_TARGETS':
        raise RuntimeError('invalid freeze status')
    checks = {
        'source_raw_sha256': sha256(RAW),
        'rules_sha256': sha256(RULES),
        'freeze_script_sha256': sha256(FREEZE_SCRIPT),
        'prediction_sha256': sha256(PRED),
    }
    for key, val in checks.items():
        if m.get(key) != val:
            raise RuntimeError(f'Frozen artifact hash mismatch: {key}')
    if m.get('target_values_read') is not False:
        raise RuntimeError('freeze manifest does not certify withheld targets')

    fr = pd.read_csv(PRED, sep='\t', compression='gzip')

    # Target reveal happens only here, after all hashes/rules are verified.
    target_cols = sorted({c for cols in TARGETS.values() for c in cols}, key=int)
    d = pd.read_csv(RAW, sep='\t', compression='gzip', usecols=ANNOT + CONTROL + target_cols, dtype={c: str for c in ANNOT})
    if d['TAIR_locus'].astype(str).tolist() != fr['gene'].astype(str).tolist():
        raise RuntimeError('gene order mismatch')

    cols = CONTROL + target_cols
    X = d[cols].apply(pd.to_numeric, errors='raise').to_numpy(float)
    lib = X.sum(axis=0)
    if np.any(lib <= 0):
        raise RuntimeError('zero library size')
    cpm = X / lib[None, :] * 1e6
    L = np.log2(cpm + 1.0)
    idx = {c: i for i, c in enumerate(cols)}
    bg = L[:, [idx[c] for c in CONTROL]].mean(axis=1)

    tech = fr['technical_universe'].astype(bool).to_numpy()
    vector_rows = []
    top_rows = []
    repl_rows = []
    comparison_rows = []

    for target_name, reps in TARGETS.items():
        target_rep = np.column_stack([L[:, idx[c]] - bg for c in reps])
        target = target_rep.mean(axis=1)

        mask = tech & np.isfinite(target)
        y = target[mask]

        models: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        comp = fr[f'pred_{target_name}'].to_numpy(float)[mask]
        models['composition'] = (np.abs(comp), comp)
        strongest = fr[f'baseline_strongest_{target_name}'].to_numpy(float)[mask]
        models['strongest_constituent_per_gene'] = (np.abs(strongest), strongest)
        for part in PARTS[target_name]:
            v = fr[part].to_numpy(float)[mask]
            models[part.replace('resp_', 'single_')] = (np.abs(v), v)

        for model_name, (score, signed) in models.items():
            rho = float(spearmanr(signed, y).statistic)
            sign_all = float(np.mean(np.sign(signed) == np.sign(y)))
            vector_rows.append({
                'target': target_name,
                'model': model_name,
                'n_genes': int(mask.sum()),
                'spearman_full_vector': rho,
                'sign_agreement_all': sign_all,
            })
            for k in KS:
                top_rows.append({'target': target_name, 'model': model_name, **topk(score, signed, y, k)})

        # Pairwise target-replicate reliability, using same frozen technical universe.
        reps_masked = target_rep[mask]
        for i in range(len(reps)):
            for j in range(i + 1, len(reps)):
                a = reps_masked[:, i]
                b = reps_masked[:, j]
                repl_rows.append({
                    'target': target_name,
                    'rep_a': reps[i],
                    'rep_b': reps[j],
                    'spearman_full_vector': float(spearmanr(a, b).statistic),
                    'sign_agreement_all': float(np.mean(np.sign(a) == np.sign(b))),
                })

        # Compact head-to-head summary against the best single baseline after reveal.
        vdf = pd.DataFrame([r for r in vector_rows if r['target'] == target_name])
        comp_rho = float(vdf.loc[vdf.model == 'composition', 'spearman_full_vector'].iloc[0])
        single_df = vdf[vdf.model.str.startswith('single_')]
        best_single_rho = float(single_df['spearman_full_vector'].max())
        best_single_name = str(single_df.loc[single_df['spearman_full_vector'].idxmax(), 'model'])
        t100 = pd.DataFrame([r for r in top_rows if r['target'] == target_name and r['k'] == 100])
        comp_ov = int(t100.loc[t100.model == 'composition', 'overlap'].iloc[0])
        single100 = t100[t100.model.str.startswith('single_')]
        best_single_ov = int(single100['overlap'].max())
        best_single100_name = str(single100.loc[single100['overlap'].idxmax(), 'model'])
        comparison_rows.append({
            'target': target_name,
            'composition_rho': comp_rho,
            'best_single_rho': best_single_rho,
            'best_single_rho_model': best_single_name,
            'delta_rho_vs_best_single': comp_rho - best_single_rho,
            'composition_top100_overlap': comp_ov,
            'best_single_top100_overlap': best_single_ov,
            'best_single_top100_model': best_single100_name,
            'delta_top100_overlap_vs_best_single': comp_ov - best_single_ov,
        })

        pd.DataFrame({
            'gene': fr['gene'],
            'technical_universe': tech,
            'target_response': target,
            'target_rep1_response': target_rep[:, 0],
            'target_rep2_response': target_rep[:, 1],
            'target_rep3_response': target_rep[:, 2],
        }).to_csv(OUT / f'TARGET_{target_name}.tsv.gz', sep='\t', index=False, compression='gzip')

    vdf = pd.DataFrame(vector_rows)
    tdf = pd.DataFrame(top_rows)
    rdf = pd.DataFrame(repl_rows)
    cdf = pd.DataFrame(comparison_rows)
    vdf.to_csv(OUT / 'blind_vector_results.tsv', sep='\t', index=False)
    tdf.to_csv(OUT / 'blind_topk_results.tsv', sep='\t', index=False)
    rdf.to_csv(OUT / 'target_replicate_reference.tsv', sep='\t', index=False)
    cdf.to_csv(OUT / 'composition_vs_best_single.tsv', sep='\t', index=False)

    comp_v = vdf[vdf.model == 'composition']
    comp_t = tdf[tdf.model == 'composition']
    summary = {
        'n_targets': len(TARGETS),
        'median_composition_rho': float(comp_v['spearman_full_vector'].median()),
        'median_composition_sign_agreement_all': float(comp_v['sign_agreement_all'].median()),
        'composition_beats_best_single_rho': int((cdf['delta_rho_vs_best_single'] > 0).sum()),
        'composition_ties_best_single_rho': int((cdf['delta_rho_vs_best_single'] == 0).sum()),
        'composition_beats_best_single_top100': int((cdf['delta_top100_overlap_vs_best_single'] > 0).sum()),
        'composition_ties_best_single_top100': int((cdf['delta_top100_overlap_vs_best_single'] == 0).sum()),
        'topk': {},
    }
    for k in KS:
        q = comp_t[comp_t.k == k]
        summary['topk'][str(k)] = {
            'median_overlap': float(q['overlap'].median()),
            'median_enrichment': float(q['enrichment'].median()),
            'median_sign_agreement_overlap': float(q['sign_agreement_overlap'].median()),
        }
    (OUT / 'BLIND_RESULT.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding='utf-8')

    print('GSE253757_COMBINATION_SCORE_SUCCESS')
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print('\nHEAD-TO-HEAD')
    print(cdf.to_string(index=False))
    print('\nCOMPOSITION TOP-K')
    print(comp_t.to_string(index=False))


if __name__ == '__main__':
    main()
