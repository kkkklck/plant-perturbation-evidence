from pathlib import Path
import gzip, hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, hypergeom

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data' / 'external_heatshock_crossstudy' / 'GSE290176'
OUT = BASE / 'crossstudy_blind'
TGT = BASE / 'GSE290176_counts_HS_Ath.txt.gz'
PRED = OUT / 'frozen_crossstudy_prediction.tsv.gz'
MANIFEST = OUT / 'FREEZE_MANIFEST.json'

CK = [f'Col0-CK{i}' for i in range(1, 7)]
HS = [f'Col0-HS{i}' for i in range(1, 4)]
MODELS = ['shared_prediction', 'pooled_donor_prediction', 'donor_5', 'donor_6', 'donor_8']
KS = [100, 250, 500, 1000]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def norm_gene(x: str) -> str:
    return str(x).strip().strip('"').upper()

with open(MANIFEST, encoding='utf-8') as f:
    manifest = json.load(f)
if manifest.get('status') != 'FROZEN_BEFORE_GSE290176_HEATSHOCK_TARGET':
    raise RuntimeError('Freeze manifest status invalid')
if sha256(PRED) != manifest['prediction_sha256']:
    raise RuntimeError('Frozen prediction hash changed')
if sha256(TGT) != manifest['target_count_file_sha256']:
    raise RuntimeError('Target count file hash changed')

pred = pd.read_csv(PRED, sep='\t')
pred['gene'] = pred['gene'].map(norm_gene)

# Target is opened only after freeze integrity checks above.
counts = pd.read_csv(TGT, sep='\t', compression='gzip')
required = ['Geneid'] + CK + HS
missing = [c for c in required if c not in counts.columns]
if missing:
    raise RuntimeError(f'Missing target columns: {missing}')
counts = counts[required].copy()
counts['gene'] = counts['Geneid'].map(norm_gene)
for c in CK + HS:
    counts[c] = pd.to_numeric(counts[c], errors='coerce').fillna(0.0)

# Collapse accidental duplicate gene IDs by summing raw counts.
counts = counts.groupby('gene', as_index=False)[CK + HS].sum()
all_cols = CK + HS
libs = counts[all_cols].sum(axis=0)
if (libs <= 0).any():
    raise RuntimeError('Zero library size in target samples')
logcpm = counts[['gene']].copy()
for c in all_cols:
    cpm = counts[c].to_numpy(dtype=float) / float(libs[c]) * 1e6
    logcpm[c] = np.log2(cpm + 1.0)

logcpm['target_response'] = logcpm[HS].mean(axis=1) - logcpm[CK].mean(axis=1)
for h in HS:
    logcpm[f'resp_{h}'] = logcpm[h] - logcpm[CK].mean(axis=1)

z = pred.merge(logcpm, on='gene', how='inner', validate='one_to_one')
if len(z) != len(pred):
    raise RuntimeError(f'Frozen gene loss after target merge: {len(pred)} -> {len(z)}')
N = len(z)
truth = z['target_response'].to_numpy(dtype=float)
truth_abs_rank = pd.Series(np.abs(truth)).rank(method='first', ascending=False).to_numpy(dtype=int)
z['target_abs_rank'] = truth_abs_rank

vector_rows = []
topk_rows = []
for model in MODELS:
    pv = z[model].to_numpy(dtype=float)
    rho = float(spearmanr(pv, truth).statistic)
    mask = (pv != 0) & (truth != 0) & np.isfinite(pv) & np.isfinite(truth)
    sign_all = float(np.mean(np.sign(pv[mask]) == np.sign(truth[mask]))) if mask.any() else np.nan
    vector_rows.append({
        'model': model,
        'n_genes': N,
        'spearman_full_vector': rho,
        'sign_agreement_all_nonzero': sign_all,
    })
    pr = pd.Series(np.abs(pv)).rank(method='first', ascending=False).to_numpy(dtype=int)
    for k in KS:
        kk = min(k, N)
        pset = set(np.where(pr <= kk)[0])
        tset = set(np.where(truth_abs_rank <= kk)[0])
        ov = sorted(pset & tset)
        overlap = len(ov)
        expected = kk * kk / N
        enrich = overlap / expected if expected > 0 else np.nan
        pval = float(hypergeom.sf(overlap - 1, N, kk, kk))
        sign = float(np.mean(np.sign(pv[ov]) == np.sign(truth[ov]))) if overlap else np.nan
        # For each true top-k gene, percentile of its predicted absolute rank; 1 is best.
        target_idx = np.array(sorted(tset), dtype=int)
        pct = 1.0 - (pr[target_idx] - 1) / max(N - 1, 1)
        topk_rows.append({
            'model': model,
            'k': kk,
            'overlap': overlap,
            'random_expected_overlap': expected,
            'enrichment': enrich,
            'hypergeom_p': pval,
            'sign_agreement_overlap': sign,
            'median_target_topk_predicted_percentile': float(np.median(pct)),
        })

vdf = pd.DataFrame(vector_rows)
tdf = pd.DataFrame(topk_rows)
vdf.to_csv(OUT / 'blind_vector_results.tsv', sep='\t', index=False)
tdf.to_csv(OUT / 'blind_topk_results.tsv', sep='\t', index=False)

# Target replicate reference: each HS replicate relative to the same mean CK profile.
rep_rows = []
rep_vecs = {h: z[f'resp_{h}'].to_numpy(dtype=float) for h in HS}
for i in range(len(HS)):
    for j in range(i + 1, len(HS)):
        a, b = HS[i], HS[j]
        va, vb = rep_vecs[a], rep_vecs[b]
        rho = float(spearmanr(va, vb).statistic)
        mask = (va != 0) & (vb != 0)
        sign = float(np.mean(np.sign(va[mask]) == np.sign(vb[mask]))) if mask.any() else np.nan
        row = {'rep_a': a, 'rep_b': b, 'spearman_full_vector': rho, 'sign_agreement_all_nonzero': sign}
        for k in [100, 500]:
            kk = min(k, N)
            ra = pd.Series(np.abs(va)).rank(method='first', ascending=False).to_numpy(dtype=int)
            rb = pd.Series(np.abs(vb)).rank(method='first', ascending=False).to_numpy(dtype=int)
            ov = len(set(np.where(ra <= kk)[0]) & set(np.where(rb <= kk)[0]))
            row[f'top{k}_overlap'] = ov
            row[f'top{k}_enrichment'] = ov / (kk * kk / N)
        rep_rows.append(row)
rdf = pd.DataFrame(rep_rows)
rdf.to_csv(OUT / 'target_replicate_reference.tsv', sep='\t', index=False)

primary_v = vdf[vdf.model == 'shared_prediction'].iloc[0].to_dict()
primary_top = tdf[tdf.model == 'shared_prediction'].set_index('k')
summary = {
    'status': 'TARGET_REVEALED_AFTER_VALID_FREEZE',
    'n_genes': N,
    'primary_model': 'shared_prediction',
    'primary_full_vector_rho': primary_v['spearman_full_vector'],
    'primary_sign_agreement_all_nonzero': primary_v['sign_agreement_all_nonzero'],
    'primary_topk': {
        str(k): {
            'overlap': int(primary_top.loc[min(k, N), 'overlap']),
            'enrichment': float(primary_top.loc[min(k, N), 'enrichment']),
            'sign_agreement_overlap': float(primary_top.loc[min(k, N), 'sign_agreement_overlap']) if pd.notna(primary_top.loc[min(k, N), 'sign_agreement_overlap']) else None,
            'hypergeom_p': float(primary_top.loc[min(k, N), 'hypergeom_p']),
        } for k in KS
    },
    'model_rho': {r['model']: float(r['spearman_full_vector']) for _, r in vdf.iterrows()},
    'target_replicate_rho_median': float(rdf['spearman_full_vector'].median()),
    'target_replicate_top100_overlap_median': float(rdf['top100_overlap'].median()),
}
with open(OUT / 'BLIND_RESULT.json', 'w', encoding='utf-8') as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)

print('GSE290176_CROSSSTUDY_HEAT_SCORE_SUCCESS')
print(json.dumps(summary, ensure_ascii=False, indent=2))
print('\nMODEL VECTOR RESULTS')
print(vdf.to_string(index=False))
print('\nPRIMARY TOP-K')
print(tdf[tdf.model == 'shared_prediction'].to_string(index=False))
print('\nTARGET REPLICATE REFERENCE')
print(rdf.to_string(index=False))
