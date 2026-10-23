from pathlib import Path
import gzip, hashlib, json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'data' / 'external_heatshock' / 'GSE121619' / 'plantvc_v1_heatshock_blind' / 'frozen_prediction.tsv.gz'
TGT = ROOT / 'data' / 'external_heatshock_crossstudy' / 'GSE290176' / 'GSE290176_counts_HS_Ath.txt.gz'
RULES = ROOT / '急性热激跨实验迁移_GSE290176_盲测规则_2026-10-20.md'
OUT = ROOT / 'data' / 'external_heatshock_crossstudy' / 'GSE290176' / 'crossstudy_blind'
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def norm_gene(x: str) -> str:
    return str(x).strip().strip('"').upper()

# Read ONLY the first field of each target line. No expression field is parsed or converted.
target_genes = []
with gzip.open(TGT, 'rb') as f:
    header = f.readline().rstrip(b'\r\n').split(b'\t')
    if not header or header[0].decode('utf-8', errors='replace').strip().strip('"') != 'Geneid':
        raise RuntimeError('Unexpected GSE290176 count-table header')
    for raw in f:
        raw = raw.rstrip(b'\r\n')
        if not raw:
            continue
        i = raw.find(b'\t')
        if i < 0:
            continue
        g = raw[:i].decode('utf-8', errors='replace').strip().strip('"').upper()
        if g:
            target_genes.append(g)

target_gene_set = set(target_genes)

src = pd.read_csv(SRC, sep='\t')
src['gene_norm'] = src['gene'].map(norm_gene)
if src['gene_norm'].duplicated().any():
    raise RuntimeError('Duplicate normalized genes in source frozen prediction')

need_cols = ['shared_prediction', 'pooled_donor_prediction', 'donor_5', 'donor_6', 'donor_8']
for c in need_cols + ['technical']:
    if c not in src.columns:
        raise RuntimeError(f'Missing source column: {c}')

keep = src['technical'].astype(bool) & src['gene_norm'].isin(target_gene_set)
frozen = src.loc[keep, ['gene_norm'] + need_cols].copy()
frozen = frozen.rename(columns={'gene_norm': 'gene'})

# Ranking is frozen here, before any target expression is opened.
for c in need_cols:
    frozen[f'{c}_abs_rank'] = frozen[c].abs().rank(method='first', ascending=False).astype(int)

frozen = frozen.sort_values('shared_prediction_abs_rank').reset_index(drop=True)
pred_path = OUT / 'frozen_crossstudy_prediction.tsv.gz'
frozen.to_csv(pred_path, sep='\t', index=False, compression='gzip')

manifest = {
    'status': 'FROZEN_BEFORE_GSE290176_HEATSHOCK_TARGET',
    'source_prediction_sha256': sha256(SRC),
    'target_count_file_sha256': sha256(TGT),
    'rules_sha256': sha256(RULES),
    'freeze_script_sha256': sha256(Path(__file__)),
    'prediction_sha256': sha256(pred_path),
    'source_primary_column': 'shared_prediction',
    'frozen_comparators': ['pooled_donor_prediction', 'donor_5', 'donor_6', 'donor_8'],
    'n_target_gene_ids_seen': len(target_gene_set),
    'n_source_technical_genes': int(src['technical'].astype(bool).sum()),
    'n_frozen_genes': int(len(frozen)),
    'target_expression_values_parsed': False,
    'target_heatshock_values_used': False,
    'target_control_values_used': False,
    'target_expression_access_at_freeze': 'Geneid field only; no numeric expression field parsed',
}
with open(OUT / 'FREEZE_MANIFEST.json', 'w', encoding='utf-8') as f:
    json.dump(manifest, f, ensure_ascii=False, indent=2)

with open(OUT / 'freeze.log', 'w', encoding='utf-8') as f:
    f.write('GSE290176_CROSSSTUDY_HEAT_FREEZE_SUCCESS\n')
    f.write(json.dumps(manifest, ensure_ascii=False, indent=2))
    f.write('\nTARGET EXPRESSION NUMERIC FIELDS NOT PARSED OR USED\n')

print('GSE290176_CROSSSTUDY_HEAT_FREEZE_SUCCESS')
print(json.dumps(manifest, ensure_ascii=False, indent=2))
print('TARGET EXPRESSION NUMERIC FIELDS NOT PARSED OR USED')
