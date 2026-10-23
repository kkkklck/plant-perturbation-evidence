from __future__ import annotations

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data' / 'external_combination' / 'GSE253757'
RAW = BASE / 'GSE253757_RawData_withAnnotations.xls.gz'
RULES = ROOT / '组合扰动_GSE253757_盲测规则_2026-10-20.md'
OUT = BASE / 'combination_blind'
OUT.mkdir(parents=True, exist_ok=True)

ANNOT = ['TAIR_locus','Description','EntrezID','EntrezIDs','Gene_symbol','Ensembl_gene_id']
PRE = [str(i) for i in range(1, 13)]
TARGET = [str(i) for i in range(13, 25)]
GROUPS = {
    'background': ['1','2','3'],
    'pat1': ['4','5','6'],
    'path1': ['7','8','9'],
    'path2': ['10','11','12'],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    # Hard firewall: only annotation + columns 1..12 are ever requested here.
    usecols = ANNOT + PRE
    if any(c in usecols for c in TARGET):
        raise RuntimeError('TARGET FIREWALL BREACH')

    d = pd.read_csv(RAW, sep='\t', compression='gzip', usecols=usecols, dtype={c: str for c in ANNOT})
    if any(c in d.columns for c in TARGET):
        raise RuntimeError('TARGET COLUMNS WERE READ')

    genes = d['TAIR_locus'].astype(str)
    X = d[PRE].apply(pd.to_numeric, errors='raise').to_numpy(float)
    lib = X.sum(axis=0)
    if np.any(lib <= 0):
        raise RuntimeError('zero library size')
    cpm = X / lib[None, :] * 1e6
    L = np.log2(cpm + 1.0)

    col_idx = {c: i for i, c in enumerate(PRE)}
    def mean_cols(cols: list[str]) -> np.ndarray:
        return L[:, [col_idx[c] for c in cols]].mean(axis=1)

    bg = mean_cols(GROUPS['background'])
    A = mean_cols(GROUPS['pat1']) - bg
    B = mean_cols(GROUPS['path1']) - bg
    C = mean_cols(GROUPS['path2']) - bg

    technical = (cpm >= 1.0).sum(axis=1) >= 3

    preds = {
        'pat1_path1': (A + B, [A, B]),
        'pat1_path2': (A + C, [A, C]),
        'path1_path2': (B + C, [B, C]),
        'pat1_path1_path2': (A + B + C, [A, B, C]),
    }

    out = pd.DataFrame({
        'gene': genes,
        'technical_universe': technical,
        'resp_pat1': A,
        'resp_path1': B,
        'resp_path2': C,
    })

    for name, (pred, parts) in preds.items():
        P = np.column_stack(parts)
        denom = np.abs(P).sum(axis=1)
        coherence = np.zeros(len(pred), dtype=float)
        np.divide(np.abs(pred), denom, out=coherence, where=denom > 0)
        strongest = P[np.arange(len(pred)), np.argmax(np.abs(P), axis=1)]
        out[f'pred_{name}'] = pred
        out[f'score_{name}'] = np.abs(pred)
        out[f'coherence_{name}'] = coherence
        out[f'baseline_strongest_{name}'] = strongest
        out[f'baseline_strongest_score_{name}'] = np.abs(strongest)

    pred_path = OUT / 'FROZEN_COMBINATION_PREDICTIONS.tsv.gz'
    out.to_csv(pred_path, sep='\t', index=False, compression='gzip')

    manifest = {
        'status': 'FROZEN_BEFORE_COMBINATION_TARGETS',
        'source_raw_sha256': sha256(RAW),
        'rules_sha256': sha256(RULES),
        'freeze_script_sha256': sha256(Path(__file__)),
        'prediction_sha256': sha256(pred_path),
        'pretarget_columns_read': PRE,
        'target_columns_not_read': TARGET,
        'n_genes_all': int(len(out)),
        'n_genes_technical': int(technical.sum()),
        'target_values_read': False,
        'primary_rule': 'signed sum of constituent single-mutant responses',
        'targets': {
            'pat1_path1': ['13','14','15'],
            'pat1_path2': ['16','17','18'],
            'path1_path2': ['19','20','21'],
            'pat1_path1_path2': ['22','23','24'],
        },
    }
    (OUT / 'FREEZE_MANIFEST.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding='utf-8')
    print('GSE253757_COMBINATION_FREEZE_SUCCESS')
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    print('COMBINATION TARGET VALUES HAVE NOT BEEN READ')


if __name__ == '__main__':
    main()
