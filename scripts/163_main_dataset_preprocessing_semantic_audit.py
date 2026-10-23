from __future__ import annotations

from pathlib import Path
import json
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'publication_controls_v1' / 'main_dataset_preprocessing_semantic_audit'
OUT.mkdir(parents=True, exist_ok=True)

rows: list[dict] = []


def add(dataset, accession, input_kind, deposited_semantics, manuscript_transform, checks, status='PASS', note=''):
    rows.append({
        'dataset': dataset,
        'accession': accession,
        'input_kind': input_kind,
        'deposited_or_local_semantics': deposited_semantics,
        'manuscript_transform': manuscript_transform,
        'checks': checks,
        'status': status,
        'note': note,
    })


def read_soft_meta(path: Path):
    processing=[]; table_header=''; vals=[]; inside=False
    with path.open('r',encoding='utf-8',errors='replace') as f:
        for line in f:
            if line.startswith('!Sample_data_processing'):
                processing.append(line.strip())
            if line.startswith('!sample_table_begin'):
                inside=True; continue
            if line.startswith('!sample_table_end'):
                break
            if inside:
                if not table_header:
                    table_header=line.strip(); continue
                parts=line.rstrip('\n').split('\t')
                if len(parts)>=2:
                    try: vals.append(float(parts[1]))
                    except Exception: pass
                if len(vals)>=1000: break
    return ' | '.join(processing), table_header, np.asarray(vals,float)

# -----------------------------------------------------------------------------
# GSE153103: supplementary file explicitly named raw_reads_counts; verify numeric content is count-like.
# -----------------------------------------------------------------------------
salt = ROOT/'data/external_salt/GSE153103/GSE153103_1_raw_reads_counts.xlsx'
h = pd.read_excel(salt, sheet_name='raw_counts', nrows=0)
cols=list(h.columns)
d = pd.read_excel(salt, sheet_name='raw_counts', usecols=cols[:min(9,len(cols))], nrows=500)
num=d.iloc[:,1:].apply(pd.to_numeric,errors='coerce').to_numpy(float)
finite=num[np.isfinite(num)]
integer_like_fraction = float(np.mean(np.isclose(finite, np.round(finite), atol=1e-9)))
nonnegative = bool(np.all(finite>=0))
condition = pd.read_excel(salt, sheet_name='condition')
condition_map = dict(zip(condition.iloc[:,0].astype(str), condition['Time'].astype(str)))
assert nonnegative and condition_map.get('RAS048')=='0h' and condition_map.get('RAS049')=='1h' and condition_map.get('RAS050')=='3h' and condition_map.get('RAS051')=='6h' and condition_map.get('RAS052')=='12h'
add('Arabidopsis salt time course','GSE153103','supplementary read-count workbook',
    'author/GEO supplementary file GSE153103_1_raw_reads_counts.xlsx; values are non-negative read-count estimates, with a minority of fractional entries',
    'library-size CPM followed by log2(CPM+1); response = treated time point minus matched 0h control within replicate',
    f'first 500 rows × {num.shape[1]} numeric columns: nonnegative={nonnegative}; exactly-integer fraction={integer_like_fraction:.4f}; condition sheet confirms RAS048/49/50/51/52 = 0/1/3/6/12h',
    note='Do not describe all entries as integer raw counts; the deposited workbook contains some fractional count estimates. CPM normalization remains well-defined.')

# -----------------------------------------------------------------------------
# GSE30098 sulfur: deposited VALUE is MAS5 linear signal, so log2 is appropriate.
# -----------------------------------------------------------------------------
sulfur = ROOT/'data/prospective_task/sulfur/source_0_3h/GSM744800.soft.txt'
proc,header,v = read_soft_meta(sulfur)
assert 'MAS 5.0' in proc and 'global scaling' in proc
assert len(v)>0 and np.all(v>0) and np.nanmedian(v)>1
add('Sulfur deficiency 0h→3h→12h','GSE30098','Affymetrix ATH1 processed sample tables',
    'MAS 5.0 probe-set Signal with global scaling; mean target intensity 250; deposited VALUE is linear positive signal',
    'log2(VALUE), then treatment minus 0h control on log2 scale',
    f'metadata contains MAS 5.0/global scaling; sampled VALUE median={float(np.median(v)):.3g}; min={float(np.min(v)):.3g}',
    note='This is the opposite of GSE46205: here the deposited Signal is not already log2.')

# -----------------------------------------------------------------------------
# GSE121619: legacy RDS expression matrix is raw UMI counts, verified exactly against metadata total_umis.
# -----------------------------------------------------------------------------
hs_log = (ROOT/'data/external_heatshock/GSE121619/plantvc_v1_heatshock_blind/freeze.log').read_text(encoding='utf-8',errors='replace')
med = re.search(r'allowed_library_size_median_relative_error=([^\n]+)',hs_log)
mx = re.search(r'allowed_library_size_max_relative_error=([^\n]+)',hs_log)
assert med and mx and float(med.group(1)) <= 1e-9 and float(mx.group(1)) <= 1e-9
add('Heat-shock held-out states','GSE121619','legacy Monocle RDS expression matrix',
    'raw UMI counts: column sums match metadata total_umis exactly in allowed cells',
    'pseudobulk per cluster/treatment; natural log1p(CPM); response = HeatShock − Control',
    f'library-size median relative error={float(med.group(1))}; max={float(mx.group(1))}',
    note='No normalization stored in the RDS expression matrix is assumed; raw-count semantics were checked numerically.')

# -----------------------------------------------------------------------------
# GSE158761: author RDS counts slot; freeze/target scripts pseudobulk raw counts and transform once.
# -----------------------------------------------------------------------------
lr_log = (ROOT/'data/external_lateral_root/GSE158761/plantvc_v1_blind/freeze.log').read_text(encoding='utf-8',errors='replace')
assert 'counts_dim=' in lr_log and 'CELL_MAPPING_OK' in lr_log
freeze_r = (ROOT/'scripts/112_gse158761_freeze.R').read_text(encoding='utf-8')
target_r = (ROOT/'scripts/114_gse158761_target.R').read_text(encoding='utf-8')
assert 'ld[["counts"]]' in freeze_r and 'log2(CPM + 1)' in freeze_r
assert 'ld[["counts"]]' in target_r and 'log2(CPM + 1)' in target_r
add('Lateral-root 0h→8h→20h','GSE158761','processed author RDS, recovered sparse counts slot',
    'sparse UMI count matrix; GEO sample documentation describes matrix.mtx as UMI counts; cells matched exactly to metadata',
    'cell-type × sample pseudobulk raw counts; CPM; log2(CPM+1); 8h/20h responses relative to matched 0h replicates',
    'freeze log reports counts matrix and exact cell mapping; both freeze and target scripts recover counts slot and transform once',
    note='No use of normalized Seurat/Monocle data layer for the response values.')

# -----------------------------------------------------------------------------
# GSE290176: GEO supplementary matrix is featureCounts raw counts.
# -----------------------------------------------------------------------------
cross = ROOT/'data/external_heatshock_crossstudy/GSE290176/GSE290176_counts_HS_Ath.txt.gz'
cc = pd.read_csv(cross,sep='\t',nrows=500)
numcols=[c for c in cc.columns if c!='Geneid']
arr=cc[numcols].apply(pd.to_numeric,errors='coerce').to_numpy(float)
finite=arr[np.isfinite(arr)]
assert np.all(finite>=0) and np.allclose(finite,np.round(finite),atol=1e-9)
add('Independent whole-seedling heat shock','GSE290176','supplementary featureCounts matrix',
    'raw gene-level read counts; GEO sample record states featureCounts generated the matrix',
    'library-size CPM followed by log2(CPM+1); target response = Col-0 heat-shock mean − Col-0 control mean',
    f'first 500 rows: nonnegative integer-like counts across {len(numcols)} sample columns',
    note='Only Col-0 control and Col-0 heat-shock columns enter the cross-study target.')

# -----------------------------------------------------------------------------
# Robertson GSE313726: Seurat v5 RNA counts layer → pseudobulk → log1p CPM.
# -----------------------------------------------------------------------------
rob_log=(ROOT/'data/external_rice/Robertson2026/pseudobulk_v0/build.log').read_text(encoding='utf-8',errors='replace')
assert 'RNA_layers=data,counts' in rob_log and 'PSEUDOBULK_COUNTS_OK' in rob_log
rob_r=(ROOT/'scripts/93_build_robertson2026_pseudobulk.R').read_text(encoding='utf-8')
assert 'LayerData(obj, assay="RNA", layer="counts")' in rob_r and 'logcpm <- log1p(cpm)' in rob_r
add('Rice drought/salinity candidate stability','GSE313726','processed Seurat v5 object, RNA counts layer',
    'RNA counts layer explicitly extracted; GEO describes integrated annotated RDS plus count matrices',
    'sample × cell-annotation pseudobulk raw counts, aggregated to tissue system where needed; CPM; natural log1p(CPM)',
    'build log: 126,139 nuclei, RNA layers data+counts, pseudobulk counts successfully created',
    note='Candidate-stability analysis uses response differences on the derived log1p(CPM) scale.')

# -----------------------------------------------------------------------------
# Immune: use author-generated Seurat/MAST FindMarkers avg_logFC tables without re-transforming.
# -----------------------------------------------------------------------------
author_de = ROOT/'data/prospective_task/immune_missing_state/scripts/scripts/05.DE_zone4_celltype6_stim.R'
txt=author_de.read_text(encoding='utf-8',errors='replace')
assert 'FindMarkers' in txt and 'test.use = "MAST"' in txt and 'logfc.threshold = 0' in txt
sample_tsv = ROOT/'data/prospective_task/immune_missing_state/source_only/Pep1/Rootcap_Columella.tsv.gz'
it = pd.read_csv(sample_tsv,sep='\t',nrows=10)
assert 'avg_logFC' in it.columns and 'pct.1' in it.columns and 'pct.2' in it.columns
add('Five-stimulus root immunity','GSE289716','author differential-expression result RDS exported through source/target firewalls',
    'Seurat FindMarkers with MAST and logfc.threshold=0; author result field avg_logFC retained directly',
    'no additional expression transform; signed avg_logFC vectors are used as source/target perturbation responses',
    'author script verified FindMarkers+MAST+zero logFC threshold; exported tables retain avg_logFC, pct.1, pct.2, adjusted P',
    note='Controls differ by treatment batch exactly as in the author analysis (Ctrl1 for flg22/Pep1; Ctrl2 for Chitin/SCOOP10/SCREW1).')

# -----------------------------------------------------------------------------
# Known corrected exception: GSE46205.
# -----------------------------------------------------------------------------
corr_summary=json.loads((ROOT/'data/publication_controls_v1/gse46205_preprocessing_correction/SUMMARY.json').read_text(encoding='utf-8'))
add('Salt stress cell-layer microarray','GSE46205','Affymetrix processed GEO sample tables',
    'GCRMA-normalized log2 signal already deposited as VALUE',
    'CORRECTED manuscript analysis uses VALUE directly; response = treated mean − matched control mean on native log2 scale',
    f"post-reveal correction audit present; historical freeze files modified={corr_summary['historical_freeze_files_modified']}",
    status='CORRECTED_POST_REVEAL',
    note='Historical scripts 131–135 applied a second log2. Original artifacts are preserved for provenance; corrected numerical results are not labeled untouched prospective scores.')

res=pd.DataFrame(rows)
res.to_csv(OUT/'main_dataset_preprocessing_semantic_audit.tsv',sep='\t',index=False)
summary={
    'status':'MAIN_DATASET_PREPROCESSING_SEMANTIC_AUDIT_COMPLETE',
    'n_datasets':int(len(res)),
    'n_pass':int((res.status=='PASS').sum()),
    'n_corrected_post_reveal':int((res.status=='CORRECTED_POST_REVEAL').sum()),
    'n_unresolved':int((~res.status.isin(['PASS','CORRECTED_POST_REVEAL'])).sum()),
    'known_exception':'GSE46205 historical double-log implementation; corrected separately and provenance downgraded.',
    'principle':'Manuscript preprocessing descriptions must follow the deposited data semantics, not infer scale from numeric magnitude alone.',
}
(OUT/'SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print('MAIN_DATASET_PREPROCESSING_SEMANTIC_AUDIT_SUCCESS')
print(json.dumps(summary,ensure_ascii=False,indent=2))
print(res[['dataset','accession','status','manuscript_transform']].to_string(index=False).replace('−','-').replace('→','->'))
