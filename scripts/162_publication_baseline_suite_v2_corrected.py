from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
V1 = ROOT/'data'/'publication_baseline_suite_v1'
CORR = ROOT/'data'/'publication_controls_v1'/'gse46205_preprocessing_correction'
OUT = ROOT/'data'/'publication_baseline_suite_v2_corrected'
OUT.mkdir(parents=True,exist_ok=True)

long = pd.read_csv(V1/'baseline_suite_long.tsv',sep='\t')
# Remove all historical GSE46205 rows affected by the double-log implementation.
long = long[long.dataset != 'GSE46205_salt_8h'].copy()

c = pd.read_csv(CORR/'corrected_8h_strong_baselines.tsv',sep='\t')
rows=[]
for _,r in c.iterrows():
    rows.append({
        'task_family':'future_time',
        'dataset':'GSE46205_salt_8h',
        'context':r['celltype'],
        'model':r['model'],
        'model_class':'strong_baseline_corrected',
        'pretarget_eligible':True,
        'n_genes':22810,
        'spearman_full_vector':r['spearman_full_vector'],
        'top100_overlap':r['top100_overlap'],
        'top100_enrichment':r['top100_enrichment'],
        'top100_direction_agreement':r['top100_direction'],
        'source_file':'data/publication_controls_v1/gse46205_preprocessing_correction/corrected_8h_strong_baselines.tsv',
        'note':'Post-reveal preprocessing correction: deposited GCRMA VALUE used directly because it is already log2 normalized. Task split and model definitions were unchanged, but corrected numerical scores are not labeled untouched prospective results.',
    })
    rows.append({
        'task_family':'future_time','dataset':'GSE46205_salt_8h','context':r['celltype'],
        'model':'random_rank_expectation','model_class':'random_expectation','pretarget_eligible':True,
        'n_genes':22810,'spearman_full_vector':0.0,'top100_overlap':10000/22810,
        'top100_enrichment':1.0,'top100_direction_agreement':np.nan,
        'source_file':'data/publication_controls_v1/gse46205_preprocessing_correction/corrected_8h_strong_baselines.tsv',
        'note':'Analytic random-ranking expectation; native-GCRMA corrected GSE46205 universe.',
    })
add=pd.DataFrame(rows)
long=pd.concat([long,add],ignore_index=True)
long.to_csv(OUT/'baseline_suite_long.tsv',sep='\t',index=False)

agg=(long.groupby(['task_family','dataset','model','model_class','pretarget_eligible'],as_index=False)
     .agg(n_contexts=('context','nunique'),median_rho=('spearman_full_vector','median'),mean_rho=('spearman_full_vector','mean'),
          median_top100_overlap=('top100_overlap','median'),mean_top100_overlap=('top100_overlap','mean'),
          median_top100_enrichment=('top100_enrichment','median'),median_top100_direction=('top100_direction_agreement','median')))
agg.to_csv(OUT/'baseline_suite_aggregate.tsv',sep='\t',index=False)

# Carry negative controls unchanged: none depends on GSE46205.
neg=pd.read_csv(V1/'negative_controls_cross_experiment.tsv',sep='\t')
neg.to_csv(OUT/'negative_controls_cross_experiment.tsv',sep='\t',index=False)

manuscript=[]
elig=agg[(agg.pretarget_eligible.astype(bool)) & (agg.model!='random_rank_expectation')].copy()
for (task,ds),g in elig.groupby(['task_family','dataset'],sort=False):
    required=int(long[(long.task_family==task)&(long.dataset==ds)].context.nunique())
    full=g[g.n_contexts==required].copy()
    best=full.sort_values(['median_rho','median_top100_overlap'],ascending=False).iloc[0]
    manuscript.append({'task_family':task,'dataset':ds,'best_model_by_median_rho':best.model,
                       'median_rho':float(best.median_rho),'median_top100_overlap':float(best.median_top100_overlap),
                       'n_contexts':int(best.n_contexts),'required_contexts':required})
pd.DataFrame(manuscript).to_csv(OUT/'manuscript_baseline_winners.tsv',sep='\t',index=False)

manifest={
    'suite_version':'publication_baseline_suite_v2_corrected',
    'supersedes':'publication_baseline_suite_v1 for manuscript use',
    'correction':'All GSE46205 rows replaced with native deposited GCRMA log2-scale calculations from script 160; other datasets copied from v1 unchanged.',
    'gse46205_prospective_label':'Corrected numerical results are post-reveal implementation audits and must not be described as untouched prospective scores.',
    'n_rows_long':int(len(long)),'n_rows_aggregate':int(len(agg)),
}
(OUT/'MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('PUBLICATION_BASELINE_SUITE_V2_CORRECTED_SUCCESS')
print(json.dumps(manifest,ensure_ascii=False,indent=2))
print(pd.DataFrame(manuscript).round(4).to_string(index=False))
