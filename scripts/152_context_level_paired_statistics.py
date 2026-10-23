from __future__ import annotations

from pathlib import Path
from itertools import product
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data'/'publication_stats_v1'/'context_level'
OUT.mkdir(parents=True,exist_ok=True)


def exact_signflip_p(diffs, alternative='greater'):
    d=np.asarray(diffs,float)
    n=len(d)
    obs=float(np.mean(d))
    vals=[]
    for signs in product([-1.0,1.0], repeat=n):
        vals.append(float(np.mean(d*np.asarray(signs))))
    vals=np.asarray(vals)
    if alternative=='greater':
        p=float(np.mean(vals>=obs-1e-12))
    else:
        p=float(np.mean(np.abs(vals)>=abs(obs)-1e-12))
    return obs,p,len(vals)

rows=[]; detail=[]

# Lateral-root 20h: same-state vs frozen broad shared across 9 cell types.
base=ROOT/'data'/'external_lateral_root'/'GSE158761'/'plantvc_v1_blind'
v=pd.read_csv(base/'blind_vector_results.tsv',sep='\t')
t=pd.read_csv(base/'blind_topk_results.tsv',sep='\t')
wide=v.pivot(index='celltype',columns='model',values='spearman_full_vector')
drho=wide['PlantVC_persistent_same_state']-wide['frozen_broad_shared']
obs,p,nperm=exact_signflip_p(drho)
rows.append({'dataset':'GSE158761_lateral_root_20h','metric':'full_vector_rho','n_contexts':len(drho),'mean_paired_difference':obs,
             'median_paired_difference':float(np.median(drho)),'wins':int(np.sum(drho>0)),'ties':int(np.sum(drho==0)),'losses':int(np.sum(drho<0)),
             'exact_signflip_p_one_sided':p,'n_permutations':nperm})
for ct,x in drho.items(): detail.append({'dataset':'GSE158761_lateral_root_20h','context':ct,'metric':'full_vector_rho','paired_difference':x})

t100=t[t.k==100].pivot(index='celltype',columns='model',values='overlap')
d=t100['PlantVC_persistent_same_state']-t100['frozen_broad_shared']
obs,p,nperm=exact_signflip_p(d)
rows.append({'dataset':'GSE158761_lateral_root_20h','metric':'top100_overlap','n_contexts':len(d),'mean_paired_difference':obs,
             'median_paired_difference':float(np.median(d)),'wins':int(np.sum(d>0)),'ties':int(np.sum(d==0)),'losses':int(np.sum(d<0)),
             'exact_signflip_p_one_sided':p,'n_permutations':nperm})
for ct,x in d.items(): detail.append({'dataset':'GSE158761_lateral_root_20h','context':ct,'metric':'top100_overlap','paired_difference':x})

# GSE46205 8h: same-state 3h vs shared-other3 3h across 4 cell layers.
base=ROOT/'data'/'prospective_router'/'GSE46205_salt'/'task_8h_blind'
v=pd.read_csv(base/'blind_8h_vector_results.tsv',sep='\t')
t=pd.read_csv(base/'blind_8h_topk_results.tsv',sep='\t')
wide=v.pivot(index='celltype',columns='model',values='spearman_full_vector')
drho=wide['same_state_3h']-wide['shared_other3_3h']
obs,p,nperm=exact_signflip_p(drho)
rows.append({'dataset':'GSE46205_salt_8h','metric':'full_vector_rho','n_contexts':len(drho),'mean_paired_difference':obs,
             'median_paired_difference':float(np.median(drho)),'wins':int(np.sum(drho>0)),'ties':int(np.sum(drho==0)),'losses':int(np.sum(drho<0)),
             'exact_signflip_p_one_sided':p,'n_permutations':nperm})
for ct,x in drho.items(): detail.append({'dataset':'GSE46205_salt_8h','context':ct,'metric':'full_vector_rho','paired_difference':x})

t100=t[t.k==100].pivot(index='celltype',columns='model',values='overlap')
d=t100['same_state_3h']-t100['shared_other3_3h']
obs,p,nperm=exact_signflip_p(d)
rows.append({'dataset':'GSE46205_salt_8h','metric':'top100_overlap','n_contexts':len(d),'mean_paired_difference':obs,
             'median_paired_difference':float(np.median(d)),'wins':int(np.sum(d>0)),'ties':int(np.sum(d==0)),'losses':int(np.sum(d<0)),
             'exact_signflip_p_one_sided':p,'n_permutations':nperm})
for ct,x in d.items(): detail.append({'dataset':'GSE46205_salt_8h','context':ct,'metric':'top100_overlap','paired_difference':x})

pd.DataFrame(rows).to_csv(OUT/'paired_context_inference.tsv',sep='\t',index=False)
pd.DataFrame(detail).to_csv(OUT/'paired_context_differences.tsv',sep='\t',index=False)
manifest={'unit_of_inference':'cell type/state context, not genes','test':'exact paired sign-flip permutation on mean difference','alternative':'same-state history > broad/shared comparator'}
(OUT/'MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('CONTEXT_LEVEL_PAIRED_STATS_SUCCESS')
print(pd.DataFrame(rows).round(5).to_string(index=False))
