from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'publication_controls_v1' / 'missing_state_strong_baselines'
OUT.mkdir(parents=True, exist_ok=True)
K=100


def ranks_abs(v):
    order=np.argsort(-np.abs(v),kind='mergesort')
    r=np.empty(len(v),dtype=int); r[order]=np.arange(1,len(v)+1); return r


def score(dataset, context, model, pred, truth):
    n=len(truth); pr=ranks_abs(pred); tr=ranks_abs(truth); kk=min(K,n)
    ov=np.intersect1d(np.where(pr<=kk)[0],np.where(tr<=kk)[0])
    exp=kk*kk/n
    return {'dataset':dataset,'context':context,'model':model,'n_genes':n,
            'spearman_full_vector':float(spearmanr(pred,truth).statistic),
            'top100_overlap':int(len(ov)),'top100_enrichment':float(len(ov)/exp),
            'top100_direction_agreement':float(np.mean(np.sign(pred[ov])==np.sign(truth[ov]))) if len(ov) else np.nan}


def common_axis(M):
    # M genes x source states. Orient rank-1 gene axis toward the signed source mean using source-only information.
    u,s,vt=np.linalg.svd(M,full_matrices=False)
    axis=u[:,0]*s[0]
    mean=M.mean(axis=1)
    if np.dot(axis,mean)<0: axis=-axis
    return axis


def source_models(M, names):
    mean=M.mean(axis=1)
    median=np.median(M,axis=1)
    axis=common_axis(M)
    rms=np.sqrt(np.mean(M*M,axis=0))
    strongest=int(np.argmax(rms))
    return {
        'signed_mean':mean,
        'signed_median':median,
        'rank1_common_axis':axis,
        'strongest_source_pretarget':M[:,strongest],
    }, names[strongest], rms

rows=[]; source_choice=[]

# -----------------------------------------------------------------------------
# GSE121619 heat shock: strict target-state holdout, fixed technical universe.
# -----------------------------------------------------------------------------
base=ROOT/'data'/'external_heatshock'/'GSE121619'/'plantvc_v1_heatshock_blind'
pred=pd.read_csv(base/'frozen_prediction.tsv.gz',sep='\t')
pred=pred[pred['technical'].astype(bool)].copy()
pred['gene']=pred['gene'].astype(str)
genes=pred['gene'].tolist()
M=np.column_stack([pred[c].to_numpy(float) for c in ['donor_5','donor_6','donor_8']])
models,strongest,rms=source_models(M,['donor_5','donor_6','donor_8'])
models['pooled_donor']=pred['pooled_donor_prediction'].to_numpy(float)
source_choice.append({'dataset':'GSE121619_heatshock','stimulus':'heatshock','selected_source':strongest,
                      'source_names':'donor_5;donor_6;donor_8','source_rms':';'.join(f'{x:.8g}' for x in rms)})
target=pd.read_csv(base/'target_profiles.tsv.gz',sep='\t')
target=target[target['technical'].astype(bool)].copy(); target['gene']=target['gene'].astype(str)
for cl in [2,3,10]:
    q=target[target.target_cluster==cl].set_index('gene').reindex(genes)
    if q['target_response'].isna().any(): raise RuntimeError(f'missing target genes cluster {cl}')
    truth=q['target_response'].to_numpy(float)
    for nm,p in models.items(): rows.append(score('GSE121619_heatshock',f'cluster_{cl}',nm,p,truth))

# -----------------------------------------------------------------------------
# Five-stimulus immune benchmark: source-only common gene universe per stimulus.
# -----------------------------------------------------------------------------
BASE=ROOT/'data'/'prospective_task'/'immune_missing_state'
SRC=BASE/'source_only'; TGT=BASE/'target_only'; BLIND=BASE/'blind_v0_1'
STIMS=['flg22','Pep1','Chitin','SCOOP10','SCREW1']
SOURCE_STATES=['Rootcap_Columella','Meristem_Epidermis','Meristem_Endodermis','Meristem_Stele',
               'Elongation_Epidermis','Elongation_Cortex','Elongation_Stele',
               'Maturation_Epidermis','Maturation_Cortex','Maturation_Endodermis']
TARGET_STATES=['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']
route=pd.read_csv(BLIND/'blind_stimulus_aggregate.tsv',sep='\t').set_index('stimulus')
for stim in STIMS:
    ser={}
    for st in SOURCE_STATES:
        d=pd.read_csv(SRC/stim/f'{st}.tsv.gz',sep='\t',usecols=['gene','avg_logFC']).drop_duplicates('gene')
        ser[st]=d.set_index('gene')['avg_logFC'].astype(float)
    genes=set(ser[SOURCE_STATES[0]].index)
    for st in SOURCE_STATES[1:]: genes &= set(ser[st].index)
    genes=sorted(genes)
    M=np.column_stack([ser[st].reindex(genes).to_numpy(float) for st in SOURCE_STATES])
    models,strongest,rms=source_models(M,SOURCE_STATES)
    source_choice.append({'dataset':'immune_missing_state','stimulus':stim,'selected_source':strongest,
                          'source_names':';'.join(SOURCE_STATES),'source_rms':';'.join(f'{x:.8g}' for x in rms)})
    allowed=bool(route.loc[stim,'gene_level_output_allowed'])
    for tgt in TARGET_STATES:
        d=pd.read_csv(TGT/stim/f'{tgt}.tsv.gz',sep='\t',usecols=['gene','avg_logFC']).drop_duplicates('gene').set_index('gene')['avg_logFC'].astype(float)
        truth=np.array([float(d[g]) if g in d.index else 0.0 for g in genes],dtype=float)
        for nm,p in models.items():
            row=score('immune_missing_state',f'{stim}|{tgt}',nm,p,truth)
            row['stimulus']=stim; row['target_state']=tgt; row['gate_allowed']=allowed
            rows.append(row)

res=pd.DataFrame(rows)
res.to_csv(OUT/'missing_state_baseline_results.tsv',sep='\t',index=False)
pd.DataFrame(source_choice).to_csv(OUT/'pretarget_strongest_source_choices.tsv',sep='\t',index=False)

# Aggregates: heatshock across 3 targets; immune across 20 targets and separately allowed/rejected.
agg=[]
for (ds,model),q in res.groupby(['dataset','model']):
    agg.append({'dataset':ds,'subset':'all','model':model,'n_contexts':len(q),
                'median_rho':float(q.spearman_full_vector.median()),'mean_rho':float(q.spearman_full_vector.mean()),
                'median_top100_overlap':float(q.top100_overlap.median()),'mean_top100_overlap':float(q.top100_overlap.mean()),
                'median_top100_enrichment':float(q.top100_enrichment.median()),
                'median_direction':float(q.top100_direction_agreement.median())})
if 'gate_allowed' in res.columns:
    iq=res[res.dataset=='immune_missing_state']
    for allowed,label in [(True,'gate_allowed'),(False,'gate_rejected')]:
        z=iq[iq.gate_allowed==allowed]
        for model,q in z.groupby('model'):
            agg.append({'dataset':'immune_missing_state','subset':label,'model':model,'n_contexts':len(q),
                        'median_rho':float(q.spearman_full_vector.median()),'mean_rho':float(q.spearman_full_vector.mean()),
                        'median_top100_overlap':float(q.top100_overlap.median()),'mean_top100_overlap':float(q.top100_overlap.mean()),
                        'median_top100_enrichment':float(q.top100_enrichment.median()),
                        'median_direction':float(q.top100_direction_agreement.median())})
agg=pd.DataFrame(agg)
agg.to_csv(OUT/'missing_state_baseline_aggregate.tsv',sep='\t',index=False)

print('MISSING_STATE_STRONG_BASELINES_SUCCESS')
print('\nAGGREGATE')
print(agg.round(4).to_string(index=False))
print('\nPRETARGET_STRONGEST_SOURCE')
print(pd.DataFrame(source_choice).to_string(index=False))
