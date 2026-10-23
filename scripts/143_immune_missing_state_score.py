from __future__ import annotations
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, hypergeom

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_task'/'immune_missing_state'
TGT=BASE/'target_only'
OUT=BASE/'blind_v0_1'
MAN=OUT/'FREEZE_MANIFEST.json'
STIMS=['flg22','Pep1','Chitin','SCOOP10','SCREW1']
TARGETS=['Rootcap_LRC','Meristem_Cortex','Elongation_Endodermis','Maturation_Stele']
KS=[100,250,500,1000]


def sha256(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            h.update(b)
    return h.hexdigest()


def read_target(p:Path)->pd.Series:
    d=pd.read_csv(p,sep='\t',usecols=['gene','avg_logFC'])
    d['gene']=d['gene'].astype(str)
    d['avg_logFC']=pd.to_numeric(d['avg_logFC'],errors='raise')
    d=d.drop_duplicates('gene',keep='first').set_index('gene')['avg_logFC']
    if not np.isfinite(d.to_numpy(float)).all():
        raise RuntimeError(f'nonfinite target values in {p}')
    return d

man=json.loads(MAN.read_text(encoding='utf-8'))
if man.get('status')!='FROZEN_BEFORE_TARGET_REVEAL': raise RuntimeError('freeze manifest status invalid')
if man.get('task_type')!='missing_state': raise RuntimeError('task type changed')
if man.get('target_states_withheld')!=TARGETS: raise RuntimeError('target list changed')
for nm,h in man['script_sha256'].items():
    p=ROOT/'scripts'/nm
    if sha256(p)!=h: raise RuntimeError(f'script changed after freeze: {nm}')
if not TGT.exists(): raise FileNotFoundError(TGT)

vector_rows=[]; topk_rows=[]; target_hashes={}
for stim in STIMS:
    route=man['routes'][stim]
    pp=ROOT/route['prediction_file']
    if sha256(pp)!=route['prediction_sha256']: raise RuntimeError(f'prediction changed {stim}')
    pred=pd.read_csv(pp,sep='\t')
    genes=pred['gene'].astype(str).tolist()
    N=len(genes)
    pv=pred.set_index('gene')['predicted_response'].reindex(genes).to_numpy(float)
    pred_rank=pred.set_index('gene')['priority_rank'].reindex(genes).to_numpy(int)
    target_hashes[stim]={}
    for target in TARGETS:
        tp=TGT/stim/f'{target}.tsv.gz'
        if not tp.exists(): raise FileNotFoundError(tp)
        target_hashes[stim][target]=sha256(tp)
        ts=read_target(tp)
        present=np.array([g in ts.index for g in genes],dtype=bool)
        tv=np.array([float(ts[g]) if g in ts.index else 0.0 for g in genes],dtype=float)
        coverage=float(present.mean())
        rho_all=float(spearmanr(pv,tv).statistic)
        rho_present=float(spearmanr(pv[present],tv[present]).statistic) if present.sum()>2 else np.nan
        mask_nonzero=(pv!=0)&(tv!=0)
        sign_all=float(np.mean(np.sign(pv[mask_nonzero])==np.sign(tv[mask_nonzero]))) if mask_nonzero.any() else np.nan
        # Truth ranking is over the frozen source-defined universe; target-omitted genes are conservatively zero.
        truth=pd.DataFrame({'gene':genes,'abs_truth':np.abs(tv)})
        truth=truth.sort_values(['abs_truth','gene'],ascending=[False,True],kind='mergesort').reset_index(drop=True)
        truth['truth_rank']=np.arange(1,N+1)
        tr=truth.set_index('gene')['truth_rank'].reindex(genes).to_numpy(int)
        vector_rows.append({
            'stimulus':stim,'target_state':target,'route':route['route'],
            'gene_level_output_allowed':bool(route['gene_level_output_allowed']),
            'n_frozen_genes':N,'target_table_coverage':coverage,
            'spearman_frozen_universe':rho_all,'spearman_target_present_genes':rho_present,
            'sign_agreement_all_nonzero':sign_all,
        })
        for k in KS:
            kk=min(k,N)
            pi=set(np.where(pred_rank<=kk)[0]); ti=set(np.where(tr<=kk)[0]); ov=sorted(pi&ti)
            n=len(ov); exp=kk*kk/N
            topk_rows.append({
                'stimulus':stim,'target_state':target,'route':route['route'],
                'gene_level_output_allowed':bool(route['gene_level_output_allowed']),
                'k':kk,'overlap':n,'random_expected_overlap':exp,'enrichment':n/exp,
                'hypergeom_p':float(hypergeom.sf(n-1,N,kk,kk)),
                'direction_agreement_overlap':float(np.mean(np.sign(pv[ov])==np.sign(tv[ov]))) if ov else np.nan,
            })

vdf=pd.DataFrame(vector_rows); tdf=pd.DataFrame(topk_rows)
vdf.to_csv(OUT/'blind_vector_results.tsv',sep='\t',index=False)
tdf.to_csv(OUT/'blind_topk_results.tsv',sep='\t',index=False)

agg=[]
for stim in STIMS:
    rv=man['routes'][stim]
    qv=vdf[vdf.stimulus==stim]
    qt=tdf[tdf.stimulus==stim]
    row={
        'stimulus':stim,'route':rv['route'],'gene_level_output_allowed':bool(rv['gene_level_output_allowed']),
        'rank1_energy':rv['source_geometry']['rank1_energy'],
        'median_pairwise_rho':rv['source_geometry']['median_pairwise_rho'],
        'effective_state_fraction':rv['source_geometry']['effective_state_fraction'],
        'median_full_rho':float(qv['spearman_frozen_universe'].median()),
        'median_present_gene_rho':float(qv['spearman_target_present_genes'].median()),
        'median_target_coverage':float(qv['target_table_coverage'].median()),
    }
    for k in KS:
        z=qt[qt.k==k]
        row[f'median_top{k}_overlap']=float(z.overlap.median())
        row[f'median_top{k}_enrichment']=float(z.enrichment.median())
        row[f'median_top{k}_direction']=float(z.direction_agreement_overlap.median())
    agg.append(row)
adf=pd.DataFrame(agg)
adf.to_csv(OUT/'blind_stimulus_aggregate.tsv',sep='\t',index=False)

accepted=adf[adf.gene_level_output_allowed]
rejected=adf[~adf.gene_level_output_allowed]
summary={
    'status':'REVEALED_AND_SCORED',
    'rule_version':man['rule_version'],
    'task_type':'missing_state',
    'stimuli_total':len(STIMS),
    'stimuli_allowed_gene_level':accepted.stimulus.tolist(),
    'stimuli_rejected_gene_level':rejected.stimulus.tolist(),
    'n_allowed':int(len(accepted)),'n_rejected':int(len(rejected)),
    'allowed_median_full_rho':float(accepted.median_full_rho.median()) if len(accepted) else None,
    'rejected_median_full_rho':float(rejected.median_full_rho.median()) if len(rejected) else None,
    'allowed_median_top100_enrichment':float(accepted.median_top100_enrichment.median()) if len(accepted) else None,
    'rejected_median_top100_enrichment':float(rejected.median_top100_enrichment.median()) if len(rejected) else None,
    'target_file_sha256':target_hashes,
    'upstream_monolithic_rds_firewall_limitation':True,
    'target_values_used_for_route_or_prediction':False,
}
(OUT/'BLIND_RESULT.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')

print('IMMUNE_MISSING_STATE_SCORE_SUCCESS')
print(json.dumps(summary,ensure_ascii=False,indent=2))
print('\nSTIMULUS_AGGREGATE')
print(adf.round(4).to_string(index=False))
print('\nTARGET_VECTOR')
print(vdf.round(4).to_string(index=False))
print('\nTOP100')
print(tdf[tdf.k==100].round(4).to_string(index=False))
