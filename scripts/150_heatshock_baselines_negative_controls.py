from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data'/'publication_baselines_v1'/'heatshock'
OUT.mkdir(parents=True,exist_ok=True)
KS=[100,250,500,1000]
SEED=20261021
B_RANDOM=10000


def common_pc1(donor_matrix, signed_mean):
    M=np.asarray(donor_matrix,float)
    u,s,vt=np.linalg.svd(M,full_matrices=False)
    pc=u[:,0]*s[0]
    if spearmanr(pc,signed_mean).statistic < 0:
        pc=-pc
    return pc


def score(dataset, context, model, pred, truth, rows_v, rows_t):
    pred=np.asarray(pred,float); truth=np.asarray(truth,float); n=len(pred)
    mask=(pred!=0)&(truth!=0)&np.isfinite(pred)&np.isfinite(truth)
    rows_v.append({'dataset':dataset,'context':str(context),'model':model,'n_genes':n,
                   'spearman_full_vector':float(spearmanr(pred,truth).statistic),
                   'sign_agreement_nonzero':float(np.mean(np.sign(pred[mask])==np.sign(truth[mask]))) if mask.any() else np.nan})
    pr=pd.Series(np.abs(pred)).rank(method='first',ascending=False).to_numpy(int)
    tr=pd.Series(np.abs(truth)).rank(method='first',ascending=False).to_numpy(int)
    for k in KS:
        kk=min(k,n); ov=sorted(set(np.where(pr<=kk)[0])&set(np.where(tr<=kk)[0])); exp=kk*kk/n
        rows_t.append({'dataset':dataset,'context':str(context),'model':model,'k':kk,'overlap':len(ov),'random_expected_overlap':exp,
                       'enrichment':len(ov)/exp,'direction_agreement_overlap':float(np.mean(np.sign(pred[ov])==np.sign(truth[ov]))) if ov else np.nan})


def random_label_top100_null(pred, truth, rng, B):
    n=len(pred); k=min(100,n)
    top_pred=np.where(pd.Series(np.abs(pred)).rank(method='first',ascending=False).to_numpy(int)<=k)[0]
    top_truth=set(np.where(pd.Series(np.abs(truth)).rank(method='first',ascending=False).to_numpy(int)<=k)[0])
    obs=len(set(top_pred)&top_truth)
    null=np.empty(B,int)
    for b in range(B):
        idx=rng.choice(n,size=k,replace=False)
        null[b]=len(set(idx.tolist())&top_truth)
    return {'observed_overlap':obs,'null_mean':float(null.mean()),'null_q975':float(np.quantile(null,.975)),
            'empirical_p':float((1+np.sum(null>=obs))/(B+1))}

rows_v=[];rows_t=[];null_rows=[]
rng=np.random.default_rng(SEED)

# Salt 3h biological negative-control signature.
salt=pd.read_csv(ROOT/'data'/'external_salt'/'GSE153103'/'plantvc_v1_6h_blind'/'FROZEN_PRETARGET_RANKING.tsv.gz',sep='\t')
salt=salt[salt.technical_universe][['gene','r3']].copy(); salt['gene']=salt.gene.astype(str).str.upper(); salt=salt.drop_duplicates('gene').set_index('gene')

# -------- GSE121619 target-state holdout --------
B=ROOT/'data'/'external_heatshock'/'GSE121619'/'plantvc_v1_heatshock_blind'
pred=pd.read_csv(B/'frozen_prediction.tsv.gz',sep='\t')
pred['gene']=pred.gene.astype(str).str.upper(); pred=pred[pred.technical].copy()
donors=pred[['donor_5','donor_6','donor_8']].to_numpy(float)
shared=pred.shared_prediction.to_numpy(float)
pred['pc1_common']=common_pc1(donors,shared)
pred['mean_abs_signedmean']=np.sign(shared)*np.mean(np.abs(donors),axis=1)
pred['salt3h_wrong_stress']=pred.gene.map(salt['r3']).fillna(0.0)

targ=pd.read_csv(B/'target_profiles.tsv.gz',sep='\t'); targ['gene']=targ.gene.astype(str).str.upper(); targ=targ[targ.technical]
for cl in [2,3,10]:
    t=targ[targ.target_cluster==cl][['gene','target_response']]
    z=pred.merge(t,on='gene',validate='one_to_one')
    models=['shared_prediction','pooled_donor_prediction','pc1_common','mean_abs_signedmean','donor_5','donor_6','donor_8','salt3h_wrong_stress']
    for m in models:
        score('GSE121619_heatshock_missing_state',cl,m,z[m],z.target_response,rows_v,rows_t)
    nr=random_label_top100_null(z.shared_prediction.to_numpy(float),z.target_response.to_numpy(float),rng,B_RANDOM)
    null_rows.append({'dataset':'GSE121619_heatshock_missing_state','context':str(cl),'negative_control':'random_gene_labels_top100','B':B_RANDOM,**nr})

# -------- GSE290176 independent cross-study target --------
B=ROOT/'data'/'external_heatshock_crossstudy'/'GSE290176'
P=B/'crossstudy_blind'/'frozen_crossstudy_prediction.tsv.gz'
pred=pd.read_csv(P,sep='\t'); pred['gene']=pred.gene.astype(str).str.upper()
donors=pred[['donor_5','donor_6','donor_8']].to_numpy(float); shared=pred.shared_prediction.to_numpy(float)
pred['pc1_common']=common_pc1(donors,shared)
pred['mean_abs_signedmean']=np.sign(shared)*np.mean(np.abs(donors),axis=1)
pred['salt3h_wrong_stress']=pred.gene.map(salt['r3']).fillna(0.0)

counts=pd.read_csv(B/'GSE290176_counts_HS_Ath.txt.gz',sep='\t')
CK=[f'Col0-CK{i}' for i in range(1,7)]; HS=[f'Col0-HS{i}' for i in range(1,4)]
counts=counts[['Geneid']+CK+HS].copy(); counts['gene']=counts.Geneid.astype(str).str.strip().str.strip('"').str.upper()
for c in CK+HS: counts[c]=pd.to_numeric(counts[c],errors='coerce').fillna(0.0)
counts=counts.groupby('gene',as_index=False)[CK+HS].sum(); libs=counts[CK+HS].sum(axis=0)
for c in CK+HS: counts[c]=np.log2(counts[c].to_numpy(float)/float(libs[c])*1e6+1)
counts['target_response']=counts[HS].mean(axis=1)-counts[CK].mean(axis=1)
z=pred.merge(counts[['gene','target_response']],on='gene',validate='one_to_one')
models=['shared_prediction','pooled_donor_prediction','pc1_common','mean_abs_signedmean','donor_5','donor_6','donor_8','salt3h_wrong_stress']
for m in models: score('GSE290176_heatshock_crossstudy','whole_seedling',m,z[m],z.target_response,rows_v,rows_t)
nr=random_label_top100_null(z.shared_prediction.to_numpy(float),z.target_response.to_numpy(float),rng,B_RANDOM)
null_rows.append({'dataset':'GSE290176_heatshock_crossstudy','context':'whole_seedling','negative_control':'random_gene_labels_top100','B':B_RANDOM,**nr})

vdf=pd.DataFrame(rows_v);tdf=pd.DataFrame(rows_t);ndf=pd.DataFrame(null_rows)
vdf.to_csv(OUT/'heatshock_baseline_vector.tsv',sep='\t',index=False)
tdf.to_csv(OUT/'heatshock_baseline_topk.tsv',sep='\t',index=False)
ndf.to_csv(OUT/'heatshock_random_label_null.tsv',sep='\t',index=False)

aggv=vdf.groupby(['dataset','model'],as_index=False).agg(n_contexts=('context','nunique'),median_rho=('spearman_full_vector','median'),median_sign=('sign_agreement_nonzero','median'))
aggt=tdf[tdf.k==100].groupby(['dataset','model'],as_index=False).agg(median_top100_overlap=('overlap','median'),median_top100_enrichment=('enrichment','median'),median_top100_direction=('direction_agreement_overlap','median'))
agg=aggv.merge(aggt,on=['dataset','model'])
agg.to_csv(OUT/'heatshock_baseline_summary.tsv',sep='\t',index=False)

manifest={'status':'POSTHOC_STRONG_BASELINE_AND_NEGATIVE_CONTROL_AUDIT','seed':SEED,'random_gene_label_draws':B_RANDOM,
          'wrong_stress_negative_control':'GSE153103 Arabidopsis salt 3h response, aligned by gene ID; missing genes set to zero',
          'target_values_used_to_define_baselines':False,
          'note':'PC1, mean-absolute donor magnitude, individual donors and salt mismatch were added after target reveal and are not relabeled as prospective.'}
(OUT/'MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('HEATSHOCK_BASELINES_NEGATIVE_CONTROLS_SUCCESS')
print('\nSUMMARY');print(agg.round(4).to_string(index=False));print('\nRANDOM LABEL NULL');print(ndf.round(4).to_string(index=False))
