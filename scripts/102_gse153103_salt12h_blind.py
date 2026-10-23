from __future__ import annotations

from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, hypergeom

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_salt" / "GSE153103"
WB = BASE / "GSE153103_1_raw_reads_counts.xlsx"
OUT = BASE / "plantvc_v1_1_12h_blind"
RULES = ROOT / "PLANTVC_V1_1_PERSISTENT_UPDATE_2026-10-20.md"
PRE = ["RAS048","RAS049","RAS050","RAS051","RAS056","RAS057","RAS058","RAS059"]
TARGET = ["RAS052","RAS060"]
KS = [100,250,500,1000]


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for c in iter(lambda:f.read(1<<20),b''): h.update(c)
    return h.hexdigest()


def read(cols):
    h=pd.read_excel(WB,sheet_name='raw_counts',nrows=0)
    g=h.columns[0]
    d=pd.read_excel(WB,sheet_name='raw_counts',usecols=[g]+cols).rename(columns={g:'gene'})
    d['gene']=d['gene'].astype(str)
    for c in cols: d[c]=pd.to_numeric(d[c],errors='raise')
    return d


def norm(d, cols):
    x=d[cols].to_numpy(float); lib=x.sum(0); c=x/lib[None,:]*1e6
    return pd.DataFrame(np.log2(c+1),columns=cols), lib


def response(L,t1,t0):
    return ((L[t1[0]]-L[t0[0]])+(L[t1[1]]-L[t0[1]]))/2


def topk(score,signed,target,k):
    n=len(score); k=min(k,n)
    pi=np.argpartition(-score,k-1)[:k]; ti=np.argpartition(-np.abs(target),k-1)[:k]
    inter=np.intersect1d(pi,ti); ov=len(inter); exp=k*k/n
    order=np.argsort(-score); ranks=np.empty(n,int); ranks[order]=np.arange(1,n+1)
    pct=1-(ranks[ti]-1)/max(n-1,1)
    return dict(k=k,overlap=ov,random_expected_overlap=exp,enrichment=ov/exp,
                hypergeom_p=float(hypergeom.sf(ov-1,n,k,k)),
                sign_agreement_overlap=float(np.mean(np.sign(signed[inter])==np.sign(target[inter]))) if ov else np.nan,
                median_target_topk_predicted_percentile=float(np.median(pct)))


def freeze():
    OUT.mkdir(parents=True,exist_ok=True)
    d=read(PRE)
    if any(c in d.columns for c in TARGET): raise RuntimeError('TARGET FIREWALL BREACH')
    L,lib=norm(d,PRE)
    r1=response(L,["RAS049","RAS057"],["RAS048","RAS056"])
    r3=response(L,["RAS050","RAS058"],["RAS048","RAS056"])
    r6=response(L,["RAS051","RAS059"],["RAS048","RAS056"])
    raw=d[PRE].to_numpy(float); cpm=raw/lib[None,:]*1e6
    technical=(cpm>=1).sum(1)>=2
    # confidence only; not used in ranking
    same13=(np.sign(r1)==np.sign(r3))&(np.sign(r1)!=0)
    same36=(np.sign(r3)==np.sign(r6))&(np.sign(r3)!=0)
    r6a=(L['RAS051']-L['RAS048']).to_numpy(); r6b=(L['RAS059']-L['RAS056']).to_numpy()
    rep6=(np.sign(r6a)==np.sign(r6b))&(np.sign(r6a)!=0)
    conf=(same13.astype(float)+same36.astype(float)+rep6.astype(float))/3
    o=pd.DataFrame({'gene':d.gene,'technical_universe':technical,'r1':r1,'r3':r3,'r6':r6,
                    'pred_v1_1':r6,'score_v1_1_abs':np.abs(r6),
                    'pred_baseline_3h':r3,'score_baseline_3h_abs':np.abs(r3),
                    'pred_baseline_1h':r1,'score_baseline_1h_abs':np.abs(r1),
                    'confidence':conf})
    f=OUT/'FROZEN_PRETARGET_RANKING.tsv.gz'; o.to_csv(f,sep='\t',index=False,compression='gzip')
    m={'status':'FROZEN_BEFORE_12H_TARGET','workbook_sha256':sha(WB),'rules_sha256':sha(RULES),
       'ranking_sha256':sha(f),'pretarget_columns_read':PRE,'target_columns_not_read':TARGET,
       'n_genes_all':len(o),'n_genes_technical':int(technical.sum()),
       'primary_rule':'latest same-system response: pred12 = r6','confidence_used_for_ranking':False}
    (OUT/'freeze_manifest.json').write_text(json.dumps(m,indent=2),encoding='utf-8')
    print('PLANTVC_V1_1_12H_FREEZE_SUCCESS'); print(json.dumps(m,indent=2)); print('12H TARGET VALUES HAVE NOT BEEN READ')


def score():
    m=json.loads((OUT/'freeze_manifest.json').read_text(encoding='utf-8'))
    f=OUT/'FROZEN_PRETARGET_RANKING.tsv.gz'
    if sha(f)!=m['ranking_sha256'] or sha(RULES)!=m['rules_sha256']: raise RuntimeError('Frozen artifact hash mismatch')
    fr=pd.read_csv(f,sep='\t')
    t=read(['RAS048','RAS052','RAS056','RAS060']); L,_=norm(t,['RAS048','RAS052','RAS056','RAS060'])
    if t.gene.tolist()!=fr.gene.astype(str).tolist(): raise RuntimeError('gene mismatch')
    a=(L['RAS052']-L['RAS048']).to_numpy(); b=(L['RAS060']-L['RAS056']).to_numpy(); r12=(a+b)/2
    models={'PlantVC_v1_1':('score_v1_1_abs','pred_v1_1'),'baseline_3h':('score_baseline_3h_abs','pred_baseline_3h'),'baseline_1h':('score_baseline_1h_abs','pred_baseline_1h')}
    vr=[]; tr=[]
    for univ,mask in [('technical',fr.technical_universe.to_numpy(bool)),('all',np.ones(len(fr),bool))]:
        y=r12[mask]
        for name,(sc,sg) in models.items():
            scorev=fr.loc[mask,sc].to_numpy(float); signed=fr.loc[mask,sg].to_numpy(float)
            vr.append({'universe':univ,'model':name,'n_genes':int(mask.sum()),'spearman_full_vector':float(spearmanr(signed,y).statistic),'sign_agreement_all':float(np.mean(np.sign(signed)==np.sign(y)))})
            for k in KS: tr.append({'universe':univ,'model':name,**topk(scorev,signed,y,k)})
    rr=[]
    for univ,mask in [('technical',fr.technical_universe.to_numpy(bool)),('all',np.ones(len(fr),bool))]:
        aa=a[mask]; bb=b[mask]
        rr.append({'universe':univ,'replicate_spearman':float(spearmanr(aa,bb).statistic),'replicate_sign_agreement_all':float(np.mean(np.sign(aa)==np.sign(bb)))})
        for k in KS: rr.append({'universe':univ,'comparison':'12h_rep1_vs_rep2',**topk(np.abs(aa),aa,bb,k)})
    pd.DataFrame(vr).to_csv(OUT/'blind_vector_results.tsv',sep='\t',index=False)
    pd.DataFrame(tr).to_csv(OUT/'blind_topk_results.tsv',sep='\t',index=False)
    pd.DataFrame(rr).to_csv(OUT/'target_replicate_reference.tsv',sep='\t',index=False)
    pd.DataFrame({'gene':fr.gene,'r12':r12,'r12_rep1':a,'r12_rep2':b}).to_csv(OUT/'target_12h_response.tsv.gz',sep='\t',index=False,compression='gzip')
    print('PLANTVC_V1_1_12H_SCORE_SUCCESS'); print(pd.DataFrame(vr).to_string(index=False)); print('\nTECHNICAL TOP-K'); print(pd.DataFrame(tr).query("universe=='technical'").to_string(index=False))


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(); p.add_argument('command',choices=['freeze','score']); a=p.parse_args()
    freeze() if a.command=='freeze' else score()
