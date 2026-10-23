from __future__ import annotations

from pathlib import Path
import gzip
import hashlib
import io
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data' / 'prospective_router' / 'Dinneny2008_salt'
CELL_MATRIX = BASE / 'GSE7641_series_matrix.txt.gz'
PRE = BASE / 'pretarget_samples'
OUT = BASE / 'router_blind'
OUT.mkdir(parents=True, exist_ok=True)
RULES = ROOT / 'Dinneny2008_自动路由前瞻盲测规则_2026-10-20.md'
ROUTER_RULES = ROOT / '来源结构路由规则_v0.1_冻结_2026-10-20.md'

CTRL = {
    'epidermis_lrc': ['GSM184889','GSM184890','GSM184891'],
    'columella': ['GSM184892','GSM184893','GSM184894'],
    'cortex': ['GSM184895','GSM184896','GSM184897'],
    'endodermis_qc': ['GSM184898','GSM184899','GSM184900'],
    'stele': ['GSM184901','GSM184902','GSM184903'],
    'protophloem': ['GSM184904','GSM184905','GSM184906'],
}
SALT = {
    'epidermis_lrc': ['GSM184907','GSM184908','GSM184909'],
    'columella': ['GSM184910','GSM184911','GSM184912'],
    'cortex': ['GSM184913','GSM184914','GSM184915'],
    'endodermis_qc': ['GSM184916','GSM184917','GSM184918'],
    'stele': ['GSM184919','GSM184920','GSM184921'],
    'protophloem': ['GSM184922','GSM184923','GSM184924'],
}
WHOLE = {
    'baseline': ['GSM184925','GSM184926'],
    '0p5h': ['GSM184927','GSM184928'],
    '1h': ['GSM184929','GSM184930'],
}
TARGET_GSMS = ['GSM184931','GSM184932']


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_series_matrix(path: Path) -> pd.DataFrame:
    lines=[]; inside=False
    with gzip.open(path,'rt',encoding='utf-8',errors='strict') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                inside=True; continue
            if line.startswith('!series_matrix_table_end'):
                break
            if inside:
                lines.append(line)
    if not lines:
        raise RuntimeError('No series matrix table')
    d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t')
    d=d.rename(columns={d.columns[0]:'probe'})
    d['probe']=d['probe'].astype(str).str.strip().str.strip('"')
    return d


def read_gsm_soft(path: Path) -> pd.DataFrame:
    lines=[];inside=False
    with path.open('r',encoding='utf-8',errors='strict') as f:
        for line in f:
            if line.startswith('!sample_table_begin'):
                inside=True;continue
            if line.startswith('!sample_table_end'):
                break
            if inside:
                lines.append(line)
    if not lines: raise RuntimeError(f'No sample table: {path}')
    d=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE'])
    d=d.rename(columns={'ID_REF':'probe','VALUE':path.stem.split('.')[0]})
    d['probe']=d.probe.astype(str)
    d.iloc[:,1]=pd.to_numeric(d.iloc[:,1],errors='raise')
    return d


def sharedness(M: np.ndarray):
    # input genes x states; transpose for SVD states x genes
    X=np.asarray(M,float).T
    u,s,vt=np.linalg.svd(X,full_matrices=False)
    energy=s*s/(s*s).sum()
    p=u[:,0]**2; p=p/p.sum()
    prs=[];loo=[]
    for i in range(X.shape[0]):
        for j in range(i+1,X.shape[0]):
            prs.append(float(spearmanr(X[i],X[j]).statistic))
        other=X[np.arange(X.shape[0])!=i].mean(axis=0)
        loo.append(float(spearmanr(X[i],other).statistic))
    return {
        'rank1_energy':float(energy[0]),
        'rank1_3_energy':float(energy[:min(3,len(energy))].sum()),
        'median_pairwise_rho':float(np.median(prs)),
        'median_loo_shared_rho':float(np.median(loo)),
        'max_rank1_loading_share':float(p.max()),
        'effective_states':float(1/np.sum(p*p)),
        'effective_state_fraction':float((1/np.sum(p*p))/X.shape[0]),
    }


def main():
    # Hard physical firewall checks.
    if (BASE/'GSE7642_series_matrix.txt.gz').exists() or (BASE/'GSE7642_family.soft.gz').exists():
        raise RuntimeError('TARGET-CONTAINING GSE7642 FULL FILE STILL PRESENT')
    for gsm in TARGET_GSMS:
        if any(BASE.rglob(f'*{gsm}*')):
            raise RuntimeError(f'TARGET SAMPLE FILE PRESENT BEFORE FREEZE: {gsm}')

    cell=read_series_matrix(CELL_MATRIX)
    expected=[x for v in CTRL.values() for x in v]+[x for v in SALT.values() for x in v]
    miss=[c for c in expected if c not in cell.columns]
    if miss: raise RuntimeError(f'Missing cell-layer columns: {miss}')
    for c in expected:
        cell[c]=pd.to_numeric(cell[c],errors='raise')

    # Whole-root pretarget files only.
    whole=None
    pre_hashes={}
    for group in WHOLE.values():
        for gsm in group:
            p=PRE/f'{gsm}.soft.txt'
            if not p.exists(): raise FileNotFoundError(p)
            pre_hashes[gsm]=sha256(p)
            q=read_gsm_soft(p)
            whole=q if whole is None else whole.merge(q,on='probe',how='inner',validate='one_to_one')

    z=cell.merge(whole,on='probe',how='inner',validate='one_to_one')
    z=z[~z.probe.str.upper().str.startswith('AFFX-')].copy()
    cols=expected+[x for v in WHOLE.values() for x in v]
    X=z[cols].to_numpy(float)
    valid=np.all(np.isfinite(X)&(X>0),axis=1)
    z=z.loc[valid].reset_index(drop=True)
    for c in cols:
        z[c]=np.log2(z[c].to_numpy(float))

    states=list(CTRL)
    R=[]
    for st in states:
        r=z[SALT[st]].mean(axis=1).to_numpy(float)-z[CTRL[st]].mean(axis=1).to_numpy(float)
        z[f'resp_{st}']=r;R.append(r)
    R=np.column_stack(R)
    geom=sharedness(R)

    baseline=z[WHOLE['baseline']].mean(axis=1).to_numpy(float)
    r05=z[WHOLE['0p5h']].mean(axis=1).to_numpy(float)-baseline
    r1=z[WHOLE['1h']].mean(axis=1).to_numpy(float)-baseline
    shared=R.mean(axis=1)
    z['whole_0p5h']=r05
    z['persistent_1h']=r1
    z['shared_1h']=shared

    # Frozen route thresholds from v0.1.
    is_shared=(geom['rank1_energy']>=0.50 and geom['median_pairwise_rho']>=0.35 and geom['effective_state_fraction']>=0.50)
    is_state=(geom['rank1_energy']<=0.40 and geom['median_pairwise_rho']<=0.15 and geom['effective_state_fraction']<=0.30)
    if is_shared:
        route='强共享型';primary='shared_1h'
    elif is_state:
        route='强状态特异型';primary='persistent_1h'
    else:
        route='路线歧义';primary=None

    # Freeze deterministic rankings for all allowed predictors, and primary if route is decisive.
    for name in ['shared_1h','persistent_1h','whole_0p5h']:
        z[f'{name}_abs_rank']=pd.Series(np.abs(z[name].to_numpy(float))).rank(method='first',ascending=False).astype(int)
    z['primary_prediction']=z[primary] if primary else np.nan
    z['primary_abs_rank']=z[f'{primary}_abs_rank'] if primary else np.nan

    keep=['probe']+[f'resp_{s}' for s in states]+['shared_1h','persistent_1h','whole_0p5h',
          'shared_1h_abs_rank','persistent_1h_abs_rank','whole_0p5h_abs_rank','primary_prediction','primary_abs_rank']
    pred_path=OUT/'FROZEN_PREDICTIONS.tsv.gz'
    z[keep].to_csv(pred_path,sep='\t',index=False,compression='gzip')

    route_info={
        'status':'FROZEN_BEFORE_4H_TARGET',
        'route':route,
        'primary_predictor':primary,
        'source_geometry':geom,
        'thresholds':{
            'shared':{'rank1_energy_min':0.50,'median_pairwise_rho_min':0.35,'effective_state_fraction_min':0.50},
            'state_specific':{'rank1_energy_max':0.40,'median_pairwise_rho_max':0.15,'effective_state_fraction_max':0.30},
        },
        'n_frozen_probes':int(len(z)),
        'target_4h_samples':['GSM184931','GSM184932'],
        'target_4h_downloaded':False,
        'target_values_read':False,
        'candidate_unit':'ATH1 probe set',
    }
    (OUT/'ROUTE_DECISION.json').write_text(json.dumps(route_info,ensure_ascii=False,indent=2),encoding='utf-8')
    manifest={
        **route_info,
        'rules_sha256':sha256(RULES),'router_rules_sha256':sha256(ROUTER_RULES),
        'cell_source_sha256':sha256(CELL_MATRIX),'pretarget_sample_sha256':pre_hashes,
        'prediction_sha256':sha256(pred_path),'freeze_script_sha256':sha256(Path(__file__)),
    }
    (OUT/'FREEZE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print('DINNENY2008_ROUTER_FREEZE_SUCCESS')
    print(json.dumps(route_info,ensure_ascii=False,indent=2))
    print('4H TARGET FILES NOT PRESENT OR READ')

if __name__=='__main__':
    main()
