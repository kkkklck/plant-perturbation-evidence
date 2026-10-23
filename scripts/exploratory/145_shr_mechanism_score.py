from __future__ import annotations
from pathlib import Path
import csv, hashlib, json, math, re, zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_mechanism'/'GSE130061_SHR_stele'
OUT=BASE/'mechanism_blind'
TARGET=BASE/'target_reveal'/'TPC2019-RA-00455R2_Supplemental_Dataset_2.xlsx'
MAN=OUT/'FREEZE_MANIFEST.json'
RULE=ROOT/'SHR机制先验_中柱非移动救援_盲测规则_2026-10-21.md'
TIERA=ROOT/'data'/'prospective_mechanism'/'SHR_direct_prior_GSE21338'/'SHR_tierA_positive_direct_targets.tsv'
TIERB=ROOT/'data'/'prospective_mechanism'/'SHR_direct_prior_GSE21338'/'SHR_tierB_chip_top25.tsv'
THIS=ROOT/'scripts'/'145_shr_mechanism_score.py'
FREEZE=ROOT/'scripts'/'144_shr_mechanism_freeze.py'


def sha(p:Path)->str:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()


def col_index(ref:str)->int:
 m=re.match(r'([A-Z]+)',ref)
 s=m.group(1); n=0
 for ch in s: n=n*26+(ord(ch)-64)
 return n-1


def read_dataset2_xlsx(path:Path)->pd.DataFrame:
 # Standard-library XLSX reader: this workbook contains one sheet and cached numeric values.
 ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
     'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
 with zipfile.ZipFile(path) as z:
  shared=[]
  try:
   root=ET.fromstring(z.read('xl/sharedStrings.xml'))
   shared=[''.join(t.text or '' for t in si.iter('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t')) for si in root]
  except KeyError:
   pass
  wb=ET.fromstring(z.read('xl/workbook.xml'))
  rel=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
  rmap={x.attrib['Id']:x.attrib['Target'] for x in rel}
  sh=list(wb.find('m:sheets',ns))[0]
  target=rmap[sh.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']].lstrip('/')
  if not target.startswith('xl/'): target='xl/'+target
  sheet=ET.fromstring(z.read(target))
  rows=[]; header=None
  for row in sheet.findall('.//m:sheetData/m:row',ns):
   rn=int(row.attrib.get('r','0'))
   vals=[None]*8
   for c in row.findall('m:c',ns):
    idx=col_index(c.attrib.get('r','A1'))
    if idx>=8: continue
    typ=c.attrib.get('t'); v=c.find('m:v',ns)
    if typ=='s':
     vals[idx]=shared[int(v.text)] if v is not None else ''
    elif typ=='inlineStr':
     vals[idx]=''.join(t.text or '' for t in c.iter('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t'))
    else:
     vals[idx]=None if v is None else v.text
   if rn==3:
    header=[str(x) if x is not None else '' for x in vals]
    continue
   if rn<4 or vals[0] is None: continue
   gene=str(vals[0]).strip().upper()
   if not re.match(r'^AT[1-5CM]G\d+',gene): continue
   try:
    num=[float(vals[i]) for i in range(1,8)]
   except (TypeError,ValueError):
    continue
   if not np.isfinite(num).all(): continue
   rows.append([gene,*num])
 if header is None: raise RuntimeError('header row 3 not found')
 expected=['AGL number','shr-2 replicate 1','shr-2 replicate 2']
 if header[:3]!=expected:
  raise RuntimeError(f'unexpected header: {header}')
 cols=['gene','shr1','shr2','rescue1','rescue2','wt1','wt2','wt3']
 d=pd.DataFrame(rows,columns=cols).drop_duplicates('gene',keep='first')
 if len(d)<1000: raise RuntimeError(f'too few genes parsed: {len(d)}')
 return d


def read_tsv(path:Path)->pd.DataFrame:
 return pd.read_csv(path,sep='\t')

man=json.loads(MAN.read_text(encoding='utf-8'))
if man.get('status')!='FROZEN_BEFORE_TARGET_VALUES': raise RuntimeError('invalid freeze status')
for p in [RULE,TIERA,TIERB,FREEZE,THIS]:
 if sha(p)!=man['hashes'][p.name]: raise RuntimeError(f'locked file changed: {p.name}')
if not TARGET.exists(): raise FileNotFoundError(TARGET)

d=read_dataset2_xlsx(TARGET)
d['shr_mean']=d[['shr1','shr2']].mean(axis=1)
d['rescue_mean']=d[['rescue1','rescue2']].mean(axis=1)
d['wt_mean']=d[['wt1','wt2','wt3']].mean(axis=1)
d['loss']=d['wt_mean']-d['shr_mean']
d['rescue_delta']=d['rescue_mean']-d['shr_mean']
d['restoration']=np.minimum(d['loss'],d['rescue_delta'])
d['rescue_closer_to_wt']=(d['rescue_mean']-d['wt_mean']).abs() < (d['shr_mean']-d['wt_mean']).abs()
d=d.sort_values('gene').reset_index(drop=True)
N=len(d)

# Deterministic truth rankings.
r=d[['gene','restoration']].sort_values(['restoration','gene'],ascending=[False,True],kind='mergesort').reset_index(drop=True)
r['restoration_rank']=np.arange(1,N+1)
r['restoration_percentile']=1-(r['restoration_rank']-1)/max(1,N-1)
d=d.merge(r[['gene','restoration_rank','restoration_percentile']],on='gene',validate='one_to_one')
a=d[['gene','rescue_delta']].copy(); a['abs_rescue']=a['rescue_delta'].abs(); a=a.sort_values(['abs_rescue','gene'],ascending=[False,True],kind='mergesort').reset_index(drop=True); a['abs_rescue_rank']=np.arange(1,N+1)
d=d.merge(a[['gene','abs_rescue_rank']],on='gene',validate='one_to_one')

K=max(1,int(math.ceil(0.10*N)))
truth_top10=set(d.loc[d.restoration_rank<=K,'gene'])
abs_top10=set(d.loc[d.abs_rescue_rank<=K,'gene'])

pa=read_tsv(TIERA); pa['gene']=pa.gene.str.upper()
pb=read_tsv(TIERB); pb['gene']=pb.gene.str.upper()
U=set(d.gene)
pa_m=pa[pa.gene.isin(U)].copy(); pb_m=pb[pb.gene.isin(U)].copy()

# Tier A primary signed mechanism test.
a_detail=pa_m.merge(d,on='gene',how='left',validate='one_to_one')
a_detail['top10_restoration']=a_detail.gene.isin(truth_top10)
a_detail['expected_direction_ok']=(a_detail.loss>0)&(a_detail.rescue_delta>0)
a_detail['toward_wt']=a_detail.rescue_closer_to_wt
A_m=len(a_detail); A_hit=int(a_detail.top10_restoration.sum()); A_exp=A_m*K/N
A_enrich=A_hit/A_exp if A_exp>0 else np.nan
A_p=float(hypergeom.sf(A_hit-1,N,K,A_m)) if A_m else np.nan
A_dir=float(a_detail.expected_direction_ok.mean()) if A_m else np.nan
A_toward=float(a_detail.toward_wt.mean()) if A_m else np.nan
A_medpct=float(a_detail.restoration_percentile.median()) if A_m else np.nan
crit=man['primary_success_rule']
A_pass=(A_m>=crit['tierA_min_measured'] and A_hit/A_m>=crit['tierA_top10pct_fraction_min'] and A_p<crit['tierA_top10pct_hypergeom_p_max'] and A_dir>=crit['tierA_expected_direction_fraction_min'])

# Tier B binding-only candidate test: no direction assumption.
b_detail=pb_m.merge(d,on='gene',how='left',validate='one_to_one')
b_detail['top10_abs_rescue']=b_detail.gene.isin(abs_top10)
B_m=len(b_detail); B_hit=int(b_detail.top10_abs_rescue.sum()); B_exp=B_m*K/N
B_enrich=B_hit/B_exp if B_exp>0 else np.nan
B_p=float(hypergeom.sf(B_hit-1,N,K,B_m)) if B_m else np.nan
B_rho=float(spearmanr(-b_detail['rank'].to_numpy(float),b_detail['rescue_delta'].abs().to_numpy(float)).statistic) if B_m>2 else np.nan

# Whole-universe descriptive rescue coherence.
whole_rho=float(spearmanr(d.loss,d.rescue_delta).statistic)
positive_loss=d.loss>0
whole_direction=float(np.mean((d.loc[positive_loss,'rescue_delta']>0))) if positive_loss.any() else np.nan

OUT.mkdir(parents=True,exist_ok=True)
a_detail[['gene','name','wt_mean','shr_mean','rescue_mean','loss','rescue_delta','restoration','restoration_rank','restoration_percentile','top10_restoration','expected_direction_ok','toward_wt']].to_csv(OUT/'tierA_gene_results.tsv',sep='\t',index=False)
b_detail[['rank','gene','name','wt_mean','shr_mean','rescue_mean','loss','rescue_delta','restoration','abs_rescue_rank','top10_abs_rescue']].to_csv(OUT/'tierB_gene_results.tsv',sep='\t',index=False)
d[['gene','wt_mean','shr_mean','rescue_mean','loss','rescue_delta','restoration','restoration_rank','restoration_percentile','abs_rescue_rank','rescue_closer_to_wt']].to_csv(OUT/'target_gene_scores.tsv.gz',sep='\t',index=False,compression='gzip')

summary={
 'status':'REVEALED_AND_SCORED',
 'target_workbook_sha256':sha(TARGET),
 'n_genes':N,
 'top10_n':K,
 'whole_loss_vs_rescue_spearman':whole_rho,
 'whole_positive_loss_rescue_same_direction_fraction':whole_direction,
 'tierA':{
  'measured':A_m,'total_prior':len(pa),'top10_hits':A_hit,'random_expected_hits':A_exp,
  'top10_enrichment':A_enrich,'hypergeom_p':A_p,'expected_direction_fraction':A_dir,
  'toward_wt_fraction':A_toward,'median_restoration_percentile':A_medpct,
  'frozen_success_rule_passed':bool(A_pass),
  'decision':'MECHANISM_SUPPORTS_LOCAL_SPECIALIST_BOOST' if A_pass else 'MECHANISM_INSUFFICIENT_TO_OVERRIDE_DEFAULT',
 },
 'tierB':{
  'measured':B_m,'total_prior':len(pb),'top10_abs_rescue_hits':B_hit,'random_expected_hits':B_exp,
  'top10_abs_rescue_enrichment':B_enrich,'hypergeom_p':B_p,
  'chip_rank_vs_abs_rescue_spearman':B_rho,
 },
 'interpretation_scope':'Tier A pass permits local priority boost only for mechanism-supported genes; never whole-transcriptome replacement.'
}
(OUT/'BLIND_RESULT.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print('SHR_MECHANISM_SCORE_SUCCESS')
print(json.dumps(summary,ensure_ascii=False,indent=2))
print('\nTIER_A')
print(a_detail[['gene','name','loss','rescue_delta','restoration','restoration_percentile','top10_restoration','expected_direction_ok','toward_wt']].round(4).to_string(index=False))
print('\nTIER_B')
print(b_detail[['rank','gene','name','rescue_delta','abs_rescue_rank','top10_abs_rescue']].round(4).to_string(index=False))
