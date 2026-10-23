from __future__ import annotations
from pathlib import Path
import hashlib, json

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_mechanism'/'GSE130061_SHR_stele'
OUT=BASE/'mechanism_blind'
TARGET=BASE/'target_reveal'/'TPC2019-RA-00455R2_Supplemental_Dataset_2.xlsx'
RULE=ROOT/'SHR机制先验_中柱非移动救援_盲测规则_2026-10-21.md'
TIERA=ROOT/'data'/'prospective_mechanism'/'SHR_direct_prior_GSE21338'/'SHR_tierA_positive_direct_targets.tsv'
TIERB=ROOT/'data'/'prospective_mechanism'/'SHR_direct_prior_GSE21338'/'SHR_tierB_chip_top25.tsv'
SCORE=ROOT/'scripts'/'145_shr_mechanism_score.py'
THIS=ROOT/'scripts'/'144_shr_mechanism_freeze.py'
OUT.mkdir(parents=True,exist_ok=True)

def sha(p:Path)->str:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

for p in [RULE,TIERA,TIERB,SCORE,THIS]:
 if not p.exists(): raise FileNotFoundError(p)
if TARGET.exists(): raise RuntimeError('target numeric workbook is present before freeze; abort')

# Check fixed prior sizes without any target information.
a=[x for x in TIERA.read_text(encoding='utf-8').splitlines()[1:] if x.strip()]
b=[x for x in TIERB.read_text(encoding='utf-8').splitlines()[1:] if x.strip()]
if len(a)!=8: raise RuntimeError(f'Tier A size changed: {len(a)}')
if len(b)!=25: raise RuntimeError(f'Tier B size changed: {len(b)}')

manifest={
 'status':'FROZEN_BEFORE_TARGET_VALUES',
 'task':'mechanism_specialist_SHR_stele_rescue',
 'target_dataset':'GSE130061 Supplemental Dataset 2',
 'target_workbook_present':False,
 'target_numeric_values_read':False,
 'primary_prior':'8 positive direct SHR targets from independent 2006 study',
 'secondary_prior':'top25 SHR ChIP-chip direct-binding candidates from independent GSE21338 / 2011 study',
 'truth_definition':{
  'loss':'mean(WT)-mean(shr-2)',
  'rescue':'mean(SHRdNLELDV)-mean(shr-2)',
  'restoration':'min(loss,rescue)',
 },
 'primary_success_rule':{
  'tierA_min_measured':6,
  'tierA_top10pct_fraction_min':0.5,
  'tierA_top10pct_hypergeom_p_max':0.05,
  'tierA_expected_direction_fraction_min':0.75,
  'scope_if_pass':'local specialist boost only; never replace whole-transcriptome default',
 },
 'secondary_rule':'Tier B enrichment in top10% absolute rescue response; descriptive only',
 'hashes':{
  RULE.name:sha(RULE),
  TIERA.name:sha(TIERA),
  TIERB.name:sha(TIERB),
  THIS.name:sha(THIS),
  SCORE.name:sha(SCORE),
 }
}
(OUT/'FREEZE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('SHR_MECHANISM_FREEZE_SUCCESS')
print(json.dumps(manifest,ensure_ascii=False,indent=2))
