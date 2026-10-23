from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, binomtest, beta

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data'/'prospective_mechanism'/'GSE130061_SHR_stele'/'mechanism_blind'
OUT=ROOT/'data'/'publication_stats_v1'/'shr_mechanism'
OUT.mkdir(parents=True,exist_ok=True)

def clopper(k,n,alpha=.05):
    lo=0.0 if k==0 else float(beta.ppf(alpha/2,k,n-k+1))
    hi=1.0 if k==n else float(beta.ppf(1-alpha/2,k+1,n-k))
    return lo,hi

scores=pd.read_csv(BASE/'target_gene_scores.tsv.gz',sep='\t')
a=pd.read_csv(BASE/'tierA_gene_results.tsv',sep='\t')
b=pd.read_csv(BASE/'tierB_gene_results.tsv',sep='\t')
N=len(scores)
# Target top 10% according to rank definition already used in frozen scoring.
top_n=int(np.floor(0.10*N))

def enrichment_test(df, hit_col, label):
    k=int(df[hit_col].sum()); n=len(df)
    # finite-population contingency table
    table=[[k,n-k],[top_n-k,(N-top_n)-(n-k)]]
    odds,p=fisher_exact(table,alternative='greater')
    lo,hi=clopper(k,n)
    bt=binomtest(k,n,p=0.10,alternative='greater')
    return {'panel':label,'selected_n':n,'top10_hits':k,'top10_random_expectation':0.10*n,'fold_enrichment':k/(0.10*n),
            'hit_fraction':k/n,'hit_fraction_CI95_low':lo,'hit_fraction_CI95_high':hi,
            'fisher_exact_oddsratio':float(odds),'fisher_exact_p_one_sided':float(p),
            'binomial_p_vs_0p10_one_sided':float(bt.pvalue)}

rows=[]
rows.append(enrichment_test(a,'top10_restoration','tierA_directional_functional'))
# Tier B uses rescue absolute top10 result; infer top10 from saved score ranking if result table has indicator, else reconstruct.
if 'top10_abs_rescue' in b.columns:
    rows.append(enrichment_test(b,'top10_abs_rescue','tierB_chip_binding'))
else:
    # Existing result table records gene and |rescue| percentile/rank; find a boolean-compatible column.
    cand=[c for c in b.columns if 'top10' in c.lower()]
    if cand:
        rows.append(enrichment_test(b,cand[0],'tierB_chip_binding'))
    else:
        # Reconstruct from target_gene_scores on the exact Tier B genes.
        gene_col='gene'
        score_gene=scores.set_index(gene_col)
        bb=b.copy()
        if 'gene' not in bb.columns:
            raise RuntimeError('Tier B gene column missing')
        # rank absolute rescue across all target genes
        if 'rescue_delta' in score_gene.columns:
            rank=score_gene['rescue_delta'].abs().rank(method='min',ascending=False)
        elif 'rescue' in score_gene.columns:
            rank=score_gene['rescue'].abs().rank(method='min',ascending=False)
        else:
            raise RuntimeError('No rescue column in target_gene_scores')
        bb['top10_reconstructed']=bb['gene'].map(rank)<=top_n
        rows.append(enrichment_test(bb,'top10_reconstructed','tierB_chip_binding'))

# Direction and toward-WT proportions for Tier A.
extra=[]
for col,label,nullp in [('expected_direction_ok','expected_direction',0.5),('toward_wt','toward_WT',0.5)]:
    k=int(a[col].sum());n=len(a);lo,hi=clopper(k,n);bt=binomtest(k,n,p=nullp,alternative='greater')
    extra.append({'property':label,'n':n,'successes':k,'fraction':k/n,'CI95_low':lo,'CI95_high':hi,
                  'binomial_null_p':nullp,'binomial_p_one_sided':float(bt.pvalue)})

pd.DataFrame(rows).to_csv(OUT/'mechanism_enrichment_exact_stats.tsv',sep='\t',index=False)
pd.DataFrame(extra).to_csv(OUT/'tierA_direction_recovery_exact_stats.tsv',sep='\t',index=False)
manifest={'independent_unit':'predefined mechanism-supported genes','target_universe_N':N,'top10_count':top_n,
          'note':'Tier A was frozen before target expression; this reanalysis adds exact confidence intervals and Fisher/binomial tests post reveal.'}
(OUT/'MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('SHR_PUBLICATION_STATS_SUCCESS')
print(pd.DataFrame(rows).round(5).to_string(index=False));print('\nDIRECTION/RECOVERY');print(pd.DataFrame(extra).round(5).to_string(index=False))
