from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[1]
ROB = ROOT / "data" / "external_rice" / "Robertson2026"
DESIGN = pd.read_csv(ROB / "sample_design.tsv", sep="\t")
OUT = ROB / "candidate_recovery_v0"
OUT.mkdir(parents=True, exist_ok=True)

sample_batch = dict(zip(DESIGN.identity.astype(str), DESIGN.batch.astype(int)))
conditions = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
control = "Well.Watered"
ks = [100, 250, 500, 1000]


def load_level(level):
    if level == "CellAnnotation":
        base = ROB / "pseudobulk_v0"; state_col = "cell_annotation"; f = "pseudobulk_log1p_cpm.tsv.gz"
    else:
        base = ROB / "tissuesystem_pseudobulk_v0"; state_col = "tissue_system"; f = "tissuesystem_log1p_cpm.tsv.gz"
    meta = pd.read_csv(base / "group_metadata.tsv", sep="\t")
    expr = pd.read_csv(base / f, sep="\t")
    genes = expr.pop("gene").astype(str).to_numpy(); X = expr.to_numpy(float)
    mask = ((X > 0).mean(axis=1) >= 0.25) & (X.mean(axis=1) >= 0.05)
    return meta, genes[mask], X[mask], state_col


def get_vectors(meta, X, state_col, cond, state):
    ctrl_by_batch = {}
    for s in meta.loc[meta.condition.eq(control), "sample"].astype(str).unique():
        ctrl_by_batch[sample_batch[s]] = s
    out=[]; labels=[]
    for ti in meta.index[(meta.condition.eq(cond)) & (meta[state_col].astype(str).eq(state))].to_numpy():
        ss=str(meta.loc[ti,"sample"]); b=sample_batch[ss]; cs=ctrl_by_batch[b]
        ci=meta.index[(meta["sample"].astype(str).eq(cs)) & (meta[state_col].astype(str).eq(state))].to_numpy()
        if len(ci)!=1: continue
        out.append(X[:,ti]-X[:,ci[0]]); labels.append((ss,cs,b))
    return out, labels

metric_rows=[]; gene_rows=[]
for level in ["CellAnnotation","TissueSystem"]:
    meta, genes, X, state_col = load_level(level); G=len(genes)
    for cond in conditions:
        for state in sorted(meta[state_col].astype(str).unique()):
            vecs, labels = get_vectors(meta,X,state_col,cond,state)
            if len(vecs)!=2: continue
            v1,v2=vecs
            o1=np.argsort(-np.abs(v1)); o2=np.argsort(-np.abs(v2))
            r1=np.empty(G,int); r1[o1]=np.arange(G)+1
            r2=np.empty(G,int); r2[o2]=np.arange(G)+1
            for k in ks:
                k=min(k,G); a=set(o1[:k]); b=set(o2[:k]); inter=sorted(a&b); ov=len(inter)
                p=float(hypergeom.sf(ov-1,G,k,k)) if ov else 1.0
                metric_rows.append(dict(level=level,condition=cond,state=state,k=k,n_genes=G,overlap=ov,expected_overlap=k*k/G,enrichment=ov/(k*k/G),hypergeom_p=p,neglog10_p=(-np.log10(p) if p>0 else np.inf)))
                if k in (100,250):
                    for gi in inter:
                        s1=np.sign(v1[gi]); s2=np.sign(v2[gi])
                        gene_rows.append(dict(level=level,condition=cond,state=state,k=k,gene=genes[gi],delta_rep1=float(v1[gi]),delta_rep2=float(v2[gi]),rank_rep1=int(r1[gi]),rank_rep2=int(r2[gi]),same_direction=bool(s1==s2 and s1!=0),direction=("up" if s1>0 and s2>0 else "down" if s1<0 and s2<0 else "discordant")))

metrics=pd.DataFrame(metric_rows); genesdf=pd.DataFrame(gene_rows)
metrics.to_csv(OUT/"candidate_overlap_significance.tsv",sep="\t",index=False)
genesdf.to_csv(OUT/"robust_intersection_genes.tsv",sep="\t",index=False)

# Condition-level recurrence across author-defined tissue systems, using k=100 and same direction.
ts=genesdf[(genesdf.level=="TissueSystem") & (genesdf.k==100) & genesdf.same_direction].copy()
rec=(ts.groupby(["condition","gene","direction"],as_index=False)
     .agg(n_tissue_systems=("state","nunique"),tissue_systems=("state",lambda x:" | ".join(sorted(set(x)))),median_abs_delta_rep1=("delta_rep1",lambda x:float(np.median(np.abs(x)))),median_abs_delta_rep2=("delta_rep2",lambda x:float(np.median(np.abs(x))))))
rec=rec.sort_values(["condition","n_tissue_systems","median_abs_delta_rep1"],ascending=[True,False,False])
rec.to_csv(OUT/"tissuesystem_recurrent_top100_genes.tsv",sep="\t",index=False)

summ=(metrics.groupby(["level","condition","k"],as_index=False)
      .agg(n_states=("state","nunique"),median_overlap=("overlap","median"),median_enrichment=("enrichment","median"),median_neglog10_p=("neglog10_p","median"),max_p=("hypergeom_p","max"),states_p_lt_0_05=("hypergeom_p",lambda x:int((x<0.05).sum())),states_p_lt_1e_6=("hypergeom_p",lambda x:int((x<1e-6).sum()))))
summ.to_csv(OUT/"candidate_overlap_significance_summary.tsv",sep="\t",index=False)
print("ROBERTSON2026_CANDIDATE_SIGNIFICANCE_SUCCESS")
print(summ.to_string(index=False))
print("\nRECURRENT_TISSUE_TOP100")
print(rec.groupby("condition").head(20).to_string(index=False))
