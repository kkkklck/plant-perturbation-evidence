from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
ROB = ROOT / "data" / "external_rice" / "Robertson2026"
DESIGN = pd.read_csv(ROB / "sample_design.tsv", sep="\t")
OUT = ROB / "candidate_recovery_v0"
OUT.mkdir(parents=True, exist_ok=True)

sample_batch = dict(zip(DESIGN.identity.astype(str), DESIGN.batch.astype(int)))
conditions = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
control = "Well.Watered"
ks = [100, 250, 500, 1000]


def load_level(level: str):
    if level == "CellAnnotation":
        base = ROB / "pseudobulk_v0"
        meta = pd.read_csv(base / "group_metadata.tsv", sep="\t")
        expr = pd.read_csv(base / "pseudobulk_log1p_cpm.tsv.gz", sep="\t")
        state_col = "cell_annotation"
    elif level == "TissueSystem":
        base = ROB / "tissuesystem_pseudobulk_v0"
        meta = pd.read_csv(base / "group_metadata.tsv", sep="\t")
        expr = pd.read_csv(base / "tissuesystem_log1p_cpm.tsv.gz", sep="\t")
        state_col = "tissue_system"
    else:
        raise ValueError(level)
    genes = expr.pop("gene").astype(str).to_numpy()
    X = expr.to_numpy(float)
    if list(expr.columns.astype(str)) != list(meta.group.astype(str)):
        raise RuntimeError(f"{level}: expression columns do not match metadata")
    mean_expr = X.mean(axis=1)
    nonzero_frac = (X > 0).mean(axis=1)
    mask = (nonzero_frac >= 0.25) & (mean_expr >= 0.05)
    return meta, genes[mask], X[mask], state_col


def analyze(level: str):
    meta, genes, X, state_col = load_level(level)
    G = len(genes)
    ctrl_by_batch = {}
    for s in meta.loc[meta.condition.eq(control), "sample"].astype(str).unique():
        ctrl_by_batch[sample_batch[s]] = s

    def get_idx(sample, state):
        z = meta.index[(meta["sample"].astype(str) == sample) & (meta[state_col].astype(str) == state)].to_numpy()
        return z[0] if len(z) == 1 else None

    rows = []
    states = sorted(meta[state_col].astype(str).unique())
    for cond in conditions:
        for state in states:
            t_idx = meta.index[(meta.condition.eq(cond)) & (meta[state_col].astype(str).eq(state))].to_numpy()
            vecs = []
            pairs = []
            for ti in t_idx:
                ss = str(meta.loc[ti, "sample"])
                b = sample_batch[ss]
                cs = ctrl_by_batch.get(b)
                if cs is None:
                    continue
                ci = get_idx(cs, state)
                if ci is None:
                    continue
                vecs.append(X[:, ti] - X[:, ci])
                pairs.append((ss, cs, b))
            if len(vecs) != 2:
                continue
            v1, v2 = vecs
            abs1 = np.abs(v1); abs2 = np.abs(v2)
            order1 = np.argsort(-abs1)
            order2 = np.argsort(-abs2)
            rank1 = np.empty(G, dtype=int); rank1[order1] = np.arange(G)
            rank2 = np.empty(G, dtype=int); rank2[order2] = np.arange(G)
            genome_rho = float(spearmanr(v1, v2).statistic)

            for k in ks:
                kk = min(k, G)
                a = set(order1[:kk].tolist()); bset = set(order2[:kk].tolist())
                ov = len(a & bset)
                expected = kk * kk / G
                enrich = ov / expected if expected > 0 else np.nan
                recall = ov / kk

                # Symmetric rank percentile: how highly does each replicate rank the other's candidates?
                p12 = 1.0 - rank2[order1[:kk]] / max(G - 1, 1)
                p21 = 1.0 - rank1[order2[:kk]] / max(G - 1, 1)
                med_pct = float(np.median(np.concatenate([p12, p21])))

                # Direction consistency among genes selected by either replicate.
                union = np.array(sorted(a | bset), dtype=int)
                sign_agree_union = float(np.mean(np.sign(v1[union]) == np.sign(v2[union]))) if len(union) else np.nan
                # Among actually overlapping top-|delta| candidates.
                inter = np.array(sorted(a & bset), dtype=int)
                sign_agree_overlap = float(np.mean(np.sign(v1[inter]) == np.sign(v2[inter]))) if len(inter) else np.nan

                # Direction-aware candidate recovery: top k/2 up + top k/2 down.
                each = max(1, kk // 2)
                up1 = set(np.argsort(-v1)[:each].tolist())
                up2 = set(np.argsort(-v2)[:each].tolist())
                dn1 = set(np.argsort(v1)[:each].tolist())
                dn2 = set(np.argsort(v2)[:each].tolist())
                same_dir_ov = len(up1 & up2) + len(dn1 & dn2)
                expected_dir = 2 * each * each / G
                dir_enrich = same_dir_ov / expected_dir if expected_dir > 0 else np.nan
                dir_recall = same_dir_ov / (2 * each)

                rows.append({
                    "level": level,
                    "condition": cond,
                    "state": state,
                    "k": kk,
                    "n_genes": G,
                    "genome_response_rho": genome_rho,
                    "abs_topk_overlap": ov,
                    "abs_topk_expected_overlap": expected,
                    "abs_topk_enrichment": enrich,
                    "abs_topk_recall": recall,
                    "crossrep_median_rank_percentile": med_pct,
                    "sign_agreement_union": sign_agree_union,
                    "sign_agreement_overlap": sign_agree_overlap,
                    "directional_topk_same_direction_overlap": same_dir_ov,
                    "directional_topk_enrichment": dir_enrich,
                    "directional_topk_recall": dir_recall,
                    "matched_pair_1": f"{pairs[0][0]}-{pairs[0][1]}:batch{pairs[0][2]}",
                    "matched_pair_2": f"{pairs[1][0]}-{pairs[1][1]}:batch{pairs[1][2]}",
                })
    return pd.DataFrame(rows)

all_rows = pd.concat([analyze("CellAnnotation"), analyze("TissueSystem")], ignore_index=True)
all_rows.to_csv(OUT / "crossrep_candidate_recovery.tsv", sep="\t", index=False)

summary = (all_rows.groupby(["level", "condition", "k"], as_index=False)
           .agg(n_states=("state", "nunique"),
                median_genome_rho=("genome_response_rho", "median"),
                median_abs_topk_enrichment=("abs_topk_enrichment", "median"),
                median_abs_topk_recall=("abs_topk_recall", "median"),
                median_crossrep_rank_percentile=("crossrep_median_rank_percentile", "median"),
                median_sign_agreement_union=("sign_agreement_union", "median"),
                median_sign_agreement_overlap=("sign_agreement_overlap", "median"),
                median_directional_topk_enrichment=("directional_topk_enrichment", "median"),
                median_directional_topk_recall=("directional_topk_recall", "median")))
summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)

print("ROBERTSON2026_CROSSREP_CANDIDATE_RECOVERY_SUCCESS")
print(summary.to_string(index=False))
