from __future__ import annotations

from pathlib import Path
import argparse, hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, hypergeom

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_salt" / "GSE153103"
WB = BASE / "GSE153103_1_raw_reads_counts.xlsx"
OUT = BASE / "plantvc_v1_6h_blind"
RULES = ROOT / "PLANTVC_V1_GSE153103_6H_FROZEN_RULES_2026-10-20.md"

PRE = ["RAS048", "RAS049", "RAS050", "RAS056", "RAS057", "RAS058"]
TARGET = ["RAS051", "RAS059"]
TIME = {
    "0h": ["RAS048", "RAS056"],
    "1h": ["RAS049", "RAS057"],
    "3h": ["RAS050", "RAS058"],
    "6h": ["RAS051", "RAS059"],
}
KS = [100, 250, 500, 1000]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_cols(cols: list[str]) -> pd.DataFrame:
    # Read header first; no numeric target values are inspected here.
    header = pd.read_excel(WB, sheet_name="raw_counts", nrows=0)
    gene_col = header.columns[0]
    missing = [c for c in cols if c not in header.columns]
    if missing:
        raise RuntimeError(f"Missing columns: {missing}")
    df = pd.read_excel(WB, sheet_name="raw_counts", usecols=[gene_col] + cols)
    df = df.rename(columns={gene_col: "gene"})
    df["gene"] = df["gene"].astype(str)
    if df["gene"].duplicated().any():
        raise RuntimeError("Duplicate gene IDs in workbook")
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="raise")
    return df


def log2cpm(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    counts = df[cols].to_numpy(float)
    libs = counts.sum(axis=0)
    if np.any(libs <= 0):
        raise RuntimeError("Non-positive library size")
    cpm = counts / libs[None, :] * 1e6
    return pd.DataFrame(np.log2(cpm + 1.0), columns=cols, index=df.index), libs


def sign_same(a, b):
    return (np.sign(a) == np.sign(b)) & (np.sign(a) != 0)


def freeze():
    OUT.mkdir(parents=True, exist_ok=True)
    # Explicit firewall: target columns are not passed to read_cols.
    df = read_cols(PRE)
    if any(c in df.columns for c in TARGET):
        raise RuntimeError("TARGET FIREWALL BREACH")
    L, libs = log2cpm(df, PRE)

    r1_rep1 = L["RAS049"] - L["RAS048"]
    r1_rep2 = L["RAS057"] - L["RAS056"]
    r3_rep1 = L["RAS050"] - L["RAS048"]
    r3_rep2 = L["RAS058"] - L["RAS056"]
    r1 = (r1_rep1 + r1_rep2) / 2.0
    r3 = (r3_rep1 + r3_rep2) / 2.0

    temporal_same = sign_same(r1.to_numpy(), r3.to_numpy())
    bonus = np.where(temporal_same, np.sqrt(np.abs(r1.to_numpy() * r3.to_numpy())), 0.0)
    magnitude = 0.5 * np.abs(r3.to_numpy()) + 0.5 * bonus
    pred = np.sign(r3.to_numpy()) * magnitude

    rep1agree = sign_same(r1_rep1.to_numpy(), r1_rep2.to_numpy()).astype(float)
    rep3agree = sign_same(r3_rep1.to_numpy(), r3_rep2.to_numpy()).astype(float)
    confidence = (rep1agree + rep3agree + temporal_same.astype(float)) / 3.0

    # Frozen technical universe: CPM >=1 in at least 2 of the six pretarget samples.
    raw = df[PRE].to_numpy(float)
    cpm = raw / libs[None, :] * 1e6
    technical = (cpm >= 1.0).sum(axis=1) >= 2

    out = pd.DataFrame({
        "gene": df["gene"],
        "technical_universe": technical,
        "r1": r1,
        "r3": r3,
        "pred_v1": pred,
        "score_v1_abs": np.abs(pred),
        "pred_baseline_3h": r3,
        "score_baseline_3h_abs": np.abs(r3),
        "pred_baseline_1h": r1,
        "score_baseline_1h_abs": np.abs(r1),
        "confidence": confidence,
        "temporal_same_direction": temporal_same,
        "replicate_agree_1h": rep1agree.astype(bool),
        "replicate_agree_3h": rep3agree.astype(bool),
    })
    for score in ["score_v1_abs", "score_baseline_3h_abs", "score_baseline_1h_abs"]:
        out[f"rank_{score}"] = out[score].rank(method="min", ascending=False).astype(int)
        out[f"rank_{score}_technical"] = np.nan
        idx = out.index[out.technical_universe]
        out.loc[idx, f"rank_{score}_technical"] = out.loc[idx, score].rank(method="min", ascending=False)

    rankfile = OUT / "FROZEN_PRETARGET_RANKING.tsv.gz"
    out.to_csv(rankfile, sep="\t", index=False, compression="gzip")

    manifest = {
        "status": "FROZEN_BEFORE_TARGET",
        "workbook_sha256": sha256(WB),
        "rules_sha256": sha256(RULES),
        "ranking_sha256": sha256(rankfile),
        "pretarget_columns_read": PRE,
        "target_columns_not_read": TARGET,
        "n_genes_all": int(len(out)),
        "n_genes_technical": int(out.technical_universe.sum()),
        "library_sizes_pretarget": {c: float(v) for c, v in zip(PRE, libs)},
        "formula": "sign(r3)*(0.5*abs(r3)+0.5*sqrt(abs(r1*r3))*I[same_sign])",
        "confidence_used_for_ranking": False,
    }
    (OUT / "freeze_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("PLANTVC_V1_GSE153103_FREEZE_SUCCESS")
    print(json.dumps(manifest, indent=2))
    print("TARGET VALUES HAVE NOT BEEN READ BY THIS COMMAND")


def topk_stats(pred_score, pred_signed, target_abs, target_signed, k):
    n = len(pred_score)
    k = min(k, n)
    pi = np.argpartition(-pred_score, k - 1)[:k]
    ti = np.argpartition(-target_abs, k - 1)[:k]
    ps, ts = set(pi.tolist()), set(ti.tolist())
    inter = sorted(ps & ts)
    overlap = len(inter)
    expected = k * k / n
    enrich = overlap / expected if expected > 0 else np.nan
    p = float(hypergeom.sf(overlap - 1, n, k, k))
    if overlap:
        dirs = np.sign(pred_signed[inter]) == np.sign(target_signed[inter])
        sign_agree = float(np.mean(dirs))
    else:
        sign_agree = np.nan
    # percentile of true top-k genes in predicted ranking (1=best)
    order = np.argsort(-pred_score)
    ranks = np.empty(n, int); ranks[order] = np.arange(1, n + 1)
    percentiles = 1.0 - (ranks[ti] - 1) / max(n - 1, 1)
    return {
        "k": k,
        "overlap": overlap,
        "random_expected_overlap": expected,
        "enrichment": enrich,
        "hypergeom_p": p,
        "sign_agreement_overlap": sign_agree,
        "median_target_topk_predicted_percentile": float(np.median(percentiles)),
    }


def score():
    manifest_path = OUT / "freeze_manifest.json"
    rankfile = OUT / "FROZEN_PRETARGET_RANKING.tsv.gz"
    if not manifest_path.exists() or not rankfile.exists():
        raise RuntimeError("Freeze artifacts missing; run freeze first")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if sha256(rankfile) != manifest["ranking_sha256"]:
        raise RuntimeError("Frozen ranking hash mismatch")
    if sha256(RULES) != manifest["rules_sha256"]:
        raise RuntimeError("Frozen rules hash mismatch")

    frozen = pd.read_csv(rankfile, sep="\t", compression="gzip")
    # Target is opened only here.
    targ = read_cols(["RAS048", "RAS051", "RAS056", "RAS059"])
    if targ["gene"].tolist() != frozen["gene"].astype(str).tolist():
        raise RuntimeError("Gene order mismatch between frozen ranking and target read")
    L, libs = log2cpm(targ, ["RAS048", "RAS051", "RAS056", "RAS059"])
    r6_rep1 = (L["RAS051"] - L["RAS048"]).to_numpy(float)
    r6_rep2 = (L["RAS059"] - L["RAS056"]).to_numpy(float)
    r6 = (r6_rep1 + r6_rep2) / 2.0

    models = {
        "PlantVC_v1": (frozen["score_v1_abs"].to_numpy(float), frozen["pred_v1"].to_numpy(float)),
        "baseline_3h": (frozen["score_baseline_3h_abs"].to_numpy(float), frozen["pred_baseline_3h"].to_numpy(float)),
        "baseline_1h": (frozen["score_baseline_1h_abs"].to_numpy(float), frozen["pred_baseline_1h"].to_numpy(float)),
    }

    rows = []
    vecrows = []
    for universe, mask in [
        ("technical", frozen["technical_universe"].to_numpy(bool)),
        ("all", np.ones(len(frozen), dtype=bool)),
    ]:
        target = r6[mask]
        ta = np.abs(target)
        for name, (score_abs, signed) in models.items():
            ps = score_abs[mask]
            sg = signed[mask]
            rho = float(spearmanr(sg, target).statistic)
            sign_all = float(np.mean(np.sign(sg) == np.sign(target)))
            vecrows.append({"universe": universe, "model": name, "n_genes": int(mask.sum()), "spearman_full_vector": rho, "sign_agreement_all": sign_all})
            for k in KS:
                z = topk_stats(ps, sg, ta, target, k)
                rows.append({"universe": universe, "model": name, **z})

    # Biological replicate ceiling/reference at 6h.
    rep_rows = []
    for universe, mask in [("technical", frozen["technical_universe"].to_numpy(bool)), ("all", np.ones(len(frozen), dtype=bool))]:
        a, b = r6_rep1[mask], r6_rep2[mask]
        rep_rows.append({
            "universe": universe,
            "n_genes": int(mask.sum()),
            "replicate_spearman": float(spearmanr(a, b).statistic),
            "replicate_sign_agreement_all": float(np.mean(np.sign(a) == np.sign(b))),
        })
        for k in KS:
            z = topk_stats(np.abs(a), a, np.abs(b), b, k)
            z.update({"universe": universe, "comparison": "6h_rep1_vs_rep2"})
            rep_rows.append(z)

    pd.DataFrame(rows).to_csv(OUT / "blind_topk_results.tsv", sep="\t", index=False)
    pd.DataFrame(vecrows).to_csv(OUT / "blind_vector_results.tsv", sep="\t", index=False)
    pd.DataFrame(rep_rows).to_csv(OUT / "target_replicate_reference.tsv", sep="\t", index=False)

    target_df = pd.DataFrame({"gene": frozen.gene, "r6": r6, "r6_rep1": r6_rep1, "r6_rep2": r6_rep2})
    target_df.to_csv(OUT / "target_6h_response.tsv.gz", sep="\t", index=False, compression="gzip")

    print("PLANTVC_V1_GSE153103_SCORE_SUCCESS")
    print(pd.DataFrame(vecrows).to_string(index=False))
    rr = pd.DataFrame(rows)
    print("\nTECHNICAL TOP-K")
    print(rr[rr.universe.eq("technical")].to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=["freeze", "score"])
    args = ap.parse_args()
    if args.command == "freeze": freeze()
    else: score()

if __name__ == "__main__":
    main()
