from __future__ import annotations

"""植物扰动证据排序器。

目标：根据已经冻结并得到支持的生物过程类型选择证据并排序候选。
这不是完整转录组模拟器，也不输出未经校准的“正确概率”。

当前支持：
- persistent：最近同状态历史；
- acute_shared：多来源急性共同响应；
- recovery_shared：多来源快速恢复共同响应；
- combination：组成单扰动有符号求和。

未支持类型直接拒判，不自动换路线。
"""

from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json
import sys
from itertools import combinations
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
VERSION = "PerturbationEvidenceRanker-v0.2-2026-10-20"
SUPPORTED = {"persistent", "acute_shared", "recovery_shared", "combination"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_path(value: str | Path) -> Path:
    p = Path(value)
    if not p.is_absolute():
        p = ROOT / p
    return p.resolve()


def read_table(path: Path) -> pd.DataFrame:
    name = path.name.lower()
    if name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return pd.read_csv(path, sep="\t")
    if name.endswith((".csv", ".csv.gz")):
        return pd.read_csv(path)
    raise ValueError(f"不支持的输入格式：{path.name}")


def parse_boolish(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    s = series.astype(str).str.strip().str.lower()
    truthy = {"true", "1", "yes", "y", "t"}
    falsy = {"false", "0", "no", "n", "f"}
    bad = ~(s.isin(truthy | falsy))
    if bad.any():
        raise ValueError(f"筛选列存在无法识别的布尔值：{sorted(s[bad].unique())[:10]}")
    return s.isin(truthy)


def apply_filter(df: pd.DataFrame, cfg: dict[str, Any]) -> np.ndarray:
    f = cfg.get("filter")
    if not f:
        return np.ones(len(df), dtype=bool)
    col = f.get("column")
    if not col or col not in df.columns or "equals" not in f:
        raise ValueError("filter 必须包含存在的 column 与 equals")
    val = f["equals"]
    if isinstance(val, bool):
        return parse_boolish(df[col]).to_numpy() == val
    return (df[col].astype(str) == str(val)).to_numpy()


def direction_label(v: np.ndarray) -> np.ndarray:
    out = np.full(len(v), "不变", dtype=object)
    out[v > 0] = "上调"
    out[v < 0] = "下调"
    return out


def safe_rho(a: np.ndarray, b: np.ndarray) -> float | None:
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3 or np.std(a[m]) == 0 or np.std(b[m]) == 0:
        return None
    return float(spearmanr(a[m], b[m]).statistic)


def abs_rank(v: np.ndarray) -> np.ndarray:
    return pd.Series(np.abs(np.asarray(v, float))).rank(method="first", ascending=False).to_numpy(int)


def top_overlap(a: np.ndarray, b: np.ndarray, k: int) -> dict[str, Any]:
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a = a[m]; b = b[m]
    n = len(a); kk = min(k, n)
    if kk == 0:
        return {"k": 0, "overlap": 0, "retention_fraction": None, "direction_agreement": None}
    ra = abs_rank(a); rb = abs_rank(b)
    ia = set(np.where(ra <= kk)[0]); ib = set(np.where(rb <= kk)[0])
    ov = sorted(ia & ib)
    sign = float(np.mean(np.sign(a[ov]) == np.sign(b[ov]))) if ov else None
    return {"k": kk, "overlap": len(ov), "retention_fraction": len(ov) / kk, "direction_agreement": sign}


def source_matrix(df: pd.DataFrame, cols: list[str]) -> np.ndarray:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"缺少证据列：{missing}")
    return np.column_stack([pd.to_numeric(df[c], errors="coerce").to_numpy(float) for c in cols])


def source_diagnostics(E: np.ndarray, cols: list[str], top_ks: list[int]) -> dict[str, Any]:
    pairs = []
    for i, j in combinations(range(E.shape[1]), 2):
        item = {"source_a": cols[i], "source_b": cols[j], "spearman": safe_rho(E[:, i], E[:, j])}
        item["topk"] = [top_overlap(E[:, i], E[:, j], k) for k in top_ks]
        pairs.append(item)
    return {"n_sources": E.shape[1], "pairwise": pairs}


def signed_consensus(E: np.ndarray, weights: np.ndarray | None = None) -> dict[str, np.ndarray]:
    if weights is None:
        weights = np.ones(E.shape[1], float)
    weights = np.asarray(weights, float)
    if len(weights) != E.shape[1] or np.any(weights < 0) or weights.sum() <= 0:
        raise ValueError("证据权重无效")
    present = np.isfinite(E)
    W = present * weights[None, :]
    denom = W.sum(axis=1)
    filled = np.where(present, E, 0.0)
    signed = np.full(E.shape[0], np.nan)
    mean_abs = np.full(E.shape[0], np.nan)
    ok = denom > 0
    signed[ok] = (filled * W).sum(axis=1)[ok] / denom[ok]
    mean_abs[ok] = (np.abs(filled) * W).sum(axis=1)[ok] / denom[ok]
    coh = np.full(E.shape[0], np.nan)
    nz = ok & (mean_abs > 0)
    coh[nz] = np.abs(signed[nz]) / mean_abs[nz]
    return {
        "prediction": signed,
        "coherence": coh,
        "support_count": present.sum(axis=1).astype(int),
        "support_fraction": present.sum(axis=1) / E.shape[1],
        "mean_abs": mean_abs,
    }


def rank_output(df: pd.DataFrame, gene_col: str, pred: np.ndarray, score: np.ndarray,
                eligible: np.ndarray, support_count: np.ndarray, support_fraction: np.ndarray,
                coherence: np.ndarray, regime: str, method: str, extra: dict[str, np.ndarray] | None = None) -> pd.DataFrame:
    out = pd.DataFrame({
        "gene": df[gene_col].astype(str),
        "eligible": eligible,
        "predicted_response": pred,
        "predicted_direction": direction_label(np.nan_to_num(pred, nan=0.0)),
        "priority_score": score,
        "evidence_support_count": support_count,
        "evidence_support_fraction": support_fraction,
        "evidence_direction_coherence": coherence,
        "regime": regime,
        "method": method,
    })
    if extra:
        for k, v in extra.items():
            out[k] = v
    q = out[out.eligible & np.isfinite(out.priority_score)].copy()
    q = q.sort_values(["priority_score", "gene"], ascending=[False, True], kind="mergesort")
    q.insert(0, "priority_rank", np.arange(1, len(q)+1))
    return q.reset_index(drop=True)


def persistent_route(df: pd.DataFrame, cfg: dict[str, Any], top_ks: list[int]):
    gene_col = cfg.get("gene_column", "gene")
    latest = cfg.get("latest_column")
    if gene_col not in df.columns or not latest or latest not in df.columns:
        raise ValueError("持续型需要 gene_column 与 latest_column")
    pred = pd.to_numeric(df[latest], errors="coerce").to_numpy(float)
    stability_cols = list(cfg.get("stability_columns", []))
    if stability_cols:
        E = source_matrix(df, stability_cols)
        diag = source_diagnostics(E, stability_cols, top_ks)
        agg = signed_consensus(E)
        support_count = agg["support_count"]
        support_fraction = agg["support_fraction"]
        coh = agg["coherence"]
    else:
        diag = {"n_sources": 1, "note": "未提供独立稳定性来源；不据此推断置信概率"}
        support_count = np.isfinite(pred).astype(int)
        support_fraction = support_count.astype(float)
        coh = np.full(len(df), np.nan)
    eligible = apply_filter(df, cfg) & np.isfinite(pred)
    ranked = rank_output(df, gene_col, pred, np.abs(pred), eligible, support_count, support_fraction, coh,
                         "persistent", f"最近同状态历史：{latest}")
    diag["interpretation"] = "当前没有经过校准的来源稳定性硬阈值；稳定性仅作诊断，不改变持续型排名。"
    return ranked, diag


def shared_route(df: pd.DataFrame, cfg: dict[str, Any], regime: str, top_ks: list[int]):
    gene_col = cfg.get("gene_column", "gene")
    donors = list(cfg.get("donor_columns", []))
    if gene_col not in df.columns or len(donors) < 2:
        raise ValueError("共享型至少需要两个 donor_columns")
    E = source_matrix(df, donors)
    w = cfg.get("donor_weights")
    weights = np.asarray(w, float) if w is not None else np.ones(len(donors), float)
    agg = signed_consensus(E, weights)
    min_support = int(cfg.get("min_support_count", len(donors)))
    eligible = apply_filter(df, cfg) & np.isfinite(agg["prediction"]) & (agg["support_count"] >= min_support)
    ranked = rank_output(df, gene_col, agg["prediction"], np.abs(agg["prediction"]), eligible,
                         agg["support_count"], agg["support_fraction"], agg["coherence"], regime,
                         "多来源有符号共同响应", {"mean_absolute_source_response": agg["mean_abs"]})
    diag = source_diagnostics(E, donors, top_ks)
    loo=[]
    if len(donors) >= 3:
        full = agg["prediction"]
        for i, name in enumerate(donors):
            keep=[j for j in range(len(donors)) if j != i]
            pred2=signed_consensus(E[:,keep],weights[keep])["prediction"]
            loo.append({"left_out":name,"topk":[top_overlap(full,pred2,k) for k in top_ks]})
    diag["leave_one_source_out"] = loo
    diag["min_support_count"] = min_support
    diag["interpretation"] = "来源两两一致性与 leave-one-source-out 排名稳定性用于诊断来源敏感性；当前不转换为统一正确概率。"
    return ranked, diag


def combination_route(df: pd.DataFrame, cfg: dict[str, Any], top_ks: list[int]):
    gene_col = cfg.get("gene_column", "gene")
    parts = list(cfg.get("constituent_columns", []))
    if gene_col not in df.columns or len(parts) < 2:
        raise ValueError("组合型至少需要两个 constituent_columns")
    E = source_matrix(df, parts)
    present = np.isfinite(E)
    min_support = int(cfg.get("min_support_count", len(parts)))
    filled = np.where(present, E, 0.0)
    pred = filled.sum(axis=1)
    support_count = present.sum(axis=1).astype(int)
    support_fraction = support_count / len(parts)
    den = np.abs(filled).sum(axis=1)
    coh = np.full(len(df), np.nan)
    np.divide(np.abs(pred), den, out=coh, where=den>0)
    eligible = apply_filter(df, cfg) & (support_count >= min_support) & np.isfinite(pred)
    absE=np.abs(E)
    strongest=np.nanmax(absE,axis=1)
    total=np.nansum(absE,axis=1)
    dominance=np.full(len(df),np.nan)
    np.divide(strongest,total,out=dominance,where=total>0)
    ranked=rank_output(df,gene_col,pred,np.abs(pred),eligible,support_count,support_fraction,coh,
                       "combination","组成单扰动有符号求和",{"constituent_dominance_gene":dominance})
    diag=source_diagnostics(E,parts,top_ks)
    rms=np.sqrt(np.nanmean(E*E,axis=0))
    diag["constituent_rms"]={parts[i]:float(rms[i]) for i in range(len(parts))}
    diag["magnitude_balance_min_over_max"]=float(np.nanmin(rms)/np.nanmax(rms)) if np.nanmax(rms)>0 else None
    diag["dominant_constituent"]=parts[int(np.nanargmax(rms))]
    diag["interpretation"] = "组合型不把高一致性解释为更高组合增益；重点报告组成互补性、幅度平衡和主导成员。"
    return ranked, diag


def write_abstain(cfg: dict[str, Any], config_path: Path, outdir: Path, reason: str):
    outdir.mkdir(parents=True,exist_ok=True)
    m={"version":VERSION,"status":"拒判","reason":reason,"regime":cfg.get("regime"),
       "config_sha256":sha256_file(config_path),"generated_utc":datetime.now(timezone.utc).isoformat()}
    (outdir/'manifest.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf-8')
    print("拒判")
    print(json.dumps(m,ensure_ascii=False,indent=2))


def run(config_path: Path) -> int:
    cfg=json.loads(config_path.read_text(encoding='utf-8'))
    regime=str(cfg.get('regime','unknown')).strip()
    outdir=resolve_path(cfg.get('output_dir','data/perturbation_evidence_ranker_output'))
    if regime not in SUPPORTED:
        write_abstain(cfg,config_path,outdir,f"类型 {regime} 当前没有已支持的排序规则，不自动换路线。")
        return 0
    inp=resolve_path(cfg['input'])
    if not inp.exists(): raise FileNotFoundError(inp)
    df=read_table(inp)
    top_ks=sorted(set(int(k) for k in cfg.get('top_k',[100,250,500,1000]) if int(k)>0))
    if regime=='persistent': ranked,diag=persistent_route(df,cfg,top_ks)
    elif regime in {'acute_shared','recovery_shared'}: ranked,diag=shared_route(df,cfg,regime,top_ks)
    else: ranked,diag=combination_route(df,cfg,top_ks)
    if len(ranked)==0:
        write_abstain(cfg,config_path,outdir,"核心证据覆盖后没有可排序基因。")
        return 0
    outdir.mkdir(parents=True,exist_ok=True)
    ranking_path=outdir/'ranking.tsv.gz'
    ranked.to_csv(ranking_path,sep='\t',index=False,compression='gzip')
    maxk=min(max(top_ks) if top_ks else 100,len(ranked))
    ranked.head(maxk).to_csv(outdir/'top_candidates.tsv',sep='\t',index=False)
    (outdir/'source_diagnostics.json').write_text(json.dumps(diag,ensure_ascii=False,indent=2),encoding='utf-8')
    rows=[]
    for k in top_ks:
        kk=min(k,len(ranked)); z=ranked.head(kk)
        rows.append({'k':kk,'median_priority_score':float(z.priority_score.median()),
                     'median_support_fraction':float(z.evidence_support_fraction.median()),
                     'median_direction_coherence':float(z.evidence_direction_coherence.median(skipna=True)) if z.evidence_direction_coherence.notna().any() else None})
    pd.DataFrame(rows).to_csv(outdir/'top_k_summary.tsv',sep='\t',index=False)
    manifest={
        'version':VERSION,'status':'有证据','task_name':cfg.get('task_name'),'regime':regime,
        'input':str(inp),'input_sha256':sha256_file(inp),'config_sha256':sha256_file(config_path),
        'ranking_sha256':sha256_file(ranking_path),'n_input_rows':int(len(df)),'n_ranked_genes':int(len(ranked)),
        'top_k':top_ks,'confidence_policy':'不输出统一正确概率；参见 source_diagnostics.json 与项目置信/拒判规则文件',
        'generated_utc':datetime.now(timezone.utc).isoformat(),
    }
    (outdir/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print('扰动证据排序完成')
    print(json.dumps(manifest,ensure_ascii=False,indent=2))
    print('\n前20候选')
    print(ranked[['priority_rank','gene','priority_score','predicted_direction','predicted_response','evidence_support_fraction','evidence_direction_coherence']].head(20).to_string(index=False))
    return 0


def main():
    ap=argparse.ArgumentParser(description='植物扰动证据排序器')
    ap.add_argument('--config',required=True)
    args=ap.parse_args()
    cp=resolve_path(args.config)
    try:
        rc=run(cp)
    except Exception as e:
        print(f'运行失败：{type(e).__name__}: {e}',file=sys.stderr)
        raise
    raise SystemExit(rc)

if __name__=='__main__':
    main()
