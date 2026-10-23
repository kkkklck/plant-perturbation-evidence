from __future__ import annotations

"""Unified PlantVC v1 candidate-prioritization entrypoint.

This is deliberately a ranking system, not a whole-transcriptome simulator.

Validated regimes currently implemented:
- persistent: use the latest reliable same-system response as the primary predictor;
- acute_shared: use the signed consensus across same-perturbation donor states;
- recovery_shared: same engine as acute_shared, with recovery/switch semantics.

Other regimes are intentionally not silently approximated. They return an
ABSTAIN manifest until a prospective rule is frozen for them.
"""

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VERSION = "PlantVC-v1.0-prototype-2026-10-20"
SUPPORTED = {"persistent", "acute_shared", "recovery_shared"}
KNOWN_BUT_UNSUPPORTED = {"mechanism", "combination", "rare_specialist", "unknown"}


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
    if name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".txt") or name.endswith(".txt.gz"):
        return pd.read_csv(path, sep="\t")
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return pd.read_csv(path)
    raise ValueError(f"Unsupported input format: {path.name}; use TSV/CSV, optionally gzipped")


def parse_boolish(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    s = series.astype(str).str.strip().str.lower()
    truthy = {"true", "1", "yes", "y", "t"}
    falsy = {"false", "0", "no", "n", "f"}
    unknown = ~(s.isin(truthy | falsy))
    if unknown.any():
        vals = sorted(s[unknown].dropna().unique().tolist())[:10]
        raise ValueError(f"Filter column contains values that are not boolean-like: {vals}")
    return s.isin(truthy)


def direction_label(x: np.ndarray) -> np.ndarray:
    out = np.full(len(x), "flat", dtype=object)
    out[x > 0] = "up"
    out[x < 0] = "down"
    return out


def weighted_consensus(M: np.ndarray, weights: np.ndarray) -> dict[str, np.ndarray]:
    """Per-gene weighted signed consensus over finite evidence columns."""
    present = np.isfinite(M)
    W = present * weights[None, :]
    wsum = W.sum(axis=1)
    filled = np.where(present, M, 0.0)

    signed_num = (filled * W).sum(axis=1)
    abs_num = (np.abs(filled) * W).sum(axis=1)

    signed_mean = np.full(M.shape[0], np.nan, dtype=float)
    mean_abs = np.full(M.shape[0], np.nan, dtype=float)
    ok = wsum > 0
    signed_mean[ok] = signed_num[ok] / wsum[ok]
    mean_abs[ok] = abs_num[ok] / wsum[ok]

    coherence = np.full(M.shape[0], np.nan, dtype=float)
    nz = ok & (mean_abs > 0)
    coherence[nz] = np.abs(signed_mean[nz]) / mean_abs[nz]

    support_count = present.sum(axis=1).astype(int)
    support_fraction = support_count / M.shape[1]

    centered = np.where(present, M - signed_mean[:, None], 0.0)
    var_num = ((centered ** 2) * W).sum(axis=1)
    sd = np.full(M.shape[0], np.nan, dtype=float)
    sd[ok] = np.sqrt(var_num[ok] / wsum[ok])

    return {
        "signed_mean": signed_mean,
        "mean_abs": mean_abs,
        "coherence": coherence,
        "support_count": support_count,
        "support_fraction": support_fraction,
        "evidence_sd": sd,
    }


def apply_filter(df: pd.DataFrame, cfg: dict[str, Any]) -> np.ndarray:
    filt = cfg.get("filter")
    if not filt:
        return np.ones(len(df), dtype=bool)
    col = filt.get("column")
    if not col or col not in df.columns:
        raise ValueError(f"Filter column not found: {col}")
    if "equals" not in filt:
        raise ValueError("filter must include an 'equals' value")
    target = filt["equals"]
    if isinstance(target, bool):
        return parse_boolish(df[col]).to_numpy() == target
    return (df[col].astype(str) == str(target)).to_numpy()


def persistent_rank(df: pd.DataFrame, cfg: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    gene_col = cfg.get("gene_column", "gene")
    latest = cfg.get("latest_column")
    if gene_col not in df.columns:
        raise ValueError(f"gene column not found: {gene_col}")
    if not latest or latest not in df.columns:
        raise ValueError(f"persistent regime requires latest_column present in input; got {latest}")

    hist = list(cfg.get("history_columns", []))
    missing = [c for c in hist if c not in df.columns]
    if missing:
        raise ValueError(f"history columns missing: {missing}")

    latest_v = pd.to_numeric(df[latest], errors="coerce").to_numpy(float)
    evidence_cols = hist + ([latest] if latest not in hist else [])
    E = np.column_stack([pd.to_numeric(df[c], errors="coerce").to_numpy(float) for c in evidence_cols])
    agg = weighted_consensus(E, np.ones(len(evidence_cols), dtype=float))

    eligible = apply_filter(df, cfg) & np.isfinite(latest_v)
    score = np.abs(latest_v)

    # Historical agreement is descriptive confidence only; it does not affect ranking.
    confidence = np.full(len(df), np.nan, dtype=float)
    if len(evidence_cols) >= 2:
        confidence = agg["coherence"] * agg["support_fraction"]

    out = pd.DataFrame({
        "gene": df[gene_col].astype(str),
        "eligible": eligible,
        "predicted_response": latest_v,
        "predicted_direction": direction_label(np.nan_to_num(latest_v, nan=0.0)),
        "priority_score": score,
        "evidence_support_count": agg["support_count"],
        "evidence_support_fraction": agg["support_fraction"],
        "evidence_direction_coherence": agg["coherence"],
        "evidence_sd": agg["evidence_sd"],
        "confidence_descriptor": confidence,
        "regime": "persistent",
        "method": f"latest same-system response: {latest}",
    })
    ranked = rank_eligible(out)
    meta = {
        "primary_evidence": latest,
        "history_columns_for_confidence_only": hist,
        "ranking_rule": f"abs({latest}) descending",
        "confidence_affects_ranking": False,
    }
    return ranked, meta


def shared_rank(df: pd.DataFrame, cfg: dict[str, Any], regime: str) -> tuple[pd.DataFrame, dict[str, Any]]:
    gene_col = cfg.get("gene_column", "gene")
    donors = list(cfg.get("donor_columns", []))
    if gene_col not in df.columns:
        raise ValueError(f"gene column not found: {gene_col}")
    if len(donors) < 2:
        raise ValueError(f"{regime} requires at least two donor_columns")
    missing = [c for c in donors if c not in df.columns]
    if missing:
        raise ValueError(f"donor columns missing: {missing}")

    wcfg = cfg.get("donor_weights")
    if wcfg is None:
        weights = np.ones(len(donors), dtype=float)
    else:
        weights = np.asarray(wcfg, dtype=float)
        if weights.shape != (len(donors),):
            raise ValueError("donor_weights length must equal donor_columns length")
        if np.any(~np.isfinite(weights)) or np.any(weights < 0) or weights.sum() <= 0:
            raise ValueError("donor_weights must be finite, nonnegative, and not all zero")

    E = np.column_stack([pd.to_numeric(df[c], errors="coerce").to_numpy(float) for c in donors])
    agg = weighted_consensus(E, weights)
    pred = agg["signed_mean"]
    eligible = apply_filter(df, cfg) & np.isfinite(pred)
    score = np.abs(pred)
    confidence = agg["coherence"] * agg["support_fraction"]

    out = pd.DataFrame({
        "gene": df[gene_col].astype(str),
        "eligible": eligible,
        "predicted_response": pred,
        "predicted_direction": direction_label(np.nan_to_num(pred, nan=0.0)),
        "priority_score": score,
        "mean_absolute_donor_response": agg["mean_abs"],
        "evidence_support_count": agg["support_count"],
        "evidence_support_fraction": agg["support_fraction"],
        "evidence_direction_coherence": agg["coherence"],
        "evidence_sd": agg["evidence_sd"],
        "confidence_descriptor": confidence,
        "regime": regime,
        "method": "weighted signed consensus across donor states",
    })
    ranked = rank_eligible(out)
    meta = {
        "primary_evidence": donors,
        "donor_weights": weights.tolist(),
        "ranking_rule": "abs(weighted signed donor consensus) descending",
        "confidence_affects_ranking": False,
        "confidence_descriptor": "direction coherence × evidence coverage; descriptive, not calibrated probability",
    }
    return ranked, meta


def rank_eligible(out: pd.DataFrame) -> pd.DataFrame:
    q = out[out["eligible"] & np.isfinite(out["priority_score"])].copy()
    # Stable deterministic tie break by gene. Confidence is deliberately NOT a tie breaker.
    q = q.sort_values(["priority_score", "gene"], ascending=[False, True], kind="mergesort")
    q.insert(0, "priority_rank", np.arange(1, len(q) + 1, dtype=int))
    return q.reset_index(drop=True)


def write_abstain(cfg: dict[str, Any], config_path: Path, outdir: Path, reason: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "version": VERSION,
        "status": "ABSTAIN",
        "reason": reason,
        "regime": cfg.get("regime"),
        "config_sha256": sha256_file(config_path),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print("PLANTVC_V1_ABSTAIN")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


def run(config_path: Path) -> int:
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    regime = str(cfg.get("regime", "unknown")).strip()
    outdir = resolve_path(cfg.get("output_dir", "data/plantvc_v1_output"))

    if regime not in SUPPORTED:
        reason = (
            f"Regime '{regime}' does not yet have a prospectively frozen PlantVC v1 ranking rule. "
            "No fallback predictor is applied."
        )
        write_abstain(cfg, config_path, outdir, reason)
        return 0

    input_path = resolve_path(cfg["input"])
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    df = read_table(input_path)

    if regime == "persistent":
        ranked, method_meta = persistent_rank(df, cfg)
    else:
        ranked, method_meta = shared_rank(df, cfg, regime)

    outdir.mkdir(parents=True, exist_ok=True)
    ranking_path = outdir / "ranking.tsv.gz"
    ranked.to_csv(ranking_path, sep="\t", index=False, compression="gzip")

    top_ks = [int(x) for x in cfg.get("top_k", [100, 250, 500, 1000]) if int(x) > 0]
    max_k = min(max(top_ks) if top_ks else min(100, len(ranked)), len(ranked))
    ranked.head(max_k).to_csv(outdir / "top_candidates.tsv", sep="\t", index=False)

    top_summary = []
    for k in top_ks:
        kk = min(k, len(ranked))
        if kk == 0:
            continue
        z = ranked.head(kk)
        top_summary.append({
            "k": kk,
            "median_priority_score": float(z["priority_score"].median()),
            "median_evidence_direction_coherence": float(z["evidence_direction_coherence"].median(skipna=True)),
            "median_evidence_support_fraction": float(z["evidence_support_fraction"].median(skipna=True)),
            "median_confidence_descriptor": float(z["confidence_descriptor"].median(skipna=True)) if z["confidence_descriptor"].notna().any() else None,
        })
    pd.DataFrame(top_summary).to_csv(outdir / "top_k_summary.tsv", sep="\t", index=False)

    manifest = {
        "version": VERSION,
        "status": "RANKED",
        "task_name": cfg.get("task_name"),
        "regime": regime,
        "input": str(input_path),
        "input_sha256": sha256_file(input_path),
        "config_sha256": sha256_file(config_path),
        "ranking_sha256": sha256_file(ranking_path),
        "n_input_rows": int(len(df)),
        "n_ranked_genes": int(len(ranked)),
        "top_k": top_ks,
        "method": method_meta,
        "important_note": "confidence_descriptor is an evidence-consistency descriptor, not a calibrated probability of correctness",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    print("PLANTVC_V1_RANK_SUCCESS")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    print("\nTOP 20")
    cols = ["priority_rank", "gene", "priority_score", "predicted_direction", "predicted_response", "evidence_direction_coherence", "evidence_support_fraction", "confidence_descriptor"]
    print(ranked[cols].head(20).to_string(index=False))
    return 0


def main() -> None:
    p = argparse.ArgumentParser(description="PlantVC v1 perturbation candidate-prioritization prototype")
    p.add_argument("--config", required=True, help="JSON task configuration")
    args = p.parse_args()
    config_path = resolve_path(args.config)
    try:
        rc = run(config_path)
    except Exception as e:
        print(f"PLANTVC_V1_ERROR: {type(e).__name__}: {e}", file=sys.stderr)
        raise
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
