#!/usr/bin/env python3
"""Continuous shared-vs-state-specific response analysis for Robertson 2026 rice stress.

No biological pass/fail threshold is used. The script reports structure across
stress types, severities, biological resolutions, and several rare-state cell
count sensitivity floors.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026"
PB = BASE / "pseudobulk"
OUT = BASE / "shared_structure_v0"
DESIGN = BASE / "sample_design.tsv"

CONDITIONS = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
CONTROL = "Well.Watered"
FLOORS = [1, 10, 20, 50, 100]
RESOLUTIONS = ["CellAnnotation", "TissueSystem"]
SEED = 20261020
NULL_PERMUTATIONS = 30


def rho(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.size < 2 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    r = spearmanr(a, b).statistic
    return float(r) if np.isfinite(r) else float("nan")


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    na = float(np.linalg.norm(a)); nb = float(np.linalg.norm(b))
    if na == 0 or nb == 0:
        return float("nan")
    return float(np.dot(a, b) / (na * nb))


def energy_spectrum(R: np.ndarray):
    gram = R @ R.T
    vals = np.linalg.eigvalsh(gram)
    vals = np.sort(np.maximum(vals, 0.0))[::-1]
    total = float(vals.sum())
    if total <= 0:
        return vals, np.full_like(vals, np.nan), np.full_like(vals, np.nan)
    frac = vals / total
    return vals, frac, np.cumsum(frac)


def null_rank1(R: np.ndarray, rng: np.random.Generator, nperm: int = NULL_PERMUTATIONS):
    out = []
    for _ in range(nperm):
        Rn = np.vstack([rng.permutation(row) for row in R])
        _, frac, _ = energy_spectrum(Rn)
        out.append(float(frac[0]))
    x = np.asarray(out, float)
    return float(np.median(x)), float(np.quantile(x, 0.95))


def load_resolution(label: str):
    p = PB / f"pseudobulk_{label}_log1pCPM.tsv.gz"
    m = PB / f"pseudobulk_{label}_meta.tsv"
    prof = pd.read_csv(p, sep="\t", compression="gzip")
    genes = prof.pop("gene").astype(str).tolist()
    meta = pd.read_csv(m, sep="\t", dtype={"replicate": str})
    # The first pseudobulk export lost column dimnames after sparse multiplication,
    # leaving deterministic 1..N headers. Recover them only when that exact pattern
    # is present and the column count matches the metadata row count/order.
    seq_headers = [str(i) for i in range(1, len(prof.columns) + 1)]
    if list(map(str, prof.columns)) == seq_headers and len(prof.columns) == len(meta):
        prof.columns = meta["pseudobulk_id"].astype(str).tolist()
    if set(prof.columns) != set(meta["pseudobulk_id"]):
        raise RuntimeError(f"{label}: profile/meta pseudobulk IDs differ")
    prof = prof.astype(float)
    return genes, prof, meta


def sample_columns(meta: pd.DataFrame, state: str, samples: list[str]) -> list[str]:
    z = meta[(meta["state"].astype(str) == state) & (meta["sample"].isin(samples))]
    d = dict(zip(z["sample"].astype(str), z["pseudobulk_id"].astype(str)))
    missing = [s for s in samples if s not in d]
    if missing:
        raise RuntimeError(f"{state}: missing samples {missing}")
    return [d[s] for s in samples]


def analyze_condition(label, prof, meta, design, condition, floor, pseudotime_map, rng):
    ctrl_samples = design.loc[design["condition_obj"] == CONTROL, "identity"].tolist()
    stress_samples = design.loc[design["condition_obj"] == condition, "identity"].tolist()
    if len(ctrl_samples) != 3 or len(stress_samples) != 2:
        raise RuntimeError(f"bad sample design for {condition}: control={ctrl_samples}, stress={stress_samples}")

    needed = ctrl_samples + stress_samples
    mm = meta[meta["sample"].isin(needed)].copy()
    coverage = mm.pivot_table(index="state", columns="sample", values="n_cells", aggfunc="sum", fill_value=0)
    for s in needed:
        if s not in coverage.columns:
            coverage[s] = 0
    states = sorted(coverage.index[(coverage[needed] >= floor).all(axis=1)].astype(str).tolist())
    if len(states) < 2:
        return None, [], [], None

    rows = []
    ctrl_by_state = {}
    stress_by_state = {}
    for st in states:
        ccols = sample_columns(meta, st, ctrl_samples)
        scols = sample_columns(meta, st, stress_samples)
        cmean = prof[ccols].to_numpy(float).mean(axis=1)
        smean = prof[scols].to_numpy(float).mean(axis=1)
        ctrl_by_state[st] = cmean
        stress_by_state[st] = smean
        rows.append(smean - cmean)
    R = np.vstack(rows)
    shared = R.mean(axis=0)
    shared_norm2 = float(np.dot(shared, shared))

    _, frac, cumulative = energy_spectrum(R)
    pairwise_vals = []
    state_rows = []
    for i, st in enumerate(states):
        y = R[i]
        loo = np.mean(np.delete(R, i, axis=0), axis=0)
        alpha = float(np.dot(y, shared) / shared_norm2) if shared_norm2 > 0 else np.nan
        proj = alpha * shared if np.isfinite(alpha) else np.full_like(y, np.nan)
        y_energy = float(np.dot(y, y))
        resid = y - proj
        resid_frac = float(np.dot(resid, resid) / y_energy) if y_energy > 0 else np.nan
        state_rows.append({
            "resolution": label,
            "condition": condition,
            "cell_floor": floor,
            "state": st,
            "rho_shared_LOO": rho(y, loo),
            "rho_shared_full": rho(y, shared),
            "cosine_shared": cosine(y, shared),
            "projection_energy_fraction": 1.0 - resid_frac if np.isfinite(resid_frac) else np.nan,
            "residual_energy_fraction": resid_frac,
            "response_l2": float(np.linalg.norm(y)),
            "mean_pseudotime": pseudotime_map.get(st, np.nan),
            "min_cells_required_samples": int(coverage.loc[st, needed].min()),
        })
    for i in range(len(states)):
        for j in range(i + 1, len(states)):
            pairwise_vals.append(rho(R[i], R[j]))

    # Replicate stability with batch-matched well-watered controls.
    rep_rows = []
    stress_design = design[design["condition_obj"] == condition].sort_values("replicate")
    ctrl_by_batch = dict(zip(
        design.loc[design["condition_obj"] == CONTROL, "batch"].astype(int),
        design.loc[design["condition_obj"] == CONTROL, "identity"].astype(str),
    ))
    if len(stress_design) == 2:
        d1, d2 = list(stress_design.itertuples(index=False))
        c1 = ctrl_by_batch[int(d1.batch)]; c2 = ctrl_by_batch[int(d2.batch)]
        for st in states:
            ids = {
                "s1": sample_columns(meta, st, [d1.identity])[0],
                "s2": sample_columns(meta, st, [d2.identity])[0],
                "c1": sample_columns(meta, st, [c1])[0],
                "c2": sample_columns(meta, st, [c2])[0],
            }
            y1 = prof[ids["s1"]].to_numpy(float) - prof[ids["c1"]].to_numpy(float)
            y2 = prof[ids["s2"]].to_numpy(float) - prof[ids["c2"]].to_numpy(float)
            rep_rows.append({
                "resolution": label,
                "condition": condition,
                "cell_floor": floor,
                "state": st,
                "rho_matched_batch_rep1_vs_rep2": rho(y1, y2),
                "sign_agreement_matched_batch": float(np.mean(np.sign(y1) == np.sign(y2))),
            })

    null_med = null_q95 = np.nan
    if floor == 1:
        null_med, null_q95 = null_rank1(R, rng)

    summary = {
        "resolution": label,
        "condition": condition,
        "cell_floor": floor,
        "n_states": len(states),
        "n_genes": R.shape[1],
        "rank1_energy_fraction": float(frac[0]),
        "rank2_cumulative_energy_fraction": float(cumulative[1]) if len(cumulative) > 1 else np.nan,
        "rank3_cumulative_energy_fraction": float(cumulative[2]) if len(cumulative) > 2 else np.nan,
        "median_pairwise_state_spearman": float(np.nanmedian(pairwise_vals)) if pairwise_vals else np.nan,
        "mean_rho_shared_LOO": float(np.nanmean([x["rho_shared_LOO"] for x in state_rows])),
        "median_rho_shared_LOO": float(np.nanmedian([x["rho_shared_LOO"] for x in state_rows])),
        "mean_projection_energy_fraction": float(np.nanmean([x["projection_energy_fraction"] for x in state_rows])),
        "median_projection_energy_fraction": float(np.nanmedian([x["projection_energy_fraction"] for x in state_rows])),
        "median_matched_batch_replicate_rho": float(np.nanmedian([x["rho_matched_batch_rep1_vs_rep2"] for x in rep_rows])) if rep_rows else np.nan,
        "mean_matched_batch_replicate_rho": float(np.nanmean([x["rho_matched_batch_rep1_vs_rep2"] for x in rep_rows])) if rep_rows else np.nan,
        "rank1_shuffle_null_median": null_med,
        "rank1_shuffle_null_q95": null_q95,
    }
    return summary, state_rows, rep_rows, shared


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    design = pd.read_csv(DESIGN, sep="\t")
    condition_map = {
        "WellWatered": "Well.Watered",
        "MildDrought": "Mild.Drought",
        "ModerateDrought": "Moderate.Drought",
        "MildSalinity": "Mild.Salinity",
        "ModerateSalinity": "Moderate.Salinity",
    }
    design["condition_obj"] = design["condition"].map(condition_map)
    if design["condition_obj"].isna().any():
        raise RuntimeError("sample design condition mapping failed")

    # Mean author pseudotime per annotation/system is descriptive context only.
    cm = pd.read_csv(PB / "cell_metadata.tsv.gz", sep="\t", compression="gzip")
    pt_maps = {
        "CellAnnotation": cm.groupby("CellAnnotation")["Pseudotime"].mean().to_dict(),
        "TissueSystem": cm.groupby("TissueSystem")["Pseudotime"].mean().to_dict(),
    }

    rng = np.random.default_rng(SEED)
    summary_rows = []
    state_rows = []
    rep_rows = []
    shared_vectors: dict[tuple[str, str], np.ndarray] = {}

    for label in RESOLUTIONS:
        genes, prof, meta = load_resolution(label)
        if len(genes) != 25985:
            raise RuntimeError(f"{label}: expected 25985 genes, got {len(genes)}")
        for floor in FLOORS:
            for cond in CONDITIONS:
                out = analyze_condition(label, prof, meta, design, cond, floor, pt_maps[label], rng)
                summary, sr, rr, shared = out
                if summary is None:
                    continue
                summary_rows.append(summary)
                state_rows.extend(sr)
                rep_rows.extend(rr)
                if floor == 1:
                    shared_vectors[(label, cond)] = shared

    summary_df = pd.DataFrame(summary_rows)
    state_df = pd.DataFrame(state_rows)
    rep_df = pd.DataFrame(rep_rows)
    summary_df.to_csv(OUT / "condition_structure_summary.tsv", sep="\t", index=False)
    state_df.to_csv(OUT / "state_structure_metrics.tsv", sep="\t", index=False)
    rep_df.to_csv(OUT / "replicate_stability.tsv", sep="\t", index=False)

    # Cross-condition similarity of shared response vectors, floor=1 only.
    cross_rows = []
    for label in RESOLUTIONS:
        for i, a in enumerate(CONDITIONS):
            for b in CONDITIONS[i+1:]:
                va = shared_vectors[(label, a)]; vb = shared_vectors[(label, b)]
                cross_rows.append({
                    "resolution": label,
                    "condition_a": a,
                    "condition_b": b,
                    "spearman_shared_vectors": rho(va, vb),
                    "cosine_shared_vectors": cosine(va, vb),
                })
    cross_df = pd.DataFrame(cross_rows)
    cross_df.to_csv(OUT / "cross_condition_shared_vector_similarity.tsv", sep="\t", index=False)

    # Developmental-position association with state-specific residuals, descriptive only.
    dev_rows = []
    d0 = state_df[(state_df["cell_floor"] == 1) & (state_df["resolution"] == "CellAnnotation")]
    for cond, g in d0.groupby("condition"):
        dev_rows.append({
            "condition": cond,
            "spearman_mean_pseudotime_vs_residual_energy": rho(g["mean_pseudotime"].to_numpy(float), g["residual_energy_fraction"].to_numpy(float)),
            "spearman_mean_pseudotime_vs_shared_rho": rho(g["mean_pseudotime"].to_numpy(float), g["rho_shared_LOO"].to_numpy(float)),
        })
    dev_df = pd.DataFrame(dev_rows)
    dev_df.to_csv(OUT / "pseudotime_structure_association.tsv", sep="\t", index=False)

    # Compact JSON at the no-extra-floor setting.
    primary = summary_df[summary_df["cell_floor"] == 1].copy()
    compact = {
        "status": "CONTINUOUS_EXTERNAL_STRUCTURE_TEST_NO_BIOLOGICAL_PASS_FAIL_GATE",
        "seed_for_descriptive_shuffle_null": SEED,
        "null_permutations": NULL_PERMUTATIONS,
        "conditions": {},
    }
    for r in primary.itertuples(index=False):
        compact["conditions"][f"{r.resolution}|{r.condition}"] = {
            "n_states": int(r.n_states),
            "rank1_energy_fraction": float(r.rank1_energy_fraction),
            "rank3_cumulative_energy_fraction": float(r.rank3_cumulative_energy_fraction),
            "median_pairwise_state_spearman": float(r.median_pairwise_state_spearman),
            "median_rho_shared_LOO": float(r.median_rho_shared_LOO),
            "median_projection_energy_fraction": float(r.median_projection_energy_fraction),
            "median_matched_batch_replicate_rho": float(r.median_matched_batch_replicate_rho),
            "rank1_shuffle_null_median": float(r.rank1_shuffle_null_median),
            "rank1_shuffle_null_q95": float(r.rank1_shuffle_null_q95),
        }
    (OUT / "summary.json").write_text(json.dumps(compact, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    # Human-readable report.
    lines = [
        "# Robertson 2026 rice stress — shared-vs-state-specific response structure",
        "",
        "No biological threshold is used as a pass/fail gate. Values below are continuous descriptive evidence.",
        "",
        "## CellAnnotation resolution (14 author cell annotations; no extra cell-count floor)",
        "",
        "| condition | states | rank1 energy | rank1-3 cumulative | median pairwise rho | median LOO-shared rho | median shared projection | median matched-batch replicate rho | shuffled rank1 median / q95 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    pc = primary[primary["resolution"] == "CellAnnotation"]
    for r in pc.itertuples(index=False):
        lines.append(
            f"| {r.condition} | {int(r.n_states)} | {r.rank1_energy_fraction:.3f} | {r.rank3_cumulative_energy_fraction:.3f} | "
            f"{r.median_pairwise_state_spearman:.3f} | {r.median_rho_shared_LOO:.3f} | {r.median_projection_energy_fraction:.3f} | "
            f"{r.median_matched_batch_replicate_rho:.3f} | {r.rank1_shuffle_null_median:.3f} / {r.rank1_shuffle_null_q95:.3f} |"
        )

    lines += ["", "## TissueSystem resolution (4 author tissue systems)", ""]
    pt = primary[primary["resolution"] == "TissueSystem"]
    for r in pt.itertuples(index=False):
        lines.append(
            f"- {r.condition}: rank1={r.rank1_energy_fraction:.3f}, rank1-3={r.rank3_cumulative_energy_fraction:.3f}, "
            f"median pairwise rho={r.median_pairwise_state_spearman:.3f}, median LOO rho={r.median_rho_shared_LOO:.3f}, "
            f"replicate rho={r.median_matched_batch_replicate_rho:.3f}"
        )

    lines += ["", "## Cross-condition similarity of shared vectors (CellAnnotation)", ""]
    for r in cross_df[cross_df["resolution"] == "CellAnnotation"].itertuples(index=False):
        lines.append(f"- {r.condition_a} vs {r.condition_b}: Spearman={r.spearman_shared_vectors:.3f}, cosine={r.cosine_shared_vectors:.3f}")

    lines += ["", "## Cell-count sensitivity", ""]
    for cond in CONDITIONS:
        g = summary_df[(summary_df["resolution"] == "CellAnnotation") & (summary_df["condition"] == cond)].sort_values("cell_floor")
        chunks = [f">={int(r.cell_floor)} cells: n={int(r.n_states)}, rank1={r.rank1_energy_fraction:.3f}, medLOO={r.median_rho_shared_LOO:.3f}" for r in g.itertuples(index=False)]
        lines.append(f"- {cond}: " + "; ".join(chunks))

    lines += ["", "## Developmental-position associations (CellAnnotation; descriptive)", ""]
    for r in dev_df.itertuples(index=False):
        lines.append(
            f"- {r.condition}: pseudotime vs residual energy rho={r.spearman_mean_pseudotime_vs_residual_energy:.3f}; "
            f"pseudotime vs shared-LOO rho={r.spearman_mean_pseudotime_vs_shared_rho:.3f}"
        )

    # most shared/divergent cell annotations at floor=1
    lines += ["", "## Most shared-like and most divergent cell annotations", ""]
    for cond in CONDITIONS:
        g = d0[d0["condition"] == cond]
        hi = g.nlargest(3, "rho_shared_LOO")
        lo = g.nsmallest(3, "rho_shared_LOO")
        lines.append(f"- {cond} shared-like: " + ", ".join(f"{r.state} ({r.rho_shared_LOO:.2f})" for r in hi.itertuples(index=False)))
        lines.append(f"  divergent: " + ", ".join(f"{r.state} ({r.rho_shared_LOO:.2f})" for r in lo.itertuples(index=False)))

    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("ROBERTSON2026_SHARED_STRUCTURE_SUCCESS")
    print(primary[["resolution","condition","n_states","rank1_energy_fraction","rank3_cumulative_energy_fraction","median_pairwise_state_spearman","median_rho_shared_LOO","median_projection_energy_fraction","median_matched_batch_replicate_rho","rank1_shuffle_null_median","rank1_shuffle_null_q95"]].to_string(index=False))
    print("\nCROSS CONDITION")
    print(cross_df.to_string(index=False))
    print("\nPSEUDOTIME")
    print(dev_df.to_string(index=False))


if __name__ == "__main__":
    main()
