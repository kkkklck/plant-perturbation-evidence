#!/usr/bin/env python3
"""Batch-paired shared-structure analysis for Robertson 2026 rice stress.

Each stress response is defined within batch:
    response_rep = stress_sample - well_watered_sample_from_same_batch
The condition response is the equal-weight mean of the two paired responses.

No biological pass/fail threshold is used. The script reports continuous
shared-vs-state-specific structure and replicate disagreement.
"""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026"
PB = BASE / "pseudobulk"
OUT = BASE / "shared_structure_v0"
DESIGN = pd.read_csv(BASE / "sample_design.tsv", sep="\t")

COND_MAP = {
    "WellWatered": "Well.Watered",
    "MildDrought": "Mild.Drought",
    "ModerateDrought": "Moderate.Drought",
    "MildSalinity": "Mild.Salinity",
    "ModerateSalinity": "Moderate.Salinity",
}
DESIGN["condition_obj"] = DESIGN["condition"].map(COND_MAP)
CONDITIONS = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
RESOLUTIONS = ["CellAnnotation", "TissueSystem"]
FLOORS = [1, 10, 20, 50, 100]


def rho(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    r = spearmanr(a, b).statistic
    return float(r) if np.isfinite(r) else np.nan


def cosine(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    return float(np.dot(a, b) / (na * nb)) if na > 0 and nb > 0 else np.nan


def structure(R):
    gram = R @ R.T
    vals = np.linalg.eigvalsh(gram)
    vals = np.sort(np.maximum(vals, 0))[::-1]
    frac = vals / vals.sum()
    cum = np.cumsum(frac)
    shared = R.mean(axis=0)
    shared_n2 = float(np.dot(shared, shared))
    pair = []
    rows = []
    for i in range(R.shape[0]):
        y = R[i]
        loo = np.mean(np.delete(R, i, axis=0), axis=0)
        alpha = float(np.dot(y, shared) / shared_n2) if shared_n2 > 0 else np.nan
        proj = alpha * shared if np.isfinite(alpha) else np.full_like(y, np.nan)
        e = float(np.dot(y, y))
        resid = y - proj
        resid_frac = float(np.dot(resid, resid) / e) if e > 0 else np.nan
        rows.append((rho(y, loo), 1.0 - resid_frac if np.isfinite(resid_frac) else np.nan))
        for j in range(i + 1, R.shape[0]):
            pair.append(rho(R[i], R[j]))
    return {
        "rank1": float(frac[0]),
        "rank2_cum": float(cum[min(1, len(cum)-1)]),
        "rank3_cum": float(cum[min(2, len(cum)-1)]),
        "median_pairwise": float(np.nanmedian(pair)),
        "median_loo": float(np.nanmedian([x[0] for x in rows])),
        "median_projection": float(np.nanmedian([x[1] for x in rows])),
        "shared": shared,
    }


def load_resolution(label):
    prof = pd.read_csv(PB / f"pseudobulk_{label}_log1pCPM.tsv.gz", sep="\t", compression="gzip")
    genes = prof.pop("gene").astype(str).tolist()
    meta = pd.read_csv(PB / f"pseudobulk_{label}_meta.tsv", sep="\t", dtype={"replicate": str})
    seq_headers = [str(i) for i in range(1, len(prof.columns)+1)]
    if list(map(str, prof.columns)) == seq_headers and len(prof.columns) == len(meta):
        prof.columns = meta["pseudobulk_id"].astype(str).tolist()
    if set(prof.columns) != set(meta["pseudobulk_id"].astype(str)):
        raise RuntimeError(f"{label}: profile/meta mismatch")
    return genes, prof.astype(float), meta


def col_for(meta, state, sample):
    z = meta[(meta.state.astype(str) == state) & (meta["sample"].astype(str) == sample)]
    if len(z) != 1:
        raise RuntimeError(f"{state}/{sample}: expected 1 pseudobulk, got {len(z)}")
    return str(z.iloc[0].pseudobulk_id)


def states_for(meta, samples, floor):
    z = meta[meta["sample"].isin(samples)]
    cov = z.pivot_table(index="state", columns="sample", values="n_cells", aggfunc="sum", fill_value=0)
    for s in samples:
        if s not in cov.columns: cov[s] = 0
    return sorted(cov.index[(cov[samples] >= floor).all(axis=1)].astype(str).tolist()), cov


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ctrl = DESIGN[DESIGN.condition_obj == "Well.Watered"]
    ctrl_by_batch = dict(zip(ctrl.batch.astype(int), ctrl.identity.astype(str)))

    summary_rows = []
    state_rows = []
    cross_vectors = {}

    for label in RESOLUTIONS:
        genes, prof, meta = load_resolution(label)
        for floor in FLOORS:
            for cond in CONDITIONS:
                sd = DESIGN[DESIGN.condition_obj == cond].sort_values("replicate")
                d1, d2 = list(sd.itertuples(index=False))
                c1 = ctrl_by_batch[int(d1.batch)]
                c2 = ctrl_by_batch[int(d2.batch)]
                samples = [d1.identity, c1, d2.identity, c2]
                states, cov = states_for(meta, samples, floor)
                if len(states) < 2:
                    continue

                R1 = []
                R2 = []
                for st in states:
                    y1 = prof[col_for(meta, st, d1.identity)].to_numpy(float) - prof[col_for(meta, st, c1)].to_numpy(float)
                    y2 = prof[col_for(meta, st, d2.identity)].to_numpy(float) - prof[col_for(meta, st, c2)].to_numpy(float)
                    R1.append(y1); R2.append(y2)
                R1 = np.vstack(R1); R2 = np.vstack(R2)
                R = (R1 + R2) / 2.0
                S = structure(R)

                # signal/noise decomposition for two paired responses:
                # mean m=(r1+r2)/2, disagreement d=(r1-r2)/2.
                # ||m||^2 versus ||d||^2 gives a transparent scale-free stability diagnostic.
                D = (R1 - R2) / 2.0
                signal_energy = np.sum(R * R, axis=1)
                disagreement_energy = np.sum(D * D, axis=1)
                ratios = signal_energy / (disagreement_energy + 1e-12)
                rep_rhos = [rho(R1[i], R2[i]) for i in range(len(states))]

                shared1 = R1.mean(axis=0)
                shared2 = R2.mean(axis=0)
                summary_rows.append({
                    "resolution": label,
                    "condition": cond,
                    "cell_floor": floor,
                    "n_states": len(states),
                    "n_genes": len(genes),
                    "rank1_energy_fraction": S["rank1"],
                    "rank2_cumulative_energy_fraction": S["rank2_cum"],
                    "rank3_cumulative_energy_fraction": S["rank3_cum"],
                    "median_pairwise_state_spearman": S["median_pairwise"],
                    "median_rho_shared_LOO": S["median_loo"],
                    "median_projection_energy_fraction": S["median_projection"],
                    "shared_vector_rep_spearman": rho(shared1, shared2),
                    "shared_vector_rep_cosine": cosine(shared1, shared2),
                    "median_state_response_rep_spearman": float(np.nanmedian(rep_rhos)),
                    "global_signal_to_disagreement_energy": float(signal_energy.sum() / (disagreement_energy.sum() + 1e-12)),
                    "median_state_signal_to_disagreement_energy": float(np.nanmedian(ratios)),
                })
                for i, st in enumerate(states):
                    state_rows.append({
                        "resolution": label,
                        "condition": cond,
                        "cell_floor": floor,
                        "state": st,
                        "min_cells_across_required_samples": int(cov.loc[st, samples].min()),
                        "replicate_response_spearman": rep_rhos[i],
                        "signal_energy": float(signal_energy[i]),
                        "disagreement_energy": float(disagreement_energy[i]),
                        "signal_to_disagreement_energy": float(ratios[i]),
                        "rho_to_condition_shared": rho(R[i], S["shared"]),
                    })
                if floor == 1:
                    cross_vectors[(label, cond)] = S["shared"]

    summary = pd.DataFrame(summary_rows)
    state = pd.DataFrame(state_rows)
    summary.to_csv(OUT / "paired_effect_structure_summary.tsv", sep="\t", index=False)
    state.to_csv(OUT / "paired_effect_state_metrics.tsv", sep="\t", index=False)

    cross = []
    for label in RESOLUTIONS:
        for i, a in enumerate(CONDITIONS):
            for b in CONDITIONS[i+1:]:
                va = cross_vectors[(label, a)]; vb = cross_vectors[(label, b)]
                cross.append({
                    "resolution": label,
                    "condition_a": a,
                    "condition_b": b,
                    "spearman_shared_vectors": rho(va, vb),
                    "cosine_shared_vectors": cosine(va, vb),
                })
    cross = pd.DataFrame(cross)
    cross.to_csv(OUT / "paired_effect_cross_condition_similarity.tsv", sep="\t", index=False)

    base = summary[summary.cell_floor == 1].copy()
    compact = {
        "status": "BATCH_PAIRED_CONTINUOUS_STRUCTURE_ANALYSIS_NO_BIOLOGICAL_PASS_FAIL_GATE",
        "conditions": {},
    }
    for r in base.itertuples(index=False):
        compact["conditions"][f"{r.resolution}|{r.condition}"] = {
            "n_states": int(r.n_states),
            "rank1_energy_fraction": float(r.rank1_energy_fraction),
            "rank3_cumulative_energy_fraction": float(r.rank3_cumulative_energy_fraction),
            "median_rho_shared_LOO": float(r.median_rho_shared_LOO),
            "shared_vector_rep_spearman": float(r.shared_vector_rep_spearman),
            "global_signal_to_disagreement_energy": float(r.global_signal_to_disagreement_energy),
            "median_state_signal_to_disagreement_energy": float(r.median_state_signal_to_disagreement_energy),
        }
    (OUT / "paired_effect_summary.json").write_text(json.dumps(compact, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Robertson 2026 batch-paired stress-effect structure",
        "",
        "Each stress replicate is contrasted with the well-watered sample from the same batch; the two paired effects are then averaged.",
        "No biological pass/fail threshold is used.",
        "",
        "## CellAnnotation (no extra cell-count floor)",
        "",
        "| condition | states | rank1 | rank1-3 | median LOO | shared-axis rep rho | state-response rep rho | global signal/disagreement | median state signal/disagreement |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in base[base.resolution == "CellAnnotation"].itertuples(index=False):
        lines.append(
            f"| {r.condition} | {int(r.n_states)} | {r.rank1_energy_fraction:.3f} | {r.rank3_cumulative_energy_fraction:.3f} | "
            f"{r.median_rho_shared_LOO:.3f} | {r.shared_vector_rep_spearman:.3f} | {r.median_state_response_rep_spearman:.3f} | "
            f"{r.global_signal_to_disagreement_energy:.3f} | {r.median_state_signal_to_disagreement_energy:.3f} |"
        )
    lines += ["", "## Cell-count sensitivity", ""]
    for cond in CONDITIONS:
        g = summary[(summary.resolution == "CellAnnotation") & (summary.condition == cond)].sort_values("cell_floor")
        chunks = [f">={int(r.cell_floor)}: n={int(r.n_states)}, rank1={r.rank1_energy_fraction:.3f}, LOO={r.median_rho_shared_LOO:.3f}, S/N={r.global_signal_to_disagreement_energy:.2f}" for r in g.itertuples(index=False)]
        lines.append(f"- {cond}: " + "; ".join(chunks))
    lines += ["", "## Paired shared-vector cross-condition similarity (CellAnnotation)", ""]
    for r in cross[cross.resolution == "CellAnnotation"].itertuples(index=False):
        lines.append(f"- {r.condition_a} vs {r.condition_b}: rho={r.spearman_shared_vectors:.3f}, cosine={r.cosine_shared_vectors:.3f}")
    (OUT / "PAIRED_EFFECT_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("ROBERTSON2026_PAIRED_EFFECT_SUCCESS")
    print(base.to_string(index=False))
    print("\nCROSS")
    print(cross.to_string(index=False))


if __name__ == "__main__":
    main()
