#!/usr/bin/env python3
"""Replicate- and control-contrast audit of Robertson 2026 shared response structure."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026"
PB = BASE / "pseudobulk"
OUT = BASE / "shared_structure_v0"
DESIGN = pd.read_csv(BASE / "sample_design.tsv", sep="\t")

CONDITION_MAP = {
    "WellWatered": "Well.Watered",
    "MildDrought": "Mild.Drought",
    "ModerateDrought": "Moderate.Drought",
    "MildSalinity": "Mild.Salinity",
    "ModerateSalinity": "Moderate.Salinity",
}
DESIGN["condition_obj"] = DESIGN["condition"].map(CONDITION_MAP)
CONDITIONS = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
RESOLUTIONS = ["CellAnnotation", "TissueSystem"]
FLOORS = [1, 50]


def rho(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if np.std(a) == 0 or np.std(b) == 0: return np.nan
    r = spearmanr(a, b).statistic
    return float(r) if np.isfinite(r) else np.nan


def cosine(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    return float(np.dot(a,b)/(na*nb)) if na > 0 and nb > 0 else np.nan


def structure(R):
    gram = R @ R.T
    vals = np.linalg.eigvalsh(gram)
    vals = np.sort(np.maximum(vals, 0))[::-1]
    frac = vals / vals.sum()
    pair = []
    loo = []
    for i in range(R.shape[0]):
        loo.append(rho(R[i], np.mean(np.delete(R, i, axis=0), axis=0)))
        for j in range(i+1, R.shape[0]):
            pair.append(rho(R[i], R[j]))
    return {
        "rank1": float(frac[0]),
        "rank3": float(np.cumsum(frac)[min(2, len(frac)-1)]),
        "median_pairwise": float(np.nanmedian(pair)),
        "median_loo": float(np.nanmedian(loo)),
    }


def load_resolution(label):
    prof = pd.read_csv(PB / f"pseudobulk_{label}_log1pCPM.tsv.gz", sep="\t", compression="gzip")
    prof.pop("gene")
    meta = pd.read_csv(PB / f"pseudobulk_{label}_meta.tsv", sep="\t", dtype={"replicate": str})
    seq_headers = [str(i) for i in range(1, len(prof.columns)+1)]
    if list(map(str, prof.columns)) == seq_headers and len(prof.columns) == len(meta):
        prof.columns = meta["pseudobulk_id"].astype(str).tolist()
    if set(prof.columns) != set(meta["pseudobulk_id"].astype(str)):
        raise RuntimeError(f"{label}: profile/meta mismatch")
    return prof.astype(float), meta


def col_for(meta, state, sample):
    z = meta[(meta.state.astype(str)==state)&(meta["sample"].astype(str)==sample)]
    if len(z) != 1: raise RuntimeError(f"{state}/{sample}: expected one pseudobulk, got {len(z)}")
    return str(z.iloc[0].pseudobulk_id)


def states_for(meta, samples, floor):
    z = meta[meta["sample"].isin(samples)]
    cov = z.pivot_table(index="state", columns="sample", values="n_cells", aggfunc="sum", fill_value=0)
    for s in samples:
        if s not in cov.columns: cov[s] = 0
    return sorted(cov.index[(cov[samples] >= floor).all(axis=1)].astype(str).tolist())


def matrix_difference(prof, meta, states, sample_a, sample_b):
    rows = []
    for st in states:
        a = prof[col_for(meta, st, sample_a)].to_numpy(float)
        b = prof[col_for(meta, st, sample_b)].to_numpy(float)
        rows.append(a-b)
    return np.vstack(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rep_rows = []
    control_rows = []
    ctrl = DESIGN[DESIGN.condition_obj == "Well.Watered"]
    ctrl_by_batch = dict(zip(ctrl.batch.astype(int), ctrl.identity.astype(str)))

    for label in RESOLUTIONS:
        prof, meta = load_resolution(label)
        for floor in FLOORS:
            # Control-vs-control structure for all three pairwise batch contrasts.
            control_samples = ctrl.sort_values("batch").identity.astype(str).tolist()
            for i in range(len(control_samples)):
                for j in range(i+1, len(control_samples)):
                    a, b = control_samples[i], control_samples[j]
                    states = states_for(meta, [a,b], floor)
                    R = matrix_difference(prof, meta, states, a, b)
                    s = structure(R)
                    control_rows.append({
                        "resolution": label, "cell_floor": floor,
                        "contrast": f"{a}_minus_{b}", "n_states": len(states), **s
                    })

            for cond in CONDITIONS:
                sd = DESIGN[DESIGN.condition_obj == cond].sort_values("replicate")
                if len(sd) != 2: raise RuntimeError(f"{cond}: expected 2 stress reps")
                d1, d2 = list(sd.itertuples(index=False))
                c1 = ctrl_by_batch[int(d1.batch)]; c2 = ctrl_by_batch[int(d2.batch)]
                needed = [d1.identity, c1, d2.identity, c2]
                states = states_for(meta, needed, floor)
                R1 = matrix_difference(prof, meta, states, d1.identity, c1)
                R2 = matrix_difference(prof, meta, states, d2.identity, c2)
                sh1 = R1.mean(axis=0); sh2 = R2.mean(axis=0)
                s1 = structure(R1); s2 = structure(R2)
                state_rep_rhos = [rho(R1[i], R2[i]) for i in range(len(states))]
                rep_rows.append({
                    "resolution": label, "condition": cond, "cell_floor": floor,
                    "n_states": len(states),
                    "shared_vector_rep_spearman": rho(sh1, sh2),
                    "shared_vector_rep_cosine": cosine(sh1, sh2),
                    "rank1_rep1": s1["rank1"], "rank1_rep2": s2["rank1"],
                    "median_pairwise_rep1": s1["median_pairwise"], "median_pairwise_rep2": s2["median_pairwise"],
                    "median_loo_rep1": s1["median_loo"], "median_loo_rep2": s2["median_loo"],
                    "median_state_response_rep_spearman": float(np.nanmedian(state_rep_rhos)),
                })

    rep = pd.DataFrame(rep_rows)
    ctl = pd.DataFrame(control_rows)
    rep.to_csv(OUT / "replicate_shared_axis_stability.tsv", sep="\t", index=False)
    ctl.to_csv(OUT / "control_control_structure.tsv", sep="\t", index=False)

    mean_ctl = ctl.groupby(["resolution","cell_floor"])[["rank1","median_pairwise","median_loo"]].median().reset_index()
    mean_ctl = mean_ctl.rename(columns={"rank1":"control_contrast_median_rank1","median_pairwise":"control_contrast_median_pairwise","median_loo":"control_contrast_median_loo"})
    out = rep.merge(mean_ctl, on=["resolution","cell_floor"], how="left")
    out.to_csv(OUT / "replicate_vs_control_audit.tsv", sep="\t", index=False)

    lines = [
        "# Robertson 2026 replicate/shared-axis audit",
        "",
        "This audit asks whether stress shared structure is reproducible across matched-batch biological replicates and how it compares with control-vs-control batch/biological contrasts.",
        "",
    ]
    for label in RESOLUTIONS:
        lines += [f"## {label}", ""]
        g = out[(out.resolution==label)&(out.cell_floor==1)]
        for r in g.itertuples(index=False):
            lines.append(
                f"- {r.condition}: shared-vector rep rho={r.shared_vector_rep_spearman:.3f}, cosine={r.shared_vector_rep_cosine:.3f}; "
                f"rank1 reps={r.rank1_rep1:.3f}/{r.rank1_rep2:.3f}; state-response rep median rho={r.median_state_response_rep_spearman:.3f}; "
                f"control-contrast median rank1={r.control_contrast_median_rank1:.3f}"
            )
        lines.append("")
    (OUT / "REPLICATE_AUDIT.md").write_text("\n".join(lines)+"\n", encoding="utf-8")

    print("ROBERTSON2026_REPLICATE_AUDIT_SUCCESS")
    print(out.to_string(index=False))
    print("\nCONTROL CONTRASTS")
    print(ctl.to_string(index=False))

if __name__ == "__main__":
    main()
