from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external_rice" / "Robertson2026" / "pseudobulk_v0"
OUT = ROOT / "data" / "external_rice" / "Robertson2026" / "shared_structure_v0"
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(BASE / "group_metadata.tsv", sep="\t")
design = pd.read_csv(ROOT / "data" / "external_rice" / "Robertson2026" / "sample_design.tsv", sep="\t")
expr = pd.read_csv(BASE / "pseudobulk_log1p_cpm.tsv.gz", sep="\t", compression="gzip")
genes = expr.pop("gene").astype(str).to_numpy()
X = expr.to_numpy(dtype=float)
if expr.columns.astype(str).tolist() != meta["group"].astype(str).tolist():
    raise RuntimeError("Expression columns do not match group_metadata order")

mean_expr = X.mean(axis=1)
nonzero_frac = (X > 0).mean(axis=1)
mask = (nonzero_frac >= 0.25) & (mean_expr >= 0.05)

conditions = ["Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity"]
control = "Well.Watered"
annotations = sorted(meta["cell_annotation"].astype(str).unique())

batch_by_sample = dict(zip(design["identity"].astype(str), design["batch"].astype(int)))
control_sample_by_batch: dict[int, str] = {}
for sample in meta.loc[meta.condition == control, "sample"].astype(str).unique():
    b = batch_by_sample[sample]
    control_sample_by_batch[b] = sample

idx_lookup: dict[tuple[str, str], int] = {}
for i, row in meta.iterrows():
    idx_lookup[(str(row["sample"]), str(row["cell_annotation"]))] = int(i)

rows = []
for ann in annotations:
    controls = sorted(control_sample_by_batch.values())
    control_vecs = {s: X[:, idx_lookup[(s, ann)]][mask] for s in controls}
    cc_rms = []
    cc_rho = []
    cc_pairs = []
    for i, a in enumerate(controls):
        for b in controls[i+1:]:
            d = control_vecs[a] - control_vecs[b]
            cc_rms.append(float(np.sqrt(np.mean(d*d))))
            cc_rho.append(float(spearmanr(control_vecs[a], control_vecs[b]).statistic))
            cc_pairs.append(f"{a}-{b}")
    baseline_noise = float(np.median(cc_rms))
    baseline_expr_rho = float(np.median(cc_rho))

    for cond in conditions:
        stress_samples = sorted(meta.loc[meta.condition == cond, "sample"].astype(str).unique())
        resp_rms = []
        resp_pairs = []
        for ss in stress_samples:
            batch = batch_by_sample[ss]
            cs = control_sample_by_batch[batch]
            ti = idx_lookup[(ss, ann)]
            ci = idx_lookup[(cs, ann)]
            r = (X[:, ti] - X[:, ci])[mask]
            resp_rms.append(float(np.sqrt(np.mean(r*r))))
            resp_pairs.append(f"{ss}-{cs}:batch{batch}")
        mean_resp_rms = float(np.mean(resp_rms))
        median_resp_rms = float(np.median(resp_rms))
        rows.append({
            "condition": cond,
            "cell_annotation": ann,
            "baseline_control_control_rms_median": baseline_noise,
            "baseline_control_expression_rho_median": baseline_expr_rho,
            "stress_response_rms_mean": mean_resp_rms,
            "stress_response_rms_median": median_resp_rms,
            "response_to_baseline_noise_ratio": mean_resp_rms / baseline_noise if baseline_noise > 0 else np.nan,
            "matched_pairs": " | ".join(resp_pairs),
            "control_pairs": " | ".join(cc_pairs),
        })

out = pd.DataFrame(rows)
out.to_csv(OUT / "signal_vs_control_noise.tsv", sep="\t", index=False)

summ = out.groupby("condition").agg(
    median_response_to_noise=("response_to_baseline_noise_ratio", "median"),
    mean_response_to_noise=("response_to_baseline_noise_ratio", "mean"),
    median_stress_response_rms=("stress_response_rms_mean", "median"),
    median_baseline_control_control_rms=("baseline_control_control_rms_median", "median"),
    median_control_expression_rho=("baseline_control_expression_rho_median", "median"),
).reset_index()
summ.to_csv(OUT / "signal_vs_control_noise_summary.tsv", sep="\t", index=False)

print("ROBERTSON2026_SIGNAL_VS_CONTROL_NOISE_SUCCESS")
print(summ.to_string(index=False))
print("\nLOWEST_RESPONSE_TO_NOISE")
print(out.sort_values("response_to_baseline_noise_ratio").groupby("condition").head(4).to_string(index=False))
print("\nHIGHEST_RESPONSE_TO_NOISE")
print(out.sort_values("response_to_baseline_noise_ratio", ascending=False).groupby("condition").head(4).to_string(index=False))
