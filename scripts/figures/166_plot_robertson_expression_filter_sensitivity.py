"""ED1 v3: source-driven expression-filter robustness; never writes source data or Word."""
from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs' / 'extended_data'
QA = OUT / 'qa'
DATA = ROOT / 'data/publication_controls_v2/robertson_wellexpressed_fullvector'
sys.path.insert(0, str(ROOT / 'scripts' / 'utilities'))
from audit_panel_alignment import require_matplotlib_panel_alignment

ORDER = ['technical_universe', 'control_mean_gt_0.5', 'control_mean_gt_1.0',
         'control_top50pct_within_technical', 'control_top25pct_within_technical']
LABELS = ['Technical\nuniverse', 'Control mean\n>0.5', 'Control mean\n>1.0',
          'Top 50%\ncontrol\nexpression', 'Top 25%\ncontrol\nexpression']
SHORT = ['Technical', 'Mean >0.5', 'Mean >1.0', 'Top 50%', 'Top 25%']
COLORS = {'Mild.Drought': '#4C718C', 'Moderate.Drought': '#88A9BE',
          'Mild.Salinity': '#98523B', 'Moderate.Salinity': '#B7674C'}
INK, GREY = '#303238', '#74797E'


def main():
    OUT.mkdir(exist_ok=True)
    QA.mkdir(exist_ok=True)
    df = pd.read_csv(DATA / 'context_level_fullvector_sensitivity.tsv', sep='\t')
    summary = pd.read_csv(DATA / 'overall_summary.tsv', sep='\t').set_index('universe').loc[ORDER]
    assert len(df) == 80 and not df.duplicated(['condition', 'tissue_system', 'universe']).any()
    assert set(df.condition) == set(COLORS)
    paired = df.pivot(index=['condition', 'tissue_system'], columns='universe',
                      values='replicate_spearman_rho').loc[:, ORDER]
    assert paired.shape == (16, 5) and paired.notna().all().all()
    assert (paired.lt(0.3).sum().to_numpy() == summary.n_rho_lt_0_3.to_numpy()).all()
    assert summary.n_rho_lt_0_3.tolist() == [16, 16, 16, 14, 14]
    assert np.allclose(paired.median().to_numpy(), summary.median_rho.to_numpy())
    quartiles = paired.quantile([0.25, 0.75])
    assert np.allclose(quartiles.loc[0.25], summary.q25_rho)
    assert np.allclose(quartiles.loc[0.75], summary.q75_rho)
    assert paired.min().min() >= -0.25 and paired.max().max() <= 0.45

    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                         'font.size': 11, 'axes.labelsize': 11, 'xtick.labelsize': 10.5,
                         'ytick.labelsize': 10.5, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
                         'axes.spines.right': False, 'axes.spines.top': False,
                         'axes.linewidth': 0.65, 'axes.edgecolor': '#686D72',
                         'text.color': INK, 'axes.labelcolor': INK,
                         'xtick.color': INK, 'ytick.color': INK, 'legend.frameon': False})
    fig = plt.figure(figsize=(11.2, 6.2), facecolor='white')
    a = fig.add_axes([0.075, 0.255, 0.47, 0.48])
    b = fig.add_axes([0.735, 0.255, 0.135, 0.48])
    fig.text(0.035, 0.947, 'Rice full-vector concordance remains limited after expression filtering',
             fontsize=14, weight='bold', va='top')
    for letter, x, title in [('a', 0.035, 'Context-level expression-filter sensitivity'),
                              ('b', 0.60, 'Overall concordance remains low')]:
        fig.text(x, 0.855, letter, weight='bold', fontsize=14)
        fig.text(x + 0.027, 0.855, title, fontsize=11.5, weight='bold')
    handles = [Line2D([0], [0], color=c, marker='o', markersize=4, lw=1,
                      label=k.replace('.', ' ').capitalize()) for k, c in COLORS.items()]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.066, 0.825), ncol=2,
               handlelength=1.5, columnspacing=1.8, labelspacing=0.45, fontsize=10.5)
    for (condition, tissue), row in paired.iterrows():
        # Exact context pairing; every source row appears once, without jitter or smoothing.
        a.plot(np.arange(5), row.to_numpy(), color=COLORS[condition], alpha=0.64,
               lw=0.85, marker='o', ms=3.2, markeredgewidth=0, zorder=3)
    a.axhline(0.3, color='#D9DDE0', lw=0.8, ls=(0, (3, 3)), zorder=0)
    a.axhline(0, color='#E8EAEC', lw=0.6, zorder=0)
    a.set(xlim=(-0.15, 4.15), ylim=(-0.25, 0.45), xticks=range(5), xticklabels=LABELS,
          yticks=[-0.2, 0, 0.2, 0.3, 0.4], ylabel='Replicate full-response Spearman ρ')
    a.tick_params(length=3, width=0.65, pad=5)
    fig.text(0.075, 0.130, '16 stress × tissue-system contexts; two biological replicates per context',
             fontsize=10.5, color=GREY)
    fig.text(0.075, 0.075, 'Dashed line: ρ = 0.3 (descriptive reference). Negative values are retained.',
             fontsize=10.5, color=GREY)

    ys = np.arange(5)[::-1]
    med = summary.median_rho.to_numpy()
    lo, hi = quartiles.loc[0.25].to_numpy(), quartiles.loc[0.75].to_numpy()
    b.hlines(ys, lo, hi, color='#939BA2', lw=2.1, zorder=2)
    b.plot(med, ys, 'o', color='#374957', ms=5.8, zorder=3)
    b.axvline(0, color='#E5E8EB', lw=0.65, zorder=0)
    b.set(xlim=(-0.11, 0.25), ylim=(-0.65, 4.65), xticks=[-0.1, 0.0, 0.2], yticks=[],
          xlabel='Median ρ (IQR)')
    b.spines['left'].set_visible(False)
    b.tick_params(axis='x', length=3, width=0.65, pad=5)
    fig.text(0.71, 0.785, 'Gene universe', fontsize=10.5, color=GREY, ha='right')
    fig.text(0.907, 0.785, 'Median', fontsize=10.5, color=GREY, ha='center')
    fig.text(0.965, 0.785, 'ρ < 0.3', fontsize=10.5, color=GREY, ha='center')
    for y, label, m, n in zip(ys, SHORT, med, summary.n_rho_lt_0_3):
        fy = 0.255 + 0.48 * (y + 0.65) / 5.3
        fig.text(0.71, fy, label, va='center', ha='right', fontsize=11)
        fig.text(0.907, fy, f'{m:.3f}', va='center', ha='center', fontsize=11)
        fig.text(0.965, fy, f'{int(n)}/16', va='center', ha='center', fontsize=11)
    fig.text(0.627, 0.130, 'Moderate filters: all 16 contexts remain below ρ = 0.3.', fontsize=10.5)
    fig.text(0.627, 0.075, 'Aggressive filters: context changes, low overall median.', fontsize=10.5)

    require_matplotlib_panel_alignment(fig, axes=[a, b], panel_ids=['a', 'b'],
        row_groups=[{'id': 'matched-plot-edges', 'panels': ['a', 'b']}],
        exemptions=[{'panels': ['a', 'b'], 'checks': ['panel-width'],
                     'reason': 'All-context trajectory hero deliberately wider than compact IQR summary.'}],
        json_out=QA / 'ExtendedData1_v3_alignment.json', strict=True)
    base = OUT / 'ExtendedData1_v3_preview'
    fig.savefig(str(base) + '.png', dpi=600, facecolor='white')
    fig.savefig(str(base) + '.pdf', facecolor='white')
    fig.savefig(str(base) + '.svg', facecolor='white')
    fig.savefig(QA / 'ExtendedData1_v3_display.png', dpi=300, facecolor='white')
    plt.close(fig)
    values = ', '.join(f'{v:.3f}' for v in med)
    legend = f'''**Extended Data Figure 1 | Rice full-vector concordance remains limited after expression filtering.**

**a,** Replicate full-response Spearman correlations for all 16 Robertson rice stress × tissue-system contexts across five gene universes, in order: the technical universe, genes with control mean >0.5 or >1.0, and the top 50% or top 25% of genes ranked by control expression within each technical universe. Each context contains two biological replicates. Lines connect the same context across filters; colours denote the four treatments. The horizontal dashed line at ρ = 0.3 is a descriptive reference, not a fitted or significance threshold. **b,** Overall median (points) and interquartile range (horizontal lines; Q25–Q75) across the 16 contexts. Median ρ values are {values}, respectively. The annotation column reports the number of contexts with ρ < 0.3. Both moderate filters leave all 16 contexts below this reference, whereas the top-50% and top-25% filters each leave 14/16 below it. Aggressive filtering changes individual contexts, including increases in some mild-drought contexts, while the overall median remains low. The results do not support a simple low-expression-noise explanation for limited full-vector concordance. No new statistical test was performed.
'''
    (OUT / 'ExtendedData1_v3_legend.md').write_text(legend, encoding='utf-8')
    report = {'source_rows': len(df), 'contexts': len(paired), 'paired_universe_order': ORDER,
              'median_source': 'overall_summary.tsv', 'summary': summary.reset_index().to_dict('records'),
              'trajectories': paired.reset_index().to_dict('records'),
              'range_retained': [float(paired.min().min()), float(paired.max().max())],
              'new_statistical_tests': 0, 'excluded_contexts': 0, 'treatment_colors': COLORS}
    (QA / 'ExtendedData1_v3_data_checks.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('ED1 rendered: all 80 source observations, 16 complete trajectories, matched summary verified.')


if __name__ == '__main__':
    main()
