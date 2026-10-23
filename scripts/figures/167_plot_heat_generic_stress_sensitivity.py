"""ED2 v3: existing source-only generic-stress sensitivity, no raw-data analysis."""
from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs' / 'extended_data'
QA = OUT / 'qa'
DATA = ROOT / 'data/publication_controls_v2/heatshock_generic_stress_sensitivity'
sys.path.insert(0, str(ROOT / 'scripts' / 'utilities'))
from audit_panel_alignment import require_matplotlib_panel_alignment

ORDER = ['none', 'exclude_union_wrong_top100', 'exclude_recurrent_wrong_top100_ge2',
         'exclude_recurrent_wrong_top250_ge2', 'exclude_panstress_top100_mean_rank',
         'exclude_panstress_top250_mean_rank']
LABELS = ['No exclusion', 'Union non-heat Top100', 'Recurrent Top100\n(≥2 non-heat programs)',
          'Recurrent Top250\n(≥2 non-heat programs)', 'Pan-stress Top100', 'Pan-stress Top250']
SMALL = ['Union non-heat Top100', 'Recurrent Top100 (≥2)', 'Recurrent Top250 (≥2)',
         'Pan-stress Top100', 'Pan-stress Top250']
HEAT, GREY, INK = '#B7674C', '#878C91', '#303238'


def main():
    OUT.mkdir(exist_ok=True)
    QA.mkdir(exist_ok=True)
    df = pd.read_csv(DATA / 'heat_transfer_after_generic_stress_exclusion.tsv', sep='\t').set_index('sensitivity').loc[ORDER]
    baseline = json.loads((DATA / 'generic_stress_candidate_baseline.json').read_text(encoding='utf-8'))
    fractions = pd.read_csv(DATA / 'original_heat_overlap_generic_fraction.tsv', sep='\t').set_index('generic_definition').loc[ORDER[1:]]
    assert baseline['common_universe_genes'] == 2870 and len(baseline['wrong_programs']) == 7
    assert baseline['generic_top100_target_overlap'] == 6
    assert baseline['true_heat_top100_target_overlap_same_common_universe'] == df.loc['none', 'top100_overlap'] == 38
    assert df.top100_overlap.tolist() == [38, 43, 42, 42, 40, 40]
    assert df.excluded_genes.tolist() == [0, 365, 182, 428, 100, 250]
    assert df.top100_direction_agreement.eq(1).all()
    assert (df.excluded_genes + df.remaining_genes).eq(baseline['common_universe_genes']).all()
    assert fractions.original_heat_overlap_n.eq(38).all()
    assert np.allclose(fractions.fraction_original_overlap_flagged_generic,
                       fractions.overlap_genes_flagged_generic / fractions.original_heat_overlap_n)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                         'font.size': 11, 'axes.labelsize': 11, 'xtick.labelsize': 10.5,
                         'ytick.labelsize': 10.5, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
                         'axes.spines.right': False, 'axes.spines.top': False,
                         'axes.linewidth': 0.65, 'axes.edgecolor': '#686D72',
                         'text.color': INK, 'axes.labelcolor': INK,
                         'xtick.color': INK, 'ytick.color': INK, 'legend.frameon': False})
    fig = plt.figure(figsize=(11.2, 6.6), facecolor='white')
    a = fig.add_axes([0.245, 0.645, 0.245, 0.125])
    b = fig.add_axes([0.77, 0.27, 0.10, 0.50])
    c = fig.add_axes([0.245, 0.17, 0.245, 0.215])
    fig.text(0.035, 0.948, 'Cross-study heat recovery persists after generic-stress exclusion',
             fontsize=14, weight='bold', va='top')
    for letter, x, y, title in [('a', 0.035, 0.855, 'Heat versus generic stress'),
                                 ('b', 0.545, 0.855, 'Recovery after generic-candidate removal'),
                                 ('c', 0.035, 0.465, 'Generic component in the original overlap')]:
        fig.text(x, y, letter, fontsize=14, weight='bold')
        fig.text(x + 0.027, y, title, fontsize=11.5, weight='bold')
    fig.text(0.062, 0.81, f"Same {baseline['common_universe_genes']:,}-gene universe", fontsize=11, color='#646B70')
    for y, key, label, col in [(1, 'generic_top100_target_overlap', 'Generic stress Top100', GREY),
                                (0, 'true_heat_top100_target_overlap_same_common_universe', 'Frozen heat Top100', HEAT)]:
        value = baseline[key]
        a.hlines(y, 0, value, color=col, lw=1.5, alpha=0.7)
        a.plot(value, y, 'o', color=col, ms=7.0)
        a.text(value + 2, y, str(value), va='center', fontsize=12, color=col, weight='bold')
    a.set(xlim=(0, 50), ylim=(-0.65, 1.65), yticks=[1, 0],
          yticklabels=['Generic stress Top100', 'Frozen heat Top100'],
          xticks=[0, 10, 20, 30, 40, 50], xlabel='Target Top100 overlap')
    a.spines['left'].set_visible(False)
    a.tick_params(axis='y', length=0, pad=9)
    a.tick_params(axis='x', length=3, width=0.65, pad=4)
    fig.text(0.062, 0.535, 'Generic score: mean absolute-rank percentile\nacross seven source-only non-heat programs.',
             fontsize=10.5, color='#646B70', linespacing=1.5)

    ys = np.arange(6)[::-1]
    b.axvline(df.loc['none', 'top100_overlap'], color='#DDD7D3', lw=0.8, ls=(0, (3, 3)), zorder=0)
    for i, (y, row) in enumerate(zip(ys, df.itertuples())):
        b.plot(row.top100_overlap, y, marker='o', ms=6.3, linestyle='none',
               markerfacecolor=HEAT if i == 0 else 'white', markeredgecolor=HEAT, markeredgewidth=1.4)
        b.text(1.04, y, str(row.top100_overlap), transform=b.get_yaxis_transform(),
               va='center', color=HEAT, fontsize=11)
        fy = 0.27 + 0.50 * (y + 0.6) / 6.2
        fig.text(0.752, fy, LABELS[i], ha='right', va='center', fontsize=10.5, linespacing=1.25)
        fig.text(0.923, fy, f'{int(row.excluded_genes)}', ha='center', va='center', fontsize=10.5)
        fig.text(0.976, fy, f'{row.spearman_full_vector:.3f}', ha='center', va='center', fontsize=10.5)
    b.set(xlim=(0, 52), ylim=(-0.6, 5.6), yticks=[], xticks=[0, 20, 40], xlabel='Heat Top100 overlap')
    b.spines['left'].set_visible(False)
    b.tick_params(axis='x', length=3, width=0.65, pad=4)
    fig.text(0.923, 0.81, 'Removed', fontsize=10.5, ha='center', color='#646B70')
    fig.text(0.976, 0.81, 'Full ρ', fontsize=10.5, ha='center', color='#646B70')
    fig.text(0.572, 0.181, 'Overlap-direction agreement: 100% in every analysis.', fontsize=10.5)
    fig.text(0.572, 0.135, 'Dashed line: no-exclusion recovery (38 genes).', fontsize=10.5, color='#646B70')

    ysmall = np.arange(5)[::-1]
    frac = fractions.fraction_original_overlap_flagged_generic.to_numpy()
    c.hlines(ysmall, 0, frac, color='#D6D9DC', lw=1.1)
    c.plot(frac, ysmall, 'o', color=GREY, ms=4.5)
    for y, v, n in zip(ysmall, frac, fractions.overlap_genes_flagged_generic):
        c.text(v + 0.014, y, f'{v:.1%}', fontsize=10.5, va='center', color='#646B70')
    c.set(xlim=(0, 0.40), ylim=(-0.6, 4.6), yticks=ysmall, yticklabels=SMALL,
          xticks=[0, 0.1, 0.2, 0.3], xlabel='Original heat overlap flagged generic')
    c.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    c.spines['left'].set_visible(False)
    c.tick_params(axis='y', length=0, pad=9)
    c.tick_params(axis='x', length=3, width=0.65, pad=4)
    fig.text(0.035, 0.037, 'Exclusions change the gene universe and rerank Top100; higher counts are not evidence of improved transfer.',
             fontsize=10.5, color='#646B70')
    # The two Top100-overlap panels are semantically comparable but deliberately
    # occupy different vertical spans; c is a smaller contextual fraction strip.
    require_matplotlib_panel_alignment(fig, axes=[a, b, c], panel_ids=['a', 'b', 'c'],
        row_groups=[{'id': 'top-boundary', 'panels': ['a', 'b']}],
        column_groups=[{'id': 'left-column', 'panels': ['a', 'c']}],
        exemptions=[{'panels': ['a', 'b'], 'checks': ['row', 'panel-width'],
                     'reason': 'Two-row baseline and six-row sensitivity have intentionally different heights and widths; top plot edges are identical.'}],
        json_out=QA / 'ExtendedData2_v3_alignment.json', strict=True)
    base = OUT / 'ExtendedData2_v3_preview'
    fig.savefig(str(base) + '.png', dpi=600, facecolor='white')
    fig.savefig(str(base) + '.pdf', facecolor='white')
    fig.savefig(str(base) + '.svg', facecolor='white')
    fig.savefig(QA / 'ExtendedData2_v3_display.png', dpi=300, facecolor='white')
    plt.close(fig)
    b0, union = df.loc['none'], df.loc['exclude_union_wrong_top100']
    pct = ', '.join(f'{v:.1%}' for v in frac)
    legend = f'''**Extended Data Figure 2 | Cross-study heat recovery persists after generic-stress exclusion.**

**a,** Candidate recovery on the {baseline['common_universe_genes']:,}-gene universe shared by the frozen heat program, the independent seedling target and seven non-heat control programs. The data-driven generic-stress ranking uses the mean absolute-rank percentile across salt 3 h, salt 6 h, Pep1, SCOOP10, SCREW1, flg22 and Chitin. Its Top100 overlaps the independent heat-target Top100 by {baseline['generic_top100_target_overlap']} genes; the frozen heat Top100 overlaps by {baseline['true_heat_top100_target_overlap_same_common_universe']} on the same universe. **b,** Frozen-heat recovery after source-only generic-stress exclusions: the union of non-heat Top100 sets, recurrent non-heat Top100 or Top250 candidates (present in at least two programs), and pan-stress Top100 or Top250 candidates ranked by mean absolute-rank percentile. Removing the union excludes {int(union.excluded_genes)} genes and leaves {int(union.top100_overlap)} heat Top100 hits. All five exclusions retain {int(df.iloc[1:].top100_overlap.min())}–{int(df.iloc[1:].top100_overlap.max())} hits. Full-response Spearman correlations range from {df.spearman_full_vector.min():.3f} to {df.spearman_full_vector.max():.3f}; overlap-direction agreement is {df.top100_direction_agreement.min():.0%} in every analysis. The dashed line denotes the no-exclusion overlap. Exclusions alter the gene universe and rerank Top100, so higher overlap counts are not interpreted as improved transfer. **c,** Fractions of the original {int(b0.top100_overlap)} heat-overlap genes flagged by the five generic definitions, in the same order: {pct}. Generic-stress contributions are present, but candidate recovery persists after exclusion. Generic sets use source-only non-heat programs, not target heat expression or a manually curated HSP list. This analysis does not establish heat-specific mechanisms. No new statistical test was performed.
'''
    (OUT / 'ExtendedData2_v3_legend.md').write_text(legend, encoding='utf-8')
    report = {'baseline': baseline, 'sensitivities': df.reset_index().to_dict('records'),
              'original_overlap_fractions': fractions.reset_index().to_dict('records'),
              'source_rows_used': {'sensitivity': 6, 'original_overlap_fractions': 5},
              'new_statistical_tests': 0, 'source_only_generic_sets': True}
    (QA / 'ExtendedData2_v3_data_checks.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('ED2 rendered: same-universe baseline, all six sensitivities and five original-overlap fractions verified.')


if __name__ == '__main__':
    main()
