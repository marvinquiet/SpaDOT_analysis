"""Render Fig. S44 and count first-place ARIs across 45 scenarios.

Usage: python render_s44.py --input-dir DIR --output-dir DIR
Input tables come from prepare_figures.py simulation. The SpaDOT branch
variance readout is separate from the WOT and raw distance analyses.
"""
from pathlib import Path
import csv
import argparse
import statistics
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input-dir', type=Path, required=True)
parser.add_argument('--output-dir', type=Path, required=True)
args = parser.parse_args()
SOURCE, ROOT = args.input_dir, args.output_dir
ROOT.mkdir(parents=True, exist_ok=True)
METHODS = ['SpaDOT', 'PRECAST', 'SEDR', 'spaGCN', 'GraphST', 'spatialPCA', 'spatiAlign', 'spaVAE']
COLORS = dict(zip(METHODS, ['#c1121f', '#68a756', '#f09232', '#4c2623', '#165451', '#d3b300', '#ded6c1', '#0082c0']))
DEGS = [100, 200, 500]
PROPS = [.9, .8, .7, .6, .5]
seeds = pd.read_csv(SOURCE / 'ARI_timepoint_mean_per_seed.csv')
summary = pd.read_csv(SOURCE / 'ARI_plot_summary.csv')
scenario_keys = ['deg', 'proportion', 'seed']
assert len(seeds) == 360
assert not seeds.duplicated(scenario_keys + ['method']).any()
assert seeds.groupby(scenario_keys).method.apply(lambda x: set(x) == set(METHODS)).all()
winner_rows, rank_rows = [], []
for key, group in seeds.groupby(scenario_keys):
    ranked = group.sort_values('mean_ari_t0_t1_t2', ascending=False)
    first = ranked.iloc[0]
    top = float(first.mean_ari_t0_t1_t2)
    assert np.isclose(ranked.mean_ari_t0_t1_t2, top, rtol=0, atol=1e-12).sum() == 1
    winner_rows.append(dict(zip(scenario_keys, key), winner=first.method, mean_ari=top,
        runner_up=ranked.iloc[1].method,
        margin=top-float(ranked.iloc[1].mean_ari_t0_t1_t2)))
    for rank, (_, r) in enumerate(ranked.iterrows(), 1):
        rank_rows.append({**r.to_dict(), 'rank': rank})
winners = pd.DataFrame(winner_rows)
assert len(winners) == 45
counts = winners.winner.value_counts().reindex(METHODS, fill_value=0)
assert counts.sum() == 45

# Independently reconstruct each seed's mean and winner from the timepoint-level scores.
independent = {}
with (SOURCE / 'included_ARI_by_timepoint.csv').open() as fh:
    for r in csv.DictReader(fh):
        k = (int(r['deg']), float(r['proportion']), int(r['seed']), r['method'])
        independent.setdefault(k, []).append(float(r['ari']))
assert len(independent) == 360 and all(len(v) == 3 for v in independent.values())
for _, r in seeds.iterrows():
    k = (r.deg, r.proportion, r.seed, r.method)
    assert abs(statistics.mean(independent[k])-r.mean_ari_t0_t1_t2) < 1e-12
for _, r in winners.iterrows():
    vals = {m: statistics.mean(independent[(r.deg, r.proportion, r.seed, m)]) for m in METHODS}
    assert max(vals, key=vals.get) == r.winner
winners.to_csv(ROOT / 'first_place_per_scenario.csv', index=False)
pd.DataFrame(rank_rows).to_csv(ROOT / 'all_method_ranks_per_scenario.csv', index=False)
count_table = pd.DataFrame({'method': METHODS, 'first_place_scenarios': counts.values,
                            'total_scenarios': 45, 'fraction': counts.values/45})
for deg in DEGS:
    count_table[f'DEG{deg}_wins'] = [int(((winners.deg == deg) & (winners.winner == m)).sum()) for m in METHODS]
count_table.to_csv(ROOT / 'first_place_counts_45_scenarios.csv', index=False)
order = sorted(METHODS, key=lambda m: (-int(counts[m]), METHODS.index(m)))
caption = ('Joint SpaDOT uses post hoc GAT block variance balancing. Other methods use their '
           'native single-slice predictions. Each scenario is ranked using the equally weighted '
           'mean ARI across time points. The upper panels summarize three simulation seeds.')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 14, 'axes.titlesize': 23,
                    'axes.labelsize': 16, 'xtick.labelsize': 12, 'ytick.labelsize': 14,
                    'pdf.fonttype': 42, 'ps.fonttype': 42})

def first_place(ax, fontsize=16):
    for y, method in enumerate(order):
        n = int(counts[method])
        color = COLORS[method]
        ax.barh(y, n, height=.62, color=color,
                edgecolor='#8a8274' if method == 'spatiAlign' else color, linewidth=1)
        ax.text(n + .65, y, f'{n}/45', va='center', fontsize=fontsize)
    ax.set_yticks(range(len(order)), labels=order)
    ax.set_ylim(len(order)-.5, -.65)
    ax.set_xlim(0, 47.4)
    ax.set_xticks([0, 9, 18, 27, 36, 45])
    ax.set_xlabel('Number of scenarios ranked first', fontsize=18, labelpad=10)
    ax.set_title('First-place counts', fontsize=22, fontweight='bold', pad=16)
    ax.tick_params(axis='y', length=0, labelsize=fontsize)
    ax.tick_params(axis='x', labelsize=14)
    ax.grid(axis='x', color='#e8e8e8', linewidth=.7)
    ax.set_axisbelow(True)
    for side in ['right', 'top', 'left']:
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color('#999999')

fig = plt.figure(figsize=(20, 24))
gs = fig.add_gridspec(5, 3, left=.10, right=.985, top=.957, bottom=.265, hspace=.33, wspace=.29)
labels = []
drawn = []
for row, prop in enumerate(PROPS):
    for col, deg in enumerate(DEGS):
        ax = fig.add_subplot(gs[row, col])
        panel = seeds[(seeds.deg == deg) & np.isclose(seeds.proportion, prop)]
        for y, method in enumerate(METHODS):
            vals = panel[panel.method == method].sort_values('seed').mean_ari_t0_t1_t2.to_numpy()
            assert len(vals) == 3
            mean, sd = float(vals.mean()), float(vals.std(ddof=1))
            expected = summary[(summary.deg == deg) & np.isclose(summary.proportion, prop) & (summary.method == method)].iloc[0]
            assert np.isclose(mean, expected['mean'], atol=1e-12, rtol=0)
            assert np.isclose(sd, expected['std'], atol=1e-12, rtol=0)
            color = COLORS[method]
            edge = '#8a8274' if method == 'spatiAlign' else color
            ax.barh(y, mean, height=.55, facecolor=color if method == 'SpaDOT' else 'white',
                    edgecolor=edge, linewidth=1.5, xerr=sd, capsize=2.5,
                    error_kw={'elinewidth': .85, 'ecolor': '#555555'}, zorder=3)
            ax.scatter(vals, y + np.array([-.075, 0, .075]), s=11, c='#333333', alpha=.65, zorder=4)
            label = ax.text(max(mean+sd, float(vals.max()), 0)+.015, y, f'{mean:.3f}',
                            ha='left', va='center', fontsize=11.5, color='#333333')
            labels.append((ax, label))
            drawn.append({'deg': deg, 'proportion': prop, 'method': method, 'mean': mean, 'std': sd})
        ax.set_xlim(-.035, 1.16)
        ax.set_xticks([0, .2, .4, .6, .8, 1.], labels=['0', '0.2', '0.4', '0.6', '0.8', '1.0'])
        ax.set_yticks(range(8), labels=METHODS)
        ax.set_ylim(7.6, -.6)
        ax.set_xlabel('Mean ARI', fontweight='bold', labelpad=5)
        ax.grid(axis='x', color='#e8e8e8', linewidth=.7)
        ax.set_axisbelow(True)
        ax.axvline(0, color='#aaaaaa', lw=.6)
        ax.tick_params(axis='y', length=0, pad=6)
        for spine in ax.spines.values():
            spine.set_color('#aaaaaa'); spine.set_linewidth(.6)
        if row == 0:
            ax.set_title(f'DEG {deg}', fontweight='bold', pad=18)
        if col == 0:
            ax.set_ylabel(f'p = {prop:.1f}', fontsize=21, fontweight='bold', labelpad=18)
bottom = fig.add_axes([.10, .055, .49, .15])
first_place(bottom, fontsize=15)
fig.legend(handles=[Patch(facecolor=COLORS[m], edgecolor='#8a8274' if m=='spatiAlign' else COLORS[m], label=m)
                    for m in METHODS], loc='center', bbox_to_anchor=(.80, .13), ncol=2,
           frameon=False, fontsize=18, columnspacing=1.8, handlelength=1.4, labelspacing=1.0,
           title='Methods', title_fontsize=19)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for ax, label in labels:
    bbox = label.get_window_extent(renderer)
    assert bbox.x1 <= ax.bbox.x1 + 1 and bbox.y0 >= ax.bbox.y0 and bbox.y1 <= ax.bbox.y1
fig.savefig(ROOT / 'S44_Simulation_ARI_with_first_place_counts.pdf', metadata={'Title': 'Simulation ARI and first-place counts across 45 scenarios', 'Subject': caption})
plt.close(fig)
pd.DataFrame(drawn).to_csv(ROOT / 'plotted_bars.csv', index=False)

fig, ax = plt.subplots(figsize=(8.8, 6.4))
first_place(ax, fontsize=16)
fig.subplots_adjust(left=.24, right=.97, bottom=.16, top=.87)
fig.savefig(ROOT / 'Simulation_first_place_counts_45scenarios.pdf', metadata={'Title': 'Simulation first-place counts across 45 scenarios', 'Subject': caption})
plt.close(fig)
