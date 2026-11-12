"""Render Fig. S42a from the timepoint scores and means produced by prepare_figures.py.

Usage: python render_s42a.py --input-dir DIR --output-dir DIR
"""
from pathlib import Path
import argparse
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
means = pd.read_csv(SOURCE/'dataset_mean_ARI.csv')
scores = pd.read_csv(SOURCE/'plotted_scores.csv')
palette = pd.read_csv(SOURCE/'method_palette.csv')
colors = dict(zip(palette.method, palette.hex))
methods = pd.read_csv(SOURCE/'method_display_order.csv').method.tolist()
assert len(means) == 24 and len(methods) == 8
panels = [('Mouse Embryo', 'Mouse Organogenesis', 3),
          ('Chicken Heart', 'Chicken Cardiogenesis', 4),
          ('Axolotl', 'Axolotl Telencephalon Neurogenesis', 3)]
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':19,
                     'pdf.fonttype':42, 'ps.fonttype':42})
fig, axes = plt.subplots(1, 3, figsize=(26, 8.4))
artists = []
for ax, (dataset, title, count) in zip(axes, panels):
    table = means[means.dataset_label == dataset].set_index('method').loc[methods]
    for method, row in table.iterrows():
        values = scores[(scores.dataset_label == dataset) & (scores.method == method)].ari.tolist()
        assert len(values) == count == row.timepoint_count
        assert abs(sum(values)/len(values)-row.mean_ari) < 1e-12
    values = table.mean_ari.to_numpy()
    ax.barh(range(8), values, height=.56,
            color=[colors[m] if m == 'SpaDOT' else 'white' for m in methods],
            edgecolor=[colors[m] for m in methods], linewidth=2.0, zorder=3)
    for y, value in enumerate(values):
        text = ax.text(max(0, value)+.015, y, f'{0 if abs(value)<.0005 else value:.3f}',
                       ha='left', va='center', fontsize=18)
        artists.append((ax, text))
    ax.set_yticks(range(8), labels=methods)
    ax.set_ylim(7.55, -.55)
    ax.set_xlim(-.035, .82)
    ax.set_xticks([0,.2,.4,.6,.8], labels=['0','0.2','0.4','0.6','0.8'])
    ax.set_xlabel('Mean ARI', fontsize=23, fontweight='bold', labelpad=10)
    title_artist = ax.set_title(title, fontsize=24, fontweight='bold', pad=22)
    ax.tick_params(axis='y', length=0, pad=9, labelsize=20)
    ax.tick_params(axis='x', labelsize=18, pad=5)
    ax.grid(color='#e9e9e9', linewidth=.75)
    ax.set_axisbelow(True)
    ax.axvline(0, color='#aaaaaa', linewidth=.75)
    for spine in ax.spines.values():
        spine.set_color('#777777'); spine.set_linewidth(.75)
fig.subplots_adjust(left=.07, right=.987, bottom=.235, top=.84, wspace=.39)
legend = fig.legend(handles=[
    Patch(facecolor=colors['SpaDOT'], edgecolor=colors['SpaDOT'], linewidth=2,
          label='SpaDOT (joint model)'),
    Patch(facecolor='white', edgecolor='#444444', linewidth=2,
          label='Single-slice methods')],
    loc='lower center', bbox_to_anchor=(.5,.025), ncol=2,
    frameon=False, fontsize=20, handlelength=2.0, handleheight=1.0,
    columnspacing=3.5, handletextpad=.8)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for ax, text in artists:
    box = text.get_window_extent(renderer)
    assert box.x1 <= ax.bbox.x1+1 and box.y0 >= ax.bbox.y0 and box.y1 <= ax.bbox.y1
for ax in axes:
    for text in [ax.title, ax.xaxis.label, *ax.get_yticklabels(), *ax.get_xticklabels()]:
        box = text.get_window_extent(renderer)
        assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.x1 and box.y1 <= fig.bbox.y1
pdf = ROOT / 'S42a_real_dataset_means_single_slice.pdf'
fig.savefig(pdf, metadata={'Title': 'Real-data mean ARI: SpaDOT and single-slice competitors'})
plt.close(fig)
