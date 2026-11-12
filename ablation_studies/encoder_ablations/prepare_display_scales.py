"""Display-only ARI restyling from verified, unchanged CSV result tables."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'spadot_grouped_ari_mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input-dir', required=True, type=Path,
                    help='Combined audited metrics and verification files')
parser.add_argument('--output-dir', required=True, type=Path)
args = parser.parse_args()
SOURCE, DESTINATION = args.input_dir, args.output_dir
CASES = ('ChickenHeart', 'AxolotlTelencephalon', 'MouseEmbryo')
TITLES = ('Chicken heart', 'Axolotl telencephalon', 'Mouse embryo')
STAGES = (('Day 4', 'Day 7', 'Day 10', 'Day 14'), ('St.44', 'St.54', 'St.57'), ('E12.5', 'E13.5', 'E15.5'))
ARMS = ('dual_native', 'gat_d1', 'gat_d2', 'gat_d3', 'gat_d4', 'gat_d6',
        'gat_d3_parammatched', 'gat_d3_expanded', 'gat_gat_matched', 'mlp_gat')
LABELS = ('SpaDOT GP+GAT', 'GAT20 depth 1', 'GAT20 depth 2', 'GAT20 depth 3',
          'GAT20 depth 4 (extension)', 'GAT20 depth 6', 'GAT20 D3 matched',
          'GAT20 D3 expanded', 'GAT10+GAT10 matched', 'MLP10+GAT10')
GROUPS = (
    (0, 1, 'SpaDOT\nreference', '#40546B', 'SpaDOT reference'),
    (1, 6, 'GAT-only\ndepth', '#6B5C91', 'GAT-only depth'),
    (6, 7, 'GAT-only\ncapacity', '#396D82', 'GAT-only capacity control'),
    (7, 8, 'GAT-only\nneighbors', '#3B7B70', 'GAT-only neighborhood control'),
    (8, 10, 'GP-branch\nreplacement', '#6F6680', 'GP-branch replacement'),
)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none'})
CMAP = plt.get_cmap('YlOrRd').copy()
CMAP.set_bad('#E2E8F0')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def text_color(rgba):
    rgb = np.asarray(rgba[:3])
    linear = np.where(rgb <= .04045, rgb / 12.92, ((rgb + .055) / 1.055) ** 2.4)
    luminance = float(linear @ np.array([.2126, .7152, .0722]))
    return '#172438' if luminance > .21 else 'white'


def rules(ax):
    for end in (0.5, 5.5, 6.5, 7.5):
        ax.axhline(end, color='white', lw=1.7)


def save(fig, stem):
    for extension in ('png', 'svg'):
        fig.savefig(DESTINATION / (stem + '.' + extension), dpi=150, facecolor='white')
    plt.close(fig)


def draw_main(primary, means, matrices, norms):
    fig = plt.figure(figsize=(24, 9.4))
    bottom, height = .235, .625
    group_ax = fig.add_axes([.012, bottom, .067, height])
    group_ax.set(xlim=(0, 1), ylim=(9.5, -.5))
    group_ax.axis('off')
    for start, end, label, color, _ in GROUPS:
        group_ax.add_patch(Rectangle((0, start - .5), 1, end - start, facecolor=color, edgecolor='white', lw=2))
        group_ax.text(.5, (start + end - 1) / 2, label, ha='center', va='center',
                      color='white', fontsize=10.5, linespacing=1.25)
    group_ax.set_title('Experiment\ncategory', fontsize=11, pad=13)
    labels_ax = fig.add_axes([.087, bottom, .183, height])
    labels_ax.set(xlim=(0, 1), ylim=(9.5, -.5))
    labels_ax.axis('off')
    for row, label in enumerate(LABELS):
        labels_ax.text(.985, row, label, ha='right', va='center', fontsize=11)
    left, total, gap = .285, .695, .027
    usable = total - 2 * gap
    widths = np.array([5, 4, 4], dtype=float) / 13 * usable
    for case, title, stages, width in zip(CASES, TITLES, STAGES, widths):
        norm, values = norms[case], matrices[case]
        cell = (width - .009) / (len(stages) + 1)
        stage_width = cell * len(stages)
        main = fig.add_axes([left, bottom, stage_width, height])
        mean = fig.add_axes([left + stage_width + .009, bottom, cell, height])
        for ax, data, ticks in ((main, values[:, :-1], stages), (mean, values[:, -1:], ('Mean',))):
            ax.imshow(np.ma.masked_invalid(data), cmap=CMAP, norm=norm, aspect='auto')
            for i in range(len(ARMS)):
                for j in range(data.shape[1]):
                    value = data[i, j]
                    ax.text(j, i, '—' if not np.isfinite(value) else f'{value:.3f}',
                            ha='center', va='center', fontsize=10.5,
                            color='#526174' if not np.isfinite(value) else text_color(CMAP(norm(value))))
            ax.set_xticks(range(data.shape[1]), ticks, fontsize=11, rotation=25, ha='right')
            ax.set_yticks([])
            rules(ax)
            for spine in ax.spines.values():
                spine.set_color('#526174')
                spine.set_linewidth(.8)
        fig.text(left + width / 2, bottom + height + .029, title, ha='center', fontsize=17)
        cb_ax = fig.add_axes([left, .13, width, .02])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cb_ax, orientation='horizontal')
        cb.set_ticks(np.linspace(norm.vmin, norm.vmax, 5))
        cb.ax.tick_params(labelsize=10)
        cb.ax.set_xticklabels([f'{x:.3f}' for x in np.linspace(norm.vmin, norm.vmax, 5)])
        cb.set_label('ARI - dataset-specific scale', fontsize=10, labelpad=5)
        left += width + gap
    fig.suptitle('Domain ARI by timepoint | training seed 2026 | executed epoch 99', fontsize=21, y=.982)
    fig.text(.54, .929, 'Yellow = lower ARI; red = higher ARI. Each dataset has its own scale; its mean uses the same scale.',
             ha='center', fontsize=12, color='#526174')
    fig.text(.53, .035, 'D4 is a later extension. Gray = missing or failed, never zero. Means require every timepoint.',
             ha='center', fontsize=11, color='#526174')
    fig.text(.53, .012, 'Category bars describe experimental design; they do not introduce group averages. Compare printed ARIs across datasets.',
             ha='center', fontsize=10, color='#526174')
    save(fig, 'primary_ARI_by_timepoint')


def draw_metric_means(means, norms):
    fig, axes = plt.subplots(2, 2, figsize=(20, 15))
    fig.subplots_adjust(left=.165, right=.985, bottom=.075, top=.92, hspace=.45, wspace=.56)
    for ax, metric in zip(axes.flat, ('ARI', 'NMI', 'PAS', 'CHAOS')):
        data = np.array([[means.loc[(means.dataset == case) & (means.arm == arm), metric].iloc[0]
                          for case in CASES] for arm in ARMS])
        if metric == 'ARI':
            for j, case in enumerate(CASES):
                ax.imshow(np.ma.masked_invalid(data[:, j:j + 1]), cmap=CMAP, norm=norms[case], aspect='auto',
                          extent=(j - .5, j + .5, 9.5, -.5))
            ax.set_xlim(-.5, 2.5)
        else:
            color = plt.get_cmap('viridis' if metric == 'NMI' else 'viridis_r').copy()
            color.set_bad('#E2E8F0')
            finite = data[np.isfinite(data)]
            norm = Normalize(0, 1 if metric == 'NMI' else float(finite.max()))
            ax.imshow(np.ma.masked_invalid(data), cmap=color, norm=norm, aspect='auto')
        for i in range(len(ARMS)):
            for j, case in enumerate(CASES):
                value = data[i, j]
                rgba = CMAP(norms[case](value)) if metric == 'ARI' else color(norm(value))
                ax.text(j, i, '—' if not np.isfinite(value) else f'{value:.3f}', ha='center', va='center',
                        fontsize=11, color='#526174' if not np.isfinite(value) else text_color(rgba))
        ax.set_xticks(range(3), TITLES, rotation=13, ha='center', fontsize=10)
        ax.set_yticks(range(len(ARMS)), LABELS, fontsize=10)
        ax.set_title(metric + (' (higher better)' if metric in ('ARI', 'NMI') else ' (lower better)'), fontsize=16)
        rules(ax)
    ari = axes[0, 0]
    pos = ari.get_position()
    for j, case in enumerate(CASES):
        cb_ax = fig.add_axes([pos.x0 + pos.width * (j / 3 + .015), pos.y0 - .083, pos.width * .30, .011])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norms[case], cmap=CMAP), cax=cb_ax, orientation='horizontal')
        cb.set_ticks([norms[case].vmin, norms[case].vmax])
        cb.ax.set_xticklabels([f'{norms[case].vmin:.3f}', f'{norms[case].vmax:.3f}'])
        cb.ax.get_xticklabels()[0].set_ha('left')
        cb.ax.get_xticklabels()[-1].set_ha('right')
        cb.ax.tick_params(labelsize=8)
    fig.text(pos.x0 + pos.width / 2, pos.y0 - .115, 'ARI: independent dataset scales; yellow to red', ha='center', fontsize=10, color='#526174')
    fig.suptitle('Equal-timepoint means | one training seed; no seed SD', fontsize=21, y=.977)
    save(fig, 'primary_metric_means')


def main():
    source_manifest = json.loads((SOURCE / 'artifact_sha256.json').read_text())
    assert all(sha(SOURCE / name) == digest for name, digest in source_manifest.items())
    if DESTINATION.exists():
        for name, digest in source_manifest.items():
            if name.endswith('.csv') or name.endswith('_spatial.png'):
                assert sha(DESTINATION / name) == digest, 'Copied scientific results were modified'
    else:
        shutil.copytree(SOURCE, DESTINATION)
    primary = pd.read_csv(SOURCE / 'primary_per_timepoint.csv', float_precision='round_trip')
    means = pd.read_csv(SOURCE / 'primary_equal_timepoint_means.csv', float_precision='round_trip')
    matrices, norms = {}, {}
    for case in CASES:
        stage_ids = sorted(primary.loc[primary.dataset == case, 'timepoint'].unique(), key=lambda x: int(str(x).removeprefix('t')))
        array = np.full((len(ARMS), len(stage_ids) + 1), np.nan)
        for i, arm in enumerate(ARMS):
            for j, tp in enumerate(stage_ids):
                r = primary.loc[(primary.dataset == case) & (primary.arm == arm) & (primary.timepoint == tp), 'ARI']
                assert len(r) <= 1
                if len(r):
                    array[i, j] = r.iloc[0]
            array[i, -1] = means.loc[(means.dataset == case) & (means.arm == arm), 'ARI'].iloc[0]
        finite = array[:, :-1][np.isfinite(array[:, :-1])]
        norms[case] = Normalize(float(finite.min()), float(finite.max()))
        matrices[case] = array
        assert np.nanmin(array[:, -1]) >= norms[case].vmin
        assert np.nanmax(array[:, -1]) <= norms[case].vmax
    draw_main(primary, means, matrices, norms)
    draw_metric_means(means, norms)
    proof = {'colormap': 'YlOrRd', 'higher_ARI_redder': True, 'normalization': 'per_dataset_observed_timepoint_min_max',
             'mean_separate_column_and_axis': True, 'mean_uses_its_dataset_timepoint_scale': True,
             'category_groups': [{'label': label, 'arms': list(ARMS[start:end]), 'color': color}
                                 for start, end, _, color, label in GROUPS],
             'scales': {case: {'minimum': norms[case].vmin, 'maximum': norms[case].vmax} for case in CASES},
             'CSV_values_unchanged': True, 'missing_values_gray': True, 'no_model_or_clustering_fit': True,
             'no_OT_solve': True, 'visual_QA_pending': True}
    (DESTINATION / 'ARI_display_verification.json').write_text(json.dumps(proof, indent=2) + '\n')
    html = (DESTINATION / 'report.html').read_text()
    if 'The left grouping bar identifies experimental categories.' not in html:
        html = html.replace('<h2>Primary20D performance</h2>', '<h2>Primary20D performance</h2><p>ARI colors run from yellow (lower) to red (higher), with a separate observed range for each dataset. Mean ARI is shown in its own separated column and uses the same scale as its dataset. Compare printed values across datasets. The left grouping bar identifies experimental categories.</p>')
    (DESTINATION / 'report.html').write_text(html)
    for name in SOURCE.glob('*.csv'):
        assert sha(name) == sha(DESTINATION / name.name)
    source_hashes = json.loads((DESTINATION / 'source_sha256.json').read_text())
    source_hashes[str(Path(__file__).resolve())] = sha(__file__)
    (DESTINATION / 'source_sha256.json').write_text(json.dumps(source_hashes, indent=2) + '\n')
    verification = json.loads((DESTINATION / 'verification.json').read_text())
    verification.update(visual_inspection_performed=False, ARI_display_updated=True, unchanged_numeric_result_tables=True,
                        dataset_specific_ARI_scales=True, grouped_experiments=True)
    verification.pop('visual_QA_record', None)
    (DESTINATION / 'verification.json').write_text(json.dumps(verification, indent=2) + '\n')
    hashes = {str(p.relative_to(DESTINATION)): sha(p) for p in sorted(DESTINATION.iterdir()) if p.is_file() and p.name != 'artifact_sha256.json'}
    (DESTINATION / 'artifact_sha256.json').write_text(json.dumps(hashes, indent=2) + '\n')
    assert all(sha(SOURCE / name) == digest for name, digest in source_manifest.items())
    print(json.dumps({'destination': str(DESTINATION), 'scales': proof['scales'], 'numeric_tables_unchanged': True}))


if __name__ == '__main__':
    main()
