"""Render PAS or CHAOS from audited saved tables using the ARI figure layout."""
from pathlib import Path
import hashlib
import json
import os
import tempfile
import argparse

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'spadot_grouped_ari_mpl'))
import s39_style as original
np, pd, plt = original.np, original.pd, original.plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch

parser = argparse.ArgumentParser(
    description=__doc__,
    epilog="Input contains primary_per_timepoint.csv and primary_equal_timepoint_means.csv from the audited parent/child ablation analysis."
)
parser.add_argument('--input-dir', required=True, type=Path)
parser.add_argument('--output-dir', required=True, type=Path)
parser.add_argument('metric', choices=('PAS', 'CHAOS'))
args = parser.parse_args()
INPUT, OUT = args.input_dir, args.output_dir
METRIC = args.metric
PDF = OUT / (METRIC + '_by_timepoint_standalone.pdf')
PNG = OUT / (METRIC + '_by_timepoint_standalone.png')
SVG = OUT / (METRIC + '_by_timepoint_standalone.svg')
CASES = ('MouseEmbryo', 'ChickenHeart', 'AxolotlTelencephalon')
NAMES = {'MouseEmbryo': 'Mouse Embryo\nOrganogenesis', 'ChickenHeart': 'Chicken Heart\nCardiogenesis',
         'AxolotlTelencephalon': 'Axolotl Telencephalon\nNeurogenesis'}
STAGES = {'MouseEmbryo': ('E12.5', 'E13.5', 'E15.5'),
          'ChickenHeart': ('Day 4', 'Day 7', 'Day 10', 'Day 14'),
          'AxolotlTelencephalon': ('St.44', 'St.54', 'St.57')}
LABELS = ('SpaDOT', 'GAT20 depth 1', 'GAT20 depth 2', 'GAT20 depth 3', 'GAT20 depth 4',
          'GAT20 depth 6', 'GAT20 D3 matched', 'GAT20 D3 expanded', 'MLP10+GAT10')
GROUPS = ((0, 1, 'Reference'), (1, 6, 'GAT-only\ndepth'),
          (6, 7, 'GAT-only\ncapacity control'), (7, 8, 'GAT-only\nneighborhood\ncontrol'),
          (8, 9, 'GP-branch\nreplacement'))
MEAN_CMAP = plt.get_cmap('Blues_r').copy()
MEAN_CMAP.set_bad('#E2E8F0')
TIMEPOINT_CMAP = plt.get_cmap('Reds_r').copy()
TIMEPOINT_CMAP.set_bad('#E2E8F0')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def render():
    primary_path, mean_path = INPUT / 'primary_per_timepoint.csv', INPUT / 'primary_equal_timepoint_means.csv'
    inputs = {str(p): sha(p) for p in (primary_path, mean_path, Path(__file__))}
    primary = pd.read_csv(primary_path, float_precision='round_trip')
    means = pd.read_csv(mean_path, float_precision='round_trip')
    matrices, norms, mean_norms, timepoint_limits = {}, {}, {}, {}
    table_rows = []
    for case in CASES:
        times = sorted(primary.loc[primary.dataset == case, 'timepoint'].unique(),
                       key=lambda x: int(str(x).removeprefix('t')))
        data = np.full((len(original.ARMS), len(times) + 1), np.nan)
        for i, arm in enumerate(original.ARMS):
            for j, timepoint in enumerate(times):
                rows = primary.loc[(primary.dataset == case) & (primary.arm == arm) &
                                   (primary.timepoint == timepoint), METRIC]
                assert len(rows) <= 1
                if len(rows):
                    data[i, j] = rows.iloc[0]
            data[i, -1] = means.loc[(means.dataset == case) & (means.arm == arm), METRIC].iloc[0]
            for j, timepoint in enumerate([*times, 'Mean']):
                table_rows.append({'dataset': case, 'arm': arm, 'timepoint': timepoint, METRIC: data[i, j]})
        finite = data[:, :-1][np.isfinite(data[:, :-1])]
        limits = {'minimum': float(finite.min()), 'maximum': float(finite.max())}
        timepoint_limits[case] = limits
        matrices[case] = data[[i for i, arm in enumerate(original.ARMS) if arm != 'gat_gat_matched']]
        norms[case] = Normalize(limits['minimum'], limits['maximum'])
        mean_values = data[:, -1][np.isfinite(data[:, -1])]
        mean_norms[case] = Normalize(float(mean_values.min()), float(mean_values.max()))

    fig = plt.figure(figsize=(24, 11.6))
    bottom, height = .235, .64
    labels = fig.add_axes([.125, bottom, .155, height])
    labels.set(xlim=(0, 1), ylim=(8.5, -.5))
    labels.axis('off')
    row_text = [labels.text(.98, row, label, ha='right', va='center', fontsize=20)
                for row, label in enumerate(LABELS)]
    # Restore the original row-name right anchor; only move the grouping marks.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    row_left = min(t.get_window_extent(renderer).x0 for t in row_text) / fig.bbox.width
    bracket_gap = 12 / (fig.get_figwidth() * 72)
    cap_x = row_left - bracket_gap
    bracket_x = cap_x - .008
    categories = fig.add_axes([0, bottom, 1, height])
    categories.set(xlim=(0, 1), ylim=(8.5, -.5))
    categories.axis('off')
    for start, end, label in GROUPS:
        top, lower = start - .43, end - .57
        categories.plot([cap_x, bracket_x, bracket_x, cap_x], [top, top, lower, lower],
                        color='#344356', lw=1.8, solid_capstyle='butt', clip_on=False)
        categories.text(bracket_x - .008, (start + end - 1) / 2, label, ha='right', va='center',
                        color='#263548', fontsize=17, linespacing=1.05)

    left, total, gap = .295, .685, .027
    widths = np.array([len(STAGES[case]) + 1 for case in CASES], dtype=float)
    widths = widths / widths.sum() * (total - 2 * gap)
    for case, width in zip(CASES, widths):
        values, norm = matrices[case], norms[case]
        cell = (width - .009) / (len(STAGES[case]) + 1.3)
        stage_width = cell * len(STAGES[case])
        mean_width = cell * 1.3
        stage_ax = fig.add_axes([left, bottom, stage_width, height])
        mean_ax = fig.add_axes([left + stage_width + .009, bottom, mean_width, height])
        for ax, data, ticks, color, scale in ((stage_ax, values[:, :-1], STAGES[case], TIMEPOINT_CMAP, norm),
                                            (mean_ax, values[:, -1:], ('Mean',), MEAN_CMAP, mean_norms[case])):
            ax.imshow(np.ma.masked_invalid(data), cmap=color, norm=scale, aspect='auto')
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    value = data[i, j]
                    ax.text(j, i, '—' if not np.isfinite(value) else f'{value:.4f}',
                            ha='center', va='center', fontsize=18,
                            color='#526174' if not np.isfinite(value) else original.text_color(color(scale(value))))
            ax.set_xticks(range(data.shape[1]), ticks, fontsize=19, rotation=25, ha='right')
            ax.set_yticks([])
            original.rules(ax)
            for spine in ax.spines.values():
                spine.set_color('#526174')
                spine.set_linewidth(.8)
        fig.text(left + width / 2, .945, NAMES[case], ha='center', va='center', fontsize=28, linespacing=1.15)
        for x, bar_width, color, scale, count in ((left, stage_width, TIMEPOINT_CMAP, norm, 3),
                                                (left + stage_width + .009, mean_width, MEAN_CMAP, mean_norms[case], 2)):
            cax = fig.add_axes([x, .128, bar_width, .021])
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=scale, cmap=color), cax=cax, orientation='horizontal')
            positions = np.linspace(scale.vmin, scale.vmax, count)
            cb.set_ticks(positions)
            cb.ax.set_xticklabels([f'{v:.3f}' for v in positions])
            cb.ax.get_xticklabels()[0].set_ha('left')
            cb.ax.get_xticklabels()[-1].set_ha('right')
            cb.ax.tick_params(labelsize=14 if count == 2 else 17)
            cb.set_label(METRIC, fontsize=18, labelpad=5)
        left += width + gap
    fig.legend(handles=[Patch(facecolor='#E2E8F0', edgecolor='#A5AFBD', label='NA')],
               loc='center right', bbox_to_anchor=(.2769, .1385), frameon=False,
               fontsize=18, handlelength=1.3, handleheight=1)
    OUT.mkdir(parents=True, exist_ok=True)
    for path in (PDF, PNG, SVG):
        fig.savefig(path, dpi=150, facecolor='white', metadata={'Title': METRIC + ' spatial continuity comparison',
                    'Subject': 'Lower is better; darker red and blue indicate lower values; seed2026 epoch99'} if path.suffix == '.pdf' else None)
    plt.close(fig)
    values_path = OUT / (METRIC + '_values.csv')
    values = pd.DataFrame(table_rows)
    values = values.loc[values.arm != 'gat_gat_matched'].copy()
    assert len(values) == 117 and values[METRIC].notna().sum() == 103
    values.to_csv(values_path, index=False)
    assert all(sha(p) == value for p, value in inputs.items())
    record = {'metric': METRIC, 'lower_is_better': True, 'darker_red_and_blue_indicate_lower_values': True,
              'source_SHA256': inputs, 'outputs_SHA256': {str(p): sha(p) for p in (PDF, PNG, SVG, values_path)},
              'dataset_order': list(CASES), 'dataset_titles': [NAMES[c].replace('\n', ' ') for c in CASES],
              'group_brackets': [{'rows': [a, b - 1], 'label': label.replace('\n', ' ')} for a, b, label in GROUPS],
              'top_title_and_subtitle_removed': True, 'footer_notes_removed': True,
              'colorbar_labels': [METRIC] * 6, 'gray_legend': 'NA', 'row_labels': list(LABELS),
              'timepoint_colormap': 'Reds_r', 'timepoint_scales': timepoint_limits,
              'mean_colormap': 'Blues_r',
              'mean_scales': {c: {'minimum': mean_norms[c].vmin, 'maximum': mean_norms[c].vmax} for c in CASES},
              'mean_normalization': 'per_dataset_observed_mean_min_max',
              'mean_colorbar_tick_decimals': 3, 'cell_decimals': 4,
              'mean_column_width_relative_to_timepoint': 1.3,
              'row_labels_original_right_alignment_restored': True,
              'row_label_right_anchor_fraction': .125 + .155 * .98,
              'row_label_leftmost_fraction': row_left,
              'bracket_right_cap_fraction': cap_x,
              'bracket_to_row_label_gap_points': bracket_gap * fig.get_figwidth() * 72,
              'font_sizes': {'category': 17, 'row_label': 20, 'cell': 18, 'timepoint': 19, 'dataset': 28, 'colorbar_tick': 17, 'mean_colorbar_tick': 14, 'colorbar_label': 18, 'NA_legend': 18},
              'NA_legend_position': {'side': 'left_of_first_colorbar', 'anchor': [.2769, .1385]},
              'saved_score_values_unchanged': True, 'omitted_display_arm': 'gat_gat_matched', 'visual_QA_pending': True}
    (OUT / (METRIC + '_display_verification.json')).write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'PDF': str(PDF), 'PNG': str(PNG), 'bytes': PDF.stat().st_size}))


if __name__ == '__main__':
    render()
