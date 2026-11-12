"""Paired Moran I differences with equal-timepoint means and timepoint SD."""
import hashlib
import json
from pathlib import Path
import argparse
from plot_helpers import DATASETS, MEAN, load_results, paired_data, sha, style_axis
from zero_crossing import draw_crossings
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, MultipleLocator
import numpy as np
import pandas as pd


def array_sha(values):
    values = np.ascontiguousarray(values)
    digest = hashlib.sha256()
    digest.update(str(values.dtype).encode())
    digest.update(str(values.shape).encode())
    digest.update(values.tobytes())
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True, help='Output from compute_moran.py')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    raw, summary = load_results(args.input_dir)
    source_hashes = {name: sha(args.input_dir / name) for name in ('per_timepoint.csv', 'summary.csv')}
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 14,
                         'axes.labelsize': 16, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, (axes, strips) = plt.subplots(2, 3, figsize=(15, 5.5), sharex='col', sharey=False,
                                      gridspec_kw={'height_ratios': (4.1, .86)})
    fig.subplots_adjust(left=.102, right=.977, top=.815, bottom=.27, wspace=.27, hspace=.12)
    box = strips[0].get_position()
    fig.text(.089, (box.y0 + box.y1) / 2,
             'First Zero\nCrossing',
             ha='right', va='center', fontsize=11.2, color='#334155')
    fig.supxlabel('Standardized spatial distance', x=.53, y=.115, fontsize=16)
    records, plotted_rows, aligned_rows = {}, [], []
    error_color = '#7c899a'
    for ax, strip, (case, title, labels) in zip(axes, strips, DATASETS):
        original, grid, stages, mean, mean_error = paired_data(raw, summary, case, labels)
        matrix = np.stack([np.interp(grid, x, delta) for _, x, delta in stages])
        assert np.array_equal(matrix.mean(0), mean)
        sd = matrix.std(0, ddof=1)
        assert np.isfinite(sd).all() and (sd >= 0).all()
        # Ten sparse display positions avoid drawing 100 overlapping error bars.
        selected = np.round(np.linspace(0, len(grid) - 1, 10)).astype(int)
        assert len(np.unique(selected)) == 10
        style_axis(ax)
        titles = {
            'MouseEmbryo': 'Mouse Embryo\nOrganogenesis',
            'ChickenHeart': 'Chicken Cardiogenesis',
            'AxolotlTelencephalon': 'Axolotl Telencephalon\nNeurogenesis',
        }
        ax.set_title(titles[case], pad=15, fontsize=18, fontweight='bold')
        line, = ax.plot(grid, mean, color=MEAN, lw=2.7, zorder=4)
        errors = ax.errorbar(grid[selected], mean[selected], yerr=sd[selected], fmt='o',
                            color=MEAN, markerfacecolor='white', markersize=4,
                            ecolor=error_color, elinewidth=1.05, capsize=3,
                            capthick=1.05, zorder=5)
        assert np.array_equal(line.get_ydata(), mean)
        segments = errors.lines[2][0].get_segments()
        for segment, x, y, s in zip(segments, grid[selected], mean[selected], sd[selected], strict=True):
            assert np.array_equal(segment[:, 0], [x, x])
            assert np.array_equal(segment[:, 1], [y - s, y + s])
        # All per-stage arrays contribute to the summary and none is drawn.
        lower, upper = mean - sd, mean + sd
        lo, hi = min(0., float(lower.min())), max(0., float(upper.max()))
        padding = .075 * (hi - lo)
        ax.set_ylim(lo - padding, hi + padding)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))
        ax.xaxis.set_major_locator(MultipleLocator(.5))
        ax.set_xlim(grid[0] - .03, grid[-1] + .03)
        crossings = draw_crossings(strip, original, grid)
        for label in strip.texts:
            label.set_fontsize(13)
        records[case] = dict(n_timepoints=len(stages), timepoints=list(labels),
                            y_limits=list(ax.get_ylim()), display_grid_indices=selected.tolist(),
                            common_grid_points=len(grid), paired_sd_ddof=1,
                            mean_sha256=array_sha(mean), SD_sha256=array_sha(sd),
                            mean_source_max_error=mean_error,
                            original_mean_Moran_zero_crossing=crossings)
        for index in selected:
            plotted_rows.append(dict(dataset=case, distance=float(grid[index]),
                                     mean=float(mean[index]), SD=float(sd[index]),
                                     n_timepoints=len(stages)))
        aligned_rows.extend(dict(dataset=case, timepoint=label, distance=float(x), delta_I=float(y))
                            for (label, _, _), values in zip(stages, matrix, strict=True)
                            for x, y in zip(grid, values, strict=True))
    axes[0].set_ylabel('Paired difference $\\Delta I$\n(SVGP − GAT)', labelpad=11)
    assert all(not a.get_shared_y_axes().joined(a, b) for i, a in enumerate(axes) for b in axes[i + 1:])
    fig.legend(handles=[Line2D([], [], color=MEAN, lw=2.7, marker='o', markerfacecolor='white',
                               markersize=4, label='Mean ±1 SD')],
               loc='lower center', bbox_to_anchor=(.532, .012), frameon=False, fontsize=14.5)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for artist in (*fig.texts, *fig.legends, *(ax.title for ax in axes), axes[0].yaxis.label):
        box = artist.get_window_extent(renderer)
        assert 0 <= box.x0 <= box.x1 <= fig.bbox.width, (getattr(artist, 'get_text', lambda: type(artist).__name__)(), box.bounds)
        assert 0 <= box.y0 <= box.y1 <= fig.bbox.height, (getattr(artist, 'get_text', lambda: type(artist).__name__)(), box.bounds)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=False)
    stem = 'Moran_I_delta_mean_SD'
    for extension, dpi in (('pdf', 300), ('png', 220)):
        fig.savefig(out / f'{stem}.{extension}', dpi=dpi, facecolor='white')
    plt.close(fig)
    pd.DataFrame(plotted_rows).to_csv(out / 'displayed_mean_SD.csv', index=False, float_format='%.17g')
    pd.DataFrame(aligned_rows).to_csv(out / 'paired_timepoint_values.csv', index=False, float_format='%.17g')
    assert all(sha(args.input_dir / name) == digest for name, digest in source_hashes.items())
    proof = dict(all_passed=True, source_sha256=source_hashes, datasets=records,
                 difference='SVGP minus GAT within each matched native distance bin',
                 pooling='Linear interpolation on shared observed distance support; equal timepoint weights',
                 errorbar='Sample SD across paired timepoint differences, ddof=1',
                 individual_timepoint_curves_drawn=False,
                 errorbar_positions='Ten equally spaced indices on the established 100-point common grid',
                 errorbar_segments_exact=True,
                 output_sha256={p.name: sha(p) for p in sorted(out.iterdir()) if p.name != 'verification.json'})
    (out / 'verification.json').write_text(json.dumps(proof, indent=2) + '\n')


if __name__ == '__main__':
    main()
