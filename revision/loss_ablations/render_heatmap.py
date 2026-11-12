"""Render verified coefficient-only loss ablations with per-dataset ARI scales."""

import argparse
from pathlib import Path
import hashlib
import json
import os

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input-dir', required=True, type=Path,
                    help='Directory with matched_loss_dataset_metrics.csv and matched_loss_timepoint_metrics.csv')
parser.add_argument('--output-dir', required=True, type=Path)
args = parser.parse_args()
ROOT, SOURCE = args.output_dir, args.input_dir
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "work" / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd


VARIANTS = [
    ("full", "Full SpaDOT"),
    ("without_ot", "Without OT penalty"),
    ("without_encoder_alignment", "Without alignment"),
    ("without_clustering_penalty", "Without clustering penalty"),
    ("dual_without_penalties", "Without auxiliary penalties"),
    ("dual_alignment_only", "With alignment only"),
]
DATASETS = [
    ("ChickenHeart", "Chicken Heart", 4),
    ("MouseEmbryo", "Mouse Embryo", 3),
    ("AxolotlTelencephalon", "Axolotl", 3),
]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text_color(rgba):
    rgb = np.asarray(rgba[:3])
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    luminance = float(linear @ np.array([0.2126, 0.7152, 0.0722]))
    return "#ffffff" if 1.05 / (luminance + 0.05) > (luminance + 0.05) / 0.05 else "#151515"


def main():
    ROOT.mkdir(exist_ok=True, parents=True)
    summary_path = SOURCE / "matched_loss_dataset_metrics.csv"
    timepoint_path = SOURCE / "matched_loss_timepoint_metrics.csv"
    summary = pd.read_csv(summary_path)
    timepoints = pd.read_csv(timepoint_path)
    assert len(summary) == 18
    assert not summary.duplicated(["dataset", "variant"]).any()
    assert set(summary["variant"]) == {key for key, _ in VARIANTS}
    assert set(summary["dataset"]) == {key for key, _, _ in DATASETS}
    assert (summary["executed_epoch"] == 99).all()
    assert (summary["clustering_status"] == "valid").all()
    assert not timepoints.duplicated(["dataset", "variant", "timepoint"]).any()
    assert (timepoints["executed_epoch"] == 99).all()
    assert np.isfinite(summary["mean_ARI"]).all()

    values = np.empty((len(VARIANTS), len(DATASETS)))
    errors = []
    for column, (dataset, _, expected_timepoints) in enumerate(DATASETS):
        for row, (variant, _) in enumerate(VARIANTS):
            record = summary.loc[(summary.dataset == dataset) & (summary.variant == variant)].iloc[0]
            per_timepoint = timepoints.loc[(timepoints.dataset == dataset) & (timepoints.variant == variant)]
            assert len(per_timepoint) == expected_timepoints
            actual = float(per_timepoint.ARI.mean())
            value = float(record.mean_ARI)
            errors.append(abs(actual - value))
            assert abs(actual - value) < 1e-12
            values[row, column] = value

    table = pd.DataFrame(
        values,
        index=pd.Index([label for _, label in VARIANTS], name="Configuration"),
        columns=[label for _, label, _ in DATASETS],
    )
    table.to_csv(ROOT / "loss-ablation-mean-ari.csv", float_format="%.16g")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "svg.fonttype": "none",
        "axes.titleweight": "medium",
        "savefig.facecolor": "white",
    })
    fig = plt.figure(figsize=(9.6, 5.7), facecolor="white")
    fig.text(0.035, 0.95, "Loss ablation", fontsize=17, weight="medium", va="top")
    fig.text(0.035, 0.895, "Mean ARI across timepoints", fontsize=11.5, color="#4d4d4d", va="top")

    greens = plt.get_cmap("Greens")
    cmap = LinearSegmentedColormap.from_list("ablation_greens", greens(np.linspace(0.09, 0.92, 256)))
    left, panel_width, gap = 0.385, 0.177, 0.025
    axes = []
    colorbars = []
    bounds = []
    for column, (_, label, _) in enumerate(DATASETS):
        x = left + column * (panel_width + gap)
        ax = fig.add_axes([x, 0.255, panel_width, 0.555])
        axes.append(ax)
        col = values[:, column]
        low, high = float(col.min()), float(col.max())
        norm = Normalize(vmin=low, vmax=high)
        image = ax.imshow(col[:, None], cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        ax.set_title(label, fontsize=12, pad=11)
        ax.set_xticks([])
        ax.set_yticks(np.arange(len(VARIANTS)))
        if column == 0:
            ax.set_yticklabels([name for _, name in VARIANTS], fontsize=11.5)
            ax.tick_params(axis="y", length=0, pad=13)
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis="y", length=0)
        ax.set_yticks(np.arange(-0.5, len(VARIANTS), 1), minor=True)
        ax.grid(which="minor", axis="y", color="white", linewidth=2.0)
        ax.tick_params(which="minor", left=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for row, value in enumerate(col):
            ax.text(0, row, f"{value:.4f}", ha="center", va="center", fontsize=13,
                    color=text_color(cmap(norm(value))))

        bar_ax = fig.add_axes([x + 0.004, 0.175, panel_width - 0.008, 0.021])
        bar = fig.colorbar(image, cax=bar_ax, orientation="horizontal", ticks=[low, high])
        colorbars.append(bar)
        bar.ax.xaxis.set_major_formatter(FormatStrFormatter("%.4f"))
        bar.ax.tick_params(labelsize=10, length=2.5, pad=4)
        tick_labels = bar.ax.get_xticklabels()
        tick_labels[0].set_ha("left")
        tick_labels[-1].set_ha("right")
        bar.outline.set_visible(False)
        bounds.append({"dataset": DATASETS[column][0], "vmin": low, "vmax": high})

    fig.text(left, 0.092, "Independent color scale for each dataset", fontsize=10, color="#4d4d4d")
    fig.text(0.035, 0.035, "Equal timepoint weights  |  Executed epoch 99  |  Seed 2026", fontsize=10, color="#4d4d4d")
    fig.canvas.draw()

    assert all(np.array_equal(np.asarray(ax.images[0].get_array())[:, 0], values[:, j]) for j, ax in enumerate(axes))
    for j, ax in enumerate(axes):
        assert ax.images[0].norm.vmin == values[:, j].min()
        assert ax.images[0].norm.vmax == values[:, j].max()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    for previous, following in zip(colorbars, colorbars[1:]):
        last_label = previous.ax.get_xticklabels()[-1].get_window_extent(renderer)
        first_label = following.ax.get_xticklabels()[0].get_window_extent(renderer)
        assert last_label.x1 + 5 < first_label.x0
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            bbox = artist.get_window_extent(renderer)
            assert bbox.x0 >= -1 and bbox.y0 >= -1 and bbox.x1 <= canvas.x1 + 1 and bbox.y1 <= canvas.y1 + 1, artist.get_text()

    png = ROOT / "loss-ablation-ari-independent-scales.png"
    svg = ROOT / "loss-ablation-ari-independent-scales.svg"
    fig.savefig(png, dpi=240)
    fig.savefig(svg)
    plt.close(fig)
    pd.DataFrame(bounds).to_csv(ROOT / "color-scale-bounds.csv", index=False)
    verification = {
        "all_passed": True,
        "metric": "ARI",
        "aggregation": "Equal-weight mean over timepoints",
        "values": int(values.size),
        "max_error_vs_timepoint_mean": max(errors),
        "source_files": {str(p): sha256(p) for p in [summary_path, timepoint_path]},
        "independent_scale": "Per-dataset min/max across all six configurations",
        "scale_bounds": bounds,
        "single_seed": 2026,
        "executed_epoch": 99,
        "chicken_without_alignment_version": "numerical_recovery_20261102",
        "output_files": {str(p): sha256(p) for p in [png, svg, ROOT / "loss-ablation-mean-ari.csv"]},
    }
    (ROOT / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    (ROOT / "README.md").write_text(
        "# Loss ablation ARI heatmap\n\n"
        "Six current-model coefficient ablations across three datasets. Each dataset has its own "
        "color scale spanning its six mean ARI values. Cell labels are actual ARIs rounded to four "
        "decimal places. Timepoints are weighted equally. Each result is from executed epoch 99, "
        "seed 2026. Independent color scales support within-dataset comparisons.\n\n"
        "The Chicken Heart alignment-off result uses the independently verified numerical OT "
        "recovery from 2026-11-02. All other entries retain their verified matched loss-ablation "
        "results. Historical architecture results are in separate reports.\n\n"
        "The plotting script validates all 18 values against the timepoint-level records. "
        "verification.json records source hashes, scale limits and output hashes.\n"
    )
    print(json.dumps({"all_passed": True, "max_mean_error": max(errors), "outputs": [str(png), str(svg)], "scales": bounds}))


if __name__ == "__main__":
    main()
