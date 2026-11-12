"""Render matched transition F1 figures from the audited WOT threshold rerun."""

import argparse
import csv
import math
from pathlib import Path


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PANELS = {
    "Equal-dataset mean": "Mean of Three Datasets",
    "ChickenHeart": "Chicken Cardiogenesis",
    "MouseEmbryo": "Mouse Organogenesis",
    "AxolotlTelencephalon": "Axolotl Telencephalon\nNeurogenesis",
}
SETTINGS = dict(epsilon=0.01, epsilon0=1, lambda1=1, lambda2=50, growth_iters=3)
SETTING_TEXT = "Waddington-OT with ε = 0.01, ε₀ = 1, λ₁ = 1, λ₂ = 50, and three growth iterations"
F1 = "f1_using_association_weights"
THRESHOLDS = [i / 100 for i in range(5, 51, 5)]


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def plot_aggregation(aggregation, prefix, filename, overall, dataset, best, output):
    metrics = (
        ("naive_accuracy", "Naive accuracy", "#2e7189", ":", 2.0, 4.5),
        ("reference_weight_retention", f"{prefix} retention", "#b76643", ":", 2.0, 4.5),
        (F1, "F1", "#333333", "-", 3.8, 6.5),
    )
    fig, axes = plt.subplots(2, 2, figsize=(17.5, 13.5), sharex=True, sharey=True)
    fig.subplots_adjust(left=.105, right=.975, bottom=.165, top=.875, wspace=.25, hspace=.85)
    curves = []
    for ax, (scope, panel_title) in zip(axes.ravel(), PANELS.items()):
        records = overall if scope == "Equal-dataset mean" else [
            row for row in dataset if row["dataset"] == scope
        ]
        data = sorted(
            (row for row in records if row["aggregation"] == aggregation and float(row["threshold"]) > 0),
            key=lambda row: float(row["threshold"]),
        )
        x = [float(row["threshold"]) for row in data]
        assert x == THRESHOLDS, (aggregation, scope, x)
        for metric, label, color, style, width, marker_size in metrics:
            y = np.array([float(row[metric]) for row in data])
            line, = ax.plot(
                x, y, color=color, linestyle=style, linewidth=width, marker="o",
                markersize=marker_size, dash_capstyle="round", label=label,
                zorder=4 if metric == F1 else 2,
            )
            assert np.array_equal(line.get_ydata(), y, equal_nan=True)
            assert line.get_linestyle() == style and line.get_linewidth() == width
            curves.append(dict(
                aggregation=aggregation, panel=scope, metric=metric, thresholds=x,
                values=[float(value) if np.isfinite(value) else None for value in y],
                linestyle=style, linewidth=width,
            ))
        selected, = [
            row for row in best if row["aggregation"] == aggregation
            and row["scope"] == scope and row["criterion"] == F1
        ]
        threshold = float(selected["best_tested_threshold"])
        peak = float(selected["best_score"])
        assert math.isclose(peak, float(data[x.index(threshold)][F1]), abs_tol=1e-12)
        assert math.isclose(peak, max(float(row[F1]) for row in data), abs_tol=1e-12)
        ax.axvline(threshold, color="#888888", linestyle="--", linewidth=1.2, zorder=1)
        ax.plot(threshold, peak, "*", markersize=15, color="#333333", zorder=5)
        ax.set_title(panel_title, fontsize=26, pad=16)
        ax.set_xlabel("Threshold", fontsize=24, labelpad=8)
        ax.set_ylabel("F1-score", fontsize=24, labelpad=10)
        ax.set_xticks(x, [f"{value:.2f}" for value in x])
        ax.tick_params(axis="x", labelbottom=True, labelsize=20, rotation=45)
        ax.tick_params(axis="y", labelleft=True, labelsize=20)
        ax.set_ylim(0, 1.04)
        ax.spines[["top", "right"]].set_visible(False)

    title = f"{prefix} Transition F1-score on Visualization Thresholds"
    fig.suptitle(title, fontsize=32, y=.985)
    legend = fig.legend(
        *axes[0, 0].get_legend_handles_labels(), loc="lower center",
        bbox_to_anchor=(.5, .02), ncol=3, frameon=False, fontsize=18,
    )
    legend.get_texts()[-1].set_fontweight("bold")
    pdf = output / "output/pdf" / filename
    png = output / f"{aggregation}_transition_f1.png"
    fig.savefig(pdf, bbox_inches="tight", metadata={
        "Title": title,
        "Subject": f"{SETTING_TEXT}. {aggregation} aggregation, thresholds 0.05 to 0.50.",
        "Creator": "Matplotlib",
    })
    fig.savefig(png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return curves, [pdf, png]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rerun", type=Path, required=True, help="Output of score_thresholds.py")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    names = ("equal_dataset_mean_results.csv", "dataset_mean_results.csv", "best_tested_thresholds.csv")
    overall, dataset, best = [read_csv(args.rerun / name) for name in names]
    (args.output / "output/pdf").mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({"font.size": 20, "pdf.fonttype": 42, "ps.fonttype": 42})
    for aggregation, prefix in (("annotation", "Annotation-based"), ("predicted", "Predicted-domain")):
        plot_aggregation(aggregation, prefix, f"{aggregation}_transition_f1.pdf",
                         overall, dataset, best, args.output)


if __name__ == "__main__":
    main()
