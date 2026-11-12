"""Score annotation and predicted-domain OT at the S41 threshold grid.

Requires NumPy. Pass the pairs directory written by
revision/transition_analysis/compute.py and a new
output directory. Predicted clusters inherit their majority annotation for
reference-edge evaluation. Dataset means weight each transition equally.
The overall curves weight the three datasets equally.
"""
import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np

# Use the reference graph shared by the transition and single-slice benchmarks.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "revision"))
from transition_analysis.reference_edges import edges

STAGES = {
    "ChickenHeart": ("Day 4", "Day 7", "Day 10", "Day 14"),
    "MouseEmbryo": ("E12.5", "E13.5", "E15.5"),
    "AxolotlTelencephalon": ("St.44", "St.54", "St.57"),
}
THRESHOLDS = [i / 100 for i in range(0, 51, 5)]
AGGREGATIONS = ("annotation", "predicted")
METRICS = ("naive_accuracy", "reference_weight_retention", "f1_using_association_weights",
           "retained_weight_fraction", "reference_edge_precision", "reference_edge_recall", "edge_level_f1")


def export(path, rows):
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_table(path):
    with np.load(path, allow_pickle=False) as table:
        mass = np.asarray(table["mass"], dtype=np.float64)
        rows, cols, row_names, col_names = (table[key].astype(str).tolist()
                                             for key in ("rows", "cols", "row_annotation", "col_annotation"))
    if (mass.shape != (len(rows), len(cols)) or len(rows) != len(set(rows)) or len(cols) != len(set(cols))
            or not np.isfinite(mass).all() or (mass < 0).any()
            or not (mass.sum(axis=1) > 0).all() or not (mass.sum(axis=0) > 0).all()):
        raise ValueError(f"Invalid domain transport table: {path}")
    return mass, rows, cols, row_names, col_names


def score_threshold(scores, reference, threshold):
    keep = (scores >= threshold) & (scores > 0)
    tp = float(scores[keep & reference].sum())
    fp = float(scores[keep & ~reference].sum())
    fn = float(scores[~keep & reference].sum())
    retained = tp + fp
    total_reference = float(scores[reference].sum())
    if total_reference <= 0:
        raise ValueError("No positive mass on reference correspondences.")
    btp, bfp, bfn = (int(mask.sum()) for mask in (keep & reference, keep & ~reference, ~keep & reference))
    return {
        "naive_accuracy": tp / retained if retained else math.nan,
        "reference_weight_retention": tp / total_reference,
        "f1_using_association_weights": 2 * tp / (retained + total_reference),
        "retained_weight_fraction": retained / float(scores.sum()),
        "reference_edge_precision": btp / (btp + bfp) if btp + bfp else math.nan,
        "reference_edge_recall": btp / (btp + bfn),
        "edge_level_f1": 2 * btp / (2 * btp + bfp + bfn),
        "retained_links": int(keep.sum()),
        "tp_weight": tp, "fp_weight": fp, "fn_weight": fn,
        "total_retained_weight": retained, "original_reference_weight": total_reference,
    }


def summarize(rows):
    return {
        "retained_links": sum(row["retained_links"] for row in rows),
        **{metric: float(np.mean([row[metric] for row in rows])) for metric in METRICS},
    }


def per_pair(pairs):
    records, all_edges = [], []
    for dataset, stages in STAGES.items():
        for i, (source_stage, target_stage) in enumerate(zip(stages, stages[1:])):
            pair = f"t{i}_t{i+1}"
            folder = pairs / dataset / pair
            annotation = load_table(folder / "annotation.npz")
            graph = edges(dataset, i, annotation[1])
            for aggregation in AGGREGATIONS:
                mass, rows, cols, row_names, col_names = (
                    annotation if aggregation == "annotation" else load_table(folder / "predicted.npz"))
                if not (set(row_names) <= set(annotation[1]) and set(col_names) <= set(annotation[2])):
                    raise ValueError(f"Unexpected predicted majority labels for {dataset} {pair}.")
                forward = mass / mass.sum(axis=1, keepdims=True)
                reverse = mass / mass.sum(axis=0, keepdims=True)
                score = np.minimum(forward, reverse)
                reference = np.array([[(s, t) in graph for t in col_names] for s in row_names], dtype=bool)
                meta = dict(aggregation=aggregation, dataset=dataset, pair=pair)
                all_edges.extend(
                    dict(**meta, source=source, target=target, source_annotation=row_names[j],
                         target_annotation=col_names[k], mass=float(mass[j, k]), forward=float(forward[j, k]),
                         reverse=float(reverse[j, k]), association_score=float(score[j, k]),
                         reference_edge=bool(reference[j, k]))
                    for j, source in enumerate(rows) for k, target in enumerate(cols))
                records.extend(
                    dict(**meta, source_stage=source_stage, target_stage=target_stage,
                         threshold=threshold, **score_threshold(score, reference, threshold))
                    for threshold in THRESHOLDS)
    return records, all_edges


def summarize_thresholds(records):
    dataset_means, overall, best = [], [], []
    for aggregation in AGGREGATIONS:
        for threshold in THRESHOLDS:
            within_threshold = []
            for dataset in STAGES:
                group = [row for row in records if (row["aggregation"], row["dataset"], row["threshold"])
                         == (aggregation, dataset, threshold)]
                within_threshold.append(dict(
                    aggregation=aggregation, dataset=dataset, threshold=threshold, pair_count=len(group),
                    undefined_accuracy_pairs=sum(not math.isfinite(row["naive_accuracy"]) for row in group),
                    **summarize(group)))
            dataset_means.extend(within_threshold)
            overall.append(dict(
                aggregation=aggregation, threshold=threshold, dataset_count=len(within_threshold),
                undefined_accuracy_pairs=sum(row["undefined_accuracy_pairs"] for row in within_threshold),
                **summarize(within_threshold)))

        scopes = {"Equal-dataset mean": [row for row in overall if row["aggregation"] == aggregation]}
        for dataset in STAGES:
            scopes[dataset] = [row for row in dataset_means if row["aggregation"] == aggregation
                               and row["dataset"] == dataset]
        for scope, rows in scopes.items():
            eligible = [row for row in rows if row["threshold"] > 0]
            for criterion in ("f1_using_association_weights", "edge_level_f1"):
                maximum = max(row[criterion] for row in eligible)
                tied = [row["threshold"] for row in eligible if abs(row[criterion] - maximum) < 1e-12]
                best.append(dict(aggregation=aggregation, scope=scope, criterion=criterion,
                                 best_tested_threshold=min(tied), best_score=maximum,
                                 tied_thresholds=";".join(f"{threshold:.2f}" for threshold in tied)))
    return dataset_means, overall, best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records, all_edges = per_pair(args.pairs)
    dataset, overall, best = summarize_thresholds(records)
    args.output.mkdir(parents=True, exist_ok=False)
    for name, rows in (("per_pair_results", records), ("association_edges", all_edges),
                       ("dataset_mean_results", dataset), ("equal_dataset_mean_results", overall),
                       ("best_tested_thresholds", best)):
        export(args.output / f"{name}.csv", rows)


if __name__ == "__main__":
    main()
