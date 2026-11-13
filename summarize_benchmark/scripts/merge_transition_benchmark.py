"""Add STIFT + shared OT to the archived naïve transition benchmark for S18.

Supply the three original OT_method_comparison.csv files under
ARCHIVE/DATASET/OT_transition_comparisons. Supply STIFT scores as CSV columns
Dataset, Transition, naive_accuracy, lambda1, lambda2, epsilon, epsilon0,
growth_iters. Use the unthresholded annotation-based shared WOT scores.
"""
import argparse
import csv
from pathlib import Path
from statistics import mean

DATASETS = {
    "Mouse Organogenesis": ("MouseEmbryo_Visium_Organogenesis", ("E12.5", "E13.5", "E15.5")),
    "Chicken Cardiogenesis": ("ChickenHeart_Visium_Cardiogenesis", ("Day 4", "Day 7", "Day 10", "Day 14")),
    "Axolotl Telencephalon Neurogenesis": (
        "AxolotlTelencephalon_Stereoseq_Developmental", ("St.44", "St.54", "St.57")),
}
METHODS = ("SpaDOT", "DeST-OT", "Moscot", "TOAST")
PARAMETERS = dict(lambda1=1, lambda2=50, epsilon=0.01, epsilon0=1, growth_iters=3)


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--stift-scores", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    scores = read(args.stift_scores)
    stift = {(r["Dataset"], r["Transition"]): float(r["naive_accuracy"]) for r in scores}
    if len(stift) != len(scores) or any(float(row[k]) != v for row in scores for k, v in PARAMETERS.items()):
        raise ValueError("Expected unique STIFT pairs with the paper WOT settings.")
    fields = ("SpaDOT (archived)", "STIFT + shared OT", *METHODS[1:])
    result = []
    used_pairs = set()
    for dataset, (folder, stages) in DATASETS.items():
        path = args.archive / folder / "OT_transition_comparisons" / "OT_method_comparison.csv"
        source = read(path)
        archive = {(r["Day"], r["method"]): float(r["naive_accuracy"]) for r in source}
        if len(archive) != len(source) or len(archive) != (len(stages) - 1) * len(METHODS):
            raise ValueError(f"Incomplete or duplicate original results: {path}")
        pair_rows = []
        for i, (a, b) in enumerate(zip(stages, stages[1:])):
            transition = f"{a} → {b}"
            pair = f"{i}->{i+1}"
            used_pairs.add((dataset, transition))
            pair_rows.append({
                "Dataset": dataset, "Transition": transition,
                "SpaDOT (archived)": archive[pair, "SpaDOT"],
                "STIFT + shared OT": stift[dataset, transition],
                **{method: archive[pair, method] for method in METHODS[1:]},
            })
        result.extend(pair_rows)
        result.append({"Dataset": dataset, "Transition": "Mean",
                       **{field: mean(row[field] for row in pair_rows) for field in fields}})
    if used_pairs != stift.keys():
        raise ValueError("STIFT scores include unexpected pairs.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=("Dataset", "Transition", *fields))
        writer.writeheader()
        writer.writerows(result)


if __name__ == "__main__":
    main()
