#!/usr/bin/env python3
"""Summarize the STIFT readouts from both manuscript simulation strategies.

Pass a CSV containing scenario, seed, fit_directory, evaluation_directory.
Directory paths may be relative to the CSV. The fit and evaluation directories
come from fit_stift.py and evaluate_domains.py, respectively. For strategy 1,
scenario IDs are 30_50, 30_30, 20_50 and 20_30. For strategy 2, use identifiers
that encode the marker-gene setting and the dominant cell-type proportion.

The output contains the individual case averages and runtimes. These can be
joined with the corresponding existing methods in the simulation benchmark.
"""
import argparse
import csv
import math
from pathlib import Path


def relative_to_manifest(value, manifest):
    return manifest.parent / value


def one_row(path):
    with path.open(newline="") as stream:
        records = list(csv.DictReader(stream))
    if len(records) != 1:
        raise ValueError(f"Expected one aggregate result in {path}.")
    return records[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    metrics = ("ARI", "AMI", "NMI", "H")
    optional = ("ASW", "PAS", "CHAOS")
    with args.cases.open(newline="") as stream:
        cases = list(csv.DictReader(stream))
    if not cases:
        raise ValueError("The case manifest is empty.")

    rows = []
    seen = set()
    optional_metrics = None
    for case in cases:
        scenario, seed = case["scenario"], int(case["seed"])
        if (scenario, seed) in seen:
            raise ValueError(f"Duplicate case {scenario}, seed {seed}.")
        seen.add((scenario, seed))
        fit = relative_to_manifest(case["fit_directory"], args.cases)
        evaluation = relative_to_manifest(case["evaluation_directory"], args.cases)
        result = one_row(evaluation / "mean_metrics.csv")
        seconds = float((fit / "STIFT_running_time.txt").read_text().strip())
        present = tuple(metric for metric in optional if metric in result)
        if optional_metrics is None:
            optional_metrics = present
        if optional_metrics != present:
            raise ValueError("All cases must have the same metric columns.")
        values = {f"mean_{key}": float(result[key]) for key in (*metrics, *present)}
        if seconds <= 0 or not all(math.isfinite(value) for value in (*values.values(), seconds)):
            raise ValueError(f"Nonfinite metrics or invalid runtime for {scenario}, seed {seed}.")
        rows.append(dict(scenario=scenario, seed=seed, method="STIFT", time=seconds, **values))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["scenario"], row["seed"])))


if __name__ == "__main__":
    main()
