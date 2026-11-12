"""Combine archived Table S8 loss-ablation metrics with the audited recovery.

Provide the original matched-loss CSV directory and the independently verified
terminal JSON from the numerical recovery study. The only replacement is the
Chicken Heart alignment-off fit. The other archived fits are retained.
"""
import argparse
import csv
import json
from pathlib import Path


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, rows):
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-dir", required=True, type=Path)
    parser.add_argument("--recovery-audit", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    audit = json.loads(args.recovery_audit.read_text())
    if not audit["all_passed"]:
        raise ValueError("The recovery audit must have passed")
    new, = [case for case in audit["cases"] if case["dataset"] == "ChickenHeart"
            and case["method"] == "without_encoder_alignment"]
    if new["executed_epoch"] != 99:
        raise ValueError("The recovery must be the executed epoch-99 result")

    mean_name = "matched_loss_dataset_metrics.csv"
    stage_name = "matched_loss_timepoint_metrics.csv"
    means = read_rows(args.archive_dir / mean_name)
    stages = read_rows(args.archive_dir / stage_name)
    row, = [record for record in means if record["dataset"] == "ChickenHeart"
            and record["variant"] == "without_encoder_alignment"]
    row.update(mean_ARI=new["mean_ARI"], mean_PAS=new["mean_PAS"], executed_epoch=99,
               clustering_status="valid", result_version="numerical_recovery_20261102")
    for record in means:
        record.setdefault("result_version", "retained_verified_result")
    if any(record["dataset"] == "ChickenHeart"
           and record["variant"] == "without_encoder_alignment" for record in stages):
        raise ValueError("Input already has a Chicken Heart alignment-off result")
    for timepoint, score in new["scores"].items():
        stages.append(dict(trial=new["trial"], dataset="ChickenHeart",
            variant="without_encoder_alignment", timepoint=timepoint,
            executed_epoch=99, spots=score["spots"], clusters=score["actual_clusters"],
            ARI=score["ARI"], PAS=score["PAS"],
            abnormal_spots=round(score["PAS"] * score["spots"]), neighbor_count=10))
    assert len(means) == 18 and len(stages) == 60
    args.output_dir.mkdir(parents=True, exist_ok=False)
    write_rows(args.output_dir / mean_name, means)
    write_rows(args.output_dir / stage_name, stages)


if __name__ == "__main__":
    main()
