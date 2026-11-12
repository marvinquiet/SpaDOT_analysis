"""Combine independently calculated ARI/AMI/H/PAS/CHAOS and median-LISI CSVs.

Each CSV uses the columns index,timepoint,value,method. Supply the output of
the corresponding dataset's calculate_python_metrics.py and calculate_lisi.R.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("python_metrics", type=Path)
    parser.add_argument("lisi_metrics", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    metrics = pd.concat([pd.read_csv(path) for path in (args.python_metrics, args.lisi_metrics)], ignore_index=True)
    groups = ["index", "method"]
    if (metrics.duplicated([*groups, "timepoint"]).any() or not np.isfinite(metrics.value).all()
            or set(metrics["index"]) != {"ARI", "AMI", "H", "PAS", "CHAOS", "median-LISI"}):
        raise ValueError("Incomplete or duplicate metric rows.")
    stages = metrics[metrics.timepoint != "avg"].groupby(groups).value.mean().sort_index()
    averages = metrics[metrics.timepoint == "avg"].set_index(groups).value.sort_index()
    pd.testing.assert_index_equal(stages.index, averages.index)
    np.testing.assert_allclose(stages, averages, rtol=0, atol=1e-12)
    metrics.to_csv(args.output, index=False, mode="x")


if __name__ == "__main__":
    main()
