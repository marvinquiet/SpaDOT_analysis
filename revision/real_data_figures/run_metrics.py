"""Recompute the six real-data benchmark metrics in the publication geometry.

Requires Python with SpaDOT 1.0.7 and the analysis dependencies, and Rscript
with the lisi package. The benchmark checkout must have the saved per-method
prediction H5ADs. Pass the matching SpaDOT latent and the STIFT reference table
for this dataset. See collect_predictions.py --help for the required schema.
"""
import argparse
from pathlib import Path
import subprocess
import sys

DATASETS = ("mouse_embryo", "chicken_heart", "axolotl")
ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=DATASETS)
    parser.add_argument("--benchmark", type=Path, default=Path("."))
    parser.add_argument("--latent", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--domain-key", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-archive", action="store_true")
    args = parser.parse_args()
    if args.dataset == "axolotl" and args.domain_key is None:
        parser.error("--domain-key is required for the Axolotl display labels")

    folder = ROOT / args.dataset
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    collect = [sys.executable, "-B", folder / "collect_predictions.py", "--benchmark", args.benchmark,
               "--latent", args.latent, "--reference", args.reference, "--output", output]
    if args.domain_key:
        collect += ["--domain-key", args.domain_key]
    if args.verify_archive:
        collect.append("--verify-archive")
    subprocess.run(collect, check=True)

    python_metrics, lisi_metrics = output / "python_metrics.csv", output / "lisi_metrics.csv"
    subprocess.run([sys.executable, "-B", folder / "calculate_python_metrics.py",
                    output / "predictions_wide.csv", args.benchmark / "summarize_benchmark/scripts/utils.py",
                    python_metrics, "--native-order", output / "native_order.csv",
                    "--audit", output / "python_metrics_audit.json"], check=True)

    lisi_inputs = {
        "mouse_embryo": ("predictions_wide.csv", "native_order.csv"),
        "chicken_heart": ("predictions_wide.csv", "lisi_order.csv"),
        "axolotl": ("predictions_lisi.csv",),
    }
    subprocess.run(["Rscript", "--vanilla", folder / "calculate_lisi.R",
                    *(output / name for name in lisi_inputs[args.dataset]), lisi_metrics], check=True)
    subprocess.run([sys.executable, "-B", ROOT / "combine_metrics.py",
                    python_metrics, lisi_metrics, output / "updated_metrics.csv"], check=True)


if __name__ == "__main__":
    main()
