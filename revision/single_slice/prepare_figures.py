"""Collect evaluate.py result JSON files for Supplementary Figs. S42–S46.

Run as 'prepare_figures.py real --results FILE_OR_DIRECTORY ... --output DIR'
for S42–S43, or 'prepare_figures.py simulation --results FILE_OR_DIRECTORY
... --rescaled S44_SCORES.csv --output DIR' for S44–S46. Directories are
searched for result.json files. The result's case identifies the real
dataset (MouseEmbryo, ChickenHeart, AxolotlTelencephalon) or simulation
(e.g. DEG100_seed1105_prop0.5). The rescaled CSV is produced by
rescale_simulation_ari.py and is joined to the original results for S44 only.
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


COLORS = {"SpaDOT": "#c1121f", "PRECAST": "#68a756", "SEDR": "#f09232", "spaGCN": "#4c2623",
          "GraphST": "#165451", "spatialPCA": "#d3b300", "spatiAlign": "#ded6c1", "spaVAE": "#0082c0"}
ALIASES = {"SpaDOT_99epoch": "SpaDOT", "SpaGCN": "spaGCN", "SpatialPCA": "spatialPCA",
           "Spatialign": "spatiAlign", "SpaVAE": "spaVAE", "GraphST_stack3000": "GraphST"}
METRICS = ("annotation_naive_accuracy", "predicted_domain_naive_accuracy", "predicted_domain_projected_naive_accuracy")
REAL_NAMES = {"MouseEmbryo": ("Mouse Embryo", "Mouse Embryo"),
              "ChickenHeart": ("Chicken Heart", "Chicken Heart"),
              "AxolotlTelencephalon": ("Axolotl", "Axolotl Telencephalon")}


def input_records(paths):
    visited = set()
    for supplied in paths:
        candidates = sorted(supplied.rglob("result.json")) if supplied.is_dir() else [supplied]
        for path in candidates:
            result = json.loads(path.read_text())
            result["method"] = ALIASES.get(result["method"], result["method"])
            key = (result["case"], result["method"], result["mode"])
            if key in visited:
                raise ValueError(f"Duplicate result: {key}")
            visited.add(key)
            yield result


def simulate_tables(ari, scores, curves, rescaled, output):
    expected = ari.loc[(ari.method == "SpaDOT") & (ari["mode"] == "multi"),
                       ["case", "timepoint", "ari"]]
    scaled = pd.read_csv(rescaled)
    old = expected.merge(scaled, on=["case", "timepoint"], how="outer", validate="one_to_one")
    if len(old) != 135 or not np.allclose(old.ari, old.original_ari, atol=1e-12, rtol=0):
        raise ValueError("The S44 rescaling does not match the original SpaDOT ARIs for all 45 cases")
    included = ari.loc[((ari.method == "SpaDOT") & (ari["mode"] == "multi")) |
                       ((ari.method != "SpaDOT") & (ari["mode"] == "single"))].copy()
    included = included.merge(scaled[["case", "timepoint", "scaled_ari"]], on=["case", "timepoint"],
                              how="left", validate="many_to_one")
    included.loc[included.method == "SpaDOT", "ari"] = included.loc[included.method == "SpaDOT", "scaled_ari"]
    included = included.drop(columns="scaled_ari")
    groups = ["deg", "proportion", "method", "mode"]
    if len(included) != 1080 or not included.groupby(groups + ["seed"]).timepoint.nunique().eq(3).all():
        raise ValueError("S44 requires all eight methods across 45 cases and three stages")
    means = included.groupby(groups + ["seed"], as_index=False).ari.mean().rename(columns={"ari": "mean_ari_t0_t1_t2"})
    summary = means.groupby(groups).mean_ari_t0_t1_t2.agg(["mean", "std", "count"]).reset_index()
    if not summary["count"].eq(3).all():
        raise ValueError("S44 requires three independent seeds per condition")
    included.to_csv(output / "included_ARI_by_timepoint.csv", index=False)
    means.to_csv(output / "ARI_timepoint_mean_per_seed.csv", index=False)
    summary.to_csv(output / "ARI_plot_summary.csv", index=False)
    scores.to_csv(output / "transition_accuracy_by_scenario_method.csv", index=False)
    folds = []
    for (deg, proportion, pair), part in scores.groupby(["deg", "proportion", "transition"]):
        for metric in METRICS:
            reference = part.loc[(part.method == "SpaDOT") & (part["mode"] == "multi")].set_index("seed")[metric].dropna()
            candidates = []
            for method in COLORS.keys() - {"SpaDOT"}:
                comp = part.loc[(part.method == method) & (part["mode"] == "single")].set_index("seed")[metric].dropna()
                paired = reference.index.intersection(comp.index)
                if len(paired) == 3:
                    candidates.append((float(comp.loc[paired].mean()), method, float(reference.loc[paired].mean())))
            if candidates:
                best, method, ref = max(candidates)
                folds.append(dict(deg=deg, proportion=proportion, transition=pair, metric=metric,
                                  SpaDOT_mean=ref, second_method=method, second_mean=best, fold=ref / best if best else None,
                                  paired_seeds=3, complete_competitors=len(candidates)))
    pd.DataFrame(folds).to_csv(output / "SpaDOT_fold_changes.csv", index=False)
    (output / "distance_CDFs.json").write_text(json.dumps(curves, allow_nan=False) + "\n")


def real_tables(ari, scores, curves, output):
    selected = ari.loc[((ari.method == "SpaDOT") & (ari["mode"] == "multi")) |
                       ((ari.method != "SpaDOT") & (ari["mode"] == "single"))].copy()
    selected["dataset_label"] = selected.case.map(lambda case: REAL_NAMES[case][0])
    selected["setting"] = selected["mode"].map({"multi": "Multi-slice", "single": "Single-slice"})
    means = selected.groupby(["dataset_label", "method", "setting"], as_index=False).agg(
        mean_ari=("ari", "mean"), timepoint_count=("timepoint", "nunique"))
    selected.to_csv(output / "plotted_scores.csv", index=False)
    means.to_csv(output / "dataset_mean_ARI.csv", index=False)
    pd.DataFrame(COLORS.items(), columns=["method", "hex"]).to_csv(output / "method_palette.csv", index=False)
    order = ["SpaDOT", "SEDR", "spaGCN", "spatialPCA", "spatiAlign", "spaVAE", "GraphST", "PRECAST"]
    pd.DataFrame({"method": order}).to_csv(output / "method_display_order.csv", index=False)
    comparison = scores.loc[((scores.method == "SpaDOT") & (scores["mode"] == "multi")) |
                            ((scores.method != "SpaDOT") & (scores["mode"] == "single"))].copy()
    comparison["Data"] = comparison.case.map(lambda case: REAL_NAMES[case][1])
    comparison["Day"] = comparison.transition.str.replace("t", "").str.replace("_", "->")
    comparison["method"] = comparison.method.replace({"GraphST": "GraphST (positional 3000)",
                                                       "spaGCN": "SpaGCN", "spatialPCA": "SpatialPCA"})
    comparison.rename(columns={"annotation_naive_accuracy": "naive_accuracy"}).to_csv(
        output / "real_transition_scores.csv", index=False)
    rows = []
    for curve in curves:
        transition_index = int(curve["transition"].split("_")[0][1:])
        for key, mode in (("single", "Single-slice"), ("multi", "Multi-slice")):
            dist = curve[key]
            rows.extend(dict(dataset=curve["case"], method=curve["method"], transition_index=transition_index,
                             transition=curve["transition"], mode=mode, distance=x, cdf=y, n_pairs=dist["n_pairs"])
                        for x, y in zip(dist["distance"], dist["cdf"]))
    pd.DataFrame(rows).to_csv(output / "raw_distance_CDF_points.csv.gz", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scope", choices=("real", "simulation"))
    parser.add_argument("--results", nargs="+", type=Path, required=True)
    parser.add_argument("--rescaled", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.scope == "simulation" and args.rescaled is None:
        parser.error("simulation requires --rescaled")
    ari_rows, score_rows, curves = [], [], []
    for result in input_records(args.results):
        meta = {key: result[key] for key in ("case", "method", "mode")}
        if args.scope == "simulation":
            match = re.fullmatch(r"DEG(\d+)_seed(\d+)_prop(0\.\d+)", result["case"])
            if not match:
                raise ValueError(f"Invalid simulation case: {result['case']}")
            deg, seed, prop = match.groups()
            meta.update(deg=int(deg), seed=int(seed), proportion=float(prop))
        ari_rows.extend({**meta, "timepoint": time, **row} for time, row in result["ari"].items())
        for pair, row in result["transitions"].items():
            score_rows.append({**meta, "transition": pair, "status": row["status"],
                               **{metric: row.get("scores", {}).get(metric) for metric in METRICS},
                               "warnings": "; ".join(row.get("warnings", []))})
        curves.extend({**meta, "transition": pair, **row} for pair, row in result["distances"].items())
    args.output.mkdir(parents=True, exist_ok=False)
    ari = pd.DataFrame(ari_rows)
    scores = pd.DataFrame(score_rows)
    if args.scope == "simulation":
        simulate_tables(ari, scores, curves, args.rescaled, args.output)
    else:
        real_tables(ari, scores, curves, args.output)


if __name__ == "__main__":
    main()
