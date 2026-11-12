#!/usr/bin/env python3
"""Score frozen STIFT embeddings with the study's clustering readouts.

Use --layout kmeans for the real-data and simulation benchmark. The
joint-louvain readout uses the tutorial's existing pooled labels.
For the separate, domain-count-matched Louvain comparison, use
--layout per-timepoint-louvain. The latter rebuilds the graph at each
timepoint and chooses a resolution by the target count alone. Use
--simulation-metrics to also call the repository's original simulation
benchmark ASW, PAS and CHAOS implementation on the same predictions.

Example:
    python evaluate_domains.py --input /path/to/fit/STIFT_latent.h5ad \
      --n-clusters 5 5 5 --output /path/to/mouse_domain_readout
"""
import argparse
import os
from pathlib import Path
import sys

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[key] = "1"


def louvain(z, spots, target):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc

    data = ad.AnnData(X=z.copy(), obs=pd.DataFrame(index=spots))
    data.obsm["STIFT"] = z
    sc.pp.neighbors(data, use_rep="STIFT", random_state=66)
    candidates = {}

    def cluster(resolution, key="_candidate"):
        sc.tl.louvain(data, resolution=resolution, key_added=key, random_state=0,
                      flavor="vtraag", directed=True, use_weights=False)
        return data.obs[key].astype(str).to_numpy()

    def candidate(resolution):
        resolution = round(float(resolution), 10)
        if resolution not in candidates:
            labels = cluster(resolution)
            candidates[resolution] = (len(set(labels)), labels)
        return candidates[resolution][0]

    def matches():
        return [value for value, (count, _) in candidates.items() if count == target]

    def nearest():
        return min(matches(), key=lambda value: (abs(value - 1.5), value))

    if candidate(1.5) != target:
        for value in np.geomspace(0.01, 4, 49):
            candidate(value)
        if not matches():
            ordered = sorted(candidates)
            for left, right in zip(ordered[:-1], ordered[1:]):
                if min(candidate(left), candidate(right)) <= target <= max(candidate(left), candidate(right)):
                    for value in np.linspace(left, right, 31):
                        candidate(value)
        if not matches():
            for value in np.geomspace(0.001, 8, 301):
                candidate(value)
        if not matches():
            raise ValueError(f"No evaluated Louvain resolution gave exactly {target} domains.")
        closest = nearest()
        ordered = sorted(candidates)
        if closest < 1.5:
            neighbor = next((value for value in ordered if value > closest), closest)
        else:
            neighbor = next((value for value in reversed(ordered) if value < closest), closest)
        for value in np.linspace(closest, neighbor, 21):
            candidate(value)
    resolution = nearest()
    labels = candidates[resolution][1]
    if not np.array_equal(labels, cluster(resolution, "_verify")):
        raise RuntimeError("Louvain partition was not reproduced.")
    return labels, resolution


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--n-clusters", type=int, nargs="+", required=True)
    parser.add_argument("--layout", choices=("kmeans", "joint-louvain", "per-timepoint-louvain"),
                        default="kmeans")
    parser.add_argument("--simulation-metrics", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    import anndata as ad
    import numpy as np
    import pandas as pd
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
    from sklearn.metrics import homogeneity_score, normalized_mutual_info_score

    data = ad.read_h5ad(args.input)
    z = np.asarray(data.obsm["STIFT"])
    if z.shape != (data.n_obs, 30) or not np.isfinite(z).all():
        raise ValueError("Expected the fitted 30-dimensional STIFT embeddings.")
    timepoints = data.obs["timepoint"].astype(str).to_numpy()
    times = sorted(set(timepoints))
    if len(times) != len(args.n_clusters) or min(args.n_clusters) < 1:
        raise ValueError("Supply one positive target cluster count per timepoint.")
    assigned = np.empty(data.n_obs, dtype=object)
    resolutions = {}
    # Freeze partitions before reading the annotations used for evaluation.
    for time, k in zip(times, args.n_clusters):
        keep = timepoints == time
        if args.layout == "kmeans":
            clusters = KMeans(n_clusters=k, random_state=2026, n_init=10,
                              init="k-means++", algorithm="lloyd").fit_predict(z[keep])
        elif args.layout == "per-timepoint-louvain":
            clusters, resolutions[time] = louvain(z[keep], data.obs_names[keep], k)
        else:
            clusters = data.obs.loc[keep, "louvain"].astype(str).to_numpy()
        assigned[keep] = (clusters.astype(str) if args.layout == "joint-louvain"
                          else [f"{time}_{label}" for label in clusters])

    obs = data.obs[["timepoint", "annotation"]].copy()
    if "original_spot" in data.obs:
        obs["original_spot"] = data.obs["original_spot"]
    else:
        obs["original_spot"] = [spot.removeprefix(f"{time}_")
                                for spot, time in zip(data.obs_names, timepoints)]
    obs["annotation"] = obs["annotation"].astype(str).str.replace("\n", " ", regex=False)
    obs[["spatial_x", "spatial_y"]] = data.obsm["spatial"]
    obs["prediction"] = assigned
    rows = []
    regional = []
    composition = []
    for time in times:
        part = obs.loc[timepoints == time]
        truth, predicted = part["annotation"], part["prediction"]
        rows.append(dict(timepoint=time, n_spots=len(part), n_clusters=predicted.nunique(),
                         resolution=resolutions.get(time), ARI=adjusted_rand_score(truth, predicted),
                         AMI=adjusted_mutual_info_score(truth, predicted),
                         NMI=normalized_mutual_info_score(truth, predicted),
                         H=homogeneity_score(truth, predicted)))
        counts = pd.crosstab(predicted, truth)
        for cluster, values in counts.iterrows():
            for label, count in values.items():
                composition.append(dict(timepoint=time, prediction=cluster, annotation=label,
                                        count=count, precision=count / values.sum(),
                                        recall=count / counts[label].sum()))
        for regions in (("Heart", "Liver"), ("MP", "DP"), ("MP", "LP"), ("DP", "LP")):
            if not set(regions) <= set(truth):
                continue
            keep = truth.isin(regions)
            a, b = regions
            regional.append(dict(timepoint=time, region_1=a, region_2=b,
                                 restricted_ARI=adjusted_rand_score(truth[keep], predicted[keep]),
                                 dominant_1=counts[a].idxmax(), dominant_2=counts[b].idxmax()))

    args.output.mkdir(parents=True, exist_ok=False)
    obs.to_csv(args.output / "predictions.csv", index_label="spot")
    metrics = pd.DataFrame(rows)
    metric_names = ["ARI", "AMI", "NMI", "H"]
    if args.simulation_metrics:
        utilities = (Path(__file__).resolve().parents[2] / "summarize_benchmark"
                     / "scripts" / "Simulation_study")
        sys.path.insert(0, str(utilities))
        from simulation_prop_utils import analyze_domain_metrics

        benchmark_input = obs.rename(columns={
            "annotation": "true_labels", "prediction": "pred_labels",
            "spatial_x": "pixel_x", "spatial_y": "pixel_y",
        })
        benchmark = analyze_domain_metrics(benchmark_input, str(args.output), "STIFT")
        for name in ("ASW", "PAS", "CHAOS"):
            metrics[name] = metrics["timepoint"].map(benchmark.loc[name])
            metric_names.append(name)
    metrics.to_csv(args.output / "metrics_by_timepoint.csv", index=False)
    metrics[metric_names].mean().to_frame().T.to_csv(
        args.output / "mean_metrics.csv", index=False
    )
    pd.DataFrame(composition).to_csv(args.output / "cluster_composition.csv", index=False)
    pd.DataFrame(regional).to_csv(args.output / "regional_scores.csv", index=False)


if __name__ == "__main__":
    main()
