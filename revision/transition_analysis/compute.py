"""Compute paper-setting spot OT and domain tables for Figures 3c–5c and S41.

Install h5py, numpy, pandas, scikit-learn and threadpoolctl. Supply the final
20-dimensional latent.h5ad at INPUT/DATASET/latent.h5ad, where DATASET is
MouseEmbryo, ChickenHeart or AxolotlTelencephalon. X must contain the embeddings
in original spot order, with obs columns timepoint (t0, t1, ...) and annotation.

python compute.py --latents INPUT --output OUTPUT
python score_thresholds.py --pairs OUTPUT/pairs --output SCORES
python plot_thresholds.py --rerun SCORES --output FIGURES

Use --domain-keys KEYS to export the display mapping with each dataset. This
expects KEYS/DATASET/domain_key.csv with stage, domain_id, paper_label, color.
The solver snapshots reproduce the original and finite-log duality-gap variants
used in the audited threshold rerun. The latter is used for Chicken Day 7–10.
"""

import argparse
import os
import sys
from pathlib import Path
from shutil import copyfile

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

import h5py
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances
from threadpoolctl import threadpool_limits

import ot_solver_original
import ot_solver_stable_gap

DATASETS = {
    "ChickenHeart": (5, 7, 7, 6),
    "MouseEmbryo": (5, 5, 5),
    "AxolotlTelencephalon": (6, 7, 8),
}
PARAMETERS = dict(epsilon=0.01, epsilon0=1, lambda1=1, lambda2=50, tau=1000, growth_iters=3)


def read_obs(file, name):
    column = file["obs"][name]
    if isinstance(column, h5py.Group):
        codes = column["codes"][:]
        if (codes < 0).any():
            raise ValueError(f"Missing values in obs[{name!r}]")
        return column["categories"].asstr()[:][codes].astype(str)
    return column.asstr()[:].astype(str)


def fit_labels(embedding, stage, clusters):
    with threadpool_limits(limits=1):
        labels = KMeans(n_clusters=clusters, random_state=2026, n_init=10,
                        init="k-means++", algorithm="lloyd").fit_predict(embedding)
    return np.array([f"{stage}_{label}" for label in labels])


def solve(embedding_a, embedding_b, solver):
    cost = pairwise_distances(embedding_a, embedding_b, metric="sqeuclidean", n_jobs=1)
    median = np.median(cost)
    if not np.isfinite(median) or median <= 0:
        raise ValueError("The median squared embedding distance must be positive and finite.")
    cost /= median
    config = dict(solver.default_config, **PARAMETERS)
    growth = np.ones(len(embedding_a))
    diagnostics = []
    for _ in range(PARAMETERS["growth_iters"]):
        diagnostic = {}

        def capture(frame, event, _):
            if event == "return" and frame.f_code is solver.optimal_transport_duality_gap.__code__:
                for key in ("duality_gap", "current_iter", "e"):
                    if key in frame.f_locals:
                        diagnostic[key] = float(frame.f_locals[key])

        previous = sys.getprofile()
        sys.setprofile(capture)
        try:
            with threadpool_limits(limits=1):
                plan = solver.optimal_transport_duality_gap(C=cost, G=growth, **config)
        finally:
            sys.setprofile(previous)
        if (not np.isfinite(plan).all() or (plan < 0).any() or plan.sum() <= 0
                or abs(diagnostic.get("duality_gap", np.inf)) > config["tolerance"]
                or diagnostic.get("e") != 5 or diagnostic.get("current_iter", np.inf) >= config["max_iter"]):
            raise RuntimeError(f"Spot OT did not converge: {diagnostic}")
        growth = plan.sum(axis=1)
        diagnostics.append(diagnostic)
    return plan, diagnostics


def aggregate(plan, source, target, source_annotation, target_annotation):
    rows = np.unique(source)
    cols = np.unique(target)
    left = (source[:, None] == rows).astype(float)
    right = (target[:, None] == cols).astype(float)
    with threadpool_limits(limits=1):
        mass = left.T @ (plan @ right)
    np.testing.assert_allclose(mass.sum(), plan.sum(), rtol=1e-12)

    def majority(labels, annotations, levels):
        names = []
        for level in levels:
            values, counts = np.unique(annotations[labels == level], return_counts=True)
            names.append(values[np.argmax(counts)])
        return np.asarray(names)

    return dict(mass=mass, rows=rows, cols=cols,
                row_annotation=majority(source, source_annotation, rows),
                col_annotation=majority(target, target_annotation, cols))


def conditionals(mass):
    row_total = mass.sum(axis=1, keepdims=True)
    col_total = mass.sum(axis=0, keepdims=True)
    forward = np.divide(mass, row_total, out=np.zeros_like(mass), where=row_total > 0)
    reverse = np.divide(mass, col_total, out=np.zeros_like(mass), where=col_total > 0)
    return forward, reverse, np.minimum(forward, reverse)


def run_dataset(dataset, clusters, args):
    with h5py.File(args.latents / dataset / "latent.h5ad") as file:
        embedding = file["X"][:]
        time = read_obs(file, "timepoint")
        annotation = read_obs(file, "annotation")
    stages = [f"t{i}" for i in range(len(clusters))]
    if sorted(set(time)) != stages or embedding.shape != (len(time), 20) or not np.isfinite(embedding).all():
        raise ValueError(f"Unexpected time points or latent shape for {dataset}.")
    stage_data = [(embedding[time == stage], annotation[time == stage]) for stage in stages]
    predicted = [fit_labels(x, stage, k) for (x, _), stage, k in zip(stage_data, stages, clusters)]
    plot_inputs = args.output / "inputs" / dataset
    plot_inputs.mkdir(parents=True)
    if args.domain_keys:
        copyfile(args.domain_keys / dataset / "domain_key.csv", plot_inputs / "domain_key.csv")
    for i, ((source_x, source_a), (target_x, target_a)) in enumerate(zip(stage_data, stage_data[1:])):
        pair = f"t{i}_t{i+1}"
        dest = args.output / "pairs" / dataset / pair
        dest.mkdir(parents=True)
        stable = dataset == "ChickenHeart" and i == 1
        solver = ot_solver_stable_gap if stable else ot_solver_original
        plan, diagnostics = solve(source_x, target_x, solver)
        if args.save_spot_plans:
            np.save(dest / "spot_growth3.npy", plan, allow_pickle=False)
        np.savez_compressed(dest / "annotation.npz", **aggregate(plan, source_a, target_a, source_a, target_a))
        table = aggregate(plan, predicted[i], predicted[i+1], source_a, target_a)
        np.savez_compressed(dest / "predicted.npz", **table)
        matrix = table["mass"]
        forward, reverse, score = conditionals(matrix)
        suffix = f"{i}_{i+1}"
        for name, values in (("transport_mass", matrix), ("joint", matrix / matrix.sum()),
                             ("forward_conditional", forward), ("reverse_conditional", reverse),
                             ("association_score", score)):
            pd.DataFrame(values, index=table["rows"], columns=table["cols"]).to_csv(
                plot_inputs / f"{name}_{suffix}.csv", index_label="source_domain")
        pd.DataFrame(diagnostics).to_csv(dest / "growth_diagnostics.csv", index=False)
        print(f"{dataset} {pair} completed", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--latents", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dataset", choices=DATASETS, help="Compute a single dataset instead of all three.")
    parser.add_argument("--domain-keys", type=Path)
    parser.add_argument("--save-spot-plans", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    for dataset, clusters in DATASETS.items():
        if args.dataset is None or dataset == args.dataset:
            run_dataset(dataset, clusters, args)


if __name__ == "__main__":
    main()
