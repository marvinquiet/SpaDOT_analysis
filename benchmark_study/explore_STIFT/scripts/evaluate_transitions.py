#!/usr/bin/env python3
"""Compare consecutive-timepoint transport using the paper's annotation score.

Modes use actual STIFT training couplings (native-stift), supplied spot
couplings in AnnData format (saved-spot-ot), or the common SpaDOT spot
OT readout on an embedding (matched-wot). The last mode defaults to
epsilon 0.01, lambda1 1 and lambda2 50. Supply alternative settings
explicitly when reproducing a sensitivity comparison. The output score
sums unthresholded bidirectional associations over the allowed edges
and divides by their sum over all pairs.

Example:
    python evaluate_transitions.py --dataset MouseEmbryo \
      --mode native-stift --input /path/to/STIFT_latent.h5ad \
      --maps /path/to/fit --output /path/to/transitions

For saved-spot-ot, --maps contains OT_0_1.h5ad, OT_1_2.h5ad, etc.
Their row/column identifiers must match the corresponding inputs.
"""
import argparse
import os
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"


def matched_wot(left, right, args):
    import numpy as np
    from sklearn.metrics import pairwise_distances
    from ot_solver import optimal_transport_duality_gap

    cost = pairwise_distances(left, right, metric="sqeuclidean", n_jobs=1)
    median = float(np.median(cost))
    if not np.isfinite(median) or median <= 0:
        raise ValueError("The median squared embedding distance must be positive.")
    cost /= median
    growth = np.ones(len(left))
    for _ in range(3):
        plan = optimal_transport_duality_gap(
            C=cost, G=growth, lambda1=args.lambda1, lambda2=args.lambda2,
            epsilon=args.epsilon, epsilon0=1, tau=1000, tolerance=1e-8,
            max_iter=1e7, batch_size=5
        )
        growth = plan.sum(axis=1)
    return plan


def saved_map(path, left, right):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    data = ad.read_h5ad(path)
    def order(ids, obs):
        expected = pd.Index(obs.index)
        actual = pd.Index(ids)
        if set(actual) != set(expected):
            if "original_spot" not in obs:
                raise ValueError(f"Spot IDs disagree with the map at {path}.")
            expected = pd.Index(obs["original_spot"].astype(str))
        if not actual.is_unique or set(actual) != set(expected):
            raise ValueError(f"Spot IDs disagree with the map at {path}.")
        return actual.get_indexer(expected)

    rows = order(data.obs_names, left)
    cols = order(data.var_names, right)
    matrix = data.X[rows][:, cols]
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)


def matrices(plan, source, target, output, prefix):
    import numpy as np
    import pandas as pd
    from benchmark_edges import onehot

    rows, cols = sorted(set(source)), sorted(set(target))
    mass = onehot(source, rows).T @ (plan @ onehot(target, cols))
    if not np.isclose(mass.sum(), plan.sum(), rtol=1e-12):
        raise ValueError("The domain aggregation did not preserve mass.")
    outgoing, incoming = mass.sum(axis=1, keepdims=True), mass.sum(axis=0, keepdims=True)
    forward = np.divide(mass, outgoing, out=np.zeros_like(mass), where=outgoing > 0)
    reverse = np.divide(mass, incoming, out=np.zeros_like(mass), where=incoming > 0)
    association = np.minimum(forward, reverse)
    for label, matrix in (("mass", mass), ("forward", forward), ("reverse", reverse),
                          ("association", association)):
        pd.DataFrame(matrix, index=rows, columns=cols).to_csv(output / f"{prefix}_{label}.csv")
    return association, rows, cols


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=("MouseEmbryo", "ChickenHeart", "Axolotl"), required=True)
    parser.add_argument("--mode", choices=("native-stift", "saved-spot-ot", "matched-wot"),
                        required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--maps", type=Path)
    parser.add_argument("--predictions", type=Path,
                        help="Optional evaluate_domains.py predictions.csv to aggregate cluster pairs.")
    parser.add_argument("--rep", default="STIFT", help="AnnData obsm key, or X, for matched-wot.")
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--lambda1", type=float, default=1)
    parser.add_argument("--lambda2", type=float, default=50)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode != "matched-wot" and args.maps is None:
        parser.error("--maps is required for native-stift or saved-spot-ot")

    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from benchmark_edges import edges, score

    data = ad.read_h5ad(args.input)
    if not data.obs_names.is_unique:
        raise ValueError("Expected unique spot identifiers.")
    timepoints = data.obs["timepoint"].astype(str)
    times = sorted(set(timepoints))
    annotation = data.obs["annotation"].astype(str)
    if args.predictions:
        predictions = pd.read_csv(args.predictions, index_col="spot", dtype=str)["prediction"]
        if not predictions.index.is_unique or set(predictions.index) != set(data.obs_names):
            raise ValueError("Prediction and input spot IDs must match.")
        predictions = predictions.reindex(data.obs_names)
    else:
        predictions = None
    if args.mode == "native-stift":
        order = pd.read_csv(args.maps / "spot_order.csv", index_col="spot", dtype=str)
        if list(order.index) != list(data.obs_names) or list(order["timepoint"]) != list(timepoints):
            raise ValueError("The embedding does not follow the native coupling spot order.")
    if args.mode == "matched-wot":
        embedding = data.X if args.rep == "X" else data.obsm[args.rep]
        embedding = embedding.toarray() if sparse.issparse(embedding) else np.asarray(embedding)
        embedding = np.asarray(embedding, dtype=np.float64)
        if not np.isfinite(embedding).all():
            raise ValueError("Expected finite embeddings.")

    args.output.mkdir(parents=True, exist_ok=False)
    rows = []
    for i, (left, right) in enumerate(zip(times[:-1], times[1:])):
        source, target = timepoints == left, timepoints == right
        if args.mode == "native-stift":
            plan = np.load(args.maps / f"Pi_{left}_{right}.npy", allow_pickle=False)
        elif args.mode == "saved-spot-ot":
            plan = saved_map(args.maps / f"OT_{i}_{i + 1}.h5ad",
                             data.obs.loc[source], data.obs.loc[target])
        else:
            plan = matched_wot(embedding[source], embedding[target], args)
        if (plan.shape != (source.sum(), target.sum()) or not np.isfinite(plan).all()
                or (plan < 0).any() or plan.sum() <= 0):
            raise ValueError(f"Invalid spot coupling for {left}->{right}.")
        labels_source = annotation[source].str.replace("\n", " ", regex=False)
        labels_target = annotation[target].str.replace("\n", " ", regex=False)
        allowed = {(a.replace("\n", " "), b.replace("\n", " "))
                   for a, b in edges(args.dataset, i, set(annotation))}
        prefix = f"{left}_{right}"
        association, source_names, target_names = matrices(
            plan, labels_source, labels_target, args.output, f"{prefix}_annotation"
        )
        rows.append(dict(source=left, target=right, mode=args.mode,
                         naive_accuracy=score(association, source_names, target_names, allowed)))
        if predictions is not None:
            matrices(plan, predictions[source], predictions[target], args.output, f"{prefix}_domain")
    pd.DataFrame(rows).to_csv(args.output / "transition_scores.csv", index=False)


if __name__ == "__main__":
    main()
