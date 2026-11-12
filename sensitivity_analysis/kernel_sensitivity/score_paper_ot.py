"""Score S47 paper-parameter transition accuracy for one trained embedding.

Supply the saved epoch-99 20D latent H5AD and an external JSON defining the
reference graph. The JSON is the "graph" object from the study case manifest,
with either a "curated" dataset name or "edges" as a list per adjacent pair.
The latent must include obs.timepoint and obs.annotation. Requires Waddington-OT,
anndata, pandas, NumPy, scipy, scikit-learn and networkx.
"""

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
from wot.ot import OTModel

from collaborator_metrics import accuracy, min_norm, norm, transition_graphs

PAPER_OT = dict(epsilon=0.01, epsilon0=1, lambda1=1.0, lambda2=50.0,
                growth_iters=3, tau=1000)


def annotation_mass(plan, annotations):
    """Sum the original spot masses by source and target annotation."""
    mass = plan.X.toarray() if sparse.issparse(plan.X) else np.asarray(plan.X)
    if not np.isfinite(mass).all() or (mass < 0).any() or not mass.sum() > 0:
        raise ValueError("OT plan must have finite nonnegative mass and a positive total")
    source = annotations.loc[plan.obs_names].to_numpy()
    target = annotations.loc[plan.var_names].to_numpy()
    rows, cols = np.unique(source), np.unique(target)
    source_membership = (source[:, None] == rows).astype(float)
    target_membership = (target[:, None] == cols).astype(float)
    values = source_membership.T @ (mass @ target_membership)
    return pd.DataFrame(values, index=rows, columns=cols)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--latent", type=Path, required=True)
    parser.add_argument("--reference-graph", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--scale", type=float, choices=(.01, .03, .1, .3, 1.), required=True)
    parser.add_argument("--seed", type=int, choices=(666, 1993, 2026), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    latent = ad.read_h5ad(args.latent)
    if latent.n_vars != 20 or not latent.obs_names.is_unique:
        parser.error("--latent must have 20 dimensions and unique spot IDs")
    values = latent.X.data if sparse.issparse(latent.X) else latent.X
    if not np.isfinite(values).all():
        parser.error("--latent contains nonfinite values")
    latent.obs["timepoint"] = latent.obs["timepoint"].astype(str)
    latent.obs["annotation"] = latent.obs["annotation"].astype(str).map(norm)
    latent.obs["day"] = latent.obs["timepoint"].astype("category").cat.codes
    latent.obs["cell_growth_rate"] = 1
    graph_spec = json.loads(args.reference_graph.read_text())
    graphs = transition_graphs(args.dataset, graph_spec, sorted(latent.obs.annotation.unique()))
    if len(graphs) != latent.obs["day"].nunique() - 1:
        parser.error("--reference-graph must define every adjacent timepoint pair")
    args.output_dir.mkdir(parents=True, exist_ok=False)

    rows = []
    with threadpool_limits(limits=1):
        model = OTModel(latent, **PAPER_OT)
        for pair, graph in enumerate(graphs):
            plan = model.compute_transport_map(pair, pair + 1)
            mass = annotation_mass(plan, latent.obs.annotation)
            association = min_norm(mass)
            if not np.isfinite(association.to_numpy()).all():
                raise ValueError(f"Nonfinite normalized association for pair {pair}")
            mass.to_csv(args.output_dir / f"annotation_transport_mass_{pair}_{pair + 1}.csv")
            association.to_csv(args.output_dir / f"annotation_conditional_min_{pair}_{pair + 1}.csv")
            rows.append(dict(dataset=args.dataset, kernel_scale=args.scale, seed=args.seed,
                             i=pair, table="WOT_annotation",
                             naive_accuracy=accuracy(association, graph)))
    pd.DataFrame(rows).to_csv(args.output_dir / "transition_metrics.csv", index=False)


if __name__ == "__main__":
    main()
