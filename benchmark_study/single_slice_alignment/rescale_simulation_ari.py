"""Compute the post hoc simulation-only SpaDOT readout for Fig. S44.

--manifest is a CSV with case, deg, seed, proportion and embedding columns.
Each embedding is the saved original 20D SpaDOT H5AD with t0/t1/t2 in
obs['timepoint'] and domain labels in obs['annotation']. Paths are relative
to --embedding-root. The first ten latent dimensions are SVGP and the last
ten are GAT. This development-selected variance adjustment is applied to
the ARI readout for all scenarios. It is not used for WOT or distances.
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits


def ari(truth, latent, k):
    fit = KMeans(n_clusters=k, random_state=2026, n_init=10, init="k-means++", algorithm="lloyd")
    with threadpool_limits(limits=1):
        return float(adjusted_rand_score(truth, fit.fit_predict(latent)))


def readout(data):
    values = np.asarray(data.X)
    if values.shape != (6000, 20) or not np.isfinite(values).all():
        raise ValueError("S44 expects 6,000 spots and 20 finite SpaDOT coordinates per simulation")
    times = data.obs.timepoint.astype(str).to_numpy()
    truth = data.obs.annotation.astype(str).to_numpy()
    stage_names = ("t0", "t1", "t2")
    masks = [times == timepoint for timepoint in stage_names]
    variance = np.array([values[mask].astype(np.float64).var(axis=0) for mask in masks])
    svgp, gat = variance[:, :10].sum(axis=1).mean(), variance[:, 10:].sum(axis=1).mean()
    if not np.isfinite([svgp, gat]).all() or min(svgp, gat) <= 0:
        raise ValueError("Both encoder blocks need positive average variance")
    multiplier = float(np.sqrt(svgp / gat))
    for timepoint, mask, k in zip(stage_names, masks, (3, 4, 5)):
        if len(np.unique(truth[mask])) != k:
            raise ValueError(f"Wrong annotation domain count for {timepoint}")
        x = values[mask]
        scaled = x.copy()
        scaled[:, 10:] *= multiplier
        yield dict(timepoint=timepoint, GAT_shared_multiplier=multiplier,
                   original_ari=ari(truth[mask], x, k), scaled_ari=ari(truth[mask], scaled, k))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--embedding-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    cases = pd.read_csv(args.manifest)
    if cases.case.duplicated().any():
        raise ValueError("Duplicate simulation case")
    rows = []
    for case in cases.to_dict("records"):
        data = ad.read_h5ad(args.embedding_root / case["embedding"])
        meta = {key: case[key] for key in ("case", "deg", "seed", "proportion")}
        rows.extend({**meta, **score} for score in readout(data))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
