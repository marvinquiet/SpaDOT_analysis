"""Compute fresh mouse domains with the released SpaDOT 1.0.7 implementation."""
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import argparse

import anndata
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score
from SpaDOT.utils import _analyze_utils

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--latent", type=Path, required=True)
parser.add_argument("--reference", type=Path, required=True, help="Verified prior plot labels keyed by spot")
parser.add_argument("--domain-key", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
stage = args.output
stage.mkdir(parents=True, exist_ok=True)
source = Path(_analyze_utils.__file__).resolve()
assert version("SpaDOT") == "1.0.7"
source_hash = sha256(source.read_bytes()).hexdigest()
assert source_hash == "6b4e43b0e5f88b3ddb8d064d6e510065a3648e2fe8f847fb54d384a5007578f2"
latent = args.latent
latent_hash = sha256(latent.read_bytes()).hexdigest()
assert latent_hash == "64c87116af1c32207692850641c7fd13de3b236ec9870cda134931063173e3df"
adata = anndata.read_h5ad(latent)
assert adata.shape == (8255, 20) and adata.obs_names.is_unique
assert np.isfinite(adata.X).all()
times = adata.obs.timepoint.astype(str)
ids = pd.Index([f"{t}_{spot}" for t, spot in zip(times, adata.obs_names)])
reference = pd.read_csv(args.reference, index_col="spot", dtype={"spot": str})
assert ids.is_unique and reference.index.is_unique and set(ids) == set(reference.index)
reference = reference.loc[ids]
np.testing.assert_array_equal(times, reference.timepoint)
np.testing.assert_array_equal(adata.obs.annotation.astype(str), reference.annotation)
np.testing.assert_array_equal(adata.obs[["pixel_x", "pixel_y"]].to_numpy(), adata.obsm["spatial"])
np.testing.assert_array_equal(adata.obsm["spatial"], reference[["spatial_x", "spatial_y"]].to_numpy())

# These are the example notebook's final KMeans settings.
predicted = _analyze_utils.KMeans_Clustering(adata, [5, 5, 5], seed=2026, n_init=10)
assert predicted.obs_names.is_unique and set(predicted.obs_names) == set(adata.obs_names)
new_labels = predicted.obs.kmeans.reindex(adata.obs_names).astype(str).to_numpy()
key = pd.read_csv(args.domain_key, dtype=str).set_index("domain_id")
plot_data = adata.obs[["timepoint", "pixel_x", "pixel_y", "annotation"]].copy()
plot_data = plot_data.rename(columns={"annotation": "true_labels"})
plot_data["prediction"] = ""
checks = []
for t, day, spots in (("t0", "E12.5", 2397), ("t1", "E13.5", 2312), ("t2", "E15.5", 3546)):
    mask = times == t
    fresh = new_labels[mask.to_numpy()]
    previous = reference.loc[mask.to_numpy(), "SpaDOT_original_KMeans"].to_numpy()
    table = pd.crosstab(fresh, previous)
    assert table.shape == (5, 5) and len(fresh) == spots
    rows, columns = linear_sum_assignment(-table.to_numpy())
    aligned = dict(zip(table.index[rows], table.columns[columns]))
    paper_labels = {new: key.loc[old, "paper_label"] for new, old in aligned.items()}
    assert set(table.columns) == set(key.index[key.stage == day])
    plot_data.loc[mask, "prediction"] = [paper_labels[x] for x in fresh]
    checks.append(dict(stage=day, spots=spots, partition_ari_vs_prior=adjusted_rand_score(previous, fresh),
                       spot_agreement_after_label_matching=float(table.to_numpy()[rows, columns].sum()/spots),
                       label_map=paper_labels))
assert (plot_data.prediction != "").all()
plot_data.to_csv(stage / "latest_spots.csv", index=False)
audit = dict(version=version("SpaDOT"), source_sha256=source_hash, latent_sha256=latent_hash,
             seed=2026, n_init=10, cluster_counts=[5, 5, 5], stages=checks)
(stage / "prediction_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
print(json.dumps(audit), flush=True)
