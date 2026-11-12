"""Align the saved Figure 3B predictions and rerun the current SpaDOT KMeans."""
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import re
import argparse

import anndata
import h5py
import numpy as np
import pandas as pd
from SpaDOT.utils import _analyze_utils

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, required=True, help="Directory for aligned predictions")
parser.add_argument("--benchmark", type=Path, default=Path("."), help="SpaDOT_analysis checkout with saved baseline predictions")
parser.add_argument("--latent", type=Path, required=True, help="Final SpaDOT 20D latent H5AD")
parser.add_argument("--reference", type=Path, required=True,
                    help="Matched STIFT CSV: spot,timepoint,annotation,spatial_x,spatial_y,SpaDOT_original_KMeans,STIFT_KMeans_per_timepoint")
parser.add_argument("--domain-key", type=Path, help="Paper domain key, required for Axolotl spatial plots")
parser.add_argument("--verify-archive", action="store_true", help="Require the original latent and reference file hashes")
args = parser.parse_args()
stage, benchmark = args.output, args.benchmark
stage.mkdir(parents=True, exist_ok=True)
methods = ["SpaDOT", "spatiAlign", "spatialPCA", "spaVAE", "SPIRAL", "SEDR", "STADIA",
           "PRECAST", "spaGCN", "STAligner", "GraphST", "STAMP", "STIFT"]
expected_counts = {"t0": 2397, "t1": 2312, "t2": 3546}
release_path = Path(_analyze_utils.__file__).resolve()
assert version("SpaDOT") == "1.0.7"
source_hash = sha256(release_path.read_bytes()).hexdigest()
if args.verify_archive:
    assert source_hash == "6b4e43b0e5f88b3ddb8d064d6e510065a3648e2fe8f847fb54d384a5007578f2"
latent_path = args.latent
latent_hash = sha256(latent_path.read_bytes()).hexdigest()
if args.verify_archive:
    assert latent_hash == "64c87116af1c32207692850641c7fd13de3b236ec9870cda134931063173e3df"
reference_path = args.reference
reference_hash = sha256(reference_path.read_bytes()).hexdigest()
if args.verify_archive:
    assert reference_hash == "c7144fb4038dda7b2e4370c98db9d1c83dbb7bedd48d7996bb79b7db067791f2"

def decode(group, field):
    node = group[field]
    if isinstance(node, h5py.Group):
        cats = node["categories"]
        values = cats.asstr()[:] if cats.dtype.kind in "OSU" else cats[:]
        codes = node["codes"][:]
        assert (codes >= 0).all(), field
        return values[codes]
    return node.asstr()[:] if node.dtype.kind in "OSU" else node[:]

def normalized_id(time, barcode, method):
    if method == "spatialPCA":
        barcode = re.sub(r"^t[012]_", "", barcode)
    if method in {"spatiAlign", "SEDR", "spaGCN"}:
        barcode = re.sub(r"(_\d+)(?:-\d+)+$", r"\1", barcode)
    return f"{time}_{barcode}"

adata = anndata.read_h5ad(latent_path)
assert adata.shape == (8255, 20) and adata.obs_names.is_unique and np.isfinite(adata.X).all()
times = adata.obs.timepoint.astype(str)
ids = pd.Index([f"{t}_{n}" for t, n in zip(times, adata.obs_names)], name="spot")
assert ids.is_unique and times.value_counts().to_dict() == expected_counts
raw_xy = adata.obs[["pixel_x", "pixel_y"]].to_numpy()
np.testing.assert_array_equal(raw_xy, adata.obsm["spatial"])
truth = adata.obs.annotation.astype(str).to_numpy()

reference = pd.read_csv(reference_path, index_col="spot", dtype={"spot": str})
assert reference.index.is_unique and set(reference.index) == set(ids)
reference = reference.loc[ids]
np.testing.assert_array_equal(reference[["spatial_x", "spatial_y"]].to_numpy(), raw_xy)
np.testing.assert_array_equal(reference.annotation, truth)
np.testing.assert_array_equal(reference.timepoint, times)
latest = _analyze_utils.KMeans_Clustering(adata, [5, 5, 5], seed=2026, n_init=10)
assert latest.obs_names.is_unique and set(latest.obs_names) == set(adata.obs_names)
labels = latest.obs.kmeans.reindex(adata.obs_names).astype(str).to_numpy()
labels = np.array([f"{t}_{c}" for t, c in zip(times, labels)])
np.testing.assert_array_equal(labels, reference.SpaDOT_original_KMeans)

table = pd.DataFrame(dict(timepoint=times.to_numpy(), pixel_x=raw_xy[:, 0],
                          pixel_y=raw_xy[:, 1], true_labels=truth), index=ids)
table["SpaDOT"] = labels
table["STIFT"] = reference.STIFT_KMeans_per_timepoint
assert not table.STIFT.isna().any()
canonical = table.reset_index().set_index(["timepoint", "pixel_x", "pixel_y"]).spot
assert canonical.index.is_unique
sources = []
native_order = pd.DataFrame({"SpaDOT": np.arange(len(ids)), "STIFT": np.arange(len(ids))}, index=ids)

for method in methods[1:-1]:
    source = benchmark / f"explore_{method}/ATAC-Visium_noremoval_results/prediction_adata.h5ad"
    with h5py.File(source) as h5:
        obs = h5["obs"]
        raw_time = decode(obs, "timepoint").astype(str)
        raw_ids = decode(obs, obs.attrs["_index"]).astype(str)
        columns = {field: decode(obs, field) for field in ("pixel_x", "pixel_y", "kmeans")}
        gt_field = "true_labels" if "true_labels" in obs else "annotation"
        columns["truth"] = decode(obs, gt_field).astype(str)
        if "pred_labels" in obs:
            columns["pred_labels"] = decode(obs, "pred_labels").astype(str)
    normalized = pd.Index([normalized_id(t, n, method) for t, n in zip(raw_time, raw_ids)])
    position = pd.MultiIndex.from_arrays([raw_time, columns["pixel_x"], columns["pixel_y"]])
    assert len(position) == 8255 and position.is_unique
    assert set(position) == set(canonical.index), (method, len(set(position) & set(canonical.index)))
    matched = pd.Index(canonical.loc[position].to_numpy())
    matches = int(np.sum(normalized.to_numpy() == matched.to_numpy()))
    assert matches == len(matched), (method, matches)
    np.testing.assert_array_equal(table.loc[matched, "true_labels"], columns["truth"])
    np.testing.assert_array_equal(table.loc[matched, "timepoint"], raw_time)
    baseline = pd.Series(columns["kmeans"].astype(str), index=matched)
    assert baseline.index.is_unique and set(baseline.index) == set(ids)
    table[method] = baseline.reindex(ids)
    native_order[method] = pd.Series(np.arange(len(matched)), index=matched).reindex(ids)
    if "pred_labels" in columns:
        check = pd.DataFrame(dict(t=raw_time, old=columns["pred_labels"], new=columns["kmeans"]))
        for _, frame in check.groupby("t"):
            assert len(frame[["old", "new"]].drop_duplicates()) == frame.old.nunique() == frame.new.nunique()
    source_digest = sha256(baseline.reindex(ids).to_csv(header=False).encode()).hexdigest()
    sources.append(dict(method=method, path=str(source), column="kmeans", spots=len(baseline),
                        all_spot_ids_gt_coords_matched=True,
                        source_order_matches_canonical=bool(matched.equals(ids)),
                        aligned_labels_sha256=source_digest))
    print(method, "8255 matched spots", flush=True)

assert set(table.columns[4:]) == set(methods) and not table.isna().any().any()
for t, frame in table.groupby("timepoint"):
    for method in ("SpaDOT", "STIFT"):
        assert frame[method].nunique() == 5
table[["timepoint", "pixel_x", "pixel_y", "true_labels", *methods]].to_csv(stage / "predictions_wide.csv")
native_order[methods].to_csv(stage / "native_order.csv")
audit = dict(spots=len(table), stages=expected_counts, methods=methods,
             spadot_version="1.0.7", spadot_source_sha256=source_hash,
             spadot_input_sha256=latent_hash, matched_stift_table_sha256=reference_hash,
             spadot_and_stift_readout=dict(method="per-stage KMeans", k=[5, 5, 5], seed=2026, n_init=10),
             latest_spadot_labels_exact_reference=True, historical_sources=sources)
(stage / "prediction_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
