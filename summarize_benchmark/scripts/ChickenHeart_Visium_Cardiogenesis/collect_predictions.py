"""Align the saved Figure 4B predictions and rerun the current SpaDOT KMeans."""
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
expected_counts = {"t0": 747, "t1": 1966, "t2": 1916, "t3": 1967}
release_path = Path(_analyze_utils.__file__).resolve()
assert version("SpaDOT") == "1.0.7"
source_hash = sha256(release_path.read_bytes()).hexdigest()
if args.verify_archive:
    assert source_hash == "6b4e43b0e5f88b3ddb8d064d6e510065a3648e2fe8f847fb54d384a5007578f2"
latent_path = args.latent
latent_hash = sha256(latent_path.read_bytes()).hexdigest()
if args.verify_archive:
    assert latent_hash == "9bc27874acc53d4cdf31d228cdd63aee04ca97cd3d94b46a66cf2e34a134ed4c"
reference_path = args.reference
reference_hash = sha256(reference_path.read_bytes()).hexdigest()
if args.verify_archive:
    assert reference_hash == "5c767251e6455b970ae4a78c6575e5498e88e5e6a9d35c3fdce775888d723b4d"

def decode(group, field):
    node = group[field]
    if isinstance(node, h5py.Group):
        cats = node["categories"]
        values = cats.asstr()[:] if cats.dtype.kind in "OSU" else cats[:]
        codes = node["codes"][:]
        assert (codes >= 0).all(), field
        return values[codes]
    return node.asstr()[:] if node.dtype.kind in "OSU" else node[:]

def barcode_id(name):
    match = re.search(r"(D\d+-[A-Za-z0-9]+_[ACGT]+-\d+)", name)
    assert match is not None, name
    return match.group(1)

def normalize_truth(labels):
    return pd.Series(labels).str.replace(r"\s+", " ", regex=True).str.strip().to_numpy()

adata = anndata.read_h5ad(latent_path)
assert adata.shape == (6596, 20) and adata.obs_names.is_unique and np.isfinite(adata.X).all()
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
np.testing.assert_array_equal(normalize_truth(reference.annotation), normalize_truth(truth))
np.testing.assert_array_equal(reference.timepoint, times)
latest = _analyze_utils.KMeans_Clustering(adata, [5, 7, 7, 6], seed=2026, n_init=10)
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

for method in ["historical_SpaDOT", *methods[1:-1]]:
    source_method = "SpaDOT" if method == "historical_SpaDOT" else method
    source = benchmark / "benchmark_study" / f"explore_{source_method}/ChickenHeart_results/prediction_adata.h5ad"
    if not source.is_file():
        source = benchmark / f"explore_{source_method}/ChickenHeart_results/prediction_adata.h5ad"
    with h5py.File(source) as h5:
        obs = h5["obs"]
        raw_time = decode(obs, "timepoint").astype(str)
        raw_ids = decode(obs, obs.attrs["_index"]).astype(str)
        columns = {field: decode(obs, field) for field in ("pixel_x", "pixel_y", "kmeans")}
        gt_field = "true_labels" if "true_labels" in obs else "annotation"
        columns["truth"] = decode(obs, gt_field).astype(str)
        if "pred_labels" in obs:
            columns["pred_labels"] = decode(obs, "pred_labels").astype(str)
    normalized = pd.Index([f"{t}_{barcode_id(n)}" for t, n in zip(raw_time, raw_ids)])
    matched = normalized
    assert len(matched) == 6596 and matched.is_unique and set(matched) == set(ids), method
    source_xy = np.column_stack([columns["pixel_x"], columns["pixel_y"]])
    latest_xy = table.loc[matched, ["pixel_x", "pixel_y"]].to_numpy()
    if method == "historical_SpaDOT":
        benchmark_xy = pd.DataFrame(source_xy, index=matched, columns=["pixel_x", "pixel_y"]).reindex(ids)
        different = np.any(source_xy != latest_xy, axis=1)
        assert np.array_equal(different, raw_time == "t0")
        latest_day4 = latest_xy[different]
        old_day4 = source_xy[different]
        design = np.column_stack([latest_day4, np.ones(len(latest_day4))])
        transform = np.linalg.lstsq(design, old_day4, rcond=None)[0]
        affine_error = float(np.max(np.abs(design @ transform - old_day4)))
        assert affine_error < 1
        print("Day 4 reference geometry aligned by all 747 IDs; latest-to-benchmark affine max residual", affine_error, flush=True)
    else:
        expected_xy = latest_xy if method == "STAMP" else benchmark_xy.loc[matched].to_numpy()
        np.testing.assert_array_equal(source_xy, expected_xy)
    np.testing.assert_array_equal(normalize_truth(table.loc[matched, "true_labels"].to_numpy()), normalize_truth(columns["truth"]))
    np.testing.assert_array_equal(table.loc[matched, "timepoint"], raw_time)
    baseline = pd.Series(columns["kmeans"].astype(str), index=matched)
    assert baseline.index.is_unique and set(baseline.index) == set(ids)
    if method == "historical_SpaDOT":
        lisi_reference_order = pd.Series(np.arange(len(matched)), index=matched).reindex(ids)
        continue
    table[method] = baseline.reindex(ids)
    native_order[method] = pd.Series(np.arange(len(matched)), index=matched).reindex(ids)
    if "pred_labels" in columns:
        check = pd.DataFrame(dict(t=raw_time, old=columns["pred_labels"], new=columns["kmeans"]))
        for _, frame in check.groupby("t"):
            assert len(frame[["old", "new"]].drop_duplicates()) == frame.old.nunique() == frame.new.nunique()
    source_digest = sha256(baseline.reindex(ids).to_csv(header=False).encode()).hexdigest()
    sources.append(dict(method=method, path=str(source), column="kmeans", spots=len(baseline),
                        all_spot_ids_gt_matched=True, coordinate_frame="latest latent" if method == "STAMP" else "historical benchmark",
                        source_order_matches_canonical=bool(matched.equals(ids)),
                        aligned_labels_sha256=source_digest))
    print(method, "6596 matched spots", flush=True)

assert set(table.columns[4:]) == set(methods) and not table.isna().any().any()
for t, frame in table.groupby("timepoint"):
    for method in ("SpaDOT", "STIFT"):
        assert frame[method].nunique() == dict(t0=5,t1=7,t2=7,t3=6)[t]
# Evaluate all methods on the manuscript's physical pixel frame for this tissue.
table[["pixel_x", "pixel_y"]] = benchmark_xy
table[["timepoint", "pixel_x", "pixel_y", "true_labels", *methods]].to_csv(stage / "predictions_wide.csv")
native_order[methods].to_csv(stage / "native_order.csv")
lisi_order = native_order.copy()
for method in methods[1:-1]:
    lisi_order[method] = lisi_reference_order
lisi_order[methods].to_csv(stage / "lisi_order.csv")
audit = dict(spots=len(table), stages=expected_counts, methods=methods,
             spadot_version="1.0.7", spadot_source_sha256=source_hash,
             spadot_input_sha256=latent_hash, matched_stift_table_sha256=reference_hash,
             spadot_and_stift_readout=dict(method="per-stage KMeans", k=[5, 7, 7, 6], seed=2026, n_init=10),
             latest_spadot_labels_exact_reference=True, historical_sources=sources,
             common_metric_geometry="historical SpaDOT pixel coordinates for every method",
             notebook_coordinate_difference=dict(stage="t0", spots=747, affine_max_abs_pixel_residual=affine_error,
                                                 affine_coefficients=transform.tolist()))
(stage / "prediction_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
