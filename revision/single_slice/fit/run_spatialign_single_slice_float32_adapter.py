"""Independent fits using the saved explore_spatiAlign native workflow.

The spatialign package is an unchanged copy of the user's installed package.
Only the dispatch and input/output paths differ from the saved joint script.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import torch
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from inputs import digest, load_case, save_json


def package_hashes(root: Path) -> dict[str, str]:
    return {
        str(path.relative_to(root)): digest(path)
        for path in sorted(root.rglob("*.py"))
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--manifest", type=Path, required=True, help="JSON produced by prepare_manifest.py")
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gpu", type=int, required=True)
    args = parser.parse_args()
    cfg = load_case(args.manifest, args.dataset)
    source = args.raw_root / cfg["relative_input"]
    if digest(source) != cfg["sha256"]:
        raise ValueError("Original raw H5AD hash mismatch")
    if args.output.exists():
        raise FileExistsError(args.output)
    if args.gpu % 2 != 0:
        raise ValueError("Native spatialign device selector accepts only even GPU IDs")
    if not torch.cuda.is_available() or args.gpu >= torch.cuda.device_count():
        raise RuntimeError("Requested CUDA device unavailable")
    assert torch.ones(1, device=f"cuda:{args.gpu}").item() == 1
    sys.path.insert(0, str(args.package_root))
    from spatialign import Spatialign  # noqa: PLC0415

    native_root = args.package_root / "spatialign"
    hashes = package_hashes(native_root)
    if not hashes:
        raise RuntimeError("Native spatialign source missing")
    data = ad.read_h5ad(source)
    if "counts" in data.layers:
        data.X = data.layers["counts"].copy()
    values = data.X.data if hasattr(data.X, "data") and hasattr(data.X, "tocsr") else np.asarray(data.X)
    if not np.array_equal(values, np.asarray(values, dtype=np.float32).astype(values.dtype)):
        raise ValueError("Float32 conversion would change raw count values")
    data.X = data.X.astype(np.float32)
    times = data.obs["timepoint"].astype(str).to_numpy()
    if set(times) != set(cfg["k"]):
        raise ValueError("Unexpected timepoint set")
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError("Nonunique source identifiers")
    if not np.array_equal(
        np.asarray(data.obsm["spatial"]),
        data.obs[["pixel_x", "pixel_y"]].to_numpy(),
    ):
        raise ValueError("Spatial coordinates changed")
    args.output.mkdir(parents=True)
    summary: dict = {
        "dataset": args.dataset,
        "method": "spatialign raw-all-genes native normalize/log + Ward",
        "source_sha256": cfg["sha256"],
        "runner_sha256": digest(Path(__file__)),
        "runner_source": Path(__file__).name,
        "input_dtype_adapter": "lossless float64 to float32 raw-count cast before native normalize/log",
        "package_hashes": hashes,
        "gpu": args.gpu,
        "native_parameters": {
            "is_norm_log": True,
            "is_scale": False,
            "is_hvg": False,
            "n_neigh": 15,
            "is_undirected": True,
            "latent_dims": 100,
            "seed": 42,
            "max_epoch": 500,
            "patient": 15,
            "clustering": "Ward agglomerative at fixed evaluation K",
        },
        "slices": {},
    }
    save_json(args.output / "progress.json", summary)
    for tp, k in cfg["k"].items():
        start = time.monotonic()
        block = args.output / tp
        block.mkdir()
        one = data[times == tp].copy()
        truth = one.obs["annotation"].astype(str).to_numpy()
        one.obs = one.obs[["timepoint", "pixel_x", "pixel_y"]].copy()
        one.write_h5ad(block / "raw_slice.h5ad")
        raw_sha = digest(block / "raw_slice.h5ad")
        model = Spatialign(
            str(block / "raw_slice.h5ad"),
            batch_key="batch",
            is_norm_log=True,
            is_scale=False,
            n_neigh=15,
            is_undirected=True,
            latent_dims=100,
            seed=42,
            gpu=args.gpu,
            save_path=str(block) + os.sep,
            is_verbose=True,
        )
        if model.device.type != "cuda":
            raise RuntimeError("Native model unexpectedly selected CPU")
        if model.dataset.n_domain != 1:
            raise RuntimeError("Fit is not independent single-slice")
        if model.dataset.merge_data.n_obs != one.n_obs:
            raise RuntimeError("Native filtering removed spots; support differs")
        model.train()
        model.alignment()
        correct = ad.read_h5ad(block / "res" / "correct_data0.h5ad")
        if not np.array_equal(correct.obs_names.to_numpy(), one.obs_names.to_numpy()):
            raise RuntimeError("Native output spot order mismatch")
        latent = np.asarray(correct.obsm["correct"])
        if latent.shape != (one.n_obs, 100) or not np.isfinite(latent).all():
            raise RuntimeError("Invalid native spatialign embedding")
        pred = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(latent)
        if len(np.unique(pred)) != k:
            raise RuntimeError("Ward did not produce K clusters")
        result = block / "result.npz"
        np.savez_compressed(
            result,
            obs_names=np.asarray(one.obs_names, dtype=str),
            predicted=pred.astype(np.int16),
            latent=latent,
            selected_genes=np.asarray(model.dataset.inner_genes, dtype=str),
        )
        summary["slices"][tp] = {
            "n_spots": one.n_obs,
            "n_genes_after_native_filter": len(model.dataset.inner_genes),
            "k": k,
            "ari": float(adjusted_rand_score(truth, pred)),
            "nmi": float(normalized_mutual_info_score(truth, pred)),
            "raw_slice_sha256": raw_sha,
            "result_sha256": digest(result),
            "checkpoint_sha256": digest(block / "ckpt" / "spatialign.bgi"),
            "corrected_h5ad_sha256": digest(block / "res" / "correct_data0.h5ad"),
            "seconds": time.monotonic() - start,
        }
        save_json(args.output / "progress.json", summary)
        print(f"COMPLETED {args.dataset} {tp}: {summary['slices'][tp]}", flush=True)
    save_json(args.output / "summary.json", summary)


if __name__ == "__main__":
    main()
