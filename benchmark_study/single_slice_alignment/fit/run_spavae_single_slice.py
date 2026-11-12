"""Independent raw-slice fits using the unmodified explore_spaVAE batch runner."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp
import torch
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from inputs import digest, load_case, save_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--manifest", type=Path, required=True, help="JSON produced by prepare_manifest.py")
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--native-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gpu", type=int, required=True)
    args = parser.parse_args()
    cfg = load_case(args.manifest, args.dataset)
    source = args.raw_root / cfg["relative_input"]
    if digest(source) != cfg["sha256"]:
        raise ValueError("Original raw source SHA mismatch")
    if args.output.exists():
        raise FileExistsError(args.output)
    native_script = args.native_source / "run_spaVAE_Batch.py"
    hashes = {p.name: digest(p) for p in sorted(args.native_source.glob("*.py"))}
    if native_script.name not in hashes:
        raise RuntimeError("Native package missing")
    if not torch.cuda.is_available() or args.gpu >= torch.cuda.device_count():
        raise RuntimeError("Requested CUDA GPU unavailable")
    assert torch.ones(1, device=f"cuda:{args.gpu}").item() == 1
    data = ad.read_h5ad(source)
    if "counts" in data.layers:
        data.X = data.layers["counts"].copy()
    times = data.obs["timepoint"].astype(str).to_numpy()
    if set(times) != set(cfg["k"]):
        raise ValueError("Unexpected timepoint set")
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError("Nonunique source identifiers")
    if not np.array_equal(data.obsm["spatial"], data.obs[["pixel_x", "pixel_y"]].to_numpy()):
        raise ValueError("Spatial coordinates differ")
    args.output.mkdir(parents=True)
    summary = {
        "dataset": args.dataset,
        "method": "spaVAE native one-batch raw-all-genes normalization + KMeans",
        "source_sha256": cfg["sha256"],
        "runner_sha256": digest(Path(__file__)),
        "native_source_hashes": hashes,
        "gpu": args.gpu,
        "native_parameters": {"maxiter": 5000, "patience": 200,
                              "GP_dim": 2, "Normal_dim": 8,
                              "batch_size": "auto", "batch_count_per_fit": 1,
                              "clustering": "KMeans random_state=0, native default n_init"},
        "slices": {},
    }
    save_json(args.output / "progress.json", summary)
    for tp, k in cfg["k"].items():
        start = time.monotonic()
        block = args.output / tp
        block.mkdir()
        part = data[times == tp].copy()
        counts = part.X.toarray() if sp.issparse(part.X) else np.asarray(part.X)
        if not np.isfinite(counts).all() or np.any(counts < 0) or np.any(counts != np.floor(counts)):
            raise ValueError("Raw counts invalid")
        selected = np.asarray(part.var_names[counts.sum(axis=0) > 0], dtype=str)
        native_h5 = block / "data.h5"
        with h5py.File(native_h5, "w") as handle:
            handle.create_dataset("X", data=counts.astype(np.float64))
            handle.create_dataset("pos", data=np.asarray(part.obsm["spatial"], dtype=np.float64))
            handle.create_dataset("batch", data=np.ones((part.n_obs, 1), dtype=np.float64))
        native_log = block / "native.log"
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
        env["OMP_NUM_THREADS"] = "4"
        env["MKL_NUM_THREADS"] = "4"
        command = [sys.executable, "-B", str(native_script),
                   "--data_file", str(native_h5),
                   "--model_file", str(block / "model.pt"),
                   "--final_latent_file", str(block / "final_latent.txt"),
                   "--denoised_counts_file", str(block / "denoised_counts.txt"),
                   "--device", "cuda"]
        with native_log.open("w") as stream:
            subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
        latent = np.loadtxt(block / "final_latent.txt", delimiter=",")
        if latent.shape != (part.n_obs, 10) or not np.isfinite(latent).all():
            raise ValueError("Native spaVAE latent invalid")
        pred = KMeans(n_clusters=k, random_state=0).fit_predict(latent)
        if len(np.unique(pred)) != k:
            raise ValueError("KMeans produced wrong number of clusters")
        result = block / "result.npz"
        np.savez_compressed(result,
                            obs_names=np.asarray(part.obs_names, dtype=str),
                            predicted=pred.astype(np.int16),
                            latent=latent,
                            selected_genes=selected)
        truth = part.obs["annotation"].astype(str).to_numpy()
        summary["slices"][tp] = {
            "n_spots": part.n_obs,
            "n_raw_genes": part.n_vars,
            "n_genes_after_native_filter": len(selected),
            "k": k,
            "ari": float(adjusted_rand_score(truth, pred)),
            "nmi": float(normalized_mutual_info_score(truth, pred)),
            "native_h5_sha256": digest(native_h5),
            "native_log_sha256": digest(native_log),
            "model_sha256": digest(block / "model.pt"),
            "latent_file_sha256": digest(block / "final_latent.txt"),
            "result_sha256": digest(result),
            "seconds": time.monotonic() - start,
        }
        save_json(args.output / "progress.json", summary)
        print(f"COMPLETED spaVAE {args.dataset} {tp}: {json.dumps(summary['slices'][tp])}", flush=True)
    save_json(args.output / "summary.json", summary)


if __name__ == "__main__":
    main()
