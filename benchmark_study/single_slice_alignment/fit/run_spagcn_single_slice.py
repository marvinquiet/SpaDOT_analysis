"""Run the saved SpaGCN single-slice workflow on original raw counts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import torch
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from inputs import digest, load_case, save_json


PACKAGE_SHA = "58cf35ee91c1caf1221da4e179bcfc9b70eeb9b313cc78e0aea0c17d4f28af05"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--manifest", type=Path, required=True, help="JSON produced by prepare_manifest.py")
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gpu", type=int, required=True)
    parser.add_argument("--package-egg", type=Path, required=True, help="Unmodified SpaGCN 1.2.7 egg used for the benchmark")
    args = parser.parse_args()
    cfg = load_case(args.manifest, args.dataset)
    source = args.raw_root / cfg["relative_input"]
    if digest(source) != cfg["sha256"] or digest(args.package_egg) != PACKAGE_SHA:
        raise ValueError("Original raw input or SpaGCN package SHA mismatch")
    if args.output.exists():
        raise FileExistsError(args.output)
    sys.path.insert(0, str(args.package_egg))
    import SpaGCN as spg  # noqa: PLC0415

    if not torch.cuda.is_available() or args.gpu >= torch.cuda.device_count():
        raise RuntimeError("CUDA unavailable")
    torch.cuda.set_device(args.gpu)
    if float(torch.ones(1, device=f"cuda:{args.gpu}").item()) != 1.0:
        raise RuntimeError("Actual CUDA check failed")
    torch.set_num_threads(1)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    data = ad.read_h5ad(source)
    if "counts" in data.layers:
        data.X = data.layers["counts"].copy()
    times = data.obs["timepoint"].astype(str).to_numpy()
    if set(times) != set(cfg["k"]):
        raise ValueError("Unexpected timepoints")
    if not np.array_equal(
        np.asarray(data.obsm["spatial"]),
        data.obs[["pixel_x", "pixel_y"]].to_numpy(),
    ):
        raise ValueError("Coordinates do not match")
    args.output.mkdir(parents=True, exist_ok=False)
    summary = {
        "dataset": args.dataset,
        "method": "SpaGCN native independent fit; all original raw genes",
        "source_sha256": cfg["sha256"],
        "package_sha256": PACKAGE_SHA,
        "script_sha256": digest(Path(__file__)),
        "gpu": args.gpu,
        "seed": 100,
        "slices": {},
    }
    save_json(args.output / "progress.json", summary)
    for tp, k in cfg["k"].items():
        started = time.monotonic()
        slice_data = data[times == tp].copy()
        matrix = slice_data.X
        values = matrix.data if sp.issparse(matrix) else np.asarray(matrix).ravel()
        if not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError("Invalid raw counts")
        # The saved explore_spaGCN single-slice workflow uses all genes,
        # scanpy's default 10,000-total normalization and log1p.
        slice_data.raw = slice_data
        sc.pp.normalize_total(slice_data)
        sc.pp.log1p(slice_data)
        xy = slice_data.obs[["pixel_x", "pixel_y"]].to_numpy(dtype=float)
        adjacency = spg.calculate_adj_matrix(
            x=xy[:, 0].tolist(),
            y=xy[:, 1].tolist(),
            beta=49,
            alpha=1,
            histology=False,
        )
        length = spg.search_l(
            0.5, adjacency, start=0.01, end=1000, tol=0.01, max_run=100
        )
        resolution = spg.search_res(
            slice_data,
            adjacency,
            length,
            k,
            start=0.7,
            step=0.1,
            tol=5e-3,
            lr=0.05,
            max_epochs=20,
            r_seed=100,
            t_seed=100,
            n_seed=100,
        )
        random.seed(100)
        np.random.seed(100)
        torch.manual_seed(100)
        clf = spg.SpaGCN()
        clf.set_l(length)
        clf.train(
            slice_data,
            adjacency,
            init_spa=True,
            init="louvain",
            res=resolution,
            tol=5e-3,
            lr=0.05,
            max_epochs=200,
        )
        pred, prob = clf.predict()
        pred = np.asarray(pred).reshape(-1).astype(str)
        prob = np.asarray(prob)
        if len(pred) != slice_data.n_obs or not np.isfinite(prob).all():
            raise ValueError("Invalid SpaGCN output")
        with torch.no_grad():
            latent_tensor, _ = clf.model.predict(clf.embed, clf.adj_exp)
        latent = latent_tensor.detach().cpu().numpy()
        if latent.shape[0] != slice_data.n_obs or not np.isfinite(latent).all():
            raise ValueError("Invalid SpaGCN latent")
        np.savez_compressed(
            args.output / f"{tp}_result.npz",
            obs_names=np.asarray(slice_data.obs_names, dtype=str),
            predicted=pred,
            latent=latent,
            probabilities=prob.astype(np.float32),
        )
        torch.save(clf.model.state_dict(), args.output / f"{tp}_model.pt")
        annotation = slice_data.obs["annotation"].astype(str).to_numpy()
        summary["slices"][tp] = {
            "n_spots": slice_data.n_obs,
            "n_genes": slice_data.n_vars,
            "k_requested": k,
            "k_predicted": len(np.unique(pred)),
            "length": float(length),
            "resolution": float(resolution),
            "ari": float(adjusted_rand_score(annotation, pred)),
            "nmi": float(normalized_mutual_info_score(annotation, pred)),
            "model_sha256": digest(args.output / f"{tp}_model.pt"),
            "result_sha256": digest(args.output / f"{tp}_result.npz"),
            "seconds": time.monotonic() - started,
        }
        save_json(args.output / "progress.json", summary)
        print(f"COMPLETED {args.dataset} {tp}: {summary['slices'][tp]}", flush=True)
    save_json(args.output / "summary.json", summary)


if __name__ == "__main__":
    main()
