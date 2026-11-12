"""Independent SEDR fits with the saved explore_SEDR preprocessing."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import torch
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from inputs import digest, load_case, save_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--manifest", type=Path, required=True, help="JSON produced by prepare_manifest.py")
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--r-script", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gpu", type=int, required=True)
    args = parser.parse_args()
    cfg = load_case(args.manifest, args.dataset)
    source = args.raw_root / cfg["relative_input"]
    if digest(source) != cfg["sha256"]:
        raise ValueError("Original raw input SHA mismatch")
    if args.output.exists():
        raise FileExistsError(args.output)
    package = args.package_root / "SEDR"
    names = [
        "__init__.py", "SEDR_model.py", "SEDR_module.py", "graph_func.py",
        "utils_func.py", "clustering_func.py",
    ]
    package_hashes = {name: digest(package / name) for name in names}
    sys.path.insert(0, str(args.package_root))
    import SEDR  # noqa: PLC0415

    if not torch.cuda.is_available() or args.gpu >= torch.cuda.device_count():
        raise RuntimeError("CUDA unavailable")
    device = torch.device(f"cuda:{args.gpu}")
    assert float(torch.ones(1, device=device).item()) == 1.0
    torch.set_num_threads(1)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    data = ad.read_h5ad(source)
    if "counts" in data.layers:
        data.X = data.layers["counts"].copy()
    times = data.obs["timepoint"].astype(str).to_numpy()
    if set(times) != set(cfg["k"]):
        raise ValueError("Unexpected timepoints")
    args.output.mkdir(parents=True)
    summary = {
        "dataset": args.dataset,
        "method": "SEDR native count-HVG2000/PCA200/no-DEC/Mclust EEE",
        "source_sha256": cfg["sha256"],
        "script_sha256": digest(Path(__file__)),
        "r_script_sha256": digest(args.r_script),
        "package_hashes": package_hashes,
        "gpu": args.gpu,
        "slices": {},
    }
    save_json(args.output / "progress.json", summary)
    for tp, k in cfg["k"].items():
        started = time.monotonic()
        slice_data = data[times == tp].copy()
        if not sp.issparse(slice_data.X):
            slice_data.X = sp.csr_matrix(np.asarray(slice_data.X))
        if not np.array_equal(
            np.asarray(slice_data.obsm["spatial"]),
            slice_data.obs[["pixel_x", "pixel_y"]].to_numpy(),
        ):
            raise ValueError("Coordinate mismatch")
        graph = SEDR.graph_construction(slice_data, 12)
        slice_data.layers["count"] = slice_data.X.copy()
        sc.pp.normalize_total(slice_data, target_sum=1e6)
        sc.pp.highly_variable_genes(
            slice_data, flavor="seurat_v3", layer="count", n_top_genes=2000
        )
        slice_data = slice_data[:, slice_data.var["highly_variable"]].copy()
        genes = np.asarray(slice_data.var_names, dtype=str)
        if len(genes) != 2000:
            raise ValueError("Native SEDR selected an unexpected gene count")
        sc.pp.scale(slice_data)
        values = slice_data.X.toarray() if sp.issparse(slice_data.X) else np.asarray(slice_data.X)
        pca = PCA(n_components=200, random_state=42).fit_transform(values)
        if not np.isfinite(pca).all():
            raise ValueError("Invalid SEDR PCA")
        sedr = SEDR.Sedr(pca, graph, mode="clustering", device=device)
        sedr.train_without_dec()
        latent, _, _, _ = sedr.process()
        latent = np.asarray(latent)
        if latent.shape[0] != slice_data.n_obs or not np.isfinite(latent).all():
            raise ValueError("Invalid SEDR embedding")
        r_input = args.output / f"{tp}_sedr.csv"
        r_output = args.output / f"{tp}_mclust.csv"
        np.savetxt(r_input, latent, delimiter=",", fmt="%.17g")
        subprocess.run(
            ["Rscript", str(args.r_script), str(r_input), str(k), str(r_output)],
            check=True,
        )
        pred = np.loadtxt(r_output, delimiter=",", dtype=np.int16).reshape(-1)
        if len(pred) != slice_data.n_obs or len(np.unique(pred)) != k:
            raise ValueError("Invalid SEDR Mclust labels")
        np.savez_compressed(
            args.output / f"{tp}_result.npz",
            obs_names=np.asarray(slice_data.obs_names, dtype=str),
            predicted=pred,
            latent=latent,
            selected_genes=genes,
        )
        torch.save(sedr.model.state_dict(), args.output / f"{tp}_model.pt")
        annotation = slice_data.obs["annotation"].astype(str).to_numpy()
        summary["slices"][tp] = {
            "n_spots": slice_data.n_obs,
            "k": k,
            "n_selected_genes": len(genes),
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
