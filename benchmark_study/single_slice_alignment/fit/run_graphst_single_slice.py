"""Independent native GraphST fits on original raw counts.

The imported GraphST package is an unchanged source snapshot. Mclust is run
in R with the same EEE/2020 settings as GraphST.utils.mclust_R.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
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
        raise ValueError("Original raw input SHA256 mismatch")
    if args.output.exists():
        raise FileExistsError(args.output)
    package = args.package_root / "GraphST"
    required = ["GraphST.py", "__init__.py", "model.py", "preprocess.py", "utils.py"]
    hashes = {name: digest(package / name) for name in required}
    if len(hashes) != 5 or not args.r_script.exists():
        raise ValueError("Missing native package or Mclust adapter")
    sys.path.insert(0, str(args.package_root))
    from GraphST.GraphST import GraphST  # noqa: PLC0415

    if not torch.cuda.is_available() or args.gpu >= torch.cuda.device_count():
        raise RuntimeError("Requested CUDA GPU is unavailable")
    device = torch.device(f"cuda:{args.gpu}")
    assert float(torch.ones(1, device=device).item()) == 1.0
    torch.set_num_threads(1)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    data = ad.read_h5ad(source)
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError("Duplicate spot or gene identifiers")
    if "counts" in data.layers:
        data.X = data.layers["counts"].copy()
    times = data.obs["timepoint"].astype(str).to_numpy()
    if set(times) != set(cfg["k"]):
        raise ValueError("Unexpected timepoints")
    if not np.array_equal(
        np.asarray(data.obsm["spatial"]),
        data.obs[["pixel_x", "pixel_y"]].to_numpy(),
    ):
        raise ValueError("Spatial coordinates differ")
    args.output.mkdir(parents=True, exist_ok=False)
    summary = {
        "dataset": args.dataset,
        "method": "GraphST native seurat_v3 HVG3000 + PCA20 + Mclust EEE",
        "source": str(source),
        "source_sha256": cfg["sha256"],
        "script_sha256": digest(Path(__file__)),
        "r_script_sha256": digest(args.r_script),
        "package_hashes": hashes,
        "gpu": args.gpu,
        "seed": 2026,
        "slices": {},
    }
    save_json(args.output / "progress.json", summary)
    for tp, k in cfg["k"].items():
        start = time.monotonic()
        slice_data = data[times == tp].copy()
        if "highly_variable" in slice_data.var:
            del slice_data.var["highly_variable"]
        if not sp.issparse(slice_data.X):
            slice_data.X = sp.csr_matrix(np.asarray(slice_data.X))
        model = GraphST(
            slice_data,
            device=device,
            epochs=600,
            random_seed=2026,
            datatype=cfg["datatype"],
        )
        trained = model.train()
        reconstruction = np.asarray(trained.obsm["emb"], dtype=np.float32)
        if reconstruction.shape != (slice_data.n_obs, 3000):
            raise ValueError(f"Unexpected GraphST reconstruction at {tp}")
        latent = PCA(n_components=20, random_state=42).fit_transform(reconstruction)
        if not np.isfinite(latent).all():
            raise ValueError("Nonfinite GraphST latent")
        r_input = args.output / f"{tp}_emb_pca.csv"
        r_output = args.output / f"{tp}_mclust.csv"
        np.savetxt(r_input, latent, delimiter=",", fmt="%.17g")
        subprocess.run(
            ["Rscript", str(args.r_script), str(r_input), str(k), str(r_output)],
            check=True,
        )
        pred = np.loadtxt(r_output, delimiter=",", dtype=np.int16).reshape(-1)
        if len(pred) != slice_data.n_obs or len(np.unique(pred)) != k:
            raise ValueError(f"Invalid Mclust result at {tp}")
        genes = np.asarray(trained.var_names[trained.var["highly_variable"]], dtype=str)
        if len(genes) != 3000:
            raise ValueError("Native GraphST did not select 3000 HVGs")
        np.savez_compressed(
            args.output / f"{tp}_result.npz",
            obs_names=np.asarray(slice_data.obs_names, dtype=str),
            predicted=pred,
            latent=latent,
            emb=reconstruction,
            selected_genes=genes,
        )
        torch.save(model.model.state_dict(), args.output / f"{tp}_model.pt")
        annotation = slice_data.obs["annotation"].astype(str).to_numpy()
        summary["slices"][tp] = {
            "n_spots": slice_data.n_obs,
            "k": k,
            "n_selected_genes": len(genes),
            "ari": float(adjusted_rand_score(annotation, pred)),
            "nmi": float(normalized_mutual_info_score(annotation, pred)),
            "model_sha256": digest(args.output / f"{tp}_model.pt"),
            "result_sha256": digest(args.output / f"{tp}_result.npz"),
            "seconds": time.monotonic() - start,
        }
        save_json(args.output / "progress.json", summary)
        print(f"COMPLETED {args.dataset} {tp}: {summary['slices'][tp]}", flush=True)
    save_json(args.output / "summary.json", summary)


if __name__ == "__main__":
    main()
