"""Lossless raw-count bridge for native R single-slice comparators."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite, mmread

from inputs import digest, load_case, save_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--manifest", type=Path, required=True, help="JSON produced by prepare_manifest.py")
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cfg = load_case(args.manifest, args.dataset)
    source = args.raw_root / cfg["relative_input"]
    if digest(source) != cfg["sha256"]:
        raise ValueError("Original source hash mismatch")
    if args.output.exists():
        raise FileExistsError(args.output)
    data = ad.read_h5ad(source)
    if "counts" in data.layers:
        data.X = data.layers["counts"].copy()
    times = data.obs["timepoint"].astype(str).to_numpy()
    if set(times) != set(cfg["k"]):
        raise ValueError("Timepoints differ from registered protocol")
    args.output.mkdir()
    report = {"dataset": args.dataset, "source_sha256": cfg["sha256"], "slices": {}}
    for tp, k in cfg["k"].items():
        part = data[times == tp].copy()
        if not np.array_equal(part.obsm["spatial"], part.obs[["pixel_x", "pixel_y"]].to_numpy()):
            raise ValueError("Coordinate mismatch")
        counts = sp.csr_matrix(part.X)
        if np.any(counts.data < 0) or not np.isfinite(counts.data).all() or np.any(counts.data != np.floor(counts.data)):
            raise ValueError("Invalid raw counts")
        block = args.output / tp
        block.mkdir()
        mtx = block / "counts_genes_by_spots.mtx"
        mmwrite(mtx, counts.T.tocoo(), field="integer")
        if (mmread(mtx).T.tocsr() != counts).nnz:
            raise ValueError("Raw count bridge round-trip mismatch")
        (block / "genes.txt").write_text("\n".join(map(str, part.var_names)) + "\n")
        (block / "spots.txt").write_text("\n".join(map(str, part.obs_names)) + "\n")
        pd.DataFrame({"spot": np.asarray(part.obs_names, dtype=str),
                      "pixel_x": part.obsm["spatial"][:, 0],
                      "pixel_y": part.obsm["spatial"][:, 1]}).to_csv(block / "coordinates.csv", index=False)
        pd.DataFrame({"spot": np.asarray(part.obs_names, dtype=str),
                      "annotation": part.obs["annotation"].astype(str).to_numpy()}).to_csv(
                          block / "posthoc_truth.csv", index=False)
        report["slices"][tp] = {
            "spots": part.n_obs, "genes": part.n_vars, "nonzero": counts.nnz,
            "k": k, "sum_counts": int(counts.sum()),
            "files_sha256": {p.name: digest(p) for p in sorted(block.iterdir())},
        }
        save_json(args.output / "progress.json", report)
    save_json(args.output / "manifest.json", report)
    print(f"EXPORTED {args.dataset}: {report['slices']}")


if __name__ == "__main__":
    main()
