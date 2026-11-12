#!/usr/bin/env python3
"""Fit STIFT using the pinned upstream ARTISTA tutorial.

Clone https://github.com/cuhklinlab/STIFT and check out
9d649a018ef0262e1e71f34b9d77439aaf9d0602. Use the original
tutorial_ARTISTA_integration.ipynb from that checkout as --tutorial.
Run with one visible CUDA device and install the upstream dependencies.
The manuscript runs used Scanpy 1.10.3, AnnData 0.10.8, NumPy 1.26.0,
SciPy 1.12.0, PyTorch 2.5.0+cu124 and torch-geometric 2.5.3.

For example:
    python fit_stift.py --stift-source /path/to/STIFT \
      --tutorial /path/to/tutorial_ARTISTA_integration.ipynb \
      --input /path/to/raw.h5ad --dataset MouseEmbryo --output /path/to/fit

Run once per dataset, or once per supplied 6,000-spot simulation case.
Strategy 1 inputs have integer X and obs.x_coor/obs.y_coor. Strategy 2
inputs have sparse integer X and obsm.spatial. Both use obs.timepoint and
obs.label or obs.annotation. Real inputs use obs.timepoint, obs.annotation
and obsm.spatial. Outputs contain the embedding, actual training couplings
and their spot order. The training input contains expression, space and
time. Ground truth is used for plotting and evaluation after the fit.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time as clock

os.environ.setdefault("MPLBACKEND", "Agg")

COMMIT = "9d649a018ef0262e1e71f34b9d77439aaf9d0602"
NOTEBOOK_SHA256 = "81c8289cd2200a322a42e2f9179ee2b7df34b5e3b8555be2588193b699393c2b"
REAL_SHAPES = {"MouseEmbryo": (8255, 22496), "ChickenHeart": (6596, 23986),
               "Axolotl": (8816, 12704)}
SIMULATIONS = ("strategy1", "strategy2")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_tutorial(source, notebook):
    head = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    modifications = subprocess.check_output(
        ["git", "-C", str(source), "status", "--porcelain", "--untracked-files=no", "--", "STIFT"],
        text=True,
    )
    if head != COMMIT or modifications:
        raise ValueError("Use the unmodified pinned upstream STIFT checkout.")
    if sha256(notebook) != NOTEBOOK_SHA256:
        raise ValueError("The ARTISTA notebook does not match the manuscript source.")
    sys.path.insert(0, str(source.resolve()))
    return json.loads(notebook.read_text())["cells"]


def fit(args):
    cells = load_tutorial(args.stift_source, args.tutorial)
    namespace = {"__name__": "__tutorial__"}

    phase_seconds = {}

    def execute(number):
        source = cells[number]["source"]
        if isinstance(source, list):
            source = "".join(source)
        print(f"Running upstream ARTISTA tutorial cell {number}", flush=True)
        start = clock.perf_counter()
        exec(compile(source, f"ARTISTA:cell_{number}", "exec"), namespace)
        if number != 3:
            namespace["torch"].cuda.synchronize()
        phase_seconds[number] = clock.perf_counter() - start

    execute(3)
    import STIFT.DESTOT as destot

    ad, np, pd, sp, torch = (namespace[key] for key in ("ad", "np", "pd", "sp", "torch"))
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise RuntimeError("Expose exactly one CUDA device to this process.")
    torch.set_num_threads(1)
    destot.get_free_gpu = lambda: 0

    raw = ad.read_h5ad(args.input)
    simulated = args.dataset in SIMULATIONS
    expected_shape = (6000, 5000) if simulated else REAL_SHAPES[args.dataset]
    if tuple(raw.shape) != expected_shape or not raw.obs_names.is_unique or not raw.var_names.is_unique:
        raise ValueError(f"Expected {expected_shape} and unique spot and gene identifiers.")
    coordinates = (raw.obs[["x_coor", "y_coor"]].to_numpy() if args.dataset == "strategy1"
                   else np.asarray(raw.obsm["spatial"]))
    if coordinates.shape != (raw.n_obs, 2) or not np.isfinite(coordinates).all():
        raise ValueError("Expected two finite spatial coordinates per spot.")
    counts = raw.layers.get("counts", raw.layers.get("count", raw.X))
    counts = sp.csr_matrix(counts)
    if (not np.isfinite(counts.data).all() or (counts.data < 0).any()
            or not np.allclose(counts.data, np.rint(counts.data), rtol=0, atol=1e-6)):
        raise ValueError("STIFT input must contain finite, nonnegative raw counts.")
    times = raw.obs["timepoint"].astype(str).to_numpy()
    sections = sorted(set(times))
    if simulated and sections != ["t0", "t1", "t2"]:
        raise ValueError("Expected t0, t1 and t2 in the simulation input.")
    label = "annotation" if "annotation" in raw.obs else "label"
    order = []
    parts = []
    for time in sections:
        indices = np.flatnonzero(times == time)
        order.extend(indices)
        original = raw.obs_names[indices]
        names = [f"{time}_{spot}" for spot in original]
        if simulated:
            obs = pd.DataFrame({"timepoint": time, "batch_names": time}, index=names)
            part = ad.AnnData(X=counts[indices], obs=obs, var=raw.var.copy())
        else:
            part = raw[indices].copy()
            part.obs_names = names
            part.obs_names.name = "unique_obs"
            part.obs["batch_names"] = time
            part.X = counts[indices]
            if "counts" in part.layers:
                del part.layers["counts"]
        part.layers["count"] = part.X.copy()
        part.obsm["spatial"] = coordinates[indices].copy()
        parts.append(part)

    args.output.mkdir(parents=True, exist_ok=False)
    namespace.update(adata_list=parts, section_ids=sections)
    execute(9)
    for i, (left, right) in enumerate(zip(sections[:-1], sections[1:])):
        coupling = np.asarray(namespace["Pi_list"][i])
        if (coupling.shape != (parts[i].n_obs, parts[i + 1].n_obs)
                or not np.isfinite(coupling).all() or (coupling < 0).any()
                or coupling.sum() <= 0):
            raise ValueError(f"Invalid native STIFT coupling for {left}->{right}.")
        np.save(args.output / f"Pi_{left}_{right}.npy", coupling, allow_pickle=False)
        np.save(args.output / f"xi_{left}_{right}.npy", namespace["xi_list"][i], allow_pickle=False)
    execute(11)
    selected = namespace["adata_list"]
    combined = namespace["adata_concat"]
    expected_genes = 5000 if simulated else 12000
    if any(part.n_vars != expected_genes for part in selected):
        raise ValueError("Unexpected per-timepoint feature count.")
    expected_order = [spot for part in parts for spot in part.obs_names]
    if list(combined.obs_names) != expected_order:
        raise ValueError("Tutorial reordered the input spots.")
    (args.output / "shared_genes_in_training_order.txt").write_text(
        "\n".join(combined.var_names) + "\n"
    )
    execute(13)
    combined = namespace["adata_concat"]
    embedding = np.asarray(combined.obsm["STIFT"])
    if embedding.shape != (raw.n_obs, 30) or not np.isfinite(embedding).all():
        raise ValueError("Expected finite 30-dimensional STIFT embeddings.")
    obs = raw.obs.iloc[order].copy()
    obs["annotation"] = obs[label].astype(str)
    obs["original_spot"] = obs.index.astype(str)
    obs.index = combined.obs_names.copy()
    obs["batch_names"] = obs["timepoint"].astype(str)
    obs[["pixel_x", "pixel_y"]] = coordinates[order]
    if not simulated:
        execute(16)
        namespace["fig"].savefig(args.output / "tutorial_UMAP.pdf", bbox_inches="tight")
        namespace["plt"].close("all")
        execute(18)
        obs["louvain"] = combined.obs["louvain"].astype(str)
    latent = ad.AnnData(X=embedding.copy(), obs=obs)
    latent.obsm["STIFT"] = embedding.copy()
    latent.obsm["spatial"] = coordinates[order].copy()
    latent.write_h5ad(args.output / "STIFT_latent.h5ad", compression="gzip")
    obs[["original_spot", "timepoint"]].to_csv(args.output / "spot_order.csv", index_label="spot")
    seconds = sum(phase_seconds[number] for number in (9, 11, 13))
    (args.output / "STIFT_running_time.txt").write_text(f"{seconds:.9f}\n")



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stift-source", type=Path, required=True)
    parser.add_argument("--tutorial", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--dataset", choices=(*REAL_SHAPES, *SIMULATIONS), required=True)
    parser.add_argument("--output", type=Path, required=True)
    fit(parser.parse_args())


if __name__ == "__main__":
    main()
