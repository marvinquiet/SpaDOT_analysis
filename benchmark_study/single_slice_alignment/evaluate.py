"""Compute stage ARI, annotation-based naïve transition accuracy and raw ECDFs.

Supply the raw H5AD (--raw) and either a native single-slice fit directory
(--fit) or an aligned saved H5AD/NPZ embedding (--embedding). Saved H5ADs
must have spot IDs in obs_names and predicted labels in --prediction-key.
An archived simulation NPZ has x, spot and pred vectors. The benchmark's
real and simulation reference graphs are built in. To provide a different
real reference graph, pass --edges CSV with transition (e.g. t0_t1), source
and target columns.
Set --joint-reference to compute distance ECDFs on the intersection of the
current and joint spot sets. Saved H5ADs with historical IDs can be matched
using the original barcode and timepoint fields. An archived simulation NPZ
can also provide its timepoint|pixel_ind keys. --output is a new JSON result.

GraphST uses its 3,000 reconstructed columns by position for OT and distances;
its ARI comes from the native PCA20/Mclust labels. OT uses the unmodified
embeddings and the same Waddington-OT settings for every method. The S44
simulation-only rescaled SpaDOT ARI is calculated separately.
"""

import argparse
import json
import os
import sys
import warnings
from pathlib import Path

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[variable] = "1"

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.spatial.distance import cdist
from sklearn.metrics import adjusted_rand_score


SIMULATION_EDGES = {
    "t0_t1": {("1", "2"), ("1", "3"), ("2", "1"), ("3", "4")},
    "t1_t2": {("1", "1"), ("2", "3"), ("3", "4"), ("4", "2"), ("4", "5")},
}
OT_PARAMS = dict(epsilon=0.01, epsilon0=1, lambda1=1, lambda2=50, growth_iters=3, tau=1000)


def dense(x):
    return x.toarray() if issparse(x) else np.asarray(x)


def canonical(value):
    value = str(value)
    return str(int(float(value))) if value.replace(".", "", 1).isdigit() else value


def canonical_ids(original, obs, reference, dataset):
    if original.isin(reference.index).all():
        if "timepoint" in obs and not np.array_equal(
                obs.timepoint.astype(str).to_numpy(), reference.loc[original, "timepoint"].astype(str).to_numpy()):
            raise ValueError("Saved and reference stage assignments disagree")
        return original
    stage_name = {"t0": "Stage44", "t1": "Stage54", "t2": "Stage57"}
    known = set(reference.index.astype(str))
    mapped = []
    for key, row in obs.iterrows():
        stage = str(row["timepoint"])
        candidates = {str(key)}
        if "barcodes" in obs:
            candidates.add(str(row["barcodes"]))
        if {"barcode_ind", "barcode"} <= set(obs):
            candidates.add(f"{row['barcode_ind']}_{row['barcode']}")
        if dataset == "AxolotlTelencephalon" and "CellID" in obs:
            candidates.add(f"{row['CellID']}_{stage_name[stage]}_{stage}")
        value = str(key)
        while "-" in value and value.rsplit("-", 1)[1].isdigit():
            value = value.rsplit("-", 1)[0]
            candidates.add(value)
        for value in list(candidates):
            if value.startswith(stage + "_"):
                candidates.add(value[len(stage) + 1:])
            if dataset == "AxolotlTelencephalon":
                candidates.add(value + "_" + stage)
        matches = [item for item in candidates & known if str(reference.at[item, "timepoint"]) == stage]
        if len(matches) != 1:
            raise ValueError(f"Expected one reference barcode for {key} at {stage}, found {matches}")
        mapped.append(matches[0])
    return pd.Index(mapped)


def load_h5ad_or_npz(path, latent_key, prediction_key, reference, dataset):
    if path.suffix == ".npz":
        # Archived benchmark NPZs contain object arrays. Use trusted outputs only.
        with np.load(path, allow_pickle=True) as data:
            ids = pd.Index(data["spot"].astype(str))
            if not ids.isin(reference.index).all() and "keys" in data:
                keys = [f"{t}|{int(i)}" for t, i in zip(reference.timepoint, reference.pixel_ind)]
                key_index = pd.Series(reference.index, index=keys)
                if not key_index.index.is_unique:
                    raise ValueError("Duplicate timepoint|pixel_ind keys in the reference")
                ids = pd.Index(key_index.loc[data["keys"].astype(str)])
            return ids, data["x"], data["pred"].astype(str)
    data = ad.read_h5ad(path)
    values = data.X if latent_key == "X" else data.obsm[latent_key]
    return (canonical_ids(data.obs_names, data.obs, reference, dataset), dense(values),
            data.obs[prediction_key].astype(str).to_numpy())


def load_fit(folder, method, stages):
    blocks = []
    for stage in stages:
        if method in {"PRECAST", "spatialPCA"}:
            latent = folder / stage / "latent.csv"
            predicted = folder / stage / "predicted.csv"
            if not latent.is_file() or not predicted.is_file():
                continue
            frame = pd.read_csv(latent, dtype={"spot": str}).set_index("spot")
            labels = pd.read_csv(predicted, dtype={"spot": str}).set_index("spot")
            blocks.append((frame.index.to_numpy(), frame.to_numpy(), labels.loc[frame.index, "predicted"].to_numpy()))
        else:
            path = (folder / f"{stage}_result.npz" if method in {"GraphST", "SEDR", "spaGCN"}
                    else folder / stage / "result.npz")
            if path.is_file():
                with np.load(path, allow_pickle=False) as data:
                    latent = data["emb"] if method == "GraphST" else data["latent"]
                    blocks.append((data["obs_names"].astype(str), latent, data["predicted"]))
    if not blocks:
        raise ValueError(f"No saved stage output in {folder}")
    return tuple(np.concatenate(items) for items in zip(*blocks))


def aligned(raw_obs, embedding):
    ids, x, pred = embedding
    ids = pd.Index(ids, dtype=str)
    if (not ids.is_unique or len(ids) != len(x) or len(ids) != len(pred)
            or not np.isfinite(x).all()):
        raise ValueError("Invalid embedding, labels or duplicate identifiers")
    obs = raw_obs.loc[ids].copy()
    obs["truth"] = obs["annotation"].map(canonical).astype(str)
    obs["pred"] = pd.Series(pred, index=ids).map(canonical).astype(str)
    obs["timepoint"] = obs["timepoint"].astype(str)
    return obs, x


def mass_by_domain(matrix, source, target):
    left, i = np.unique(source, return_inverse=True)
    right, j = np.unique(target, return_inverse=True)
    rows = np.eye(len(left))[i]
    cols = np.eye(len(right))[j]
    mass = rows.T @ (matrix @ cols)
    if not np.isclose(mass.sum(), matrix.sum(), rtol=1e-12):
        raise ValueError("Aggregated transport mass does not match spot-level mass")
    return left, right, mass


def association_accuracy(mass, rows, columns, edges):
    rsum = mass.sum(axis=1, keepdims=True)
    csum = mass.sum(axis=0, keepdims=True)
    forward = np.divide(mass, rsum, out=np.zeros_like(mass), where=rsum > 0)
    reverse = np.divide(mass, csum, out=np.zeros_like(mass), where=csum > 0)
    q = np.minimum(forward, reverse)
    correct = np.array([[(a, b) in edges for b in columns] for a in rows])
    return float(q[correct].sum() / q.sum()), q


def score(matrix, source, target, edges):
    sy, ty = source["truth"].to_numpy(), target["truth"].to_numpy()
    sp, tp = source["pred"].to_numpy(), target["pred"].to_numpy()
    rows, columns, mass = mass_by_domain(matrix, sy, ty)
    domain_rows, domain_columns, domain_mass = mass_by_domain(matrix, sp, tp)
    fs = np.array([[np.mean(sy[sp == p] == label) for label in rows] for p in domain_rows])
    ft = np.array([[np.mean(ty[tp == p] == label) for label in columns] for p in domain_columns])
    annotation, q = association_accuracy(mass, rows, columns, edges)
    predicted, _ = association_accuracy(domain_mass, rows[fs.argmax(1)], columns[ft.argmax(1)], edges)
    projected, _ = association_accuracy(fs.T @ domain_mass @ ft, rows, columns, edges)
    return {
        "annotation_naive_accuracy": annotation,
        "predicted_domain_naive_accuracy": predicted,
        "predicted_domain_projected_naive_accuracy": projected,
    }, {"rows": rows.tolist(), "columns": columns.tolist(), "mass": mass.tolist(), "score": q.tolist()}


def ecdf(x, y, simulation):
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    values = np.concatenate([cdist(x[i:i + 128], y).ravel() for i in range(0, len(x), 128)])
    values.sort()
    n = len(values)
    if simulation:
        ix = np.unique(np.rint(np.linspace(0, n - 1, 1201)).astype(int))
        return {"distance": values[ix].tolist(), "cdf": ((ix + 1) / n).tolist(), "n_pairs": n}
    ix = np.unique(np.concatenate([np.linspace(0, n - 1, 400).astype(int),
                                   (np.geomspace(1, n, 180).astype(int) - 1),
                                   n - np.geomspace(1, n, 180).astype(int)]))
    points = np.unique(values[ix])
    return {"distance": points.tolist(), "cdf": (np.searchsorted(values, points, side="right") / n).tolist(),
            "n_pairs": n}


def evaluate(args):
    import wot

    raw = ad.read_h5ad(args.raw, backed="r")
    raw_obs = raw.obs.copy()
    raw.file.close()
    stages = sorted(raw_obs["timepoint"].astype(str).unique())
    embedding = (load_fit(args.fit, args.method, stages) if args.fit else
                 load_h5ad_or_npz(args.embedding, args.latent_key, args.prediction_key, raw_obs, args.case))
    obs, x = aligned(raw_obs, embedding)
    reference = (aligned(raw_obs, load_h5ad_or_npz(args.joint_reference, args.joint_latent_key,
                                                    args.joint_prediction_key, raw_obs, args.case))
                 if args.joint_reference else None)
    if reference is not None and not args.simulation and x.shape[1] != reference[1].shape[1]:
        raise ValueError("Real-data distance comparisons require the same dimensions in both modes")
    if args.simulation:
        edges = SIMULATION_EDGES
    elif args.edges is not None:
        edge_frame = pd.read_csv(args.edges, dtype=str)
        edges = {pair: set(zip(block.source, block.target)) for pair, block in edge_frame.groupby("transition")}
    else:
        if args.case not in {"MouseEmbryo", "ChickenHeart", "AxolotlTelencephalon"}:
            raise ValueError("Provide --edges for a case outside the three real benchmark datasets")
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "revision"))
        from transition_analysis.reference_edges import edges as reference_edges

        regions = set(raw_obs.annotation.map(canonical))
        edges = {f"{left}_{right}": reference_edges(args.case, i, regions)
                 for i, (left, right) in enumerate(zip(stages, stages[1:]))}
    record = {"case": args.case, "method": args.method, "mode": args.mode, "ari": {},
              "transitions": {}, "distances": {}, "params": OT_PARAMS, "local_pca_requested": 30}
    for stage in stages:
        subset = obs.loc[obs.timepoint.eq(stage)]
        if not subset.empty:
            record["ari"][stage] = {"ari": float(adjusted_rand_score(subset.truth, subset.pred)),
                                    "spots": len(subset), "clusters": subset.pred.nunique()}
    days = obs.timepoint.map({stage: i for i, stage in enumerate(stages)})
    ot = ad.AnnData(x.copy(), obs=pd.DataFrame({"day": days, "cell_growth_rate": 1.0}, index=obs.index))
    for i, (left, right) in enumerate(zip(stages, stages[1:])):
        name = f"{left}_{right}"
        if name not in edges:
            raise ValueError(f"Missing reference graph for {name}")
        a = obs.timepoint.eq(left).to_numpy()
        b = obs.timepoint.eq(right).to_numpy()
        if not a.any() or not b.any():
            record["transitions"][name] = {"status": "missing", "error": "Missing fitted stage"}
            continue
        if reference is not None:
            ref_obs, ref_x = reference
            source = obs.index[a].intersection(ref_obs.index[ref_obs.timepoint.eq(left)], sort=False)
            target = obs.index[b].intersection(ref_obs.index[ref_obs.timepoint.eq(right)], sort=False)
            if len(source) and len(target):
                record["distances"][name] = {
                    "single": ecdf(x[obs.index.get_indexer(source)], x[obs.index.get_indexer(target)], args.simulation),
                    "multi": ecdf(ref_x[ref_obs.index.get_indexer(source)], ref_x[ref_obs.index.get_indexer(target)], args.simulation),
                    "matched_source": len(source), "matched_target": len(target),
                }
        try:
            with warnings.catch_warnings(record=True) as messages:
                warnings.simplefilter("always")
                if args.simulation:
                    pair = ot[a | b].copy()
                    pair.obs["day"] = pair.obs["day"] - i
                    model = wot.ot.OTModel(pair, local_pca=30, **OT_PARAMS)
                    transport = model.compute_transport_map(0, 1)
                else:
                    model = wot.ot.OTModel(ot, local_pca=30, **OT_PARAMS)
                    transport = model.compute_transport_map(i, i + 1)
            ri = transport.obs_names.get_indexer(obs.index[a])
            ci = transport.var_names.get_indexer(obs.index[b])
            if (ri < 0).any() or (ci < 0).any():
                raise ValueError("Transport IDs do not cover the two stages")
            matrix = dense(transport.X)[np.ix_(ri, ci)]
            if not np.isfinite(matrix).all() or (matrix < 0).any() or matrix.sum() <= 0:
                raise ValueError("Invalid transport mass")
            scores, detail = score(matrix, obs.loc[a], obs.loc[b], edges[name])
            record["transitions"][name] = {"status": "complete", "scores": scores, "details": detail,
                                           "warnings": [str(item.message) for item in messages]}
        except (ValueError, FloatingPointError, OverflowError, RuntimeError) as error:
            record["transitions"][name] = {"status": "failed", "error": str(error)}
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--mode", choices=("single", "multi"), required=True)
    parser.add_argument("--simulation", action="store_true")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--fit", type=Path)
    source.add_argument("--embedding", type=Path)
    parser.add_argument("--latent-key", default="X")
    parser.add_argument("--prediction-key", default="kmeans")
    parser.add_argument("--edges", type=Path)
    parser.add_argument("--joint-reference", type=Path)
    parser.add_argument("--joint-latent-key", default="X")
    parser.add_argument("--joint-prediction-key", default="kmeans")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, allow_nan=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
