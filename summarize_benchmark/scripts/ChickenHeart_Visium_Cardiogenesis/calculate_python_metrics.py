"""Recompute Chicken Heart domain metrics with the original benchmark formulas."""

import argparse
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.special import xlogy
from scipy.stats import hypergeom
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score, homogeneity_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

METHODS = (
    "SpaDOT", "spatiAlign", "spatialPCA", "spaVAE", "SPIRAL", "SEDR",
    "STADIA", "PRECAST", "spaGCN", "STAligner", "GraphST", "STAMP", "STIFT",
)
SPOT_COUNTS = {"t0": 747, "t1": 1966, "t2": 1916, "t3": 1967}
COORDINATES = ["pixel_x", "pixel_y"]
LABELS = ["true_labels", *METHODS]
K = 10


def original_continuity_functions(path):
    source = path.read_text()
    names = {"_compute_PAS", "_compute_CHAOS"}
    nodes = [node for node in ast.parse(source).body
             if isinstance(node, ast.FunctionDef) and node.name in names]
    assert len(nodes) == 2 and {node.name for node in nodes} == names
    scope = {"np": np, "NearestNeighbors": NearestNeighbors, "StandardScaler": StandardScaler}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), scope)
    return scope["_compute_PAS"], scope["_compute_CHAOS"]


def annotation_metrics_independent(truth, predicted):
    """Check the scores from the contingency counts and hypergeometric null."""
    _, row_ids = np.unique(truth, return_inverse=True)
    _, col_ids = np.unique(predicted, return_inverse=True)
    table = np.zeros((row_ids.max() + 1, col_ids.max() + 1), dtype=np.int64)
    np.add.at(table, (row_ids, col_ids), 1)
    n = int(table.sum())
    rows, cols = table.sum(1), table.sum(0)
    pair_count = lambda counts: float((counts * (counts - 1) // 2).sum())
    row_pairs, col_pairs = pair_count(rows), pair_count(cols)
    expected_pairs = row_pairs * col_pairs / (n * (n - 1) / 2)
    ari = (pair_count(table) - expected_pairs) / ((row_pairs + col_pairs) / 2 - expected_pairs)

    p_true, p_pred = rows / n, cols / n
    entropy_true = -xlogy(p_true, p_true).sum()
    entropy_pred = -xlogy(p_pred, p_pred).sum()
    expected_counts = np.outer(rows, cols) / n
    mi = (xlogy(table, table / expected_counts) / n).sum()
    emi = 0.0
    for a in rows:
        for b in cols:
            overlap = np.arange(max(1, a + b - n), min(a, b) + 1)
            contribution = overlap / n * np.log(overlap * n / (a * b))
            emi += np.dot(hypergeom.pmf(overlap, n, a, b), contribution)
    ami = (mi - emi) / ((entropy_true + entropy_pred) / 2 - emi)
    return np.array([ari, ami, mi / entropy_true])


def independent_neighborhoods(coordinates):
    """Return nearer spots and any candidates tied at the tenth-neighbor boundary."""
    tree = cKDTree(coordinates)
    distances, indices = tree.query(coordinates, k=K + 2)
    assert np.array_equal(indices[:, 0], np.arange(len(coordinates)))
    nearest = indices[:, 1:K + 1]
    tied_rows = np.flatnonzero(np.isclose(distances[:, K], distances[:, K + 1],
                                         rtol=1e-12, atol=1e-12))
    ties = {}
    for i in tied_rows:
        radius = distances[i, K]
        candidates = np.array(tree.query_ball_point(coordinates[i], radius * (1 + 1e-12)))
        candidates = candidates[candidates != i]
        candidate_distances = np.linalg.norm(coordinates[candidates] - coordinates[i], axis=1)
        boundary = np.isclose(candidate_distances, radius, rtol=1e-12, atol=1e-12)
        strictly_nearer = candidates[~boundary]
        tied = candidates[boundary]
        assert (candidate_distances[~boundary] < radius).all()
        assert len(strictly_nearer) < K < len(strictly_nearer) + len(tied)
        ties[i] = strictly_nearer, tied
    return nearest, ties


def pas_bounds(labels, nearest, ties):
    """Tie handling can vary across trees, so check the admissible exact PAS range."""
    mismatches = (labels[nearest] != labels[:, None]).sum(1)
    lower = mismatches.copy()
    upper = mismatches.copy()
    for i, (near, tied) in ties.items():
        need = K - len(near)
        different = int((labels[tied] != labels[i]).sum())
        already = int((labels[near] != labels[i]).sum())
        lower[i] = already + max(0, need - (len(tied) - different))
        upper[i] = already + min(need, different)
    return (lower > K / 2).mean(), (upper > K / 2).mean()


def chaos_independent(standardized, labels):
    total, count = 0.0, 0
    for label in np.unique(labels):
        points = standardized[labels == label]
        if len(points) > 2:
            distances = cKDTree(points).query(points, k=2)[0][:, 1]
            total += distances.sum()
            count += len(points)
    assert count > 0
    return total / count


def read_labels(path):
    data = pd.read_csv(path, dtype=str, keep_default_na=False)
    expected = {"spot", "timepoint", *COORDINATES, *LABELS}
    assert expected <= set(data), f"Missing columns: {sorted(expected - set(data))}"
    assert data.spot.is_unique and data.spot.str.strip().ne("").all()
    assert data.groupby("timepoint").size().to_dict() == SPOT_COUNTS
    for label in LABELS:
        assert not data[label].str.strip().str.lower().isin(
            {"", "na", "nan", "none", "null", "<na>"}).any(), f"Missing label in {label}"
    data[COORDINATES] = data[COORDINATES].apply(pd.to_numeric, errors="raise")
    assert np.isfinite(data[COORDINATES].to_numpy()).all()
    assert not data.duplicated(["timepoint", *COORDINATES]).any()
    return data


def read_native_order(path, data):
    order = pd.read_csv(path, dtype=str, keep_default_na=False).set_index("spot")
    assert order.index.is_unique and set(order.index) == set(data.spot)
    assert set(METHODS) <= set(order.columns)
    expected = np.arange(len(data))
    order = order[list(METHODS)].copy()
    for method in METHODS:
        assert order[method].str.fullmatch(r"[0-9]+").all(), method
        order[method] = pd.to_numeric(order[method], errors="raise")
        np.testing.assert_array_equal(np.sort(order[method].to_numpy()), expected, err_msg=method)
    return order


def compute(data, compute_pas, compute_chaos, native_order=None):
    rows, audits = [], []
    for stage in SPOT_COUNTS:
        frame = data.loc[data.timepoint == stage]
        xy = frame[COORDINATES].to_numpy(dtype=float)
        canonical_truth = frame.true_labels.to_numpy()
        nearest, ties = independent_neighborhoods(xy)
        standardized = (xy - xy.mean(axis=0)) / xy.std(axis=0)
        for method in METHODS:
            if native_order is None:
                ordered = frame
            else:
                positions = native_order.loc[frame.spot, method].to_numpy()
                ordered = frame.iloc[np.argsort(positions)]
            native_xy = ordered[COORDINATES].to_numpy(dtype=float)
            truth, predicted = ordered.true_labels.to_numpy(), ordered[method].to_numpy()
            canonical_labels = frame[method].to_numpy()
            metrics = {
                "ARI": adjusted_rand_score(truth, predicted),
                "AMI": adjusted_mutual_info_score(truth, predicted),
                "H": homogeneity_score(truth, predicted),
                "PAS": compute_pas(native_xy, predicted, k=K),
                "CHAOS": compute_chaos(native_xy, predicted),
            }
            original = np.array([metrics[name] for name in ("ARI", "AMI", "H")])
            independent = annotation_metrics_independent(canonical_truth, canonical_labels)
            np.testing.assert_allclose(original, independent, rtol=0, atol=1e-10)
            chaos_difference = abs(metrics["CHAOS"] - chaos_independent(standardized, canonical_labels))
            assert chaos_difference <= 1e-12, (stage, method, chaos_difference)
            lower, upper = pas_bounds(canonical_labels, nearest, ties)
            assert lower - 1e-12 <= metrics["PAS"] <= upper + 1e-12
            assert all(np.isfinite(value) for value in metrics.values())
            audits.append(dict(method=method, timepoint=stage, spots=len(frame),
                               clusters=len(set(predicted)), pas_tied_spots=len(ties),
                               pas_lower=float(lower), pas_upper=float(upper),
                               max_annotation_difference=float(np.max(np.abs(original - independent))),
                               chaos_difference=float(chaos_difference)))
            rows.extend(dict(index=metric, timepoint=stage, value=value, method=method)
                        for metric, value in metrics.items())
    per_stage = pd.DataFrame(rows, columns=["index", "timepoint", "value", "method"])
    assert len(per_stage) == len(SPOT_COUNTS) * len(METHODS) * 5
    averages = per_stage.groupby(["index", "method"], sort=False, as_index=False).value.mean()
    averages["timepoint"] = "avg"
    combined = pd.concat([per_stage, averages.loc[:, per_stage.columns]], ignore_index=True)
    assert len(combined) == 325
    assert not combined.duplicated(["index", "timepoint", "method"]).any()
    return combined, audits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", type=Path)
    parser.add_argument("original_utils", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--audit", type=Path)
    parser.add_argument("--native-order", type=Path, help="Per-method original row positions keyed by spot")
    args = parser.parse_args()
    data = read_labels(args.labels)
    native_order = read_native_order(args.native_order, data) if args.native_order else None
    with threadpool_limits(limits=1):
        metrics, audits = compute(data, *original_continuity_functions(args.original_utils), native_order)
    metrics.to_csv(args.output, index=False)
    if args.audit:
        audit = dict(all_passed=True, labels_sha256=hashlib.sha256(args.labels.read_bytes()).hexdigest(),
                     original_utils_sha256=hashlib.sha256(args.original_utils.read_bytes()).hexdigest(),
                     native_order_sha256=(hashlib.sha256(args.native_order.read_bytes()).hexdigest()
                                          if args.native_order else None),
                     spot_counts=SPOT_COUNTS, methods=list(METHODS), checks=audits)
        args.audit.write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
