"""Compute the S40a distance-binned Moran correlograms.

Provide the frozen 20D posterior-mean latent H5AD for each dataset, with 10 GP
dimensions followed by 10 GAT dimensions, obs.timepoint and obsm['spatial'].
The paper used the final seed-2026 checkpoint at epoch 99. Requires Python
with anndata, NumPy and pandas. No model fitting, clustering or OT is performed.
"""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import latent_structure
from moran_helpers import weight_statistics

BINS = 10
CASES = ('MouseEmbryo', 'ChickenHeart', 'AxolotlTelencephalon')
TIMEPOINTS = {'MouseEmbryo': ('t0', 't1', 't2'), 'ChickenHeart': ('t0', 't1', 't2', 't3'),
              'AxolotlTelencephalon': ('t0', 't1', 't2')}
SPOTS = {'MouseEmbryo': (2397, 2312, 3546), 'ChickenHeart': (747, 1966, 1916, 1967),
         'AxolotlTelencephalon': (1477, 2929, 4410)}
BRANCHES = ('GP10', 'GAT10')
NORMALIZATIONS = ('trace_normalized', 'ridge_whitened')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def summarize_timepoints(table, grid_points=100):
    """Align observed distance coordinates, never bin ordinals or extrapolation."""
    assert grid_points >= 2
    assert not table.duplicated(['dataset', 'timepoint', 'branch', 'normalization', 'bin']).any()
    summaries, aligned, supports = [], [], {}
    for case in CASES:
        case_rows = table[table.dataset == case]
        assert set(case_rows.timepoint) == set(TIMEPOINTS[case])
        # Geometry determines support, identically for both branches and norms.
        bounds = []
        for tp in TIMEPOINTS[case]:
            geometry = case_rows[(case_rows.timepoint == tp) &
                                 (case_rows.branch == 'GP10') &
                                 (case_rows.normalization == 'trace_normalized')].sort_values('distance_mean')
            x = geometry.distance_mean.to_numpy(dtype=float)
            assert len(x) >= 2 and np.isfinite(x).all() and (np.diff(x) > 0).all()
            bounds.append((float(x[0]), float(x[-1])))
            for branch in BRANCHES:
                for norm in NORMALIZATIONS:
                    part = case_rows[(case_rows.timepoint == tp) & (case_rows.branch == branch) &
                                     (case_rows.normalization == norm)].sort_values('distance_mean')
                    assert np.array_equal(part.distance_mean.to_numpy(), x)
        lower = max(x[0] for x in bounds)
        upper = min(x[1] for x in bounds)
        assert lower < upper, ('No common observed distance support', case)
        grid = np.linspace(lower, upper, grid_points)
        supports[case] = dict(lower=lower, upper=upper, grid_points=grid_points,
                              timepoints=list(TIMEPOINTS[case]), observed_bounds=[list(x) for x in bounds],
                              interpolation='linear on bin-mean spatial distances; no extrapolation')
        for branch in BRANCHES:
            for norm in NORMALIZATIONS:
                values = []
                for tp in TIMEPOINTS[case]:
                    part = case_rows[(case_rows.timepoint == tp) & (case_rows.branch == branch) &
                                     (case_rows.normalization == norm)].sort_values('distance_mean')
                    x = part.distance_mean.to_numpy(dtype=float)
                    y = part.Moran_I_trace.to_numpy(dtype=float)
                    assert np.isfinite(y).all() and grid[0] >= x[0] and grid[-1] <= x[-1]
                    value = np.interp(grid, x, y)
                    values.append(value)
                    aligned.extend(dict(dataset=case, timepoint=tp, branch=branch, normalization=norm,
                                        distance=float(d), Moran_I=float(v)) for d, v in zip(grid, value))
                matrix = np.stack(values)
                assert matrix.shape == (len(TIMEPOINTS[case]), grid_points)
                mean, sd = matrix.mean(0), matrix.std(0, ddof=1)
                summaries.extend(dict(dataset=case, branch=branch, normalization=norm, distance=float(d),
                                      mean=float(mu), SD=float(s), n_timepoints=len(values),
                                      lower_SD=float(mu-s), upper_SD=float(mu+s))
                                 for d, mu, s in zip(grid, mean, sd))
    return pd.DataFrame(summaries), pd.DataFrame(aligned), supports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mouse-latent', required=True, type=Path)
    parser.add_argument('--chicken-latent', required=True, type=Path)
    parser.add_argument('--axolotl-latent', required=True, type=Path)
    parser.add_argument('--destination', required=True, type=Path)
    args = parser.parse_args()
    output = args.destination
    if output.exists():
        parser.error('--destination must be a new directory')
    latent_paths = dict(zip(CASES, (args.mouse_latent, args.chicken_latent, args.axolotl_latent)))
    protected, references = {}, {}

    rows, dimensions, provenance = [], [], {}
    max_formula_error = 0.0
    counts = []
    for case, path in latent_paths.items():
        protected[path] = sha(path)
        references[case] = {'latent': str(path), 'latent_sha256': protected[path]}
        latent = ad.read_h5ad(path)
        labels = latent.obs.timepoint.astype(str).to_numpy()
        spatial = np.asarray(latent.obsm['spatial'], dtype=np.float64)
        embedding = np.asarray(latent.X, dtype=np.float64)
        assert latent.n_vars == 20 and latent.obs_names.is_unique
        assert np.isfinite(spatial).all() and np.isfinite(embedding).all()
        assert set(labels) == set(TIMEPOINTS[case])
        assert tuple(int((labels == tp).sum()) for tp in TIMEPOINTS[case]) == SPOTS[case]
        provenance[case] = {}
        for tp_index, tp in enumerate(TIMEPOINTS[case]):
            seed = 20261104 + tp_index
            use = labels == tp
            coords = spatial[use]
            coords = (coords - coords.mean(0)) / coords.std(0)
            i, j, total = latent_structure.unordered_pairs(len(coords), 100_000, seed)
            distance = np.linalg.norm(coords[i] - coords[j], axis=1)
            assert len(i) == 100_000 and total == len(coords) * (len(coords) - 1) // 2
            index_sha = latent_structure.array_sha(np.column_stack([i, j]))
            distance_sha = latent_structure.array_sha(distance)
            edges = np.quantile(distance, np.linspace(0, 1, BINS + 1))
            assert (np.diff(edges) > 0).all()
            groups = np.clip(np.searchsorted(edges, distance, side='right') - 1, 0, BINS - 1)
            count = np.bincount(groups, minlength=BINS)
            assert len(count) == BINS and count.min() > 0 and count.sum() == 100_000
            counts.extend(count.tolist())
            provenance[case][tp] = dict(
                sampling_seed=seed, pair_index_sha256=index_sha,
                pair_distance_sha256=distance_sha, sampled_unordered_pairs=len(i),
                distance_edges=edges.tolist(), pairs_per_bin=count.tolist())
            for branch, columns in latent_structure.BRANCHES.items():
                raw = embedding[use, columns]
                normalized, _ = latent_structure.normalize_latent(raw)
                centered = raw - raw.mean(axis=0)
                for normalization, values in normalized.items():
                    for k in range(BINS):
                        select = groups == k
                        left, right = i[select], j[select]
                        stats = weight_statistics(values, left, right)
                        if normalization == 'trace_normalized':
                            numerator = np.einsum('ij,ij->i', centered[left], centered[right]).mean()
                            denominator = np.einsum('ij,ij->i', centered, centered).mean()
                            error = abs(stats['Moran_I_trace'] - numerator / denominator)
                            assert error < 1e-12, (case, tp, branch, k)
                            max_formula_error = max(max_formula_error, error)
                        arrays = {key: stats.pop(key) for key in
                                  ('dimension_I', 'dimension_variance_weights', 'valid_dimensions')}
                        base = dict(dataset=case, timepoint=tp, branch=branch,
                                    normalization=normalization, bin=k,
                                    distance_left=float(edges[k]), distance_right=float(edges[k + 1]),
                                    distance_mean=float(distance[select].mean()))
                        rows.append(dict(base, **stats))
                        for dim in range(10):
                            valid = arrays['valid_dimensions'][dim]
                            dimensions.append(dict(base, dimension=dim,
                                Moran_I=float(arrays['dimension_I'][dim]) if valid else None,
                                variance_weight=float(arrays['dimension_variance_weights'][dim]),
                                status='defined' if valid else 'undefined_constant_dimension'))
    table, per_dimension = pd.DataFrame(rows), pd.DataFrame(dimensions)
    assert len(table) == 400 and len(per_dimension) == 4000
    summary, aligned, supports = summarize_timepoints(table)
    assert len(summary) == 1200 and len(aligned) == 4000
    for path, digest in protected.items():
        assert sha(path) == digest, path
    output.mkdir()
    for name, frame in dict(per_timepoint=table, per_dimension=per_dimension,
                            summary=summary, interpolated_timepoints=aligned).items():
        frame.to_csv(output / (name + '.csv'), index=False, float_format='%.17g')
    write(output / 'pair_provenance.json', provenance)
    write(output / 'distance_support.json', supports)
    write(output / 'references.json', references)
    write(output / 'verification.json', dict(
        all_passed=True, distance_bins=BINS, datasets=list(CASES),
        n_timepoints={case: len(TIMEPOINTS[case]) for case in CASES},
        paired_samples_between_branches=True,
        sampled_pairs_per_timepoint=100_000, pairs_per_bin_range=[min(counts), max(counts)],
        max_primary_direct_formula_error=max_formula_error,
        aggregation='Equal timepoint mean after linear interpolation on common observed bin-mean distance support',
        shade='Mean +/- one sample SD across timepoints', extrapolation=False,
        input_sha256={str(path): digest for path, digest in protected.items()},
        script_sha256=sha(__file__), no_training=True))
    write(output / 'artifact_sha256.json',
                {path.name: sha(path) for path in sorted(output.iterdir())})


if __name__ == '__main__':
    main()
