"""Render original mean curves and paired, stage-wise Moran I differences."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATASETS = (
    ('MouseEmbryo', 'Mouse Embryo Organogenesis', ('E12.5', 'E13.5', 'E15.5')),
    ('ChickenHeart', 'Chicken Heart Cardiogenesis', ('Day 4', 'Day 7', 'Day 10', 'Day 14')),
    ('AxolotlTelencephalon', 'Axolotl Telencephalon Neurogenesis', ('St.44', 'St.54', 'St.57')),
)
BRANCHES = ('GP10', 'GAT10')
BRANCH_COLORS = ('#2458a6', '#d46a19')
STAGE_COLORS = ('#6a9dba', '#b57d50', '#629b83', '#9875b0')
ZERO = '#667384'
MEAN = '#18222f'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def load_results(source):
    manifest = read_json(source / 'artifact_sha256.json')
    for filename, digest in manifest.items():
        assert sha(source / filename) == digest, filename
    proof = read_json(source / 'verification.json')
    assert proof['all_passed'] and proof['distance_bins'] == 10
    tables = []
    for file in ('per_timepoint.csv', 'summary.csv'):
        table = pd.read_csv(source / file, float_precision='round_trip')
        table = table[table.normalization.eq('trace_normalized')]
        assert set(table.dataset) == {row[0] for row in DATASETS}
        tables.append(table)
    return tuple(tables)


def paired_data(raw, summary, case, labels):
    original = {}
    for branch in BRANCHES:
        part = summary[(summary.dataset == case) & (summary.branch == branch)].sort_values('distance')
        assert len(part) == 100 and part.n_timepoints.eq(len(labels)).all()
        original[branch] = part
    grid = original['GP10'].distance.to_numpy()
    assert np.array_equal(grid, original['GAT10'].distance)
    assert np.all(np.diff(grid) > 0)

    stages = []
    all_support = []
    case_raw = raw[raw.dataset == case]
    assert set(case_raw.timepoint) == {f't{k}' for k in range(len(labels))}
    for k, label in enumerate(labels):
        stage = case_raw[case_raw.timepoint == f't{k}']
        gp, gat = (stage[stage.branch == branch].sort_values('bin') for branch in BRANCHES)
        for block in (gp, gat):
            assert np.array_equal(block.bin, np.arange(10))
            assert np.isfinite(block.Moran_I_trace).all()
        for field in ('distance_left', 'distance_right', 'distance_mean', 'sampled_unordered_pairs'):
            assert np.array_equal(gp[field], gat[field]), (case, label, field)
        distance = gp.distance_mean.to_numpy()
        delta = gp.Moran_I_trace.to_numpy() - gat.Moran_I_trace.to_numpy()
        assert np.all(np.diff(distance) > 0)
        assert grid[0] >= distance[0] and grid[-1] <= distance[-1]
        all_support.append((distance[0], distance[-1]))
        stages.append((label, distance, delta))
    assert grid[0] == max(a for a, _ in all_support)
    assert grid[-1] == min(b for _, b in all_support)
    stage_delta = np.stack([np.interp(grid, distance, delta) for _, distance, delta in stages])
    mean_delta = stage_delta.mean(axis=0)
    expected = original['GP10']['mean'].to_numpy() - original['GAT10']['mean'].to_numpy()
    error = float(np.max(np.abs(mean_delta - expected)))
    assert error < 1e-12
    return original, grid, stages, mean_delta, error


def style_axis(ax):
    ax.spines[['right', 'top']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color('#687484')
    ax.tick_params(axis='both', length=4, color='#687484', labelcolor='#202b38')
    ax.grid(axis='y', color='#e1e7ef', lw=.75)
    ax.set_axisbelow(True)
    ax.axhline(0, color=ZERO, lw=1.05, ls=(0, (4, 3)), zorder=2)
