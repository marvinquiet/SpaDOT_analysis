"""Shared arm order and cell styles for the S39 plots."""
import tempfile
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "spadot_figures"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ARMS = ('dual_native', 'gat_d1', 'gat_d2', 'gat_d3', 'gat_d4', 'gat_d6',
        'gat_d3_parammatched', 'gat_d3_expanded', 'gat_gat_matched', 'mlp_gat')
def text_color(rgba):
    rgb = np.asarray(rgba[:3])
    linear = np.where(rgb <= .04045, rgb / 12.92, ((rgb + .055) / 1.055) ** 2.4)
    luminance = float(linear @ np.array([.2126, .7152, .0722]))
    return '#172438' if luminance > .21 else 'white'


def rules(ax):
    for end in (0.5, 5.5, 6.5, 7.5):
        ax.axhline(end, color='white', lw=1.7)
