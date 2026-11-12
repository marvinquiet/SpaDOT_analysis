"""Draw the first zero-crossing comparison for the Moran curves."""
import numpy as np
from plot_helpers import BRANCH_COLORS

def first_zero_crossing(x, y):
    hits = np.flatnonzero((y[:-1] > 0) & (y[1:] <= 0))
    assert hits.size > 0
    k = hits[0]
    return float(x[k] - y[k] * (x[k + 1] - x[k]) / (y[k + 1] - y[k]))


def draw_crossings(ax, original, grid):
    ax.set_ylim(0, 1.1)
    ax.set_yticks([])
    ax.spines[['top', 'left', 'right']].set_visible(False)
    ax.spines['bottom'].set_color('#687484')
    ax.tick_params(axis='x', length=4, color='#687484', labelcolor='#202b38')
    crossings = {}
    for branch, label, color, y in zip(
        ('GP10', 'GAT10'), ('SVGP', 'GAT'), BRANCH_COLORS, (.83, .32),
    ):
        crossing = first_zero_crossing(grid, original[branch]['mean'].to_numpy())
        crossings[branch] = crossing
        ax.axhline(y, color='#e4eaf1', linewidth=.8, zorder=0)
        ax.plot(crossing, y, marker='o', ms=5, color=color, zorder=2)
        ax.annotate(f'{label}  {crossing:.2f}', (crossing, y), xytext=(7, 0),
                    textcoords='offset points', va='center', color=color,
                    fontsize=10, fontweight='bold',
                    bbox=dict(facecolor='white', edgecolor='none', pad=.45))
    return crossings
