"""Visualize dissections as colored polygon plots."""
from __future__ import annotations
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
from .explore import random_dissection
from .geometry import area2


# Warm wood palette (light to dark)
WOOD_COLORS = [
    '#e8c99b', '#deb887', '#d4a574', '#cd853f', '#c9a06c',
    '#c4956a', '#bc8f5e', '#b8860b', '#a0785a', '#a0522d',
    '#9c7c5c', '#8b7355', '#8b6914', '#8b4513', '#704030',
]


def plot_dissection(
    pieces,
    width=12, height=12,
    title=None,
    ax=None,
    show_labels=True,
):
    """Plot a dissection with colored pieces."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    else:
        fig = ax.figure

    n = len(pieces)

    # Sort pieces by area (largest first) for color assignment
    indexed = sorted(enumerate(pieces), key=lambda x: -area2(x[1]))
    color_map = {}
    for rank, (i, _) in enumerate(indexed):
        color_map[i] = WOOD_COLORS[rank % len(WOOD_COLORS)]

    for i, poly in enumerate(pieces):
        # Flip y for display (y=0 at bottom)
        verts = [(x, height - y) for x, y in poly]
        verts.append(verts[0])  # close polygon
        xs, ys = zip(*verts)

        ax.fill(xs, ys, color=color_map[i], edgecolor='#3d2b1e',
                linewidth=1.5, zorder=2)

        if show_labels:
            # Label at centroid
            cx = sum(x for x, y in poly) / len(poly)
            cy = height - sum(y for x, y in poly) / len(poly)
            a = area2(poly) // 2
            ax.text(cx, cy, f'{a}', ha='center', va='center',
                    fontsize=8, fontweight='bold', color='#2a1510',
                    zorder=3)

    ax.set_xlim(-0.3, width + 0.3)
    ax.set_ylim(-0.3, height + 0.3)
    ax.set_aspect('equal')
    ax.set_facecolor('#241a0e')

    # Grid
    for i in range(width + 1):
        ax.axvline(i, color='#3d2b1e', linewidth=0.3, alpha=0.3, zorder=1)
    for i in range(height + 1):
        ax.axhline(i, color='#3d2b1e', linewidth=0.3, alpha=0.3, zorder=1)

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color('#3d2b1e')
        spine.set_linewidth(2)

    if title:
        ax.set_title(title, fontsize=12, fontweight='bold', color='#e8a862',
                     pad=8)

    return fig, ax


def champion_gallery(champions: dict, width=12, height=12, outpath='/tmp/champions.png'):
    """Plot a gallery of champion prime dissections.

    champions: dict mapping n -> seed
    """
    ns = sorted(champions.keys())
    ncols = min(4, len(ns))
    nrows = (len(ns) + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5.5 * nrows))
    fig.patch.set_facecolor('#1c1410')

    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes.reshape(1, -1)
    elif ncols == 1:
        axes = axes.reshape(-1, 1)

    for idx, n in enumerate(ns):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]

        seed = champions[n]
        pieces = random_dissection(width, height, n, seed=seed)

        if pieces:
            from .experiment import count_tilings_geometric
            geo = count_tilings_geometric(pieces, width, height, res=2)
            plot_dissection(pieces, width, height,
                          title=f'n={n}  |  {geo} tilings',
                          ax=ax, show_labels=True)
        else:
            ax.set_visible(False)

    # Hide unused axes
    for idx in range(len(ns), nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle('Champion Prime Dissections', fontsize=16,
                fontweight='bold', color='#e8a862', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outpath, dpi=150, bbox_inches='tight',
               facecolor=fig.get_facecolor())
    print(f'Saved to {outpath}')
    return outpath


if __name__ == '__main__':
    # Champion seeds from the search
    champions = {
        3: 31,
        4: 116,
        5: 343,
        6: 133,
    }
    champion_gallery(champions)
