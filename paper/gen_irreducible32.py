"""Generate figure showing all 32 irreducible (fault-free) Stomachion tilings.

Usage: cd ~/git/wooden && python paper/gen_irreducible32.py
"""
from __future__ import annotations
import sys
import os
import time

# Ensure the package root is on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from solver.stomachion import PIECE_VERTS, PIECE_NAMES
from solver.solver import build_solver
from solver.geometry import orientations, translate


# ---------------------------------------------------------------------------
# Palette
# ---------------------------------------------------------------------------
COLORS = [
    '#e8c99b', '#a0522d', '#deb887', '#6b4226', '#cd853f', '#8b7355', '#d4a574',
    '#704030', '#c9a06c', '#8b4513', '#d2a679', '#966600', '#c4956a', '#bc8f5e',
]
STROKES = [
    '#a07840', '#602a10', '#9a7848', '#3d2218', '#875218', '#504028', '#8a5e30',
    '#3c1c14', '#846028', '#502008', '#8b6238', '#5e4000', '#7d5528', '#7a5628',
]


# ---------------------------------------------------------------------------
# Reconstruct piece vertices from placement
# ---------------------------------------------------------------------------
def placement_vertices(pl, piece_verts):
    """Reconstruct the absolute vertices of a placed piece."""
    verts = piece_verts[pl.piece_idx]
    orients = orientations(verts)
    oriented = orients[pl.orient_idx]
    placed = translate(oriented, pl.tx, pl.ty)
    return placed


# ---------------------------------------------------------------------------
# Vertex-based canonical form under D4
# ---------------------------------------------------------------------------
def solution_vertex_key(sol, piece_verts):
    """Compute a canonical vertex-based key for a solution.

    Returns a sorted tuple of per-piece sorted-vertex-tuples.
    This ignores piece labels (just the set of regions matters).
    """
    regions = []
    for pl in sol:
        verts = placement_vertices(pl, piece_verts)
        regions.append(tuple(sorted(verts)))
    return tuple(sorted(regions))


def d4_transform_vertex(x, y, board_size, rot, flip):
    """Apply D4 transform to a vertex coordinate.

    flip: reflect x -> board_size - x (applied first)
    rot: number of 90-degree CCW rotations of (x,y) -> (board_size-y, x)
    """
    if flip:
        x = board_size - x
    for _ in range(rot):
        x, y = board_size - y, x
    return (x, y)


def d4_transform_region(region, board_size, rot, flip):
    """Transform a sorted vertex tuple under D4."""
    transformed = tuple(sorted(
        d4_transform_vertex(x, y, board_size, rot, flip)
        for x, y in region
    ))
    return transformed


def canonical_under_d4_vertices(vertex_key, board_size):
    """Return the lex-minimum D4 transform of a vertex-based tiling key."""
    best = vertex_key
    for flip in (False, True):
        for rot in range(4):
            transformed = tuple(sorted(
                d4_transform_region(region, board_size, rot, flip)
                for region in vertex_key
            ))
            if transformed < best:
                best = transformed
    return best


# ---------------------------------------------------------------------------
# Fault-free (irreducible) filter
# ---------------------------------------------------------------------------
def is_fault_free(sol, piece_verts, board_size=12) -> bool:
    """Check if a solution is fault-free (no straight line divides all pieces to one side).

    A fault line is an axis-aligned line or a diagonal that no piece crosses.
    We check:
      - Horizontal lines y = k for k = 1..11
      - Vertical lines x = k for k = 1..11
      - Diagonal y = x
      - Anti-diagonal y + x = 12
    """
    # Get vertices for each placed piece
    piece_vertex_sets = []
    for pl in sol:
        verts = placement_vertices(pl, piece_verts)
        piece_vertex_sets.append(verts)

    # For each potential fault line, check if any piece crosses it.
    # A piece crosses a line if it has vertices on both sides of it.
    # Note: some vertices may lie exactly on the line, so we need to check
    # for vertices strictly on each side.

    # However, using just vertices is insufficient for convex pieces with
    # edges crossing the line while vertices are on the same side. But for
    # Stomachion pieces (all convex polygons on integer lattice), if a polygon
    # has all vertices on one side (or on) a line, then the entire polygon is
    # on that side. So vertex-based checking is sufficient.

    # Check horizontal fault lines y = k (k = 1..11)
    for k in range(1, board_size):
        all_one_side = True
        for verts in piece_vertex_sets:
            has_above = any(y < k for x, y in verts)  # y < k means above in y-down
            has_below = any(y > k for x, y in verts)  # y > k means below in y-down
            if has_above and has_below:
                all_one_side = False
                break
        if all_one_side:
            return False

    # Check vertical fault lines x = k (k = 1..11)
    for k in range(1, board_size):
        all_one_side = True
        for verts in piece_vertex_sets:
            has_left = any(x < k for x, y in verts)
            has_right = any(x > k for x, y in verts)
            if has_left and has_right:
                all_one_side = False
                break
        if all_one_side:
            return False

    # Check diagonal y = x
    # A piece crosses if it has vertices with y > x AND vertices with y < x
    all_one_side = True
    for verts in piece_vertex_sets:
        has_above = any(y > x for x, y in verts)
        has_below = any(y < x for x, y in verts)
        if has_above and has_below:
            all_one_side = False
            break
    if all_one_side:
        return False

    # Check anti-diagonal y + x = 12
    # A piece crosses if it has vertices with y + x > 12 AND y + x < 12
    all_one_side = True
    for verts in piece_vertex_sets:
        has_above = any(y + x > board_size for x, y in verts)
        has_below = any(y + x < board_size for x, y in verts)
        if has_above and has_below:
            all_one_side = False
            break
    if all_one_side:
        return False

    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    board_size = 12
    res = 2

    print("=" * 60)
    print("Generating all 32 irreducible Stomachion tilings")
    print("=" * 60)
    sys.stdout.flush()

    # Step 1: Build and solve
    dlx, placements = build_solver(PIECE_VERTS, board_size, res, verbose=True)
    sys.stdout.flush()

    print("  Solving (this may take ~10 minutes)...")
    sys.stdout.flush()
    t0 = time.time()
    count = dlx.solve(max_solutions=0)
    t1 = time.time()
    print(f"  Raw solutions: {count} ({t1 - t0:.1f}s)")
    sys.stdout.flush()

    # Extract solutions as lists of Placement objects
    raw_solutions = []
    for sol_row_ids in dlx.solutions:
        raw_solutions.append([placements[rid] for rid in sol_row_ids])

    # Step 2: Deduplicate under D4 using vertex-based canonical forms
    print("  Deduplicating under D4...")
    sys.stdout.flush()
    t2 = time.time()

    seen = {}  # canonical_key -> first solution
    for sol in raw_solutions:
        vkey = solution_vertex_key(sol, PIECE_VERTS)
        canon = canonical_under_d4_vertices(vkey, board_size)
        if canon not in seen:
            seen[canon] = sol

    t3 = time.time()
    print(f"  Geometrically distinct: {len(seen)} ({t3 - t2:.1f}s)")
    sys.stdout.flush()

    # Step 3: Filter for fault-free
    print("  Filtering for fault-free tilings...")
    sys.stdout.flush()
    irreducible = []
    for canon, sol in seen.items():
        if is_fault_free(sol, PIECE_VERTS, board_size):
            irreducible.append((canon, sol))

    print(f"  Fault-free (irreducible): {len(irreducible)}")
    sys.stdout.flush()

    if len(irreducible) != 32:
        print(f"  WARNING: Expected 32, got {len(irreducible)}")
        sys.stdout.flush()

    # Step 4: Sort for consistent ordering
    irreducible.sort(key=lambda x: x[0])

    # Step 5: Render
    n_tilings = len(irreducible)
    ncols = 4
    nrows = (n_tilings + ncols - 1) // ncols
    print(f"  Rendering {n_tilings} tilings in {nrows}x{ncols} grid...")
    sys.stdout.flush()

    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 4 * nrows))
    fig.patch.set_facecolor('#1c1410')

    if nrows == 1:
        axes = axes.reshape(1, -1)

    for idx, (canon, sol) in enumerate(irreducible):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]

        # Draw each piece
        for pl in sol:
            verts = placement_vertices(pl, PIECE_VERTS)
            piece_idx = pl.piece_idx

            # Flip y for display (y=0 at bottom)
            display_verts = [(x, board_size - y) for x, y in verts]
            display_verts.append(display_verts[0])  # close polygon
            xs, ys = zip(*display_verts)

            ax.fill(xs, ys, color=COLORS[piece_idx], edgecolor=STROKES[piece_idx],
                    linewidth=1.2, zorder=2)

            # Label at centroid
            cx = sum(x for x, y in verts) / len(verts)
            cy = board_size - sum(y for x, y in verts) / len(verts)
            label = PIECE_NAMES[piece_idx]

            # Choose label color based on background brightness
            bg_color = COLORS[piece_idx]
            r_val = int(bg_color[1:3], 16)
            g_val = int(bg_color[3:5], 16)
            b_val = int(bg_color[5:7], 16)
            brightness = 0.299 * r_val + 0.587 * g_val + 0.114 * b_val
            label_color = '#1a0e08' if brightness > 140 else '#f0d8b0'

            ax.text(cx, cy, label, ha='center', va='center',
                    fontsize=6, fontweight='bold', color=label_color, zorder=3)

        ax.set_xlim(-0.15, board_size + 0.15)
        ax.set_ylim(-0.15, board_size + 0.15)
        ax.set_aspect('equal')
        ax.set_facecolor('#241a0e')
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color('#3d2b1e')
            spine.set_linewidth(1.5)
        ax.set_title(f'#{idx + 1}', fontsize=8, fontweight='bold',
                     color='#e8a862', pad=3)

    # Hide any unused panels
    for idx in range(n_tilings, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle('All 32 irreducible Stomachion tilings', fontsize=20,
                 fontweight='bold', color='#e8a862', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.99], h_pad=0.8, w_pad=0.5)

    outpath = os.path.join(os.path.dirname(__file__), 'irreducible32.png')
    fig.savefig(outpath, dpi=200, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved to {outpath}")


if __name__ == '__main__':
    main()
