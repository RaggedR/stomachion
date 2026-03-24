"""Explore alternative dissections using the coproduct structure.

Strategy: choose a fault line that splits the 12×12 square into two
independent sub-regions. Design piece sets for each sub-region.
Total solutions = product of sub-region solutions.

To beat 536, we need the product to exceed 536.

For a horizontal fault line y=h splitting into:
  - Top: 12×h rectangle, k pieces
  - Bottom: 12×(12-h) rectangle, (14-k) pieces
we need Sol(top) × Sol(bottom) > 536.

This module generates random dissections of rectangles by repeated
lattice-point cuts, counts their tilings, and searches for maxima.
"""
from __future__ import annotations
import random
import time
import math
from .geometry import (
    Polygon, shoelace2, area2, rasterize, orientations,
    normalize, translate, bounding_box,
)
from .pieces import Placement
from .dlx import DLX


# ═══════════════════════════════════════════════════════
# Polygon cutting — split convex polygons along lattice lines
# ═══════════════════════════════════════════════════════

def boundary_lattice_points(poly: Polygon) -> list[tuple[int, int]]:
    """Return all lattice points on the boundary of an integer-coord polygon."""
    pts = set()
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        dx, dy = x2 - x1, y2 - y1
        g = math.gcd(abs(dx), abs(dy)) if (dx or dy) else 0
        if g == 0:
            pts.add((x1, y1))
            continue
        sx, sy = dx // g, dy // g
        for t in range(g + 1):
            pts.add((x1 + t * sx, y1 + t * sy))
    return sorted(pts)


def _edge_index(poly: Polygon, p: tuple[int, int]) -> int | None:
    """Find which edge of poly contains point p.
    Returns edge index i (p is on edge poly[i]→poly[i+1])."""
    n = len(poly)
    px, py = p
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        # Cross product = 0 means collinear
        cross = (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1)
        if cross != 0:
            continue
        # Check bounds
        if (min(x1, x2) <= px <= max(x1, x2) and
            min(y1, y2) <= py <= max(y1, y2)):
            return i
    return None


def _clean_poly(verts: list[tuple[int, int]]) -> Polygon:
    """Remove consecutive duplicates and collinear points."""
    if not verts:
        return ()
    # Remove consecutive duplicates
    clean = [verts[0]]
    for v in verts[1:]:
        if v != clean[-1]:
            clean.append(v)
    if len(clean) > 1 and clean[-1] == clean[0]:
        clean.pop()
    # Remove collinear points
    if len(clean) < 3:
        return tuple(clean)
    result = []
    n = len(clean)
    for i in range(n):
        p = clean[(i - 1) % n]
        q = clean[i]
        r = clean[(i + 1) % n]
        cross = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
        if cross != 0:
            result.append(q)
    return tuple(result) if len(result) >= 3 else tuple(clean)


def cut_polygon(poly: Polygon, p1: tuple[int, int], p2: tuple[int, int]):
    """Cut polygon along line from p1 to p2 (both boundary lattice points).
    Returns (poly_a, poly_b) or None if invalid cut."""
    e1 = _edge_index(poly, p1)
    e2 = _edge_index(poly, p2)
    if e1 is None or e2 is None or e1 == e2:
        return None

    n = len(poly)

    # Ensure e1 < e2
    if e1 > e2:
        e1, e2 = e2, e1
        p1, p2 = p2, p1

    # Poly A: p1, vertices (e1+1)..e2, p2
    a = [p1]
    for i in range(e1 + 1, e2 + 1):
        if poly[i] != p1:
            a.append(poly[i])
    if poly[e2] != p2:
        a.append(p2)

    # Poly B: p2, vertices (e2+1)..n-1, 0..e1, p1
    b = [p2]
    for i in range(e2 + 1, n):
        b.append(poly[i])
    for i in range(0, e1 + 1):
        b.append(poly[i])
    if poly[e1] != p1 and b[-1] != p1:
        b.append(p1)

    a = _clean_poly(a)
    b = _clean_poly(b)

    if len(a) < 3 or len(b) < 3:
        return None
    if area2(a) == 0 or area2(b) == 0:
        return None

    return a, b


# ═══════════════════════════════════════════════════════
# Random dissection generator
# ═══════════════════════════════════════════════════════

def random_dissection(
    width: int,
    height: int,
    num_pieces: int,
    seed: int | None = None,
) -> list[Polygon] | None:
    """Generate a random dissection of a rectangle into num_pieces pieces.

    Uses repeated random cuts: choose a piece, choose two boundary lattice
    points, cut. All pieces remain convex with integer vertices.

    Returns list of piece polygons, or None if generation failed.
    """
    if seed is not None:
        random.seed(seed)

    rect = ((0, 0), (width, 0), (width, height), (0, height))
    pieces = [rect]

    for _ in range(num_pieces - 1):
        # Pick a random piece to cut (prefer larger pieces)
        areas = [area2(p) for p in pieces]
        total = sum(areas)
        if total == 0:
            return None

        # Weighted random choice by area
        r = random.randint(1, total)
        cumul = 0
        chosen = 0
        for i, a in enumerate(areas):
            cumul += a
            if cumul >= r:
                chosen = i
                break

        poly = pieces[chosen]
        bpts = boundary_lattice_points(poly)
        if len(bpts) < 4:
            return None

        # Try random pairs of boundary points
        found = False
        attempts = 0
        while attempts < 50 and not found:
            p1, p2 = random.sample(bpts, 2)
            result = cut_polygon(poly, p1, p2)
            if result is not None:
                a, b = result
                if area2(a) >= 2 and area2(b) >= 2:  # min area threshold
                    pieces[chosen] = a
                    pieces.append(b)
                    found = True
            attempts += 1

        if not found:
            return None

    return pieces


# ═══════════════════════════════════════════════════════
# Sub-region tiling solver
# ═══════════════════════════════════════════════════════

def count_tilings(
    pieces: list[Polygon],
    width: int,
    height: int,
    res: int = 2,
) -> int:
    """Count distinct tilings of a width×height rectangle using given pieces.

    Returns the raw count (includes D4 symmetries and piece permutations).
    """
    grid_w = width * res
    grid_h = height * res
    num_cells = grid_w * grid_h
    num_pieces = len(pieces)

    # Enumerate all placements for each piece within the rectangle
    all_placements = []
    for pi, pverts in enumerate(pieces):
        pverts_norm = normalize(pverts)
        orients = orientations(pverts_norm)

        for oi, oriented in enumerate(orients):
            mnx, mny, mxx, mxy = bounding_box(oriented)
            ow = mxx - mnx
            oh = mxy - mny

            # Rasterize once at origin
            base_cells = _rasterize_rect(oriented, grid_w, grid_h, res)
            if not base_cells:
                continue

            # Convert to (row, col) offsets
            offsets = tuple((c // grid_w, c % grid_w) for c in base_cells)

            for tx in range(0, width - ow + 1):
                for ty in range(0, height - oh + 1):
                    dr = ty * res
                    dc = tx * res
                    cells = frozenset(
                        (r + dr) * grid_w + (c + dc)
                        for r, c in offsets
                        if 0 <= r + dr < grid_h and 0 <= c + dc < grid_w
                    )
                    if len(cells) == len(offsets):  # all cells within bounds
                        all_placements.append(Placement(
                            piece_idx=pi, orient_idx=oi,
                            tx=tx, ty=ty, cells=cells,
                        ))

    if not all_placements:
        return 0

    # Build DLX
    num_primary = num_pieces + num_cells
    dlx = DLX(num_primary)

    for pl_idx, pl in enumerate(all_placements):
        cols = [pl.piece_idx]
        cols.extend(num_pieces + c for c in pl.cells)
        dlx.add_row(cols, row_id=pl_idx)

    return dlx.solve()


def _rasterize_rect(
    verts: Polygon,
    grid_w: int,
    grid_h: int,
    res: int,
) -> frozenset[int]:
    """Rasterize polygon onto a grid_w × grid_h fine grid."""
    inv = 1.0 / res
    half = 0.5 * inv

    mnx, mny, mxx, mxy = bounding_box(verts)
    c_lo = max(int(mnx * res), 0)
    c_hi = min(int(mxx * res), grid_w - 1)
    r_lo = max(int(mny * res), 0)
    r_hi = min(int(mxy * res), grid_h - 1)

    n = len(verts)
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]

    cells = []
    for r in range(r_lo, r_hi + 1):
        py = r * inv + half
        for c in range(c_lo, c_hi + 1):
            px = c * inv + half
            inside = False
            j = n - 1
            for i in range(n):
                yi, yj = ys[i], ys[j]
                if (yi > py) != (yj > py):
                    xi, xj = xs[i], xs[j]
                    if px < xi + (py - yi) * (xj - xi) / (yj - yi):
                        inside = not inside
                j = i
            if inside:
                cells.append(r * grid_w + c)
    return frozenset(cells)


# ═══════════════════════════════════════════════════════
# Coproduct exploration — fault-line-first search
# ═══════════════════════════════════════════════════════

def explore_fault_line(
    fault_y: int = 6,
    pieces_top: int = 7,
    pieces_bot: int = 7,
    trials: int = 50,
    res: int = 2,
    verbose: bool = True,
):
    """Explore dissections split by horizontal fault line y=fault_y.

    For each trial:
      1. Generate random dissection of top rectangle (12×fault_y, pieces_top pieces)
      2. Generate random dissection of bottom rectangle (12×(12-fault_y), pieces_bot pieces)
      3. Count tilings of each half
      4. Product = total solution count

    Returns the best (product, top_dissection, bottom_dissection) found.
    """
    width = 12
    h_top = fault_y
    h_bot = 12 - fault_y

    best_product = 0
    best_top = None
    best_bot = None
    best_counts = (0, 0)

    if verbose:
        print(f"Exploring fault line y={fault_y}: "
              f"top {width}×{h_top} ({pieces_top}pc) × "
              f"bot {width}×{h_bot} ({pieces_bot}pc)")
        print(f"  Target: product > 536")
        print()

    for trial in range(trials):
        t0 = time.time()

        top = random_dissection(width, h_top, pieces_top, seed=trial * 2)
        bot = random_dissection(width, h_bot, pieces_bot, seed=trial * 2 + 1)

        if top is None or bot is None:
            if verbose:
                print(f"  Trial {trial+1:3d}: generation failed")
            continue

        # Verify areas
        top_area = sum(area2(p) for p in top)
        bot_area = sum(area2(p) for p in bot)
        if top_area != 2 * width * h_top or bot_area != 2 * width * h_bot:
            if verbose:
                print(f"  Trial {trial+1:3d}: area mismatch")
            continue

        ct = count_tilings(top, width, h_top, res)
        cb = count_tilings(bot, width, h_bot, res)
        product = ct * cb

        dt = time.time() - t0

        if verbose:
            marker = " ★" if product > best_product and product > 0 else ""
            print(f"  Trial {trial+1:3d}: top={ct:5d} × bot={cb:5d} = {product:8d}  ({dt:.1f}s){marker}")

        if product > best_product:
            best_product = product
            best_top = top
            best_bot = bot
            best_counts = (ct, cb)

    if verbose:
        print()
        print(f"Best: {best_counts[0]} × {best_counts[1]} = {best_product}")
        if best_product > 536:
            print(f"  ★ BEATS ARCHIMEDES (536)!")
        else:
            print(f"  Still below 536")

    return best_product, best_top, best_bot


def explore_all_splits(trials_per_split: int = 30, res: int = 2):
    """Try all horizontal fault lines y=1..11 with various piece splits."""
    print("=" * 70)
    print("Coproduct Exploration: Can we beat 536?")
    print("=" * 70)
    print()

    results = []

    for fault_y in [6, 4, 3, 8, 5, 7, 2, 9, 10]:
        h_top = fault_y
        h_bot = 12 - fault_y

        # Try a few piece splits for each fault line
        for k in range(max(2, h_top // 2), min(13, h_top * 3 + 1)):
            pieces_top = k
            pieces_bot = 14 - k

            if pieces_top < 2 or pieces_bot < 2:
                continue

            # Area feasibility: each piece needs area ≥ 1
            if 12 * h_top < pieces_top or 12 * h_bot < pieces_bot:
                continue

            product, _, _ = explore_fault_line(
                fault_y=fault_y,
                pieces_top=pieces_top,
                pieces_bot=pieces_bot,
                trials=trials_per_split,
                res=res,
                verbose=True,
            )
            results.append((product, fault_y, pieces_top, pieces_bot))

    # Summary
    results.sort(reverse=True)
    print()
    print("=" * 70)
    print("Top 10 results:")
    print("=" * 70)
    for prod, fy, pt, pb in results[:10]:
        print(f"  y={fy:2d}, {pt:2d}+{pb:2d} pieces: {prod:8d}")

    return results


if __name__ == "__main__":
    # Quick exploration: y=6 split with 7+7 pieces
    explore_fault_line(fault_y=6, pieces_top=7, pieces_bot=7, trials=50)
