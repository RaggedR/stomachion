"""Experiment: distribution of geometric tiling counts across random dissections.

For each n (number of pieces), generates many random dissections of a 12×12
square, counts geometric solutions (up to D4 + piece permutations), and
reports the distribution.
"""
from __future__ import annotations
import sys
import time
from collections import Counter
from .explore import random_dissection, area2, _rasterize_rect
from .geometry import Polygon, orientations, normalize, bounding_box, translate
from .pieces import Placement
from .dlx import DLX
from .symmetry import count_geometric_solutions, canonical_tiling, canonical_under_d4


def count_tilings_geometric(
    pieces: list[Polygon],
    width: int,
    height: int,
    res: int = 2,
) -> int:
    """Count geometrically distinct tilings of a width×height rectangle.

    Mods out D4 of the rectangle and piece label permutations.
    Returns the number of geometric solutions.
    """
    grid_w = width * res
    grid_h = height * res
    num_cells = grid_w * grid_h
    num_pieces = len(pieces)

    # Enumerate placements
    all_placements = []
    for pi, pverts in enumerate(pieces):
        pverts_norm = normalize(pverts)
        orients = orientations(pverts_norm)

        for oi, oriented in enumerate(orients):
            mnx, mny, mxx, mxy = bounding_box(oriented)
            ow = mxx - mnx
            oh = mxy - mny

            base_cells = _rasterize_rect(oriented, grid_w, grid_h, res)
            if not base_cells:
                continue

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
                    if len(cells) == len(offsets):
                        all_placements.append(Placement(
                            piece_idx=pi, orient_idx=oi,
                            tx=tx, ty=ty, cells=cells,
                        ))

    if not all_placements:
        return 0

    # Build and solve DLX
    num_primary = num_pieces + num_cells
    dlx = DLX(num_primary)
    for pl_idx, pl in enumerate(all_placements):
        cols = [pl.piece_idx]
        cols.extend(num_pieces + c for c in pl.cells)
        dlx.add_row(cols, row_id=pl_idx)

    raw_count = dlx.solve()
    if raw_count == 0:
        return 0

    # Extract solutions and count geometric
    raw_solutions = [[all_placements[rid] for rid in sol] for sol in dlx.solutions]

    # For rectangles, use the appropriate grid side
    # D4 only applies to squares; for rectangles, use the rectangle's symmetry group
    if width == height:
        grid_side = width * res
        seen = set()
        for sol in raw_solutions:
            key = canonical_tiling(sol)
            canon = canonical_under_d4(key, grid_side)
            seen.add(canon)
        return len(seen)
    else:
        # Rectangle symmetry: Klein 4-group (id, reflect-x, reflect-y, rotate-180)
        grid_side_w = width * res
        grid_side_h = height * res
        seen = set()
        for sol in raw_solutions:
            key = canonical_tiling(sol)
            # Generate all 4 Klein-4 transforms
            variants = set()
            for flip_h in (False, True):
                for flip_v in (False, True):
                    transformed = tuple(sorted(
                        _klein4_transform(cells, grid_side_w, grid_side_h, flip_h, flip_v)
                        for cells in key
                    ))
                    variants.add(transformed)
            canon = min(variants)
            seen.add(canon)
        return len(seen)


def _klein4_transform(
    cells: frozenset[int], w: int, h: int, flip_h: bool, flip_v: bool
) -> frozenset[int]:
    """Apply Klein-4 transform to cells in a w×h grid."""
    result = set()
    for idx in cells:
        r, c = divmod(idx, w)
        if flip_h:
            c = w - 1 - c
        if flip_v:
            r = h - 1 - r
        result.add(r * w + c)
    return frozenset(result)


def run_experiment(n: int, num_seeds: int = 500, res: int = 2):
    """Run the experiment for n-piece dissections of a 12×12 square."""
    width = height = 12
    counts = Counter()
    total_ok = 0

    print(f"n={n}: generating {num_seeds} random dissections of {width}×{height} square")
    print(f"  Counting GEOMETRIC solutions (up to D4 + piece permutations)")
    print()

    t0 = time.time()
    for seed in range(num_seeds):
        d = random_dissection(width, height, n, seed=seed)
        if d is None:
            continue
        if sum(area2(p) for p in d) != 2 * width * height:
            continue
        total_ok += 1

        geo = count_tilings_geometric(d, width, height, res)
        counts[geo] += 1

        if total_ok % 50 == 0:
            elapsed = time.time() - t0
            print(f"  ...{total_ok} done ({elapsed:.1f}s)", flush=True)

    dt = time.time() - t0
    print(f"\n  {total_ok} valid dissections ({dt:.1f}s)")
    print()

    # Histogram
    print(f"  {'geo_count':>10s}  {'freq':>5s}  {'%':>6s}  histogram")
    print(f"  {'-'*55}")
    for k in sorted(counts.keys()):
        pct = 100 * counts[k] / total_ok
        bar = '#' * min(int(pct * 0.8), 60)
        print(f"  {k:10d}  {counts[k]:5d}  {pct:5.1f}%  {bar}")

    print()
    vals = []
    for k, freq in counts.items():
        vals.extend([k] * freq)
    if vals:
        import statistics
        print(f"  Min:    {min(vals)}")
        print(f"  Max:    {max(vals)}")
        print(f"  Median: {statistics.median(vals)}")
        print(f"  Mean:   {statistics.mean(vals):.2f}")
        if len(vals) > 1:
            print(f"  Stdev:  {statistics.stdev(vals):.2f}")

    return counts


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    run_experiment(n, seeds)
