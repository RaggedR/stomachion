"""Champion search for maximal prime (irreducible) dissections.

A dissection is "prime" (irreducible) if it has no forced fault lines —
no axis-aligned straight line from boundary to boundary that all pieces
respect.  Equivalently, every candidate cut line y=k or x=k (for integer k
strictly between the boundaries) is crossed by at least one piece.

For each n from 3 to 14, we generate many random dissections, filter for
irreducible ones, count geometric tilings for each, and track the champion
(highest geometric count).
"""
from __future__ import annotations
import sys
import time
from collections import Counter

from .explore import random_dissection, area2
from .experiment import count_tilings_geometric, _klein4_transform
from .geometry import Polygon, orientations, normalize, bounding_box
from .pieces import Placement
from .dlx import DLX
from .explore import _rasterize_rect
from .symmetry import canonical_tiling, canonical_under_d4


# ═══════════════════════════════════════════════════════
# Irreducibility check
# ═══════════════════════════════════════════════════════

def is_irreducible(
    pieces: list[Polygon],
    width: int,
    height: int,
) -> bool:
    """Check whether a dissection is irreducible (has no forced fault lines).

    A fault line is an axis-aligned line y=k (k=1..height-1) or x=k
    (k=1..width-1) such that every piece lies entirely on one side.

    For convex polygons with integer vertices, a piece crosses the line
    y=k iff it has vertices with y < k AND vertices with y > k.  (Vertices
    exactly ON the line don't count — such a piece could be entirely on
    either side.)

    Actually, for a polygon that has vertices strictly on both sides of y=k,
    it must cross the line.  A piece with min_y < k < max_y (strict) crosses.
    But a piece with min_y == k or max_y == k could have all its interior on
    one side.  For integer-vertex convex polygons: min_y < k < max_y (in vertex
    y-coordinates) is necessary and sufficient for crossing.

    We use the simpler and fully correct check: a piece crosses y=k iff
    min(vertex y-coords) < k AND max(vertex y-coords) > k.  This works for
    any polygon (convex or not) as long as it has positive area and integer
    vertices.

    Returns True if irreducible (no fault lines), False if reducible.
    """
    # Precompute y-ranges and x-ranges for each piece
    y_ranges = []
    x_ranges = []
    for poly in pieces:
        ys = [v[1] for v in poly]
        xs = [v[0] for v in poly]
        y_ranges.append((min(ys), max(ys)))
        x_ranges.append((min(xs), max(xs)))

    # Check horizontal fault lines y=k
    for k in range(1, height):
        # Does any piece cross this line?
        crossed = False
        for y_min, y_max in y_ranges:
            if y_min < k < y_max:
                crossed = True
                break
        if not crossed:
            return False

    # Check vertical fault lines x=k
    for k in range(1, width):
        crossed = False
        for x_min, x_max in x_ranges:
            if x_min < k < x_max:
                crossed = True
                break
        if not crossed:
            return False

    return True


# ═══════════════════════════════════════════════════════
# Geometric counting with early termination
# ═══════════════════════════════════════════════════════

def count_tilings_geometric_bounded(
    pieces: list[Polygon],
    width: int,
    height: int,
    res: int = 2,
    max_solutions: int = 0,
) -> int:
    """Count geometrically distinct tilings, with optional early DLX cutoff.

    Like count_tilings_geometric from experiment.py, but passes max_solutions
    to DLX.solve() so we don't explode memory on high-count dissections.

    If max_solutions > 0, the DLX solver stops after finding that many raw
    solutions.  The returned geometric count is then a lower bound (but
    still correct modulo symmetry for the solutions found).

    Returns: geometric tiling count (exact if max_solutions==0 or if
             total raw solutions < max_solutions).
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

    raw_count = dlx.solve(max_solutions=max_solutions)
    if raw_count == 0:
        return 0

    # Extract solutions and count geometric
    raw_solutions = [
        [all_placements[rid] for rid in sol] for sol in dlx.solutions
    ]

    if width == height:
        grid_side = width * res
        seen = set()
        for sol in raw_solutions:
            key = canonical_tiling(sol)
            canon = canonical_under_d4(key, grid_side)
            seen.add(canon)
        return len(seen)
    else:
        grid_side_w = width * res
        grid_side_h = height * res
        seen = set()
        for sol in raw_solutions:
            key = canonical_tiling(sol)
            variants = set()
            for flip_h in (False, True):
                for flip_v in (False, True):
                    transformed = tuple(sorted(
                        _klein4_transform(
                            cells, grid_side_w, grid_side_h, flip_h, flip_v
                        )
                        for cells in key
                    ))
                    variants.add(transformed)
            canon = min(variants)
            seen.add(canon)
        return len(seen)


# ═══════════════════════════════════════════════════════
# Single-n champion finder
# ═══════════════════════════════════════════════════════

def find_champion(
    n: int,
    num_seeds: int = 2000,
    width: int = 12,
    height: int = 12,
    res: int = 2,
    max_solutions: int = 10000,
    verbose: bool = True,
) -> tuple[int, list[Polygon] | None, dict]:
    """Find the irreducible dissection with the most geometric tilings.

    For a given n (number of pieces), generates num_seeds random dissections,
    filters for valid and irreducible ones, counts geometric tilings for each,
    and returns the champion.

    Returns:
        (max_geo_count, champion_pieces, stats_dict)

    stats_dict contains:
        total_valid:     number of valid dissections generated
        total_irreducible: number passing irreducibility filter
        geo_counts:      Counter of geometric counts seen
        champion_seed:   seed that produced the champion
    """
    total_valid = 0
    total_irreducible = 0
    geo_counts = Counter()
    max_geo = 0
    champion_pieces = None
    champion_seed = -1

    t0 = time.time()

    for seed in range(num_seeds):
        d = random_dissection(width, height, n, seed=seed)
        if d is None:
            continue
        if sum(area2(p) for p in d) != 2 * width * height:
            continue
        total_valid += 1

        if not is_irreducible(d, width, height):
            continue
        total_irreducible += 1

        geo = count_tilings_geometric_bounded(
            d, width, height, res, max_solutions=max_solutions
        )
        geo_counts[geo] += 1

        if geo > max_geo:
            max_geo = geo
            champion_pieces = d
            champion_seed = seed
            if verbose:
                elapsed = time.time() - t0
                print(
                    f"  n={n} seed={seed:5d}: NEW CHAMPION geo={geo} "
                    f"({elapsed:.1f}s)",
                    flush=True,
                )

        if verbose and (seed + 1) % 100 == 0:
            elapsed = time.time() - t0
            print(
                f"  n={n} progress: {seed+1}/{num_seeds} seeds, "
                f"{total_valid} valid, {total_irreducible} irred, "
                f"best={max_geo} ({elapsed:.1f}s)",
                flush=True,
            )

    stats = {
        "total_valid": total_valid,
        "total_irreducible": total_irreducible,
        "geo_counts": geo_counts,
        "champion_seed": champion_seed,
    }
    return max_geo, champion_pieces, stats


# ═══════════════════════════════════════════════════════
# Full champion search across n values
# ═══════════════════════════════════════════════════════

def champion_search(
    n_min: int = 3,
    n_max: int = 14,
    seeds_per_n: int = 2000,
    width: int = 12,
    height: int = 12,
    res: int = 2,
    max_solutions: int = 10000,
):
    """Run the champion search for each n in [n_min, n_max].

    Prints progress and a final summary table.
    """
    print("=" * 72)
    print(f"Prime Dissection Champion Search")
    print(f"  Board: {width}x{height}, res={res}")
    print(f"  n range: {n_min}..{n_max}")
    print(f"  Seeds per n: {seeds_per_n}")
    print(f"  Max raw solutions (DLX cutoff): {max_solutions}")
    print("=" * 72)
    print()

    results = {}

    for n in range(n_min, n_max + 1):
        print(f"--- n = {n} ---")
        t0 = time.time()

        max_geo, champion, stats = find_champion(
            n,
            num_seeds=seeds_per_n,
            width=width,
            height=height,
            res=res,
            max_solutions=max_solutions,
            verbose=True,
        )

        dt = time.time() - t0
        results[n] = {
            "max_geo": max_geo,
            "champion": champion,
            "stats": stats,
            "time": dt,
        }

        tv = stats["total_valid"]
        ti = stats["total_irreducible"]
        cs = stats["champion_seed"]
        print(
            f"  n={n} DONE: {tv} valid, {ti} irred, "
            f"champion_geo={max_geo} (seed={cs}), {dt:.1f}s"
        )

        # Mini histogram of the irreducible geometric counts
        gc = stats["geo_counts"]
        if gc:
            print(f"  Distribution (irreducible only):")
            for k in sorted(gc.keys()):
                print(f"    geo={k:6d}: {gc[k]:4d}")
        print()

    # Summary table
    print()
    print("=" * 72)
    print("SUMMARY: Prime Dissection Champions")
    print("=" * 72)
    print()
    print(
        f"{'n':>3s} | {'total':>6s} | {'irred':>6s} | "
        f"{'%irred':>6s} | {'champ_geo':>10s} | {'seed':>6s} | {'time':>7s}"
    )
    print(
        f"{'---':>3s}-+-{'------':>6s}-+-{'------':>6s}-+-"
        f"{'------':>6s}-+-{'----------':>10s}-+-{'------':>6s}-+-{'-------':>7s}"
    )

    for n in range(n_min, n_max + 1):
        r = results[n]
        s = r["stats"]
        tv = s["total_valid"]
        ti = s["total_irreducible"]
        pct = 100 * ti / tv if tv > 0 else 0
        mg = r["max_geo"]
        cs = s["champion_seed"]
        dt = r["time"]
        print(
            f"{n:3d} | {tv:6d} | {ti:6d} | "
            f"{pct:5.1f}% | {mg:10d} | {cs:6d} | {dt:6.1f}s"
        )

    return results


# ═══════════════════════════════════════════════════════
# CLI entry point
# ═══════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Search for champion prime (irreducible) dissections"
    )
    parser.add_argument(
        "--n-min", type=int, default=3,
        help="Minimum number of pieces (default: 3)",
    )
    parser.add_argument(
        "--n-max", type=int, default=14,
        help="Maximum number of pieces (default: 14)",
    )
    parser.add_argument(
        "--seeds", type=int, default=2000,
        help="Number of random seeds per n (default: 2000)",
    )
    parser.add_argument(
        "--width", type=int, default=12,
        help="Board width (default: 12)",
    )
    parser.add_argument(
        "--height", type=int, default=12,
        help="Board height (default: 12)",
    )
    parser.add_argument(
        "--res", type=int, default=2,
        help="Rasterization resolution (default: 2)",
    )
    parser.add_argument(
        "--max-solutions", type=int, default=10000,
        help="DLX early cutoff for raw solutions (default: 10000)",
    )

    args = parser.parse_args()

    champion_search(
        n_min=args.n_min,
        n_max=args.n_max,
        seeds_per_n=args.seeds,
        width=args.width,
        height=args.height,
        res=args.res,
        max_solutions=args.max_solutions,
    )
