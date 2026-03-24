"""Experiment: distribution of total AND irreducible (fault-free) geometric
tilings across random dissections.

For each n from 3 to 10, generates 300 random dissections of a 12x12 square,
solves each with DLX, and computes:
  - Total geometric tilings (mod D4, vertex-based canonical)
  - Fault-free (irreducible) geometric tilings: solutions where NO axis-aligned
    or diagonal fault line exists, deduped under D4

A "fault line" in a solution means ALL pieces in that solution sit entirely
on one side of the line.  We check:
  - Horizontal: y = k*res for k = 1..11
  - Vertical:   x = k*res for k = 1..11
  - Diagonal:   y = x  and  y = 12 - x  (checked via grid-cell centers)
"""
from __future__ import annotations

import sys
import time
import statistics
from collections import Counter

from .explore import random_dissection, area2, _rasterize_rect
from .geometry import Polygon, orientations, normalize, bounding_box
from .pieces import Placement
from .dlx import DLX
from .symmetry import canonical_tiling, canonical_under_d4


# ===================================================================
# Fault-line checking for a single raw solution
# ===================================================================

def solution_has_fault(sol_placements, res=2, gw=24):
    """Check whether a solution has any fault line.

    Checks axis-aligned lines y=k*res, x=k*res for k=1..11,
    plus diagonals y=x and y=12-x (in geometric coords).

    Returns True if ANY fault line exists (i.e. solution is reducible).
    """
    # Axis-aligned fault lines
    for k in range(1, 12):
        rb = k * res
        y_ok = x_ok = True
        for pl in sol_placements:
            if y_ok:
                above = below = False
                for c in pl.cells:
                    if c // gw < rb:
                        above = True
                    else:
                        below = True
                    if above and below:
                        break
                if above and below:
                    y_ok = False
            if x_ok:
                left = right = False
                for c in pl.cells:
                    if c % gw < rb:
                        left = True
                    else:
                        right = True
                    if left and right:
                        break
                if left and right:
                    x_ok = False
            if not y_ok and not x_ok:
                break
        if y_ok or x_ok:
            return True

    # Diagonal fault lines: y=x and y=12-x
    for diag in ('main', 'anti'):
        ok = True
        for pl in sol_placements:
            a = b = False
            for c in pl.cells:
                r, cc = divmod(c, gw)
                gy = (r + 0.5) / res
                gx = (cc + 0.5) / res
                if diag == 'main':
                    v = gy - gx
                else:
                    v = gy - (12 - gx)
                if v > 0:
                    a = True
                elif v < 0:
                    b = True
                if a and b:
                    break
            if a and b:
                ok = False
                break
        if ok:
            return True

    return False


# ===================================================================
# Core: build placements, solve, extract both counts
# ===================================================================

def analyse_dissection(
    pieces: list[Polygon],
    width: int,
    height: int,
    res: int = 2,
    max_solutions: int = 5000,
) -> tuple[int, int, bool]:
    """Solve a dissection and return (total_geo, irred_geo, capped).

    total_geo:  number of geometrically distinct tilings (mod D4)
    irred_geo:  number of fault-free geometrically distinct tilings (mod D4)
    capped:     True if DLX hit max_solutions limit
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
        return 0, 0, False

    # Build and solve DLX
    num_primary = num_pieces + num_cells
    dlx = DLX(num_primary)
    for pl_idx, pl in enumerate(all_placements):
        cols = [pl.piece_idx]
        cols.extend(num_pieces + c for c in pl.cells)
        dlx.add_row(cols, row_id=pl_idx)

    raw_count = dlx.solve(max_solutions=max_solutions)
    capped = (max_solutions > 0 and raw_count >= max_solutions)

    if raw_count == 0:
        return 0, 0, False

    # Extract solutions as lists of Placement objects
    raw_solutions = [
        [all_placements[rid] for rid in sol] for sol in dlx.solutions
    ]

    grid_side = width * res  # square board

    # --- Total geometric count (mod D4) ---
    seen_total = set()
    for sol in raw_solutions:
        key = canonical_tiling(sol)
        canon = canonical_under_d4(key, grid_side)
        seen_total.add(canon)
    total_geo = len(seen_total)

    # --- Fault-free geometric count (mod D4) ---
    seen_irred = set()
    for sol in raw_solutions:
        if not solution_has_fault(sol, res=res, gw=grid_w):
            key = canonical_tiling(sol)
            canon = canonical_under_d4(key, grid_side)
            seen_irred.add(canon)
    irred_geo = len(seen_irred)

    return total_geo, irred_geo, capped


# ===================================================================
# Run experiment for one value of n
# ===================================================================

def run_one_n(
    n: int,
    num_seeds: int = 300,
    width: int = 12,
    height: int = 12,
    res: int = 2,
    max_solutions: int = 5000,
):
    """Run the experiment for n-piece dissections.

    Returns (total_counts, irred_counts, num_valid, num_capped).
    """
    total_counts = []   # geometric count per valid dissection
    irred_counts = []   # fault-free geometric count per valid dissection
    num_valid = 0
    num_capped = 0

    t0 = time.time()

    for seed in range(num_seeds):
        d = random_dissection(width, height, n, seed=seed)
        if d is None:
            continue
        if sum(area2(p) for p in d) != 2 * width * height:
            continue
        num_valid += 1

        total_geo, irred_geo, capped = analyse_dissection(
            d, width, height, res, max_solutions=max_solutions
        )
        total_counts.append(total_geo)
        irred_counts.append(irred_geo)
        if capped:
            num_capped += 1

        if num_valid % 50 == 0:
            elapsed = time.time() - t0
            print(
                f"  n={n}: {num_valid} done / {seed+1} seeds "
                f"({elapsed:.1f}s)",
                flush=True,
            )

    dt = time.time() - t0
    print(
        f"  n={n}: FINISHED — {num_valid} valid from {num_seeds} seeds, "
        f"{num_capped} capped, {dt:.1f}s",
        flush=True,
    )

    return total_counts, irred_counts, num_valid, num_capped


# ===================================================================
# Print distribution for one n
# ===================================================================

def print_distribution(label: str, counts: list[int]):
    """Print a frequency table for a list of counts."""
    if not counts:
        print(f"  {label}: no data")
        return

    freq = Counter(counts)
    total = len(counts)

    print(f"  {label} distribution ({total} dissections):")
    print(f"    {'count':>8s}  {'freq':>5s}  {'%':>6s}  histogram")
    print(f"    {'-'*50}")
    for k in sorted(freq.keys()):
        pct = 100 * freq[k] / total
        bar = '#' * min(int(pct * 0.8), 60)
        print(f"    {k:8d}  {freq[k]:5d}  {pct:5.1f}%  {bar}")

    print(f"    Min={min(counts)}, Max={max(counts)}, "
          f"Median={statistics.median(counts):.1f}, "
          f"Mean={statistics.mean(counts):.2f}")
    if len(counts) > 1:
        print(f"    Stdev={statistics.stdev(counts):.2f}")
    print()


# ===================================================================
# Main
# ===================================================================

def main():
    width = height = 12
    res = 2
    num_seeds = 300
    max_solutions = 5000

    print("=" * 78)
    print("Irreducible Tiling Experiment")
    print(f"  Board: {width}x{height}, res={res}")
    print(f"  Seeds per n: {num_seeds}")
    print(f"  Max raw solutions (DLX cutoff): {max_solutions}")
    print(f"  n range: 3..10")
    print("=" * 78)
    print(flush=True)

    # Collect results for summary table
    summary = {}

    for n in range(3, 11):
        print(f"--- n = {n} ---", flush=True)
        t0 = time.time()

        total_counts, irred_counts, num_valid, num_capped = run_one_n(
            n, num_seeds=num_seeds, width=width, height=height,
            res=res, max_solutions=max_solutions,
        )

        dt = time.time() - t0
        print()

        print_distribution(f"n={n} TOTAL geometric", total_counts)
        print_distribution(f"n={n} IRREDUCIBLE (fault-free) geometric", irred_counts)

        # Store summary stats
        if total_counts:
            t_med = statistics.median(total_counts)
            t_mean = statistics.mean(total_counts)
            t_max = max(total_counts)
        else:
            t_med = t_mean = t_max = 0

        if irred_counts:
            i_med = statistics.median(irred_counts)
            i_mean = statistics.mean(irred_counts)
            i_max = max(irred_counts)
        else:
            i_med = i_mean = i_max = 0

        pct_irred = (100 * i_mean / t_mean) if t_mean > 0 else 0

        summary[n] = {
            'num_valid': num_valid,
            'num_capped': num_capped,
            'total_med': t_med,
            'total_mean': t_mean,
            'total_max': t_max,
            'irred_med': i_med,
            'irred_mean': i_mean,
            'irred_max': i_max,
            'pct_irred': pct_irred,
            'time': dt,
        }

    # ===== Summary table =====
    print()
    print("=" * 100)
    print("SUMMARY TABLE")
    print("=" * 100)
    print()
    print(
        f"{'n':>3s} | {'valid':>5s} | {'cap':>3s} | "
        f"{'total_med':>10s} | {'total_mean':>10s} | {'total_max':>10s} | "
        f"{'irred_med':>10s} | {'irred_mean':>10s} | {'irred_max':>10s} | "
        f"{'%irred_mean':>11s}"
    )
    print(
        f"{'---':>3s}-+-{'-----':>5s}-+-{'---':>3s}-+-"
        f"{'----------':>10s}-+-{'----------':>10s}-+-{'----------':>10s}-+-"
        f"{'----------':>10s}-+-{'----------':>10s}-+-{'----------':>10s}-+-"
        f"{'-----------':>11s}"
    )

    for n in range(3, 11):
        s = summary[n]
        print(
            f"{n:3d} | {s['num_valid']:5d} | {s['num_capped']:3d} | "
            f"{s['total_med']:10.1f} | {s['total_mean']:10.2f} | {s['total_max']:10d} | "
            f"{s['irred_med']:10.1f} | {s['irred_mean']:10.2f} | {s['irred_max']:10d} | "
            f"{s['pct_irred']:10.1f}%"
        )

    print()
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
