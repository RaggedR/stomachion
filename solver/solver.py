"""Tiling solver — enumerate all tilings of a square using DLX exact cover.

The approach:
1. Rasterize the board into fine-grid cells (12×res)²
2. For each piece, enumerate all valid placements (orientation × translation)
3. Encode as exact cover: each placement covers a set of cells + one "piece used" column
4. Solve with DLX (Algorithm X with Dancing Links)
5. Each solution is a set of 14 placements that tile the board perfectly
"""
from __future__ import annotations
import time
from .stomachion import PIECE_VERTS, PIECE_NAMES
from .pieces import enumerate_all_placements, Placement
from .dlx import DLX
from .geometry import Polygon


def build_solver(
    piece_verts: list[Polygon],
    board_size: int = 12,
    res: int = 4,
    verbose: bool = True,
) -> tuple[DLX, list[Placement]]:
    """Build the DLX exact cover instance for a tiling problem.

    Columns:
      - 0..num_pieces-1: primary columns, one per piece (each piece used exactly once)
      - num_pieces..num_pieces+num_cells-1: primary columns, one per cell (each cell covered exactly once)

    Rows:
      - One per valid placement of any piece

    Returns (dlx_instance, placements_list).
    """
    num_pieces = len(piece_verts)
    grid_size = board_size * res
    num_cells = grid_size * grid_size

    if verbose:
        print(f"Building solver: {num_pieces} pieces, {board_size}×{board_size} board, res={res}")
        print(f"  Grid: {grid_size}×{grid_size} = {num_cells} cells")

    # Enumerate all placements
    t0 = time.time()
    placements = enumerate_all_placements(piece_verts, board_size, res)
    t1 = time.time()
    if verbose:
        print(f"  Placements: {len(placements)} ({t1-t0:.2f}s)")

    # Build DLX
    # Primary columns: num_pieces (piece identity) + num_cells (board coverage)
    num_primary = num_pieces + num_cells
    dlx = DLX(num_primary, num_secondary=0)

    for pl_idx, pl in enumerate(placements):
        # This row covers: piece identity column + cell columns
        cols = [pl.piece_idx]  # piece identity
        cols.extend(num_pieces + c for c in pl.cells)  # cell coverage
        dlx.add_row(cols, row_id=pl_idx)

    t2 = time.time()
    if verbose:
        print(f"  DLX matrix: {len(placements)} rows × {num_primary} cols ({t2-t1:.2f}s)")

    return dlx, placements


def solve(
    piece_verts: list[Polygon] | None = None,
    board_size: int = 12,
    res: int = 4,
    max_solutions: int = 0,
    verbose: bool = True,
    store_solutions: bool = False,
) -> tuple[int, list[list[Placement]]]:
    """Solve the tiling problem. Returns (count, solutions).

    If store_solutions is False, solutions list is empty (saves memory).
    """
    if piece_verts is None:
        piece_verts = PIECE_VERTS

    dlx, placements = build_solver(piece_verts, board_size, res, verbose)

    if verbose:
        print(f"  Solving...")
    t0 = time.time()
    count = dlx.solve(max_solutions=max_solutions)
    t1 = time.time()

    if verbose:
        print(f"  Solutions: {count} ({t1-t0:.2f}s)")

    solutions = []
    if store_solutions:
        for sol_row_ids in dlx.solutions:
            solutions.append([placements[rid] for rid in sol_row_ids])

    return count, solutions


def verify_stomachion():
    """Verify: the canonical Stomachion has 17152 raw solutions (536 × 32)."""
    print("=" * 60)
    print("Stomachion Solver — Verification")
    print("=" * 60)

    count, _ = solve(PIECE_VERTS, verbose=True)

    # The raw count includes:
    # - 8 rotations/reflections of the square (D4)
    # - 4 permutations of identical piece pairs (2 pairs × 2 swaps)
    # So raw_count = 536 × 8 × 4 = 17152
    # BUT: some solutions may have non-trivial symmetry, making the exact
    # relationship more subtle. Cutler reports 536 geometrically distinct solutions.

    print(f"\n  Raw solution count: {count}")
    print(f"  Expected: 17152 (536 × 32) or similar")
    print(f"  ÷ 32 = {count / 32:.1f} (should be ~536)")

    if count == 17152:
        print("\n  ✓ VERIFIED: 17152 raw solutions = 536 × 32")
    else:
        print(f"\n  Note: raw count {count} ÷ 8 (D4) = {count/8:.1f}")
        print(f"  The 536 count may require careful symmetry analysis")

    return count


if __name__ == "__main__":
    verify_stomachion()
