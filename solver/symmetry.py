"""Symmetry reduction for tiling solutions.

Given raw DLX solutions, compute the number of geometrically distinct
tilings by modding out:
  1. D4 symmetry of the square (8 elements: 4 rotations × 2 reflections)
  2. Permutations of congruent pieces

A "geometric solution" is a set of regions (polygons) that tile the square.
Two raw solutions are geometrically equivalent if they produce the same
set of regions, regardless of which piece label occupies which region.
"""
from __future__ import annotations


def _cells_to_key(cells) -> tuple[int, ...]:
    """Convert a cell set to a hashable, totally-ordered key."""
    return tuple(sorted(cells))


def canonical_tiling(placements) -> tuple[tuple[int, ...], ...]:
    """Convert a list of Placements into a canonical geometric key.

    The key represents the SET of regions (as sorted cell tuples), sorted
    lexicographically. Independent of piece labels.
    """
    return tuple(sorted(_cells_to_key(pl.cells) for pl in placements))


def d4_transform_cells(cells, grid_side: int, rot: int, flip: bool) -> tuple[int, ...]:
    """Apply a D4 transform to a set of cell indices.

    rot: 0-3 (number of 90° CW rotations)
    flip: horizontal reflection before rotation
    grid_side: side length of the fine grid (e.g. 24 for res=2)

    Returns a sorted tuple of transformed cell indices.
    """
    result = []
    for idx in cells:
        r, c = divmod(idx, grid_side)
        if flip:
            c = grid_side - 1 - c
        for _ in range(rot):
            r, c = c, grid_side - 1 - r
        result.append(r * grid_side + c)
    result.sort()
    return tuple(result)


def d4_transform_tiling(
    tiling_key: tuple[tuple[int, ...], ...],
    grid_side: int,
    rot: int,
    flip: bool,
) -> tuple[tuple[int, ...], ...]:
    """Apply a D4 transform to a canonical tiling key."""
    return tuple(sorted(
        d4_transform_cells(cells, grid_side, rot, flip)
        for cells in tiling_key
    ))


def canonical_under_d4(
    tiling_key: tuple[tuple[int, ...], ...],
    grid_side: int,
) -> tuple[tuple[int, ...], ...]:
    """Return the lexicographically smallest D4 transform of a tiling key."""
    best = tiling_key
    for flip in (False, True):
        for rot in range(4):
            t = d4_transform_tiling(tiling_key, grid_side, rot, flip)
            if t < best:
                best = t
    return best


def count_geometric_solutions(
    raw_solutions: list[list],
    board_size: int = 12,
    res: int = 2,
) -> int:
    """Count geometrically distinct tilings from raw DLX solutions.

    Mods out:
      - Piece label permutations (by using cell sets, not piece IDs)
      - D4 symmetry of the square (by canonicalizing under all 8 transforms)
    """
    grid_side = board_size * res
    seen = set()

    for sol in raw_solutions:
        key = canonical_tiling(sol)
        canon = canonical_under_d4(key, grid_side)
        seen.add(canon)

    return len(seen)
