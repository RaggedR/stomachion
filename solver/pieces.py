"""Piece orientation and placement enumeration.

Generates all valid placements of puzzle pieces within a square board.
Each placement is a (piece_index, orientation, translation) with its
rasterized cell set, ready for exact cover encoding.
"""
from __future__ import annotations
from dataclasses import dataclass
from .geometry import (
    Polygon, orientations, normalize, translate, rasterize,
    bounding_box, area2,
)


@dataclass(frozen=True)
class Placement:
    """A specific placement of a piece on the board."""
    piece_idx: int       # which piece
    orient_idx: int      # which orientation (0..7)
    tx: int              # translation x
    ty: int              # translation y
    cells: frozenset     # rasterized cell indices at given resolution

    def __repr__(self):
        return f"Placement(piece={self.piece_idx}, orient={self.orient_idx}, t=({self.tx},{self.ty}), ncells={len(self.cells)})"


def get_orientations(verts: Polygon) -> list[Polygon]:
    """Return distinct orientations of a piece (up to 8 under D4)."""
    return orientations(verts)


def enumerate_placements(
    piece_idx: int,
    verts: Polygon,
    board_size: int = 12,
    res: int = 4,
) -> list[Placement]:
    """Generate all valid placements of a piece within the board.

    For each orientation of the piece, try all integer translations
    that keep the piece fully within [0, board_size] × [0, board_size].
    Rasterize each valid placement to a cell set.

    Returns list of Placement objects.
    """
    orients = get_orientations(verts)
    placements = []
    grid_side = board_size * res

    for oi, oriented in enumerate(orients):
        minx, miny, maxx, maxy = bounding_box(oriented)
        w = maxx - minx
        h = maxy - miny

        # Rasterize once at origin position, then offset for each translation.
        # This avoids re-running PIP for every (tx, ty).
        base_cells = rasterize(oriented, res)
        if not base_cells:
            continue

        # Convert flat indices to (row, col) pairs.
        # The piece is normalized to origin, so these are absolute positions
        # for the piece at (0, 0). Translation just adds (ty*res, tx*res).
        offsets = tuple(
            (c // grid_side, c % grid_side)
            for c in base_cells
        )

        for tx in range(0, board_size - w + 1):
            for ty in range(0, board_size - h + 1):
                # Shift to absolute grid position
                dr = ty * res
                dc = tx * res
                cells = frozenset(
                    (r + dr) * grid_side + (c + dc)
                    for r, c in offsets
                )
                placements.append(Placement(
                    piece_idx=piece_idx,
                    orient_idx=oi,
                    tx=tx,
                    ty=ty,
                    cells=cells,
                ))

    return placements


def enumerate_all_placements(
    piece_verts: list[Polygon],
    board_size: int = 12,
    res: int = 4,
) -> list[Placement]:
    """Enumerate placements for all pieces.

    Returns a flat list of all placements for all pieces.
    """
    all_placements = []
    for i, verts in enumerate(piece_verts):
        pls = enumerate_placements(i, verts, board_size, res)
        all_placements.extend(pls)
    return all_placements


def placement_summary(piece_verts: list[Polygon], board_size: int = 12, res: int = 4):
    """Print a summary of placement counts per piece."""
    total = 0
    for i, verts in enumerate(piece_verts):
        orients = get_orientations(verts)
        pls = enumerate_placements(i, verts, board_size, res)
        print(f"  Piece {i:2d}: {len(orients)} orientations, {len(pls):5d} placements")
        total += len(pls)
    print(f"  Total: {total} placements")
    return total
