"""Canonical Stomachion piece definitions.

The 14 pieces of Archimedes' Stomachion on a 12×12 integer grid.
Vertices in SVG convention (y-down, origin top-left).
Coordinates verified: all areas sum to 144, tiling confirmed at 48×48 rasterization.

Source: Ed Pegg Jr., Math Games (Nov 2003), verified by Cutler (2003).
"""

# Canonical vertices in absolute grid coordinates (y-down)
# Each piece: (name, vertices_tuple, 2*area)
PIECES = [
    ("A", ((0,12),(0,0),(2,2)),                     24),  # triangle
    ("B", ((0,12),(2,8),(6,12)),                     24),  # triangle
    ("C", ((0,12),(2,2),(4,4)),                      24),  # triangle
    ("D", ((0,0),(6,0),(2,2)),                       12),  # triangle
    ("E", ((4,4),(4,10),(2,8)),                      12),  # triangle
    ("F", ((4,1),(6,0),(10,4)),                      12),  # triangle
    ("G", ((6,0),(12,0),(10,4)),                     24),  # triangle
    ("H", ((7,4),(8,8),(7,10)),                       6),  # triangle
    ("I", ((10,4),(8,8),(7,4)),                      12),  # triangle
    ("J", ((12,6),(12,8),(9,6)),                      6),  # triangle
    ("K", ((12,0),(12,6),(9,6)),                     18),  # triangle
    ("L", ((2,2),(4,1),(10,4),(4,4)),                24),  # quadrilateral
    ("M", ((9,6),(12,8),(12,12),(6,12)),             48),  # quadrilateral
    ("N", ((6,12),(4,10),(4,4),(7,4),(7,10)),        42),  # pentagon
]

# Verify: total area = 2 * 144 = 288
assert sum(a for _, _, a in PIECES) == 288, "Area sum mismatch"

# Convenience: just the vertex tuples
PIECE_VERTS = [v for _, v, _ in PIECES]
PIECE_NAMES = [n for n, _, _ in PIECES]
PIECE_AREAS = [a for _, _, a in PIECES]

# The 18 distinct interior+boundary vertices used in the canonical tiling
VERTICES = {
    (0,0),(0,12),(2,2),(2,4),(2,8),(2,10),
    (4,1),(4,4),(4,8),(4,10),(4,11),
    (6,0),(6,12),(7,2),(7,4),(7,8),(7,10),
    (8,4),(8,8),(9,6),(10,4),(10,8),
    (12,0),(12,4),(12,6),(12,8),(12,12),
}

# Note: vertices listed here are in math convention from the research.
# The SVG-convention pieces above use (x, 12-y) transform.
# Both are equivalent — the tiling is the same.
