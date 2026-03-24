"""Polygon geometry for integer-lattice dissection puzzles.

Core operations for the Stomachion solver: area computation, point-in-polygon
testing, rasterization onto a fine grid, and D4 symmetry transforms. All
coordinates live on a 12x12 integer lattice; arithmetic is exact (integer)
wherever possible.
"""

from __future__ import annotations

Polygon = tuple[tuple[int, int], ...]


def shoelace2(verts: Polygon) -> int:
    """Return 2 x signed area (always integer for integer coords).

    Positive if vertices are CCW, negative if CW.
    Uses the shoelace formula: sum of (x_i * y_{i+1} - x_{i+1} * y_i).
    """
    n = len(verts)
    s = 0
    for i in range(n):
        x0, y0 = verts[i]
        x1, y1 = verts[(i + 1) % n]
        s += x0 * y1 - x1 * y0
    return s


def area2(verts: Polygon) -> int:
    """Return 2 x area (absolute value)."""
    return abs(shoelace2(verts))


def point_in_polygon(px: float, py: float, verts: Polygon) -> bool:
    """Ray-casting point-in-polygon test with consistent boundary assignment.

    Casts a ray in the +x direction from (px, py) and counts crossings.
    The half-open interval convention on edge y-ranges ensures each boundary
    point is assigned to exactly one polygon when tiling.
    """
    n = len(verts)
    inside = False
    j = n - 1
    for i in range(n):
        xi, yi = verts[i]
        xj, yj = verts[j]
        # Edge crosses the ray's y-level if exactly one endpoint is above py.
        # Half-open interval: include lower endpoint, exclude upper.
        if (yi > py) != (yj > py):
            # x-coordinate where edge crosses y = py
            x_cross = xi + (py - yi) * (xj - xi) / (yj - yi)
            if px < x_cross:
                inside = not inside
        j = i
    return inside


def rasterize(verts: Polygon, res: int = 4) -> frozenset[int]:
    """Rasterize polygon to fine-grid cell indices.

    Grid is (12*res) x (12*res). Cell (r, c) has its center at
    ((c + 0.5) / res, (r + 0.5) / res) in polygon-coordinate space.
    Returns frozenset of flat indices r * (12*res) + c.

    Uses a bounding-box pre-filter to avoid testing every cell.
    """
    side = 12 * res
    inv = 1.0 / res
    half = 0.5 * inv

    # Compute bounding box in grid-cell coordinates to limit iteration.
    min_x, min_y, max_x, max_y = bounding_box(verts)
    c_lo = max(int(min_x * res), 0)
    c_hi = min(int(max_x * res), side - 1)
    r_lo = max(int(min_y * res), 0)
    r_hi = min(int(max_y * res), side - 1)

    # Pre-extract vertex coords for the inner loop to avoid repeated
    # tuple unpacking and attribute lookups.
    n = len(verts)
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]

    cells: list[int] = []
    for r in range(r_lo, r_hi + 1):
        py = r * inv + half
        row_offset = r * side
        for c in range(c_lo, c_hi + 1):
            px = c * inv + half
            # Inline ray-cast for performance (avoids function-call overhead
            # on every cell).
            inside = False
            j = n - 1
            for i in range(n):
                yi = ys[i]
                yj = ys[j]
                if (yi > py) != (yj > py):
                    xi = xs[i]
                    xj = xs[j]
                    if px < xi + (py - yi) * (xj - xi) / (yj - yi):
                        inside = not inside
                j = i
            if inside:
                cells.append(row_offset + c)
    return frozenset(cells)


def normalize(verts: Polygon) -> Polygon:
    """Translate so min_x = min_y = 0."""
    min_x = min(v[0] for v in verts)
    min_y = min(v[1] for v in verts)
    return tuple((x - min_x, y - min_y) for x, y in verts)


def rotate_cw(verts: Polygon) -> Polygon:
    """Rotate 90 degrees clockwise: (x, y) -> (max_y - y, x). Normalized."""
    max_y = max(v[1] for v in verts)
    return normalize(tuple((max_y - y, x) for x, y in verts))


def reflect_h(verts: Polygon) -> Polygon:
    """Reflect horizontally: (x, y) -> (max_x - x, y). Normalized."""
    max_x = max(v[0] for v in verts)
    return normalize(tuple((max_x - x, y) for x, y in verts))


def translate(verts: Polygon, tx: int, ty: int) -> Polygon:
    """Translate by (tx, ty)."""
    return tuple((x + tx, y + ty) for x, y in verts)


def _canonical_key(verts: Polygon) -> Polygon:
    """Return a canonical vertex ordering for shape comparison.

    Sorts vertices lexicographically so that two polygons with the same
    vertex set but different starting points / winding orders compare equal.
    """
    return tuple(sorted(verts))


def orientations(verts: Polygon) -> list[Polygon]:
    """Return all distinct orientations (up to 8) under the D4 group.

    Generates 4 rotations of the original and 4 rotations of the horizontal
    reflection, normalizes each, and deduplicates using canonical vertex
    ordering.
    """
    seen: set[Polygon] = set()
    result: list[Polygon] = []

    current = normalize(verts)
    for _ in range(4):
        key = _canonical_key(current)
        if key not in seen:
            seen.add(key)
            result.append(current)
        current = rotate_cw(current)

    current = reflect_h(normalize(verts))
    for _ in range(4):
        key = _canonical_key(current)
        if key not in seen:
            seen.add(key)
            result.append(current)
        current = rotate_cw(current)

    return result


def bounding_box(verts: Polygon) -> tuple[int, int, int, int]:
    """Return (min_x, min_y, max_x, max_y)."""
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    return min(xs), min(ys), max(xs), max(ys)
