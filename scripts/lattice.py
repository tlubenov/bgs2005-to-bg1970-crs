"""
Triangulation of a regular (but ragged-edged) point lattice, in plain numpy.

The control points in source_data/ are not a scatter -- they are a clean square
lattice (1 km or 2 km step) clipped to each zone's territory, so they need no
general Delaunay routine. Snapping them back to their integer lattice indices
and cutting each occupied cell into triangles is both cheaper and *better* than
Delaunay for this job:

  * Delaunay triangulates the convex hull. The zones are not convex (their
    boundary is the staircase edge of the map-sheet grid), so Delaunay spans the
    concavities with long, thin triangles that interpolate across territory that
    was never sampled -- which is exactly where the previous version of this
    project saw its worst hold-out errors. Cell-wise triangulation covers the
    sampled region and nothing else, so a query outside the data is reported as
    outside the data instead of being silently extrapolated.
  * Every triangle is half of a lattice cell, so no slivers, and the
    interpolation error is bounded by the curvature of the transformation over
    one cell rather than over an arbitrarily long hull edge.

Cells with all four corners present yield two triangles; cells with exactly
three present yield one, which recovers the diagonal steps of the boundary.
"""
from __future__ import annotations

import numpy as np


def infer_lattice(xy: np.ndarray) -> tuple[float, float, float, float]:
    """Return (x0, y0, dx, dy): the lattice origin and step implied by `xy`.
    The step is the smallest positive difference between distinct sorted
    coordinate values along each axis."""
    out = []
    for axis in (0, 1):
        v = np.unique(xy[:, axis])
        if v.size < 2:
            raise ValueError(f"axis {axis} has fewer than two distinct coordinates")
        d = np.diff(v)
        out.append((v[0], d.min()))
    (x0, dx), (y0, dy) = out
    return x0, y0, dx, dy


def lattice_indices(xy: np.ndarray, x0: float, y0: float, dx: float, dy: float,
                    tol: float = 1e-6) -> tuple[np.ndarray, np.ndarray]:
    """Map each point to its integer (i, j) lattice cell index, verifying that it
    really sits on the lattice."""
    fi = (xy[:, 0] - x0) / dx
    fj = (xy[:, 1] - y0) / dy
    i, j = np.rint(fi).astype(np.int64), np.rint(fj).astype(np.int64)
    off = np.maximum(np.abs(fi - i), np.abs(fj - j))
    if off.max() > tol:
        raise ValueError(
            f"{int((off > tol).sum())} point(s) are not on the inferred "
            f"{dx:g}x{dy:g} lattice (worst offset {off.max()*max(dx,dy):.4f} m)"
        )
    return i, j


def triangulate_lattice(xy: np.ndarray) -> np.ndarray:
    """Triangulate a regular point lattice with a ragged outline.

    Returns an (M, 3) int array of vertex indices into `xy`, counter-clockwise,
    covering every lattice cell that has at least three of its four corners
    sampled. Raises if `xy` is not on a regular lattice."""
    x0, y0, dx, dy = infer_lattice(xy)
    i, j = lattice_indices(xy, x0, y0, dx, dy)

    ni, nj = i.max() + 2, j.max() + 2
    grid = np.full((ni, nj), -1, dtype=np.int64)
    grid[i, j] = np.arange(len(xy), dtype=np.int64)

    # Corners of every candidate cell (i, j): a=(i,j) b=(i+1,j) c=(i,j+1) d=(i+1,j+1)
    a = grid[:-1, :-1]
    b = grid[1:, :-1]
    c = grid[:-1, 1:]
    d = grid[1:, 1:]

    tris: list[np.ndarray] = []

    full = (a >= 0) & (b >= 0) & (c >= 0) & (d >= 0)
    if full.any():
        av, bv, cv, dv = a[full], b[full], c[full], d[full]
        tris.append(np.column_stack([av, bv, dv]))
        tris.append(np.column_stack([av, dv, cv]))

    # Exactly three corners present -> one triangle, wound counter-clockwise.
    present = (a >= 0).astype(np.int8) + (b >= 0) + (c >= 0) + (d >= 0)
    three = present == 3
    for missing, corners in (
        (d, (a, b, c)),   # a,b,c present: CCW is a -> b -> c
        (c, (a, b, d)),   # a,b,d present: CCW is a -> b -> d
        (b, (a, d, c)),   # a,c,d present: CCW is a -> d -> c
        (a, (b, d, c)),   # b,c,d present: CCW is b -> d -> c
    ):
        m = three & (missing < 0)
        if m.any():
            tris.append(np.column_stack([corners[0][m], corners[1][m], corners[2][m]]))

    if not tris:
        raise ValueError("no lattice cell had three or more sampled corners")
    return np.vstack(tris)


if __name__ == "__main__":
    # Self-test on an L-shaped lattice: the concave corner must stay uncovered,
    # and every emitted triangle must be counter-clockwise with the right area.
    step = 1000.0
    have = lambda ii, jj: not (ii >= 3 and jj >= 3)
    pts = [(ii * step, jj * step) for ii in range(6) for jj in range(6) if have(ii, jj)]
    xy = np.array(pts, dtype=np.float64)
    tri = triangulate_lattice(xy)

    # Independently count what the cell census says we should get.
    n_full = n_three = 0
    for ii in range(5):
        for jj in range(5):
            k = sum(have(*c) for c in ((ii, jj), (ii + 1, jj), (ii, jj + 1), (ii + 1, jj + 1)))
            n_full += k == 4
            n_three += k == 3
    expect_tris = 2 * n_full + n_three
    expect_area = (n_full + 0.5 * n_three) * step * step

    p3 = xy[tri]
    cross = ((p3[:, 1, 0] - p3[:, 0, 0]) * (p3[:, 2, 1] - p3[:, 0, 1])
             - (p3[:, 2, 0] - p3[:, 0, 0]) * (p3[:, 1, 1] - p3[:, 0, 1]))
    assert (cross > 0).all(), "some triangles are clockwise"
    assert np.allclose(np.abs(cross) / 2, step * step / 2), "unexpected triangle areas"
    total = np.abs(cross).sum() / 2

    print(f"  {len(xy)} vertices -> {len(tri)} triangles (expected {expect_tris}: "
          f"{n_full} full cells + {n_three} three-corner cells), all CCW")
    print(f"  covered area {total/1e6:.1f} km² (expected {expect_area/1e6:.1f} km²); "
          f"the 9-cell concave notch is correctly left uncovered")
    assert len(tri) == expect_tris and np.isclose(total, expect_area)
    print("  OK -- triangulate_lattice covers the sampled region exactly, concavity excluded.")
