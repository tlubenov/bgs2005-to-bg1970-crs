"""Cell-wise triangulation of a ragged lattice."""
import numpy as np

from _util import raises
from bgks1970.lattice import infer_lattice, triangulate_lattice

STEP = 1000.0


def _l_shape():
    have = lambda i, j: not (i >= 3 and j >= 3)
    pts = [(i * STEP, j * STEP) for i in range(6) for j in range(6) if have(i, j)]
    return np.array(pts, dtype=np.float64), have


def test_infers_origin_and_step():
    xy, _ = _l_shape()
    x0, y0, dx, dy = infer_lattice(xy)
    assert (x0, y0, dx, dy) == (0.0, 0.0, STEP, STEP)


def test_covers_the_sampled_region_and_nothing_else():
    """The concave notch must stay uncovered -- that is the whole point of not
    using Delaunay, which would span it with long hull triangles."""
    xy, have = _l_shape()
    tri = triangulate_lattice(xy)

    n_full = n_three = 0
    for i in range(5):
        for j in range(5):
            k = sum(have(*c) for c in ((i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)))
            n_full += k == 4
            n_three += k == 3
    assert len(tri) == 2 * n_full + n_three

    p = xy[tri]
    cross = ((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
             - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))
    assert (cross > 0).all(), "triangles must be counter-clockwise"
    assert np.allclose(np.abs(cross) / 2, STEP * STEP / 2)
    assert np.isclose(np.abs(cross).sum() / 2, (n_full + 0.5 * n_three) * STEP * STEP)


def test_rejects_a_scatter():
    rng = np.random.default_rng(0)
    xy = rng.uniform(0, 10_000, size=(200, 2))
    with raises(ValueError):
        triangulate_lattice(xy)
