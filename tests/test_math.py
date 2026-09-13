"""The numeric core, checked by its defining properties rather than golden values."""
import numpy as np

from bgks1970.lcc import (
    GRS80_E2, _m, _t, cone_constant, lcc_2sp_forward, lcc_forward_cone,
    scale_factor, similarity_2d, similarity_2d_inverse,
)
from bgks1970.optim import least_squares_lm

LAT0 = np.radians(42.6678756833333)
LON0 = np.radians(25.5)
LAT = np.radians(np.array([41.5, 43.0, 44.0]))
LON = np.radians(np.array([23.0, 25.5, 28.0]))


def test_tangent_scale_factor_is_exactly_one():
    """A tangent LCC has scale factor 1 at its standard parallel, by definition."""
    for deg in (36.0, 42.6678756833333, 43.0, 53.0):
        phi = np.radians(deg)
        assert abs(scale_factor(phi, np.sin(phi)) - 1) < 1e-13, deg


def test_2sp_scale_factor_is_one_at_both_parallels():
    l1, l2 = np.radians(42.0), np.radians(43.3333333333333)
    n = cone_constant(l1, l2)
    F = _m(l1, GRS80_E2) / (n * _t(l1, np.sqrt(GRS80_E2)) ** n)
    for phi in (l1, l2):
        assert abs(scale_factor(phi, n, F=F) - 1) < 1e-13


def test_tangent_is_the_continuous_limit_of_2sp():
    """Shrinking the parallel separation must shrink the difference quadratically."""
    p = np.radians(43.0)
    xa, ya = lcc_2sp_forward(LAT, LON, LAT0, p, p, LON0)
    prev, ratios = None, []
    for d in (1e-2, 1e-3, 1e-4, 1e-5):
        xb, yb = lcc_2sp_forward(LAT, LON, LAT0, p - d, p + d, LON0)
        diff = np.hypot(xa - xb, ya - yb).max()
        if prev is not None:
            ratios.append(prev / diff)
        prev = diff
    assert all(90 < r < 110 for r in ratios), ratios


def test_cone_form_matches_2sp_up_to_a_uniform_scale():
    l1, l2 = np.radians(42.0), np.radians(43.3333333333333)
    x2, y2 = lcc_2sp_forward(LAT, LON, LAT0, l1, l2, LON0)
    xc, yc = lcc_forward_cone(LAT, LON, LAT0, LON0, cone_constant(l1, l2))
    ratio = np.hypot(x2, y2) / np.hypot(xc, yc)
    assert np.allclose(ratio, ratio[0], rtol=1e-12)


def test_similarity_round_trips():
    x, y = lcc_2sp_forward(LAT, LON, LAT0, np.radians(42.0), np.radians(43.3), LON0)
    X, Y = similarity_2d(x, y, 0.004, 1.00012, 8.5e6, 4.6e6)
    xr, yr = similarity_2d_inverse(X, Y, 0.004, 1.00012, 8.5e6, 4.6e6)
    assert np.hypot(xr - x, yr - y).max() < 1e-6


def test_lm_recovers_a_known_model():
    """Levenberg-Marquardt must recover noise-free synthetic data exactly."""
    truth = np.array([24.5, 41.9, 43.7, 0.004, 1.00013, 8_500_000.0, 4_600_000.0])
    rng = np.random.default_rng(7)
    lat = np.radians(rng.uniform(41.3, 44.2, 2000))
    lon = np.radians(rng.uniform(22.4, 28.6, 2000))

    def model(p, la, lo):
        lon0, lat1, lat2, th, s, tx, ty = p
        x, y = lcc_2sp_forward(la, lo, LAT0, np.radians(lat1), np.radians(lat2),
                               np.radians(lon0))
        return similarity_2d(x, y, th, s, tx, ty)

    E, N = model(truth, lat, lon)

    def resid(p, la, lo, Eo, No):
        Em, Nm = model(p, la, lo)
        return np.concatenate([Em - Eo, Nm - No])

    x0 = np.array([25.5, 42.0, 43.3, 0.0, 1.0, 8_400_000.0, 4_500_000.0])
    res = least_squares_lm(resid, x0, args=(lat, lon, E, N))
    Ef, Nf = model(res.x, lat, lon)
    assert np.hypot(Ef - E, Nf - N).max() < 1e-6
