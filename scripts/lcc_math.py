"""
Vectorized ellipsoidal Lambert Conformal Conic forward projection (Snyder 1987,
"Map Projections: A Working Manual", eqs. 15-1 to 15-9), plus a 2D similarity
(rotate + uniform scale + translate) used to compose the empirical model fitted
by fit_analytic_lcc.py.

Two entry points, because the fit needs both:

  lcc_2sp_forward   the familiar two-standard-parallel form, parameterized by
                    (lat1, lat2). Handles lat1 == lat2 as the tangent limit
                    rather than dividing by zero.
  lcc_forward_cone  the same projection parameterized directly by its cone
                    constant n. This is the form fit_analytic_lcc.py actually
                    fits, because n is the only part of the conic's shape that
                    the data can identify -- see that module's docstring.
"""
from __future__ import annotations

import numpy as np

GRS80_A = 6378137.0
GRS80_INVF = 298.257222101
GRS80_E2 = (2 * (1 / GRS80_INVF)) - (1 / GRS80_INVF) ** 2


def _m(phi, e2):
    return np.cos(phi) / np.sqrt(1 - e2 * np.sin(phi) ** 2)


def _t(phi, e):
    return np.tan(np.pi / 4 - phi / 2) / (
        ((1 - e * np.sin(phi)) / (1 + e * np.sin(phi))) ** (e / 2)
    )


def cone_constant(lat1, lat2, e2=GRS80_E2):
    """Cone constant n of an LCC with standard parallels lat1/lat2 (radians).
    For lat1 == lat2 this is the tangent case, n = sin(lat1)."""
    if np.isclose(lat1, lat2, rtol=0, atol=1e-12):
        return np.sin(lat1)
    e = np.sqrt(e2)
    return (np.log(_m(lat1, e2)) - np.log(_m(lat2, e2))) / (np.log(_t(lat1, e)) - np.log(_t(lat2, e)))


def lcc_forward_cone(lat, lon, lat0, lon0, n, a=GRS80_A, e2=GRS80_E2):
    """LCC forward projection parameterized by the cone constant n directly.

    lat, lon, lat0, lon0 in radians. The projection's scale factor F is taken from
    the tangent parallel asin(n), which is the standard 1SP-equivalent normalization;
    any other choice differs only by a uniform scale, which the similarity applied
    afterwards absorbs exactly. Returns (x, y) relative to the natural origin
    (false easting/northing = 0)."""
    e = np.sqrt(e2)
    phi_t = np.arcsin(np.clip(n, -1.0 + 1e-15, 1.0 - 1e-15))
    F = _m(phi_t, e2) / (n * _t(phi_t, e) ** n)

    rho = a * F * _t(lat, e) ** n
    rho0 = a * F * _t(lat0, e) ** n

    x = rho * np.sin(n * (lon - lon0))
    y = rho0 - rho * np.cos(n * (lon - lon0))
    return x, y


def lcc_2sp_forward(lat, lon, lat0, lat1, lat2, lon0, a=GRS80_A, e2=GRS80_E2):
    """lat, lon, lat0, lat1, lat2, lon0 in radians. Returns (x, y) with false
    easting/northing = 0 (i.e. relative to the projection's natural origin)."""
    e = np.sqrt(e2)
    n = cone_constant(lat1, lat2, e2)
    F = _m(lat1, e2) / (n * _t(lat1, e) ** n)

    rho = a * F * _t(lat, e) ** n
    rho0 = a * F * _t(lat0, e) ** n

    x = rho * np.sin(n * (lon - lon0))
    y = rho0 - rho * np.cos(n * (lon - lon0))
    return x, y


def similarity_2d(x, y, theta, s, tx, ty):
    """Rotate by theta (radians), scale by s, then translate by (tx, ty)."""
    c, sn = np.cos(theta), np.sin(theta)
    X = tx + s * (x * c - y * sn)
    Y = ty + s * (x * sn + y * c)
    return X, Y


def similarity_2d_inverse(X, Y, theta, s, tx, ty):
    c, sn = np.cos(theta), np.sin(theta)
    x0, y0 = (X - tx) / s, (Y - ty) / s
    x = x0 * c + y0 * sn
    y = -x0 * sn + y0 * c
    return x, y


def scale_factor(lat, n, a=GRS80_A, e2=GRS80_E2, F=None):
    """Point scale factor k of an LCC with cone constant n at latitude `lat` (radians).
    k = n·rho / (a·m). Equals exactly 1 at the projection's standard parallel(s)."""
    e = np.sqrt(e2)
    if F is None:
        phi_t = np.arcsin(np.clip(n, -1.0 + 1e-15, 1.0 - 1e-15))
        F = _m(phi_t, e2) / (n * _t(phi_t, e) ** n)
    return n * (a * F * _t(lat, e) ** n) / (a * _m(lat, e2))


if __name__ == "__main__":
    lat0, lon0 = np.radians(42.6678756833333), np.radians(25.5)
    lat, lon = np.radians(np.array([41.5, 43.0, 44.0])), np.radians(np.array([23.0, 25.5, 28.0]))

    # 1) Defining property of the tangent branch: scale factor is exactly 1 at asin(n).
    for deg in (36.0, 42.6678756833333, 43.0, 53.0):
        phi = np.radians(deg)
        k = scale_factor(phi, np.sin(phi))
        print(f"  tangent LCC at phi={deg:8.4f} deg: k(phi) - 1 = {k-1:.3e}")
        assert abs(k - 1) < 1e-13

    # 2) Defining property of the 2SP branch: scale factor is exactly 1 at BOTH parallels.
    l1, l2 = np.radians(42.0), np.radians(43.3333333333333)
    n = cone_constant(l1, l2)
    e = np.sqrt(GRS80_E2)
    F = _m(l1, GRS80_E2) / (n * _t(l1, e) ** n)
    for name, phi in (("lat1", l1), ("lat2", l2)):
        k = scale_factor(phi, n, F=F)
        print(f"  2SP LCC at {name}={np.degrees(phi):8.4f} deg: k - 1 = {k-1:.3e}")
        assert abs(k - 1) < 1e-13

    # 3) The tangent branch is the continuous limit of the 2SP branch. Shrinking the
    #    parallel separation must shrink the positional difference quadratically,
    #    until the 2SP log-difference loses precision to cancellation.
    p_ = np.radians(43.0)
    xa, ya = lcc_2sp_forward(lat, lon, lat0, p_, p_, lon0)
    print("  tangent limit vs 2SP as the parallels close in (expect ~100x per decade):")
    prev, ratios = None, []
    for d in (1e-2, 1e-3, 1e-4, 1e-5):
        xb, yb = lcc_2sp_forward(lat, lon, lat0, p_ - d, p_ + d, lon0)
        diff = np.hypot(xa - xb, ya - yb).max()
        if prev is not None:
            ratios.append(prev / diff)
        print(f"    separation {2*d:.0e} rad -> max diff {diff:.3e} m"
              + (f"   (shrank {ratios[-1]:.0f}x)" if ratios else ""))
        prev = diff
    assert all(90 < r < 110 for r in ratios), f"not quadratic convergence: {ratios}"

    # 4) lcc_forward_cone must agree with lcc_2sp_forward up to a uniform scale.
    x2, y2 = lcc_2sp_forward(lat, lon, lat0, l1, l2, lon0)
    xc, yc = lcc_forward_cone(lat, lon, lat0, lon0, n)
    ratio = np.hypot(x2, y2) / np.hypot(xc, yc)
    print(f"  cone-form vs 2SP form: scale ratio {ratio.min():.12f} .. {ratio.max():.12f} (must be constant)")
    assert np.allclose(ratio, ratio[0], rtol=1e-12)

    # 5) similarity_2d_inverse must undo similarity_2d.
    X, Y = similarity_2d(x2, y2, 0.004, 1.00012, 8.5e6, 4.6e6)
    xr, yr = similarity_2d_inverse(X, Y, 0.004, 1.00012, 8.5e6, 4.6e6)
    print(f"  similarity round-trip: max residual {np.hypot(xr-x2, yr-y2).max():.3e} m")
    assert np.hypot(xr - x2, yr - y2).max() < 1e-6
    print("  OK -- lcc_math self-checks pass.")
