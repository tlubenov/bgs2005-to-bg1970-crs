#!/usr/bin/env python3
"""
Fit the best-possible closed-form model for bg_k{n}_{r}_sel (EPSG:7801,
BGS2005/CCS2005) -> bg_k{n}_{r}_sel_tr ("Coordinate system 1970" zone K{n},
parameters unpublished), for each of Bulgaria's four zones (K3/K5/K7/K9).

Why this model shape
---------------------
Bulgaria's "1970 system" defining parameters were never published, and each zone
is described as a conformal conic projection on a *rotated* plane grid -- i.e.
projection plus a small in-plane rotation, not a plain Lambert Conformal Conic.
The exact BGS2005 -> Pulkovo 1942(58) datum shift used internally is unknown too
(the Krasovsky 1940 ellipsoid is taken from every _sel_tr.prj). Rather than guess
a 7-parameter geocentric datum shift -- ill-conditioned over an area as small as
one Bulgarian zone, and largely redundant with a final planar correction anyway --
this script fits an empirical model straight against the real correspondences:

    lon, lat (BGS2005/GRS80, from the _sel layer via EPSG:7801 -> EPSG:7798)
        --> Lambert Conformal Conic, cone constant n
        --> 2D similarity (rotate theta, scale s, translate tx, ty)
        --> E, N  (fitted to match the _sel_tr layer)

What the data can and cannot identify  (IMPORTANT)
--------------------------------------------------
It is tempting to fit the textbook parameter set (lon0, lat1, lat2) plus the
similarity, and report the fitted lon0/lat1/lat2 as "the recovered zone
parameters". That would be wrong, and this project used to do it. Those seven
parameters are NOT identifiable from planar correspondences, because:

  * changing lon0 rotates the projection plane about the cone apex, which the
    similarity's own rotation theta absorbs exactly (only theta - n·lon0 is fixed);
  * changing (lat1, lat2) at constant cone constant n changes only the projection's
    overall scale factor F, which the similarity's scale s absorbs exactly.

Run `--identifiability` to see this measured: starting points scattered across
lon0 = 23°..27.5° converge to parameter sets differing by *degrees*, yet agree on
every coordinate to under two micrometres, and agree on n to nine significant
figures. So the conic shape recovered from this data is exactly one number, the
cone constant n, plus a planar similarity.

This script therefore fits the identifiable 5-parameter form directly:

    free:  n, theta, s, tx, ty
    fixed: lon0 = 25.5 deg (a stated convention -- CCS2005's own central meridian;
           any other choice gives the same coordinates with a different theta)
           lat1 = lat2 = asin(n), the tangent form of the same cone
           lat0 = 42.6678756833333 deg (redundant with ty, so pinned)

and reports n as the recovered quantity. Fewer parameters, no flat directions,
and nothing claimed that the data does not actually determine.

How accuracy is measured
------------------------
The model is fitted on the 1 km set only and scored on the 2 km set, which shares
no ids and no coordinates with it (see bg_points.py). Nothing in the 2 km set
influences the fit, so its residuals are a clean out-of-sample number.

The result will NOT be exact -- it cannot reproduce the local network-adjustment
structure baked into the real historical 1970 system -- so the residuals are
reported in full. For survey-grade work use build_tin_grid.py instead.

Usage:
    python3 fit_analytic_lcc.py                    # fits all four zones
    python3 fit_analytic_lcc.py --zone k5          # fits just one
    python3 fit_analytic_lcc.py --identifiability  # reproduce the argument above
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from pyproj import Transformer

from bg_points import (
    FIT_STEP, QUADRANT, SRC_GEOGRAPHIC_EPSG, SRC_PROJECTED_EPSG, TEST_STEP, ZONES,
    load_matched_points,
)
from lcc_math import (
    GRS80_A, GRS80_E2, cone_constant, lcc_2sp_forward, lcc_forward_cone, similarity_2d,
)
from optim import least_squares_lm

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output"

# Fixed by convention -- see "What the data can and cannot identify" above.
LAT0_FIXED_DEG = 42.6678756833333
LON0_FIXED_DEG = 25.5

# Starting guess per zone. K3/K9 are the "west" zones (their _sel_tr eastings carry
# an 8-prefix), K5/K7 the "east" zones (9-prefix); only tx/ty need to reflect that.
# n starts at sin(lat0), i.e. a cone tangent at the CCS2005 origin parallel.
N0 = np.sin(np.radians(LAT0_FIXED_DEG))
X0_BY_ZONE = {
    "k3": np.array([N0, 0.0, 1.0, 8_500_000.0, 4_700_000.0]),
    "k9": np.array([N0, 0.0, 1.0, 8_500_000.0, 4_550_000.0]),
    "k5": np.array([N0, 0.0, 1.0, 9_500_000.0, 4_600_000.0]),
    "k7": np.array([N0, 0.0, 1.0, 9_500_000.0, 4_750_000.0]),
}


def model_xy(params, lat_rad, lon_rad):
    """The identifiable 5-parameter model: LCC(cone constant n) then similarity."""
    n, theta, s, tx, ty = params
    x, y = lcc_forward_cone(
        lat_rad, lon_rad,
        lat0=np.radians(LAT0_FIXED_DEG),
        lon0=np.radians(LON0_FIXED_DEG),
        n=n, a=GRS80_A, e2=GRS80_E2,
    )
    return similarity_2d(x, y, theta, s, tx, ty)


def residuals(params, lat_rad, lon_rad, E_obs, N_obs):
    E, N = model_xy(params, lat_rad, lon_rad)
    return np.concatenate([E - E_obs, N - N_obs])


def to_geographic(src_xy):
    to_geog = Transformer.from_crs(SRC_PROJECTED_EPSG, SRC_GEOGRAPHIC_EPSG, always_xy=True)
    lon_deg, lat_deg = to_geog.transform(src_xy[:, 0], src_xy[:, 1])
    return np.radians(np.asarray(lat_deg)), np.radians(np.asarray(lon_deg))


def report_stats(label, E_model, N_model, E_obs, N_obs):
    d = np.hypot(E_model - E_obs, N_model - N_obs)
    print(
        f"  {label:<24s} n={len(d):>6d}  RMSE={np.sqrt(np.mean(d**2)):7.3f} m  "
        f"mean={d.mean():7.3f}  median={np.median(d):7.3f}  "
        f"p95={np.percentile(d,95):7.3f}  max={d.max():8.3f}"
    )
    return d


def pipeline_string(params) -> str:
    """PROJ pipeline for the fitted model. Input lon,lat in degrees (BGS2005/GRS80),
    output E,N in the zone's KS1970 plane."""
    n, theta, s, tx, ty = params
    phi = np.degrees(np.arcsin(n))
    c, sn = np.cos(theta), np.sin(theta)
    return (
        f"+proj=pipeline"
        f" +step +proj=lcc +lat_0={LAT0_FIXED_DEG:.10f} +lat_1={phi:.10f} +lat_2={phi:.10f}"
        f" +lon_0={LON0_FIXED_DEG:.10f} +x_0=0 +y_0=0 +ellps=GRS80"
        f" +step +proj=affine +xoff={tx:.6f} +yoff={ty:.6f}"
        f" +s11={s*c:.12f} +s12={-s*sn:.12f} +s21={s*sn:.12f} +s22={s*c:.12f}"
    )


def verify_pipeline(params, lat_rad, lon_rad) -> float:
    """Confirm PROJ's own lcc implementation reproduces this module's model, so the
    published pipeline string really is the thing that was fitted. Returns max
    disagreement in metres."""
    tr = Transformer.from_pipeline(pipeline_string(params))
    E_proj, N_proj = tr.transform(np.degrees(lon_rad), np.degrees(lat_rad))
    E_py, N_py = model_xy(params, lat_rad, lon_rad)
    return float(np.hypot(np.asarray(E_proj) - E_py, np.asarray(N_proj) - N_py).max())


def fit_zone(zone: str) -> dict:
    print(f"\n===== zone {zone} ({QUADRANT[zone]}) =====")
    _, fit_src, fit_dst = load_matched_points(zone, FIT_STEP)
    _, val_src, val_dst = load_matched_points(zone, TEST_STEP)
    print(f"  fit set  ({FIT_STEP}) : {len(fit_src)} matched points")
    print(f"  test set ({TEST_STEP}) : {len(val_src)} matched points, disjoint from the fit set")

    print(f"Converting _sel eastings/northings ({SRC_PROJECTED_EPSG}) "
          f"-> geodetic lon/lat ({SRC_GEOGRAPHIC_EPSG}) ...")
    fit_lat, fit_lon = to_geographic(fit_src)
    val_lat, val_lon = to_geographic(val_src)

    print(f"Fitting (Levenberg-Marquardt, 5 identifiable unknowns) on the {FIT_STEP} set ...")
    res = least_squares_lm(
        residuals, X0_BY_ZONE[zone], args=(fit_lat, fit_lon, fit_dst[:, 0], fit_dst[:, 1]),
    )
    print(f"  {'converged' if res.converged else 'STOPPED'} in {res.n_iter} iterations ({res.message})")

    n, theta, s, tx, ty = res.x
    phi_deg = np.degrees(np.arcsin(n))
    print("Fitted parameters:")
    print(f"  cone constant n              : {n:.9f}")
    print(f"  equivalent tangent parallel  : {phi_deg:.6f} deg  (= asin n)")
    print(f"  central meridian (convention): {LON0_FIXED_DEG:.6f} deg")
    print(f"  latitude of origin (pinned)  : {LAT0_FIXED_DEG:.8f} deg")
    print(f"  post-projection rotation     : {theta:.10f} rad ({np.degrees(theta)*3600:.3f} arcsec)")
    print(f"  post-projection scale        : {s:.10f}  ({(s-1)*1e6:+.2f} ppm)")
    print(f"  post-projection translate    : tx={tx:.3f} m, ty={ty:.3f} m")

    print("Accuracy (Euclidean residual in meters):")
    Ef, Nf = model_xy(res.x, fit_lat, fit_lon)
    report_stats(f"fit set ({FIT_STEP})", Ef, Nf, fit_dst[:, 0], fit_dst[:, 1])
    Ev, Nv = model_xy(res.x, val_lat, val_lon)
    d_test = report_stats(f"out-of-sample ({TEST_STEP})", Ev, Nv, val_dst[:, 0], val_dst[:, 1])

    if np.sqrt(np.mean(d_test**2)) > 1.0:
        print("  NOTE: out-of-sample RMSE is above 1 m -- treat these analytic parameters as "
              "a coarse approximation and use build_tin_grid.py for survey-grade accuracy.")

    pipeline = pipeline_string(res.x)
    dev = verify_pipeline(res.x, val_lat, val_lon)
    print(f"PROJ cross-check: emitted pipeline vs this module's model = {dev*1000:.4f} mm max")
    assert dev < 1e-3, "emitted PROJ pipeline does not reproduce the fitted model"
    print(f"  {pipeline}")

    result = {
        "zone": zone,
        "quadrant": QUADRANT[zone],
        "model": "lcc(cone constant n, lon0/lat0 fixed by convention) + similarity2d(theta,s,tx,ty)",
        "identifiable_parameters": ["n", "theta_rad", "scale", "tx", "ty"],
        "note": (
            "lon0 and the standard parallels are NOT recoverable from planar correspondences: "
            "lon0 trades off exactly against theta, and the standard parallels against scale. "
            "Only the cone constant n is determined. See fit_analytic_lcc.py --identifiability."
        ),
        "cone_constant_n": n,
        "tangent_parallel_deg": phi_deg,
        "lon0_convention_deg": LON0_FIXED_DEG,
        "lat0_fixed_deg": LAT0_FIXED_DEG,
        "theta_rad": theta,
        "scale": s,
        "tx": tx,
        "ty": ty,
        "ellps": "GRS80",
        "input": "lon,lat degrees in BGS2005 geographic (EPSG:7798)",
        "output": f"E,N meters in KS1970 zone K{zone[1:]} (unpublished 'Coordinate system 1970')",
        "proj_pipeline": pipeline,
        "fit_step": FIT_STEP,
        "test_step": TEST_STEP,
        "accuracy_m": {
            "out_of_sample_rmse": float(np.sqrt(np.mean(d_test**2))),
            "out_of_sample_p95": float(np.percentile(d_test, 95)),
            "out_of_sample_max": float(d_test.max()),
            "n_fit": int(len(fit_src)),
            "n_test": int(len(val_src)),
        },
    }

    OUTPUT_DIR.mkdir(exist_ok=True)
    out_path = OUTPUT_DIR / f"lcc_affine_fit_{zone}.json"
    out_path.write_text(json.dumps(result, indent=2))
    print(f"Saved fitted parameters to {out_path}")
    return result


# --------------------------------------------------------------------------------
# Identifiability demonstration
# --------------------------------------------------------------------------------

def _resid7(p, lat, lon, E, N):
    """Residuals of the OVER-parameterized textbook form, used only to demonstrate
    that (lon0, lat1, lat2) are not identifiable."""
    lon0, lat1, lat2, theta, s, tx, ty = p
    x, y = lcc_2sp_forward(lat, lon, np.radians(LAT0_FIXED_DEG),
                           np.radians(lat1), np.radians(lat2), np.radians(lon0),
                           a=GRS80_A, e2=GRS80_E2)
    Em, Nm = similarity_2d(x, y, theta, s, tx, ty)
    return np.concatenate([Em - E, Nm - N])


def _model7(p, lat, lon):
    lon0, lat1, lat2, theta, s, tx, ty = p
    x, y = lcc_2sp_forward(lat, lon, np.radians(LAT0_FIXED_DEG),
                           np.radians(lat1), np.radians(lat2), np.radians(lon0),
                           a=GRS80_A, e2=GRS80_E2)
    return similarity_2d(x, y, theta, s, tx, ty)


def _resid7(p, lat, lon, E, N):
    """Residuals of the OVER-parameterized textbook form, used only to demonstrate
    that (lon0, lat1, lat2) are not identifiable."""
    Em, Nm = _model7(p, lat, lon)
    return np.concatenate([Em - E, Nm - N])


def identifiability_report(zone: str):
    print(f"\n===== identifiability check, zone {zone} =====")
    print("Fitting the over-parameterized 7-parameter form (lon0, lat1, lat2, theta, s, tx, ty)")
    print("from widely scattered starting points, then comparing what they agree on.\n")
    _, fit_src, fit_dst = load_matched_points(zone, FIT_STEP)
    _, val_src, val_dst = load_matched_points(zone, TEST_STEP)
    lat, lon = to_geographic(fit_src)
    vlat, vlon = to_geographic(val_src)
    tx0, ty0 = X0_BY_ZONE[zone][3], X0_BY_ZONE[zone][4]

    starts = [(23.0, 41.5, 43.5), (25.5, 42.6, 42.7), (27.5, 41.0, 44.5)]
    print(f"  {'start lon0':>10s} | {'fitted lon0':>12s} {'lat1':>10s} {'lat2':>10s} "
          f"{'theta(arcsec)':>14s} {'scale-1(ppm)':>13s} | {'cone n':>14s} | {'RMSE(m)':>8s}")
    coords = []
    for lon0, l1, l2 in starts:
        r = least_squares_lm(_resid7, np.array([lon0, l1, l2, 0.0, 1.0, tx0, ty0]),
                             args=(lat, lon, fit_dst[:, 0], fit_dst[:, 1]), max_iter=400)
        p = r.x
        Ev, Nv = _model7(p, vlat, vlon)
        coords.append((np.asarray(Ev), np.asarray(Nv)))
        rmse = np.sqrt(np.mean((Ev - val_dst[:, 0]) ** 2 + (Nv - val_dst[:, 1]) ** 2))
        nn = cone_constant(np.radians(p[1]), np.radians(p[2]))
        print(f"  {lon0:10.1f} | {p[0]:12.6f} {p[1]:10.5f} {p[2]:10.5f} "
              f"{np.degrees(p[3])*3600:14.1f} {(p[4]-1)*1e6:13.2f} | {nn:14.9f} | {rmse:8.4f}")

    print("\n  Pairwise coordinate disagreement over the out-of-sample points:")
    for i in range(1, len(coords)):
        d = np.hypot(coords[i][0] - coords[0][0], coords[i][1] - coords[0][1]).max()
        print(f"    solution 1 vs solution {i+1}: {d*1e6:.3f} micrometres")
    print("\n  => lon0/lat1/lat2 differ by degrees between solutions; the coordinates they")
    print("     produce, and the cone constant n, do not. Only n is recoverable.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zone", choices=ZONES, help="only this zone (default: all four)")
    ap.add_argument("--identifiability", action="store_true",
                    help="demonstrate that lon0/lat1/lat2 are not recoverable, then exit")
    args = ap.parse_args()
    zones = [args.zone] if args.zone else ZONES

    if args.identifiability:
        for z in zones:
            identifiability_report(z)
        return

    results = [fit_zone(z) for z in zones]

    print(f"\n===== summary (out-of-sample accuracy on the {TEST_STEP} set) =====")
    print(f"  {'zone':<5s} {'cone n':>13s} {'tangent phi':>12s} {'rot(\")':>9s} "
          f"{'scale(ppm)':>11s} {'n_test':>7s} {'RMSE(m)':>8s} {'p95(m)':>7s} {'max(m)':>7s}")
    for r in results:
        a = r["accuracy_m"]
        print(f"  {r['zone']:<5s} {r['cone_constant_n']:13.9f} {r['tangent_parallel_deg']:12.6f} "
              f"{np.degrees(r['theta_rad'])*3600:9.1f} {(r['scale']-1)*1e6:11.2f} "
              f"{a['n_test']:7d} {a['out_of_sample_rmse']:8.3f} {a['out_of_sample_p95']:7.3f} "
              f"{a['out_of_sample_max']:7.3f}")


if __name__ == "__main__":
    main()
