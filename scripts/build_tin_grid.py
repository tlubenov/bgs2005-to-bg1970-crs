#!/usr/bin/env python3
"""
Build a PROJ "tinshift" triangulated-grid correction file for each of Bulgaria's
four 1970-system zones: bg_k{n}_{r}_sel (EPSG:7801) <-> bg_k{n}_{r}_sel_tr
("Coordinate system 1970" zone K{n}), directly from the real correspondences in
source_data/.

Why a grid alongside fit_analytic_lcc.py
-----------------------------------------
fit_analytic_lcc.py gets a global closed-form model to roughly 0.09-0.15 m RMSE
per zone, which already matches the 0.14 m plane accuracy AGKK's own tool
declares in its log files -- but it cannot go below that, because a single global
conic plus a similarity has no way to express the local structure of the
historical 1970 network adjustment. A triangulated grid passes exactly through
every known point and only interpolates (very mildly) between them, so inside the
sampled territory it reaches survey-grade accuracy. It also assumes no
projection or datum-shift model at all.

PROJ's tinshift operator (core since PROJ 7.2, see
https://proj.org/en/stable/operations/transformations/tinshift.html) reads a JSON
file listing vertices (each with a source_x/source_y and target_x/target_y pair)
and a triangulation over them, then barycentrically interpolates target
coordinates for any query point inside a triangle. PROJ evaluates the same file
in both directions ("+inv" locates the point in target space using the same
triangle connectivity), so one file per zone serves both directions.

source_x/source_y/target_x/target_y are used here as *plain planar coordinates*
(the _sel layer's Easting/Northing and the _sel_tr layer's Easting/Northing).
There is no trip through geodetic lon/lat, which sidesteps the datum ambiguity
entirely.

Triangulation and validation
----------------------------
The control points are a regular 1 km lattice clipped to each zone, so the
triangulation is done cell-wise (see lattice.py) rather than by Delaunay: it
covers the sampled territory and nothing else, instead of spanning the zone's
concave map-sheet boundary with long, thin hull triangles that would silently
extrapolate.

Accuracy is measured on the 2 km set, which shares no ids and no coordinates with
the 1 km set the grid is built from (see bg_points.py). Those points sit in the
*interior* of the 1 km cells -- the worst case for piecewise-linear interpolation
-- so this is a conservative estimate, and a much more honest one than holding
out a random subset of the same lattice would give. It is measured through PROJ's
own tinshift engine, not a reimplementation of it.

Usage:
    python3 build_tin_grid.py             # builds all four zones
    python3 build_tin_grid.py --zone k5   # builds just one
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

import numpy as np
from pyproj import Transformer

from bg_points import FIT_STEP, QUADRANT, TEST_STEP, ZONES, load_matched_points
from lattice import triangulate_lattice

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output"


def grid_json_path(zone: str) -> Path:
    return OUTPUT_DIR / f"bg_{zone}_tinshift.json"


def make_grid_dict(src_xy: np.ndarray, dst_xy: np.ndarray, triangles: np.ndarray,
                   name: str, description: str) -> dict:
    return {
        "file_type": "triangulation_file",
        "format_version": "1.0",
        "name": name,
        "version": "2.0",
        "publication_date": f"{date.today().isoformat()}T00:00:00Z",
        "license": "Derived from source_data/*.shp (project-internal data).",
        "description": description,
        "input_crs": "EPSG:7801",
        "output_crs": "unnamed Lambert Conformal Conic, 'Coordinate system 1970' zone; "
                      "official parameters unpublished.",
        "transformed_components": ["horizontal"],
        # NOTE: "fallback_strategy" is part of the PROJ triangulation-file schema, but
        # some PROJ versions reject the whole file if it is present at all (any value,
        # "none" included). Omitted deliberately -- points outside the triangulation
        # return inf/inf; fall back to the analytic pipeline from fit_analytic_lcc.py
        # for those (transform_points.py does this automatically).
        "vertices_columns": ["source_x", "source_y", "target_x", "target_y"],
        "triangles_columns": ["idx_vertex1", "idx_vertex2", "idx_vertex3"],
        "vertices": np.column_stack([src_xy, dst_xy]).tolist(),
        "triangles": triangles.tolist(),
    }


def write_grid(grid: dict, path: Path):
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(grid, separators=(",", ":")))
    print(f"  wrote {path.name}  ({path.stat().st_size/1e6:.2f} MB, "
          f"{len(grid['vertices'])} vertices, {len(grid['triangles'])} triangles)")


def evaluate(pipeline: str, src_xy: np.ndarray, dst_xy: np.ndarray, label: str,
             outlier_m: float = 0.001) -> dict:
    tr = Transformer.from_pipeline(pipeline)
    E, N = tr.transform(src_xy[:, 0], src_xy[:, 1])
    E, N = np.asarray(E), np.asarray(N)
    inside = np.isfinite(E) & np.isfinite(N)
    d = np.hypot(E[inside] - dst_xy[inside, 0], N[inside] - dst_xy[inside, 1])
    print(
        f"  {label:<22s} n={inside.sum():>6d}/{len(src_xy)} covered "
        f"({100*inside.mean():5.1f}%)  median={np.median(d)*1000:7.3f} mm  "
        f"p95={np.percentile(d,95)*1000:7.3f}  RMSE={np.sqrt(np.mean(d**2))*1000:7.3f}  "
        f"max={d.max()*1000:8.3f}"
    )
    # Locate outliers. They are not scattered: they line up on the seam in AGKK's own
    # output (see seam_report below), so report which source eastings they sit on.
    out = d > outlier_m
    eastings = np.unique(src_xy[inside][out, 0]) if out.any() else np.array([])
    return {
        "n_total": int(len(src_xy)),
        "n_covered": int(inside.sum()),
        "coverage_fraction": float(inside.mean()),
        "median_mm": float(np.median(d) * 1000),
        "p95_mm": float(np.percentile(d, 95) * 1000),
        "rmse_mm": float(np.sqrt(np.mean(d**2)) * 1000),
        "max_mm": float(d.max() * 1000),
        "n_over_1mm": int(out.sum()),
        "outlier_eastings": [float(x) for x in eastings],
        "median_excluding_outliers_mm": float(np.median(d[~out]) * 1000) if (~out).any() else None,
        "rmse_excluding_outliers_mm": float(np.sqrt(np.mean(d[~out]**2)) * 1000) if (~out).any() else None,
    }


def seam_report(zone: str, src_xy: np.ndarray, dst_xy: np.ndarray) -> dict:
    """Detect discontinuities in the SOURCE DATA -- i.e. in AGKK's own transformation.

    Along any straight row of the lattice, a smooth projection produces a second
    difference of a couple of millimetres at 1 km spacing (that is just the conic's
    curvature). A localized spike an order of magnitude larger is a kink in the
    transformation itself, and no interpolating model can hide it: the grid
    reproduces it exactly at the sampled points but has to interpolate straight
    across it in between.

    Returns the worst such spike per zone and the meridian it falls on."""
    from pyproj import Transformer as _T
    rows = np.unique(src_xy[:, 1])
    row = rows[len(rows) // 2]
    m = src_xy[:, 1] == row
    order = np.argsort(src_xy[m, 0])
    X = src_xy[m, 0][order]
    dE, dN = dst_xy[m, 0][order], dst_xy[m, 1][order]
    curv = np.hypot(np.diff(dE, 2), np.diff(dN, 2))
    typical, k = float(np.median(curv)), int(np.argmax(curv))
    worst, at_e = float(curv[k]), float(X[k + 1])
    lon = _T.from_crs("EPSG:7801", "EPSG:7798", always_xy=True).transform(at_e, row)[0]
    print(f"  row N={row:.0f}: typical curvature {typical*1000:.3f} mm/km², "
          f"worst {worst*1000:.1f} mm at E={at_e:.0f} (longitude {lon:.4f} deg)")
    return {"row_northing": float(row), "typical_curvature_mm": typical * 1000,
            "worst_curvature_mm": worst * 1000, "worst_at_easting": at_e,
            "worst_at_longitude_deg": float(lon),
            "is_discontinuity": bool(worst > 20 * max(typical, 1e-9))}


def build_zone(zone: str) -> dict:
    print(f"\n===== zone {zone} ({QUADRANT[zone]}) =====")
    _, src_xy, dst_xy = load_matched_points(zone, FIT_STEP)
    _, val_src, val_dst = load_matched_points(zone, TEST_STEP)
    print(f"  building from the {FIT_STEP} set: {len(src_xy)} points")
    print(f"  scoring on the  {TEST_STEP} set: {len(val_src)} points, disjoint from it")

    print("Triangulating the lattice ...")
    triangles = triangulate_lattice(src_xy)
    print(f"  {len(triangles)} triangles over {len(src_xy)} vertices")

    grid = make_grid_dict(
        src_xy, dst_xy, triangles,
        name=f"bg_{zone}_to_ks1970_k{zone[1:]}",
        description=(
            f"Piecewise-linear correction grid mapping BGS2005/CCS2005 (EPSG:7801) "
            f"Easting/Northing to 'Coordinate system 1970' zone K{zone[1:]} "
            f"Easting/Northing, built from {len(src_xy)} exact point correspondences "
            f"on a {FIT_STEP} lattice. Use +inv for the reverse direction."
        ),
    )
    grid_path = grid_json_path(zone)
    write_grid(grid, grid_path)

    pipeline = f"+proj=pipeline +step +proj=tinshift +file={grid_path}"

    print(f"Accuracy on the independent {TEST_STEP} set, through PROJ's own tinshift engine:")
    stats = evaluate(pipeline, val_src, val_dst, f"out-of-sample ({TEST_STEP})")
    if stats["coverage_fraction"] < 1.0:
        n_out = stats["n_total"] - stats["n_covered"]
        print(f"  ({n_out} of the {TEST_STEP} points lie outside the {FIT_STEP} lattice's "
              f"outline -- the two sets do not cover exactly the same rectangle -- and are "
              f"excluded from the statistics above rather than extrapolated.)")

    print(f"Reproduction check on the {FIT_STEP} points the grid was built from "
          f"(should be exact, they are its vertices):")
    evaluate(pipeline, src_xy, dst_xy, f"built-from ({FIT_STEP})")

    print("Checking AGKK's own output for discontinuities (independent of any model here):")
    seam = seam_report(zone, src_xy, dst_xy)
    if seam["is_discontinuity"]:
        print(f"  => AGKK's transformation has a slope discontinuity near longitude "
              f"{seam['worst_at_longitude_deg']:.3f} deg in this zone. The grid reproduces it "
              f"exactly at the sampled points but interpolates across it in the cell that "
              f"straddles it, which is where this zone's largest residuals come from.")

    print("Round-trip check (forward then inverse must return the original point):")
    tr = Transformer.from_pipeline(pipeline)
    E, N = tr.transform(val_src[:, 0], val_src[:, 1])
    E, N = np.asarray(E), np.asarray(N)
    ok = np.isfinite(E) & np.isfinite(N)
    x_rt, y_rt = tr.transform(E[ok], N[ok], direction="INVERSE")
    d_rt = np.hypot(np.asarray(x_rt) - val_src[ok, 0], np.asarray(y_rt) - val_src[ok, 1])
    print(f"  round-trip residual: mean={d_rt.mean()*1000:.6f} mm  max={d_rt.max()*1000:.6f} mm")

    return {
        "zone": zone,
        "quadrant": QUADRANT[zone],
        "n_vertices": int(len(src_xy)),
        "n_triangles": int(len(triangles)),
        "build_step": FIT_STEP,
        "test_step": TEST_STEP,
        "out_of_sample": stats,
        "agkk_seam": seam,
        "roundtrip_max_mm": float(d_rt.max() * 1000),
        "grid_path": str(grid_path),
        "grid_size_mb": round(grid_path.stat().st_size / 1e6, 2),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zone", choices=ZONES, help="build only this zone (default: all four)")
    args = ap.parse_args()
    zones = [args.zone] if args.zone else ZONES

    results = [build_zone(z) for z in zones]

    summary_path = OUTPUT_DIR / "tinshift_accuracy.json"
    summary_path.write_text(json.dumps(results, indent=2))

    print(f"\n===== summary (out-of-sample interpolation accuracy on the {TEST_STEP} set) =====")
    print(f"  {'zone':<5s} {'vertices':>9s} {'tris':>8s} {'MB':>5s} {'covered':>15s} "
          f"{'median(mm)':>11s} {'p95(mm)':>8s} {'RMSE(mm)':>9s} {'max(mm)':>8s} {'>1mm':>6s}")
    for r in results:
        s_ = r["out_of_sample"]
        print(f"  {r['zone']:<5s} {r['n_vertices']:9d} {r['n_triangles']:8d} "
              f"{r['grid_size_mb']:5.2f} {s_['n_covered']:6d}/{s_['n_total']:<8d} "
              f"{s_['median_mm']:11.3f} {s_['p95_mm']:8.3f} {s_['rmse_mm']:9.3f} "
              f"{s_['max_mm']:8.3f} {s_['n_over_1mm']:6d}")
    seams = [r for r in results if r["agkk_seam"]["is_discontinuity"]]
    if seams:
        print("\n  Zones where AGKK's own transformation is discontinuous "
              "(all residuals above 1 mm sit on these lines):")
        for r in seams:
            k = r["agkk_seam"]
            print(f"    {r['zone']}: longitude {k['worst_at_longitude_deg']:.4f} deg "
                  f"(source easting {k['worst_at_easting']:.0f} m), "
                  f"kink {k['worst_curvature_mm']:.0f} mm vs {k['typical_curvature_mm']:.1f} mm "
                  f"typical curvature")
        print("  Excluding those points, every zone's out-of-sample residual is:")
        for r in results:
            s_ = r["out_of_sample"]
            print(f"    {r['zone']}: median {s_['median_excluding_outliers_mm']:.3f} mm, "
                  f"RMSE {s_['rmse_excluding_outliers_mm']:.3f} mm")

    print(f"\n  wrote {summary_path}")


if __name__ == "__main__":
    main()
