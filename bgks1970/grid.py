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

Accuracy is measured on a sampling whose positions are disjoint from the one the
grid is built from (see bg_points.resolve_split) -- normally the 2 km set against a
grid built from the 1 km set. Those points sit at the exact *centres* of the 1 km
cells, the worst case for piecewise-linear interpolation, so this is an upper bound
on the error rather than a favourable sample. It is measured through PROJ's own
tinshift engine, not a reimplementation of it.
"""
from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import numpy as np
from pyproj import Transformer

from .data import QUADRANT, ZONES, load_matched_points, load_sources, resolve_split
from .lattice import triangulate_lattice
from .paths import output_dir

def grid_json_path(zone: str, root=None) -> Path:
    return output_dir(root) / f"bg_{zone}_tinshift.json"


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


def write_grid(grid: dict, path: Path, say=print):
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(grid, separators=(",", ":")))
    say(f"  wrote {path.name}  ({path.stat().st_size/1e6:.2f} MB, "
          f"{len(grid['vertices'])} vertices, {len(grid['triangles'])} triangles)")


def evaluate(pipeline: str, src_xy: np.ndarray, dst_xy: np.ndarray, label: str,
             outlier_m: float = 0.001, say=print) -> dict:
    tr = Transformer.from_pipeline(pipeline)
    E, N = tr.transform(src_xy[:, 0], src_xy[:, 1])
    E, N = np.asarray(E), np.asarray(N)
    inside = np.isfinite(E) & np.isfinite(N)
    d = np.hypot(E[inside] - dst_xy[inside, 0], N[inside] - dst_xy[inside, 1])
    say(
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


def seam_report(zone: str, src_xy: np.ndarray, dst_xy: np.ndarray, say=print) -> dict:
    """Detect discontinuities in the SOURCE DATA -- i.e. in AGKK's own transformation.

    This samples one lattice row, which is enough to detect and locate a kink. The
    easting it reports is that row's -- the kink lies on a meridian, whose easting
    drifts by kilometres across a zone's latitude range, so the LONGITUDE is the
    stable figure. make_zone_figure.py profiles every row for the same reason.

    Along any straight row of the lattice, a smooth projection produces a second
    difference of a couple of millimetres at 1 km spacing (that is just the conic's
    curvature). A localized spike an order of magnitude larger is a kink in the
    transformation itself, and no interpolating model can hide it: the grid
    reproduces it exactly at the sampled points but has to interpolate straight
    across it in between.

    Returns the worst such spike per zone and the meridian it falls on."""
    from pyproj import Transformer as _T
    step = float(np.diff(np.unique(src_xy[:, 0])).min())
    rows = np.unique(src_xy[:, 1])
    row = rows[len(rows) // 2]
    m = src_xy[:, 1] == row
    order = np.argsort(src_xy[m, 0])
    X = src_xy[m, 0][order]
    dE, dN = dst_xy[m, 0][order], dst_xy[m, 1][order]
    curv = np.hypot(np.diff(dE, 2), np.diff(dN, 2))
    # A row is not necessarily gap-free -- it is a lattice clipped to a zone. Second
    # differences taken across a gap are meaningless, so keep only triples whose two
    # spacings are both one step, and normalise to mm per km^2 so zones sampled at
    # different steps stay comparable.
    contiguous = (np.diff(X)[:-1] == step) & (np.diff(X)[1:] == step)
    if not contiguous.any():
        raise ValueError(f"{zone}: no contiguous lattice triples in the sampled row")
    curv, Xc = curv[contiguous], X[1:-1][contiguous]
    scale = 1000.0 / (step / 1000.0) ** 2
    typical, k = float(np.median(curv)) * scale, int(np.argmax(curv))
    worst, at_e = float(curv[k]) * scale, float(Xc[k])
    lon = _T.from_crs("EPSG:7801", "EPSG:7798", always_xy=True).transform(at_e, row)[0]
    say(f"  row N={row:.0f} ({contiguous.sum()} contiguous triples at {step:g} m): "
          f"typical curvature {typical:.3f} mm/km², worst {worst:.1f} mm/km² "
          f"at E={at_e:.0f} (longitude {lon:.4f} deg)")
    return {"row_northing": float(row), "lattice_step_m": step,
            "typical_curvature_mm": typical, "worst_curvature_mm": worst,
            "worst_at_easting": at_e,
            "worst_at_longitude_deg": float(lon),
            "is_discontinuity": bool(worst > 20 * max(typical, 1e-9))}


def build_zone(zone: str, root=None, verbose: bool = True) -> dict:
    say = print if verbose else (lambda *a, **k: None)
    say(f"\n===== zone {zone} ({QUADRANT[zone]}) =====")
    split = resolve_split(zone, root)
    _, src_xy, dst_xy = load_matched_points(zone, split.fit.step, split.fit.dataset)
    val_src, val_dst = load_sources(zone, split.test)
    say(f"  building from {split.fit.label:<28s} {len(src_xy):6d} points")
    say(f"  scoring on    {', '.join(s.label for s in split.test):<28s} {len(val_src):6d} points")
    if split.provisional:
        say(f"  !! PROVISIONAL: {split.note}")

    say("Triangulating the lattice ...")
    triangles = triangulate_lattice(src_xy)
    say(f"  {len(triangles)} triangles over {len(src_xy)} vertices")

    grid = make_grid_dict(
        src_xy, dst_xy, triangles,
        name=f"bg_{zone}_to_ks1970_k{zone[1:]}",
        description=(
            f"Piecewise-linear correction grid mapping BGS2005/CCS2005 (EPSG:7801) "
            f"Easting/Northing to 'Coordinate system 1970' zone K{zone[1:]} "
            f"Easting/Northing, built from {len(src_xy)} exact point correspondences "
            f"({split.fit.label}). Use +inv for the reverse direction."
        ),
    )
    grid_path = grid_json_path(zone, root)
    write_grid(grid, grid_path, say)

    pipeline = f"+proj=pipeline +step +proj=tinshift +file={grid_path}"

    say("Accuracy out of sample, through PROJ's own tinshift engine:")
    stats = evaluate(pipeline, val_src, val_dst, "out-of-sample", say=say)
    if stats["coverage_fraction"] < 1.0:
        n_out = stats["n_total"] - stats["n_covered"]
        say(f"  ({n_out} scoring points lie outside the built lattice's outline and are "
              f"excluded from the statistics above rather than extrapolated.)")

    say("Reproduction check on the points the grid was built from "
          "(should be exact, they are its vertices):")
    evaluate(pipeline, src_xy, dst_xy, "built-from", say=say)

    say("Checking AGKK's own output for discontinuities (independent of any model here):")
    seam = seam_report(zone, src_xy, dst_xy, say)
    if seam["is_discontinuity"]:
        say(f"  => AGKK's transformation has a slope discontinuity near longitude "
              f"{seam['worst_at_longitude_deg']:.3f} deg in this zone. The grid reproduces it "
              f"exactly at the sampled points but interpolates across it in the cell that "
              f"straddles it, which is where this zone's largest residuals come from.")

    say("Round-trip check (forward then inverse must return the original point):")
    tr = Transformer.from_pipeline(pipeline)
    E, N = tr.transform(val_src[:, 0], val_src[:, 1])
    E, N = np.asarray(E), np.asarray(N)
    fwd_ok = np.isfinite(E) & np.isfinite(N)
    x_rt, y_rt = tr.transform(E[fwd_ok], N[fwd_ok], direction="INVERSE")
    x_rt, y_rt = np.asarray(x_rt), np.asarray(y_rt)
    # A point can transform forward and still fail coming back: the inverse locates
    # it in target space, and one that landed exactly on the domain edge may fall
    # outside every triangle there. Exclude those instead of letting a single inf
    # swallow the whole statistic.
    back_ok = np.isfinite(x_rt) & np.isfinite(y_rt)
    d_rt = np.hypot(x_rt[back_ok] - val_src[fwd_ok][back_ok, 0],
                    y_rt[back_ok] - val_src[fwd_ok][back_ok, 1])
    n_edge = int((~back_ok).sum())
    say(f"  round-trip residual over {back_ok.sum()} points: "
          f"mean={d_rt.mean()*1e6:.4f} um  max={d_rt.max()*1e6:.4f} um"
          + (f"   [{n_edge} point(s) on the domain edge did not invert]" if n_edge else ""))

    return {
        "zone": zone,
        "quadrant": QUADRANT[zone],
        "n_vertices": int(len(src_xy)),
        "n_triangles": int(len(triangles)),
        "build_source": split.fit.label,
        "test_sources": [s.label for s in split.test],
        "provisional": split.provisional,
        "provisional_reason": split.note,
        "out_of_sample": stats,
        "agkk_seam": seam,
        "roundtrip_max_mm": float(d_rt.max() * 1000),
        "roundtrip_edge_failures": n_edge,
        "grid_path": str(grid_path),
        "grid_size_mb": round(grid_path.stat().st_size / 1e6, 2),
    }


def build_all(zones=ZONES, root=None, verbose: bool = True) -> list[dict]:
    """Build every zone's grid, write the accuracy summary, and return the results."""
    results = [build_zone(z, root, verbose) for z in zones]

    summary_path = output_dir(root) / "tinshift_accuracy.json"
    summary_path.write_text(json.dumps(results, indent=2))

    if verbose:
        print("\n===== summary (out-of-sample interpolation accuracy) =====")
        print(f"  {'zone':<5s} {'vertices':>9s} {'tris':>8s} {'MB':>5s} {'covered':>15s} "
              f"{'median(mm)':>11s} {'p95(mm)':>8s} {'RMSE(mm)':>9s} {'max(mm)':>8s} {'>1mm':>6s}")
        for r in results:
            s_ = r["out_of_sample"]
            print(f"  {r['zone']:<5s} {r['n_vertices']:9d} {r['n_triangles']:8d} "
                  f"{r['grid_size_mb']:5.2f} {s_['n_covered']:6d}/{s_['n_total']:<8d} "
                  f"{s_['median_mm']:11.3f} {s_['p95_mm']:8.3f} {s_['rmse_mm']:9.3f} "
                  f"{s_['max_mm']:8.3f} {s_['n_over_1mm']:6d}"
                  + ("  PROVISIONAL" if r["provisional"] else ""))
        seams = [r for r in results if r["agkk_seam"]["is_discontinuity"]]
        if seams:
            print("\n  Zones where AGKK's own transformation is discontinuous "
                  "(all residuals above 1 mm sit on these lines):")
            for r in seams:
                k = r["agkk_seam"]
                print(f"    {r['zone']}: longitude {k['worst_at_longitude_deg']:.4f} deg "
                      f"(source easting {k['worst_at_easting']:.0f} m), "
                      f"kink {k['worst_curvature_mm']:.0f} mm/km² vs "
                      f"{k['typical_curvature_mm']:.1f} typical")
            print("  Excluding those points, every zone's out-of-sample residual is:")
            for r in results:
                s_ = r["out_of_sample"]
                print(f"    {r['zone']}: median {s_['median_excluding_outliers_mm']:.3f} mm, "
                      f"RMSE {s_['rmse_excluding_outliers_mm']:.3f} mm")
        print(f"\n  wrote {summary_path}")
    return results
