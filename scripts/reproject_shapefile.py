#!/usr/bin/env python3
"""
Reproject a whole point/line/polygon vector file between BGS2005/CCS2005 (EPSG:7801)
and "Coordinate system 1970" zone K{n}, producing a ready-to-open file for QGIS -- no
custom CRS registration needed, since the output coordinates are already in the target
system.

Why this exists instead of `ogr2ogr -ct "+proj=pipeline ... +proj=tinshift ..."`:
`ogr2ogr`'s vector-layer reprojection path does its own internal validity sampling of
the coordinate operation and permanently disables the whole transform object after a
handful of early failures (points right at the grid's domain edge) -- this is a GDAL
robustness rough edge for custom/unregistered operations, not a problem with the
pipeline itself: the lower-level `gdaltransform` CLI and pyproj both apply the exact
same pipeline correctly (verified). This script uses pyproj directly (same engine as
transform_points.py) and GDAL/OGR only for reading and writing vector geometry, so it
sidesteps that GDAL quirk. The pyproj Transformer is built ONCE (loading the grid JSON
a single time) and every coordinate in the file is transformed in one batched call --
building it per-point/per-feature was tried first and was far too slow (each build
re-parses the multi-MB grid file).

Usage:
    python3 reproject_shapefile.py --zone k3 --direction forward \
        --in source_data/bg_k3_1km_sel.shp --out /path/to/bg_k3_as_ks1970.gpkg

    python3 reproject_shapefile.py --zone k5 --direction reverse \
        --in some_ks1970_k5_layer.shp --out reprojected.gpkg

Output CRS: since the true destination projection is unofficial/unpublished, the
output file is written WITHOUT a coordinate reference system (the coordinates are
correct, but the file's own .prj/CRS metadata field is left blank) -- this matches how
source_data/bg_k{n}_{r}_sel_tr.shp itself was delivered, with a .prj naming the datum
and projection type but carrying no numeric parameter. If you want a CRS attached for
on-the-fly reprojection of *other* layers against it, register the custom CRS from
output/bg_k{n}_ks1970_fitted.wkt (build_qgis_crs_wkt.py) in QGIS -- but note that CRS
carries the analytic model, which is ~0.1 m, not the grid used to write this file.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from osgeo import ogr

from bg_points import ZONES
from transform_points import build_analytic_pipeline, build_grid_pipeline
from pyproj import Transformer

DRIVER_BY_EXT = {
    ".shp": "ESRI Shapefile",
    ".gpkg": "GPKG",
    ".geojson": "GeoJSON",
    ".json": "GeoJSON",
}


def collect_points(geom: "ogr.Geometry", out_coords: list):
    """Depth-first walk of a geometry, appending every vertex's (x, y) to
    out_coords and returning a matching nested "shape" descriptor used to
    scatter transformed coordinates back in scatter_points()."""
    n_geoms = geom.GetGeometryCount()
    if n_geoms > 0:
        return [collect_points(geom.GetGeometryRef(i), out_coords) for i in range(n_geoms)]
    n_points = geom.GetPointCount()
    start = len(out_coords)
    for i in range(n_points):
        out_coords.append((geom.GetX(i), geom.GetY(i)))
    return (start, n_points)


def scatter_points(geom: "ogr.Geometry", shape, X: np.ndarray, Y: np.ndarray):
    if isinstance(shape, list):
        for i, sub in enumerate(shape):
            scatter_points(geom.GetGeometryRef(i), sub, X, Y)
        return
    start, n_points = shape
    for i in range(n_points):
        geom.SetPoint_2D(i, float(X[start + i]), float(Y[start + i]))


def batch_transform(x: np.ndarray, y: np.ndarray, zone: str, direction: str, method: str, fallback: bool):
    pyproj_dir = "FORWARD" if direction == "forward" else "INVERSE"
    if method == "grid":
        tr = Transformer.from_pipeline(build_grid_pipeline(zone))
        X, Y = tr.transform(x, y, direction=pyproj_dir)
        X, Y = np.asarray(X, dtype=np.float64), np.asarray(Y, dtype=np.float64)
        bad = ~(np.isfinite(X) & np.isfinite(Y))
        if bad.any():
            print(f"[reproject_shapefile] {bad.sum()} point(s) outside grid coverage", file=sys.stderr)
            if fallback:
                print("  falling back to the analytic model for those.", file=sys.stderr)
                tr_a = Transformer.from_pipeline(build_analytic_pipeline(zone))
                Xa, Ya = tr_a.transform(x[bad], y[bad], direction=pyproj_dir)
                X[bad], Y[bad] = Xa, Ya
        return X, Y
    if method == "analytic":
        tr = Transformer.from_pipeline(build_analytic_pipeline(zone))
        X, Y = tr.transform(x, y, direction=pyproj_dir)
        return np.asarray(X, dtype=np.float64), np.asarray(Y, dtype=np.float64)
    raise ValueError(f"unknown method {method!r}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zone", choices=ZONES, required=True)
    ap.add_argument("--direction", choices=["forward", "reverse"], default="forward",
                     help="forward: CCS2005 -> KS1970 K{n}. reverse: KS1970 K{n} -> CCS2005.")
    ap.add_argument("--method", choices=["grid", "analytic"], default="grid")
    ap.add_argument("--no-fallback", action="store_true")
    ap.add_argument("--in", dest="in_path", type=Path, required=True)
    ap.add_argument("--out", dest="out_path", type=Path, required=True)
    args = ap.parse_args()

    driver_name = DRIVER_BY_EXT.get(args.out_path.suffix.lower())
    if driver_name is None:
        ap.error(f"unrecognized output extension {args.out_path.suffix!r}; use .shp, .gpkg, or .geojson")

    src_ds = ogr.Open(str(args.in_path))
    if src_ds is None:
        ap.error(f"could not open {args.in_path}")
    src_layer = src_ds.GetLayer()
    n = src_layer.GetFeatureCount()

    print(f"Reading {n} features and collecting vertices ...")
    geoms, shapes = [], []
    coords: list = []
    for feat in src_layer:
        geom = feat.GetGeometryRef().Clone()
        geoms.append(geom)
        shapes.append(collect_points(geom, coords))
    x = np.array([c[0] for c in coords], dtype=np.float64)
    y = np.array([c[1] for c in coords], dtype=np.float64)
    print(f"  {len(coords)} vertices total")

    print(f"Transforming ({args.direction}, zone {args.zone}, method {args.method}), Transformer built once ...")
    X, Y = batch_transform(x, y, args.zone, args.direction, args.method, not args.no_fallback)

    for geom, shape in zip(geoms, shapes):
        scatter_points(geom, shape, X, Y)

    if args.out_path.exists():
        args.out_path.unlink()
    driver = ogr.GetDriverByName(driver_name)
    out_ds = driver.CreateDataSource(str(args.out_path))
    out_layer = out_ds.CreateLayer(args.out_path.stem, geom_type=src_layer.GetGeomType())

    src_defn = src_layer.GetLayerDefn()
    for i in range(src_defn.GetFieldCount()):
        out_layer.CreateField(src_defn.GetFieldDefn(i))
    out_defn = out_layer.GetLayerDefn()

    src_layer.ResetReading()
    for feat, geom in zip(src_layer, geoms):
        out_feat = ogr.Feature(out_defn)
        out_feat.SetGeometry(geom)
        for i in range(src_defn.GetFieldCount()):
            out_feat.SetField(src_defn.GetFieldDefn(i).GetName(), feat.GetField(i))
        out_layer.CreateFeature(out_feat)

    out_ds = None
    print(f"Wrote {n} features to {args.out_path} (no CRS assigned -- see module docstring)")


if __name__ == "__main__":
    main()
