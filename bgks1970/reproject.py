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
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from osgeo import ogr
from pyproj import Transformer

from .transform import build_analytic_pipeline, build_grid_pipeline

ogr.UseExceptions()

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


def batch_transform(x: np.ndarray, y: np.ndarray, zone: str, direction: str,
                    method: str, fallback: bool, root=None):
    pyproj_dir = "FORWARD" if direction == "forward" else "INVERSE"
    if method == "grid":
        tr = Transformer.from_pipeline(build_grid_pipeline(zone, root))
        X, Y = tr.transform(x, y, direction=pyproj_dir)
        X, Y = np.asarray(X, dtype=np.float64), np.asarray(Y, dtype=np.float64)
        bad = ~(np.isfinite(X) & np.isfinite(Y))
        if bad.any():
            print(f"[reproject_shapefile] {bad.sum()} point(s) outside grid coverage", file=sys.stderr)
            if fallback:
                print("  falling back to the analytic model for those.", file=sys.stderr)
                tr_a = Transformer.from_pipeline(build_analytic_pipeline(zone, root))
                Xa, Ya = tr_a.transform(x[bad], y[bad], direction=pyproj_dir)
                X[bad], Y[bad] = Xa, Ya
        return X, Y
    if method == "analytic":
        tr = Transformer.from_pipeline(build_analytic_pipeline(zone, root))
        X, Y = tr.transform(x, y, direction=pyproj_dir)
        return np.asarray(X, dtype=np.float64), np.asarray(Y, dtype=np.float64)
    raise ValueError(f"unknown method {method!r}")


def reproject_file(in_path: Path, out_path: Path, zone: str, direction: str = "forward",
                   method: str = "grid", fallback: bool = True, root=None,
                   verbose: bool = True) -> int:
    """Reproject every vertex of a vector file. Returns the feature count written.

    The pyproj Transformer is built once and every coordinate goes through one
    batched call -- building it per feature re-parses the multi-MB grid each time."""
    say = print if verbose else (lambda *a, **k: None)
    driver_name = DRIVER_BY_EXT.get(Path(out_path).suffix.lower())
    if driver_name is None:
        raise ValueError(f"unrecognized output extension {Path(out_path).suffix!r}; "
                         f"use one of {sorted(DRIVER_BY_EXT)}")

    src_ds = ogr.Open(str(in_path))
    if src_ds is None:
        raise FileNotFoundError(in_path)
    src_layer = src_ds.GetLayer()
    n = src_layer.GetFeatureCount()

    say(f"Reading {n} features and collecting vertices ...")
    geoms, shapes, coords = [], [], []
    for feat in src_layer:
        geom = feat.GetGeometryRef().Clone()
        geoms.append(geom)
        shapes.append(collect_points(geom, coords))
    x = np.array([c[0] for c in coords], dtype=np.float64)
    y = np.array([c[1] for c in coords], dtype=np.float64)
    say(f"  {len(coords)} vertices total")

    say(f"Transforming ({direction}, zone {zone}, method {method}) ...")
    X, Y = batch_transform(x, y, zone, direction, method, fallback, root)
    for geom, shape in zip(geoms, shapes):
        scatter_points(geom, shape, X, Y)

    out_path = Path(out_path)
    if out_path.exists():
        out_path.unlink()
    out_ds = ogr.GetDriverByName(driver_name).CreateDataSource(str(out_path))
    out_layer = out_ds.CreateLayer(out_path.stem, geom_type=src_layer.GetGeomType())
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
    src_ds = None
    say(f"Wrote {n} features to {out_path} (no CRS assigned -- see module docstring)")
    return n
