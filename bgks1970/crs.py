"""
Generate a QGIS-importable WKT2 "Custom CRS" for each zone's fitted analytic model
(fit_analytic_lcc.py), so the whole conic + rotation/scale/translation model
can be used as a single ordinary CRS -- selectable in QGIS's Custom CRS manager, usable
in "Set Layer/Project CRS", "Reproject Layer", on-the-fly reprojection, etc. -- without
having to paste a raw pipeline string anywhere.

This works via PROJ's "PROJ-based operation method": a WKT2 DERIVEDPROJCRS whose
DERIVINGCONVERSION carries an arbitrary +proj=pipeline string as its METHOD. Verified here against
the installed PROJ/pyproj: CRS.from_wkt() parses it, and
Transformer.from_crs("EPSG:7801", crs) both transforms and round-trips correctly.

IMPORTANT: this WKT encodes the ANALYTIC model only (0.09-0.15 m RMSE, see README),
not the sub-millimetre tinshift grid -- a "CRS" is fundamentally a formula, it cannot
carry a triangulated correction grid. For survey-grade accuracy use the grid via
build_tin_grid.py + transform_points.py / reproject_shapefile.py instead; this WKT is
for convenient approximate on-the-fly use directly inside QGIS/any WKT-aware GIS tool.

The CRS is named "unofficial" on purpose: it is an empirical reconstruction, not a
published definition, and its conic parameters are a convention rather than recovered
values (see fit_analytic_lcc.py on identifiability).
"""
from __future__ import annotations

import numpy as np
from pyproj import CRS, Transformer

from .data import ZONES, load_sources, resolve_split
from .paths import output_dir
from .transform import build_analytic_pipeline

BASE_PROJCRS_WKT = '''BASEPROJCRS["BGS2005 / CCS2005",
   BASEGEOGCRS["BGS2005",
     DATUM["Bulgaria Geodetic System 2005",
       ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],
     PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]]],
   CONVERSION["CCS2005",
     METHOD["Lambert Conic Conformal (2SP)",ID["EPSG",9802]],
     PARAMETER["Latitude of false origin",42.6678756833333,ANGLEUNIT["degree",0.0174532925199433]],
     PARAMETER["Longitude of false origin",25.5,ANGLEUNIT["degree",0.0174532925199433]],
     PARAMETER["Latitude of 1st standard parallel",42,ANGLEUNIT["degree",0.0174532925199433]],
     PARAMETER["Latitude of 2nd standard parallel",43.3333333333333,ANGLEUNIT["degree",0.0174532925199433]],
     PARAMETER["Easting at false origin",500000,LENGTHUNIT["metre",1]],
     PARAMETER["Northing at false origin",4725824.3591,LENGTHUNIT["metre",1]]],
   CS[Cartesian,2],AXIS["easting",east,ORDER[1],LENGTHUNIT["metre",1]],AXIS["northing",north,ORDER[2],LENGTHUNIT["metre",1]],
   ID["EPSG",7801]]'''


def build_wkt(zone: str, root=None) -> str:
    full_pipeline = build_analytic_pipeline(zone, root)
    n = zone[1:]
    return (
        f'DERIVEDPROJCRS["BG 1970 K{n} (fitted, unofficial)",\n'
        f' {BASE_PROJCRS_WKT},\n'
        f' DERIVINGCONVERSION["CCS2005 -> KS1970 K{n}, fitted pipeline",\n'
        f'   METHOD["PROJ-based operation method: {full_pipeline}"]],\n'
        f' CS[Cartesian,2],\n'
        f' AXIS["easting",east,ORDER[1],LENGTHUNIT["metre",1]],\n'
        f' AXIS["northing",north,ORDER[2],LENGTHUNIT["metre",1]]]'
    )


def build_zone(zone: str, root=None, verbose: bool = True):
    say = print if verbose else (lambda *a, **k: None)
    wkt = build_wkt(zone, root)

    crs = CRS.from_wkt(wkt)  # raises if malformed
    tr = Transformer.from_crs("EPSG:7801", crs, always_xy=True)

    # Score against the out-of-sample set, so this number is comparable to the one
    # fit_analytic_lcc.py reports rather than being a training-set figure.
    split = resolve_split(zone, root)
    src_xy, dst_xy = load_sources(zone, split.test)
    X, Y = tr.transform(src_xy[:, 0], src_xy[:, 1])
    d = np.hypot(np.asarray(X) - dst_xy[:, 0], np.asarray(Y) - dst_xy[:, 1])
    say(f"{zone}: WKT parses OK; out-of-sample RMSE over {len(src_xy)} points "
          f"({', '.join(s.label for s in split.test)}) = {np.sqrt(np.mean(d**2)):.3f} m "
          f"(matches the analytic model's own accuracy)")

    # Round-trip through the CRS must return the input.
    xb, yb = tr.transform(np.asarray(X), np.asarray(Y), direction="INVERSE")
    rt = np.hypot(np.asarray(xb) - src_xy[:, 0], np.asarray(yb) - src_xy[:, 1]).max()
    say(f"  round-trip through the CRS: max {rt*1000:.6f} mm")
    assert rt < 1e-3

    out_path = output_dir(root) / f"bg_{zone}_ks1970_fitted.wkt"
    out_path.write_text(wkt + "\n")
    say(f"  wrote {out_path}")


def build_all(zones=ZONES, root=None, verbose: bool = True) -> None:
    for zone in zones:
        build_zone(zone, root, verbose)
