"""
Shared helpers for loading matched control points for one of Bulgaria's four
"1970 system" Lambert zones (K3/K5/K7/K9).

Data layout in source_data/ -- per zone n in {3,5,7,9} and per sampling step
r in {1km, 2km}:

  bg_k{n}_{r}_sel.shp      Bulgaria Geodetic System 2005 / CCS2005 (EPSG:7801)
                            -> fully known Lambert Conformal Conic (2SP) definition.
  bg_k{n}_{r}_sel_tr.shp   "Coordinate system 1970", zone K{n} (Lambert Conformal
                            Conic, Pulkovo 1942(58) / Krasovsky 1940 ellipsoid),
                            produced by AGKK's official BGSTrans tool
                            (see bg_k{n}_{r}_sel_tr.log: "Plane transform
                            accuracy - 0.14 meters"). Its .prj names the datum and
                            the projection *type* but carries no numeric parameter.

Each `_sel` / `_sel_tr` pair shares an integer "id" field, so the same feature in
both files is the same physical point before/after AGKK's transformation. That
gives tens of thousands of exact correspondences per zone to recover the unknown
projection from -- no manually surveyed control points needed.

Why two sampling steps
----------------------
The 1 km and 2 km sets are *disjoint*: they share no ids, and their lattices are
offset from each other by 120 m, so no 2 km point coincides with a 1 km point
(verified -- zero coordinate collisions in all four zones). Every model here is
therefore fitted/built on the 1 km set alone and scored on the 2 km set, which
is a genuinely independent sample rather than a random split of one lattice.
For an interpolating model (build_tin_grid.py) that distinction matters a lot:
a random hold-out of a single lattice leaves test points sitting on the edges of
the triangles built from their own neighbours, which flatters the result. The
2 km points instead land in the *interior* of the 1 km cells, which is where
piecewise-linear interpolation is at its worst.

bg_ext/bg_zones1970.shp holds the four official zone extents (one polygon per
zone, field CLIST = 3/5/7/9) in WGS 84 / UTM zone 35N; load_zone_polygon() reads
them, reprojected to whichever CRS you ask for.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from osgeo import ogr

ogr.UseExceptions()

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "source_data"
EXT_DIR = ROOT / "bg_ext"

ZONES = ["k3", "k5", "k7", "k9"]

#: Sampling step of each shapefile set. "1km" is the dense set every model is
#: built from; "2km" is the independent set every model is scored against.
STEPS = ["1km", "2km"]
FIT_STEP = "1km"
TEST_STEP = "2km"

#: Compass quadrant each zone covers (visible directly in the destination
#: extents: K3/K9 eastings carry an "8xxxxxxx" prefix, K5/K7 a "9xxxxxxx" one).
QUADRANT = {"k3": "northwest", "k5": "southeast", "k7": "northeast", "k9": "southwest"}

SRC_PROJECTED_EPSG = "EPSG:7801"  # BGS2005 / CCS2005, fully defined
SRC_GEOGRAPHIC_EPSG = "EPSG:7798"  # BGS2005 geographic (GRS80), used for the analytic fit
ZONES_EXT_EPSG = "EPSG:32635"  # WGS 84 / UTM zone 35N -- the CRS bg_zones1970.shp is stored in

# Full known definition of every bg_k{n}_{r}_sel layer's CRS, EPSG:7801 (BGS2005 /
# CCS2005). Identical across all four zones -- only the destination projection is
# zone-specific and unknown.
SRC_LCC_PROJ4 = (
    "+proj=lcc +lat_0=42.6678756833333 +lat_1=42 +lat_2=43.3333333333333"
    " +lon_0=25.5 +x_0=500000 +y_0=4725824.3591 +ellps=GRS80"
)


def src_shp(zone: str, step: str = FIT_STEP) -> Path:
    return DATA_DIR / f"bg_{zone}_{step}_sel.shp"


def dst_shp(zone: str, step: str = FIT_STEP) -> Path:
    return DATA_DIR / f"bg_{zone}_{step}_sel_tr.shp"


def read_id_xy(shp_path: Path) -> dict[int, tuple[float, float]]:
    """Read {id: (x, y)} from a point shapefile. x/y are read in the file's native
    (Easting, Northing) storage order, which is what OGR's GetX()/GetY() return
    regardless of the CRS's authority-defined axis order."""
    ds = ogr.Open(str(shp_path))
    if ds is None:
        raise FileNotFoundError(shp_path)
    layer = ds.GetLayer()
    out: dict[int, tuple[float, float]] = {}
    for feat in layer:
        fid = feat.GetField("id")
        geom = feat.GetGeometryRef()
        out[fid] = (geom.GetX(), geom.GetY())
    return out


def load_matched_points(zone: str, step: str = FIT_STEP) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns (ids, src_xy, dst_xy) as aligned arrays for one zone ('k3'/'k5'/'k7'/'k9')
    at one sampling step ('1km'/'2km'). src_xy is in EPSG:7801 (the _sel layer), dst_xy in
    the unknown KS1970 zone system (the _sel_tr layer). Only ids present in both are kept."""
    if step not in STEPS:
        raise ValueError(f"step must be one of {STEPS}, got {step!r}")
    src = read_id_xy(src_shp(zone, step))
    dst = read_id_xy(dst_shp(zone, step))

    common = sorted(set(src) & set(dst))
    missing_src = set(dst) - set(src)
    missing_dst = set(src) - set(dst)
    if missing_src or missing_dst:
        print(
            f"[bg_points:{zone}/{step}] WARNING: {len(missing_src)} ids only in _sel_tr, "
            f"{len(missing_dst)} ids only in _sel -- dropped from the matched set."
        )

    ids = np.array(common, dtype=np.int64)
    src_xy = np.array([src[i] for i in common], dtype=np.float64)
    dst_xy = np.array([dst[i] for i in common], dtype=np.float64)
    return ids, src_xy, dst_xy


def load_zone_polygon(zone: str, to_epsg: str = SRC_PROJECTED_EPSG) -> list[np.ndarray]:
    """Read zone K{n}'s official extent from bg_ext/bg_zones1970.shp as a list of
    rings (Nx2 arrays of x/y), reprojected from UTM 35N into `to_epsg`."""
    from pyproj import Transformer

    ds = ogr.Open(str(EXT_DIR / "bg_zones1970.shp"))
    if ds is None:
        raise FileNotFoundError(EXT_DIR / "bg_zones1970.shp")
    layer = ds.GetLayer()
    layer.SetAttributeFilter(f"CLIST = {int(zone[1:])}")
    feat = next(iter(layer), None)
    if feat is None:
        raise ValueError(f"no polygon with CLIST={zone[1:]} in bg_zones1970.shp")

    tr = Transformer.from_crs(ZONES_EXT_EPSG, to_epsg, always_xy=True)
    geom = feat.GetGeometryRef()
    rings = []
    for i in range(geom.GetGeometryCount()):
        ring = geom.GetGeometryRef(i)
        pts = np.array(ring.GetPoints(), dtype=np.float64)[:, :2]
        x, y = tr.transform(pts[:, 0], pts[:, 1])
        rings.append(np.column_stack([x, y]))
    return rings


if __name__ == "__main__":
    for zone in ZONES:
        print(f"--- {zone} ({QUADRANT[zone]}) ---")
        for step in STEPS:
            ids, src_xy, dst_xy = load_matched_points(zone, step)
            print(f"  {step}: {len(ids)} matched points")
            print(f"    src E {src_xy[:,0].min():12.3f} .. {src_xy[:,0].max():12.3f}   "
                  f"N {src_xy[:,1].min():12.3f} .. {src_xy[:,1].max():12.3f}")
            print(f"    dst E {dst_xy[:,0].min():12.3f} .. {dst_xy[:,0].max():12.3f}   "
                  f"N {dst_xy[:,1].min():12.3f} .. {dst_xy[:,1].max():12.3f}")
        rings = load_zone_polygon(zone)
        print(f"  bg_ext extent: {len(rings)} ring(s), {sum(len(r) for r in rings)} vertices")
