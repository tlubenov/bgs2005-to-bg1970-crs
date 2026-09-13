"""
Matched control points for Bulgaria's four "1970 system" Lambert zones.

Per zone n in {3,5,7,9} and per sampling step r in {1km, 2km}, each dataset holds:

  bg_k{n}_{r}_sel.shp      Bulgaria Geodetic System 2005 / CCS2005 (EPSG:7801)
                            -> fully known Lambert Conformal Conic (2SP) definition.
  bg_k{n}_{r}_sel_tr.shp   "Coordinate system 1970", zone K{n} (Lambert Conformal
                            Conic, Pulkovo 1942(58) / Krasovsky 1940 ellipsoid),
                            produced by AGKK's official BGSTrans tool
                            (see bg_k{n}_{r}_sel_tr.log: "Plane transform
                            accuracy - 0.14 meters"). Its .prj names the datum and
                            the projection *type* but carries no numeric parameter.

Each `_sel` / `_sel_tr` pair shares an integer "id", so the same feature in both
files is the same physical point before and after AGKK's transformation. The id is
a per-file row key only, and the join is always made within one pair, so the fact
that a zone's 1 km and 2 km files reuse some id values is harmless.

Why two sampling steps
----------------------
A zone's 1 km and 2 km lattices are offset from each other by exactly half a
kilometre, so a 2 km point sits at the *centre* of a 1 km cell and none coincides
with a 1 km node. Every model is therefore fitted on the 1 km set and scored on the
2 km set, which is a genuinely independent sample rather than a random split of one
lattice. For an interpolating model that distinction matters: a random hold-out of
one lattice leaves test points on the edges of triangles built from their own
neighbours, which flatters the result, whereas cell centres are the worst case for
piecewise-linear interpolation. The reported errors are an upper bound.

Incomplete exports
------------------
BGSTrans writes the transformed shapefile last, so an interrupted export leaves a
`.dbf` and `.shx` with no `.shp` and a truncated `.log`. resolve_split() detects
that and falls back to the best complete combination for the zone, marking the
result provisional; once the export is re-run the same call picks the full data
back up with no code change.
"""
from __future__ import annotations

from typing import Iterable, NamedTuple

import numpy as np
from osgeo import ogr

from .datasets import Dataset
from .paths import ext_dir

ogr.UseExceptions()

ZONES = ("k3", "k5", "k7", "k9")
STEPS = ("1km", "2km")

PRIMARY_DATASET = "source_data_buff20km"
FALLBACK_DATASET = "source_data"

#: Compass quadrant each zone covers (visible in the destination extents:
#: K3/K9 eastings carry an "8xxxxxxx" prefix, K5/K7 a "9xxxxxxx" one).
QUADRANT = {"k3": "northwest", "k5": "southeast", "k7": "northeast", "k9": "southwest"}

SRC_PROJECTED_EPSG = "EPSG:7801"    # BGS2005 / CCS2005, fully defined
SRC_GEOGRAPHIC_EPSG = "EPSG:7798"   # BGS2005 geographic (GRS80), used for the fit
ZONES_EXT_EPSG = "EPSG:32635"       # WGS 84 / UTM 35N -- bg_zones1970.shp's CRS

# Full known definition of every _sel layer's CRS, EPSG:7801. Identical across all
# four zones -- only the destination projection is zone-specific and unknown.
SRC_LCC_PROJ4 = (
    "+proj=lcc +lat_0=42.6678756833333 +lat_1=42 +lat_2=43.3333333333333"
    " +lon_0=25.5 +x_0=500000 +y_0=4725824.3591 +ellps=GRS80"
)


def src_stem(zone: str, step: str) -> str:
    return f"bg_{zone}_{step}_sel"


def dst_stem(zone: str, step: str) -> str:
    return f"bg_{zone}_{step}_sel_tr"


class Source(NamedTuple):
    """One usable (dataset, sampling step) combination for a zone."""
    dataset: Dataset
    step: str

    @property
    def label(self) -> str:
        return f"{self.dataset.label}/{self.step}"


class Split(NamedTuple):
    """How a zone's data divides into what a model is built from and what it is
    scored on. `note` is set only when the preferred arrangement was unavailable."""
    fit: Source
    test: tuple[Source, ...]
    note: str | None

    @property
    def provisional(self) -> bool:
        return self.note is not None


def pair_available(dataset: Dataset, zone: str, step: str) -> bool:
    """True when both halves of a pair are complete -- every shapefile part, not
    just the .shp."""
    return dataset.has_layer(src_stem(zone, step)) and dataset.has_layer(dst_stem(zone, step))


def resolve_split(zone: str, root=None, dataset_name: str = PRIMARY_DATASET) -> Split:
    """Pick what to fit on and what to score against, for one zone.

    Normally the primary dataset's 1 km set against its 2 km set. When a pair is
    incomplete, fall back to the densest complete pair available and score it
    against the other dataset, whose positions are disjoint from it -- and say so,
    so a reduced result is never mistaken for a full one."""
    primary = Dataset.find(dataset_name, root).require()
    if pair_available(primary, zone, "1km") and pair_available(primary, zone, "2km"):
        return Split(Source(primary, "1km"), (Source(primary, "2km"),), None)

    missing = [s for s in STEPS if not pair_available(primary, zone, s)]
    fallback_step = next((s for s in STEPS if pair_available(primary, zone, s)), None)

    others: tuple[Source, ...] = ()
    secondary = Dataset.find(FALLBACK_DATASET, root)
    if secondary.exists and secondary.name != primary.name:
        others = tuple(Source(secondary, s) for s in STEPS
                       if pair_available(secondary, zone, s))

    if fallback_step is None or not others:
        raise FileNotFoundError(
            f"zone {zone}: no usable data. Incomplete pairs in {primary.label}: {missing}"
        )
    note = (
        f"{primary.label} has no complete {'/'.join(missing)} pair for {zone} "
        f"(the transformed .shp was never written -- check its .log). "
        f"Fitted on {primary.label}/{fallback_step} instead and scored against "
        f"{', '.join(s.label for s in others)}, whose positions are disjoint from it. "
        f"Re-run the export and this picks the full data back up automatically."
    )
    return Split(Source(primary, fallback_step), others, note)


def read_id_xy(ogr_path: str) -> dict[int, tuple[float, float]]:
    """Read {id: (x, y)} from a point layer. x/y come back in the file's native
    (Easting, Northing) storage order, which is what OGR returns regardless of the
    CRS's authority-defined axis order."""
    ds = ogr.Open(ogr_path)
    if ds is None:
        raise FileNotFoundError(ogr_path)
    layer = ds.GetLayer()
    out: dict[int, tuple[float, float]] = {}
    for feat in layer:
        geom = feat.GetGeometryRef()
        out[feat.GetField("id")] = (geom.GetX(), geom.GetY())
    count = layer.GetFeatureCount()
    del ds
    if len(out) != count:
        raise ValueError(f"{ogr_path}: 'id' is not unique, the join would drop points")
    return out


def load_matched_points(zone: str, step: str = "1km", dataset: Dataset | None = None,
                        root=None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(ids, src_xy, dst_xy) for one zone at one sampling step. src_xy is EPSG:7801,
    dst_xy the unknown KS1970 zone system. Only ids present in both are kept."""
    if step not in STEPS:
        raise ValueError(f"step must be one of {STEPS}, got {step!r}")
    ds = (dataset or Dataset.find(PRIMARY_DATASET, root)).require()

    src = read_id_xy(ds.layer_path(src_stem(zone, step)))
    dst = read_id_xy(ds.layer_path(dst_stem(zone, step)))

    common = sorted(set(src) & set(dst))
    only_src, only_dst = set(src) - set(dst), set(dst) - set(src)
    if only_src or only_dst:
        print(f"[data:{zone}/{step}] WARNING: {len(only_dst)} ids only in _sel_tr, "
              f"{len(only_src)} ids only in _sel -- dropped from the matched set.")

    ids = np.array(common, dtype=np.int64)
    return (ids,
            np.array([src[i] for i in common], dtype=np.float64),
            np.array([dst[i] for i in common], dtype=np.float64))


def load_sources(zone: str, sources: Iterable[Source]) -> tuple[np.ndarray, np.ndarray]:
    """Concatenate the matched points of several sources into one scoring set."""
    src_parts, dst_parts = [], []
    for s in sources:
        _, sx, dx = load_matched_points(zone, s.step, s.dataset)
        src_parts.append(sx)
        dst_parts.append(dx)
    return np.vstack(src_parts), np.vstack(dst_parts)


def _rings_from_geometry(geom, transformer) -> list[np.ndarray]:
    """Flatten a Polygon or MultiPolygon into a list of Nx2 rings, reprojected."""
    rings: list[np.ndarray] = []

    def walk(g):
        if g.GetGeometryName() in ("MULTIPOLYGON", "GEOMETRYCOLLECTION"):
            for i in range(g.GetGeometryCount()):
                walk(g.GetGeometryRef(i))
            return
        for i in range(g.GetGeometryCount()):
            pts = np.array(g.GetGeometryRef(i).GetPoints(), dtype=np.float64)[:, :2]
            if transformer is not None:
                x, y = transformer.transform(pts[:, 0], pts[:, 1])
                pts = np.column_stack([x, y])
            rings.append(pts)

    walk(geom)
    return rings


def load_zone_polygon(zone: str, to_epsg: str = SRC_PROJECTED_EPSG,
                      buffered: bool = False, root=None) -> list[np.ndarray]:
    """Zone K{n}'s extent as a list of rings in `to_epsg`.

    `buffered=False` is the official zone boundary; `buffered=True` is the 20 km
    buffered outline the primary dataset actually samples."""
    from pyproj import Transformer

    base = ext_dir(root)
    if buffered:
        path, native = base / "bg_zones1970_buff20km.geojson", SRC_PROJECTED_EPSG
    else:
        path, native = base / "bg_zones1970.shp", ZONES_EXT_EPSG
    if not path.exists():
        raise FileNotFoundError(path)

    ds = ogr.Open(str(path))
    layer = ds.GetLayer()
    layer.SetAttributeFilter(f"CLIST = {int(zone[1:])}")
    feat = next(iter(layer), None)
    if feat is None:
        raise ValueError(f"no polygon with CLIST={zone[1:]} in {path.name}")
    tr = None if native == to_epsg else Transformer.from_crs(native, to_epsg, always_xy=True)
    rings = _rings_from_geometry(feat.GetGeometryRef(), tr)
    del ds
    return rings


def load_sampled_extent(to_epsg: str = SRC_PROJECTED_EPSG, root=None) -> list[np.ndarray]:
    """Outline of the whole region the primary dataset samples: the union of the
    four 20 km-buffered zone extents.

    A union rather than four outlines because buffering each zone outward by 20 km
    makes neighbouring buffers overlap heavily -- four outlines would be clutter,
    while their union shows the one fact worth showing: how far past the zone mosaic
    the data reaches."""
    from pyproj import Transformer

    path = ext_dir(root) / "bg_zones1970_buff20km.geojson"
    ds = ogr.Open(str(path))
    layer = ds.GetLayer()
    union = None
    for feat in layer:
        g = feat.GetGeometryRef().Clone()
        union = g if union is None else union.Union(g)
    tr = None if to_epsg == SRC_PROJECTED_EPSG else Transformer.from_crs(
        SRC_PROJECTED_EPSG, to_epsg, always_xy=True)
    rings = _rings_from_geometry(union, tr)
    del ds
    return rings
