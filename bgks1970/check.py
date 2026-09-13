"""
Audit the input shapefiles before anything is fitted from them.

BGSTrans writes the transformed shapefile LAST. If a run is interrupted it leaves a
`.dbf` and `.shx` behind with no `.shp`, and a truncated `.log` -- which looks, to a
casual listing, like a complete export. This makes that state impossible to miss,
and reports exactly how far the run got, because the orphaned `.dbf` still holds the
ids of every record that was processed.

Checked per layer:

  * all four shapefile parts present (the `.shp` is the one that goes missing)
  * the AGKK `.log` carries its accuracy report rather than just a timestamp
  * `id` is unique -- the _sel <-> _sel_tr join is by id, and a repeat would
    silently drop points
  * no null or duplicated geometry
  * the source points really are on a regular lattice, and at what step (the
    transformed side is warped by construction, so it is not checked for this)
  * `_sel` and `_sel_tr` hold the same number of records
"""
from __future__ import annotations

import struct

import numpy as np

from .data import PRIMARY_DATASET, STEPS, ZONES, dst_stem, read_id_xy, src_stem
from .datasets import Dataset


def dbf_record_count(raw: bytes) -> int:
    """Record count from a .dbf header -- readable even when its .shp is missing."""
    return struct.unpack("<I", raw[4:8])[0]


def dbf_ids(raw: bytes) -> list[int]:
    """Every id in a .dbf, in file order. Works on the orphan left by a failed run."""
    n = struct.unpack("<I", raw[4:8])[0]
    header = struct.unpack("<H", raw[8:10])[0]
    reclen = struct.unpack("<H", raw[10:12])[0]
    out = []
    for i in range(n):
        off = header + i * reclen
        out.append(int(raw[off + 1:off + reclen].decode("latin1").strip()))
    return out


def _layer_stats(ogr_path: str, problems: list[str], lattice: bool) -> dict:
    from osgeo import ogr
    ds = ogr.Open(ogr_path)
    layer = ds.GetLayer()
    ids, xy = [], []
    for feat in layer:
        geom = feat.GetGeometryRef()
        ids.append(feat.GetField("id"))
        xy.append((np.nan, np.nan) if geom is None or geom.IsEmpty()
                  else (geom.GetX(), geom.GetY()))
    del ds
    xy = np.array(xy, dtype=np.float64)
    name = ogr_path.rsplit("/", 1)[-1]
    info: dict = {"n": len(ids)}

    if len(set(ids)) != len(ids):
        problems.append(f"{name}: 'id' is not unique -- the join would drop points")
    n_null = int(np.isnan(xy).any(axis=1).sum())
    if n_null:
        problems.append(f"{name}: {n_null} null/empty geometries")

    good = xy[~np.isnan(xy).any(axis=1)]
    n_dupe = len(good) - len({tuple(p) for p in np.round(good, 3)})
    if n_dupe:
        problems.append(f"{name}: {n_dupe} duplicated positions")

    ux, uy = np.unique(good[:, 0]), np.unique(good[:, 1])
    if lattice and ux.size > 1 and uy.size > 1:
        sx, sy = float(np.diff(ux).min()), float(np.diff(uy).min())
        info["step"] = sx
        on = (np.allclose((good[:, 0] - ux.min()) % sx, 0, atol=1e-6)
              and np.allclose((good[:, 1] - uy.min()) % sy, 0, atol=1e-6))
        if not on:
            problems.append(f"{name}: points are not on a regular {sx:g} m lattice")
    return info


def _log_state(dataset: Dataset, stem: str, problems: list[str]) -> str:
    """AGKK's log carries a full report on success and little else on failure."""
    name = f"{stem}.log"
    if not dataset.has(name):
        return "no .log"
    text = dataset.read_text(name)
    if "accuracy" in text.lower():
        return "report present"
    problems.append(f"{name}: truncated ({len(text)} bytes, no accuracy report) -- "
                    f"that run did not finish")
    return f"TRUNCATED ({len(text)} B)"


def audit(dataset_name: str = PRIMARY_DATASET, root=None, verbose: bool = True) -> list[str]:
    """Audit one dataset. Returns the list of problems (empty when healthy)."""
    ds = Dataset.find(dataset_name, root).require()
    problems: list[str] = []
    if verbose:
        print(f"Auditing {ds.label} in {ds.root}\n")

    for zone in ZONES:
        for step in STEPS:
            src, dst = src_stem(zone, step), dst_stem(zone, step)
            label = f"{zone}/{step}"
            log_state = _log_state(ds, dst, problems)

            if not ds.has_layer(src):
                problems.append(f"{src}.shp: source layer incomplete "
                                f"{ds.missing_parts(src)}")
                if verbose:
                    print(f"  {label:<10s} SOURCE INCOMPLETE {ds.missing_parts(src)}")
                continue
            s_info = _layer_stats(ds.layer_path(src), problems, lattice=True)

            missing = ds.missing_parts(dst)
            if missing:
                # The .shp is gone, but the orphaned .dbf still says how far it got.
                got = dbf_record_count(ds.read_bytes(f"{dst}.dbf")) if ds.has(f"{dst}.dbf") else 0
                short = s_info["n"] - got
                if verbose:
                    print(f"  {label:<10s} src {s_info['n']:6d} @ {s_info.get('step',0):.0f} m"
                          f"  |  TRANSFORMED MISSING {missing}")
                    print(f"  {'':<10s}   orphaned .dbf holds {got} of {s_info['n']} "
                          f"({short} short)  |  log: {log_state}")
                if ds.has(f"{dst}.dbf") and short > 0:
                    done = set(dbf_ids(ds.read_bytes(f"{dst}.dbf")))
                    src_ids = list(read_id_xy(ds.layer_path(src)))
                    dropped = [i for i in src_ids if i not in done]
                    tail = src_ids[-short:] == dropped
                    if verbose:
                        print(f"  {'':<10s}   ids not processed: {dropped[:8]}"
                              f"{' ...' if len(dropped) > 8 else ''}")
                        print(f"  {'':<10s}   " + ("they are the LAST records in the file -- "
                              "a clean truncation, not rejected points" if tail else
                              "they are scattered through the file -- specific points "
                              "were rejected"))
                problems.append(f"{dst}.shp: not written ({short} records short of {s_info['n']})")
                continue

            d_info = _layer_stats(ds.layer_path(dst), problems, lattice=False)
            same = d_info["n"] == s_info["n"]
            if not same:
                problems.append(f"{dst}.shp: {d_info['n']} records vs {s_info['n']} in _sel")
            if verbose:
                print(f"  {label:<10s} src {s_info['n']:6d} @ {s_info.get('step',0):.0f} m"
                      f"  |  tr {d_info['n']:6d}  {'OK' if same else 'COUNT MISMATCH'}"
                      f"  |  log: {log_state}")

    if verbose:
        print()
        if problems:
            print(f"{len(problems)} problem(s):")
            for p in problems:
                print(f"  - {p}")
        else:
            print("All layers complete and consistent.")
    return problems
