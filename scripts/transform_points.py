#!/usr/bin/env python3
"""
Transform points between BGS2005/CCS2005 (EPSG:7801) and "Coordinate system 1970"
zone K{n}, forward or reverse, for any of Bulgaria's four zones (K3/K5/K7/K9),
using either:

  --method grid      output/bg_k{n}_tinshift.json      (build_tin_grid.py)
                      Sub-millimetre inside the sampled territory of the zone;
                      returns inf/inf outside it.
  --method analytic  output/lcc_affine_fit_k{n}.json   (fit_analytic_lcc.py)
                      Closed-form, works everywhere, 0.09-0.15 m RMSE typical.

By default (--method grid, which is the recommended/most accurate option) points
that fall outside the grid's coverage are automatically retried with the analytic
model instead (disable with --no-fallback).

Examples:
    # single point, CCS2005 -> KS1970 K3
    python3 transform_points.py --zone k3 --xy 300000 4800000

    # single point, KS1970 K5 -> CCS2005
    python3 transform_points.py --zone k5 --direction reverse --xy 9500000 4600000

    # batch CSV (columns named x,y by default)
    python3 transform_points.py --zone k9 --csv in.csv --out out.csv

    # force the closed-form analytic model instead of the grid
    python3 transform_points.py --zone k7 --method analytic --xy 500000 4800000
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from pyproj import Transformer

from bg_points import SRC_LCC_PROJ4, ZONES

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output"


def grid_json_path(zone: str) -> Path:
    return OUTPUT_DIR / f"bg_{zone}_tinshift.json"


def analytic_json_path(zone: str) -> Path:
    return OUTPUT_DIR / f"lcc_affine_fit_{zone}.json"


def build_analytic_pipeline(zone: str) -> str:
    """Full CCS2005 Easting/Northing -> KS1970 Easting/Northing pipeline: undo the known
    EPSG:7801 projection, then apply the fitted conic + similarity from
    fit_analytic_lcc.py. `phi` below is the fitted cone's tangent parallel asin(n);
    lon0 is a convention, not a recovered value -- see fit_analytic_lcc.py."""
    p = json.loads(analytic_json_path(zone).read_text())
    lat0_fixed, lon0 = p["lat0_fixed_deg"], p["lon0_convention_deg"]
    phi = p["tangent_parallel_deg"]
    theta, s, tx, ty = p["theta_rad"], p["scale"], p["tx"], p["ty"]
    c, sn = np.cos(theta), np.sin(theta)
    return (
        f"+proj=pipeline"
        f" +step +inv {SRC_LCC_PROJ4}"
        f" +step +proj=lcc +lat_0={lat0_fixed:.10f} +lat_1={phi:.10f} +lat_2={phi:.10f}"
        f" +lon_0={lon0:.10f} +x_0=0 +y_0=0 +ellps=GRS80"
        f" +step +proj=affine +xoff={tx:.6f} +yoff={ty:.6f}"
        f" +s11={s*c:.12f} +s12={-s*sn:.12f} +s21={s*sn:.12f} +s22={s*c:.12f}"
    )


def build_grid_pipeline(zone: str) -> str:
    return f"+proj=pipeline +step +proj=tinshift +file={grid_json_path(zone)}"


def transform(x, y, zone: str, direction: str, method: str, fallback: bool):
    x = np.atleast_1d(np.asarray(x, dtype=np.float64))
    y = np.atleast_1d(np.asarray(y, dtype=np.float64))
    pyproj_dir = "FORWARD" if direction == "forward" else "INVERSE"

    if method == "grid":
        tr = Transformer.from_pipeline(build_grid_pipeline(zone))
        X, Y = tr.transform(x, y, direction=pyproj_dir)
        X, Y = np.asarray(X, dtype=np.float64), np.asarray(Y, dtype=np.float64)
        bad = ~(np.isfinite(X) & np.isfinite(Y))
        if bad.any():
            if fallback:
                print(
                    f"[transform_points] {bad.sum()} point(s) outside grid coverage, "
                    "falling back to the analytic model for those.",
                    file=sys.stderr,
                )
                tr_a = Transformer.from_pipeline(build_analytic_pipeline(zone))
                Xa, Ya = tr_a.transform(x[bad], y[bad], direction=pyproj_dir)
                X[bad], Y[bad] = Xa, Ya
            else:
                print(f"[transform_points] WARNING: {bad.sum()} point(s) outside grid coverage (inf).",
                      file=sys.stderr)
        return X, Y

    if method == "analytic":
        tr = Transformer.from_pipeline(build_analytic_pipeline(zone))
        X, Y = tr.transform(x, y, direction=pyproj_dir)
        return np.asarray(X, dtype=np.float64), np.asarray(Y, dtype=np.float64)

    raise ValueError(f"unknown method {method!r}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zone", choices=ZONES, required=True, help="which 1970-system zone (k3/k5/k7/k9)")
    ap.add_argument("--method", choices=["grid", "analytic"], default="grid")
    ap.add_argument("--direction", choices=["forward", "reverse"], default="forward",
                     help="forward: CCS2005 -> KS1970 K{n}. reverse: KS1970 K{n} -> CCS2005.")
    ap.add_argument("--no-fallback", action="store_true",
                     help="disable automatic analytic fallback for out-of-coverage points (grid method only)")
    ap.add_argument("--xy", nargs=2, type=float, metavar=("X", "Y"), help="transform a single point")
    ap.add_argument("--csv", type=Path, help="batch-transform points from a CSV file")
    ap.add_argument("--out", type=Path, help="output CSV path (required with --csv)")
    ap.add_argument("--x-field", default="x")
    ap.add_argument("--y-field", default="y")
    args = ap.parse_args()

    if args.xy:
        X, Y = transform([args.xy[0]], [args.xy[1]], args.zone, args.direction, args.method, not args.no_fallback)
        print(f"{X[0]:.4f} {Y[0]:.4f}")
        return

    if args.csv:
        if not args.out:
            ap.error("--csv requires --out")
        with args.csv.open(newline="") as f:
            rows = list(csv.DictReader(f))
        x = [float(r[args.x_field]) for r in rows]
        y = [float(r[args.y_field]) for r in rows]
        X, Y = transform(x, y, args.zone, args.direction, args.method, not args.no_fallback)
        fieldnames = list(rows[0].keys()) + ["x_out", "y_out"] if rows else ["x_out", "y_out"]
        with args.out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r, xo, yo in zip(rows, X, Y):
                r["x_out"], r["y_out"] = f"{xo:.4f}", f"{yo:.4f}"
                w.writerow(r)
        print(f"wrote {len(rows)} rows to {args.out}")
        return

    ap.error("provide --xy X Y or --csv IN --out OUT")


if __name__ == "__main__":
    main()
