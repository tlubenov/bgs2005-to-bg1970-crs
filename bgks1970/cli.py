"""
Command-line interface: `bgks1970 <command>`.

`reconstruct` is the one-command entry point -- it audits the input archives, fits
the analytic model, builds the grids, writes the QGIS CRS definitions, regenerates
the paper's figures and renders its PDF editions. The individual stages are also
exposed as their own commands for when only one needs redoing.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from . import __version__
from .data import ZONES
from .pipeline import STAGE_ORDER


def _add_common(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("--root", help="project directory holding the data archives "
                                   "(default: cwd, else the installed repo)")
    ap.add_argument("-q", "--quiet", action="store_true", help="suppress progress output")


def _zone_arg(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("--zone", choices=ZONES, help="only this zone (default: all four)")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="bgks1970",
        description="Recover Bulgaria's unpublished KS1970 Lambert zones (K3/K5/K7/K9) "
                    "from BGS2005/CCS2005.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Run `bgks1970 reconstruct` to rebuild every artifact in one go.",
    )
    ap.add_argument("--version", action="version", version=f"bgks1970 {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("reconstruct", help="run the whole pipeline (one command)")
    _add_common(p); _zone_arg(p)
    p.add_argument("--stages", default=",".join(STAGE_ORDER),
                   help=f"comma-separated subset of: {','.join(STAGE_ORDER)}")
    p.add_argument("--keep-going", action="store_true",
                   help="continue past a failing stage instead of stopping")
    p.add_argument("--browser", help="path to Chrome/Chromium for the PDF stage")

    p = sub.add_parser("check", help="audit the input archives")
    _add_common(p)
    p.add_argument("--dataset", default=None, help="dataset name (default: the primary one)")

    p = sub.add_parser("fit", help="fit the closed-form model per zone")
    _add_common(p); _zone_arg(p)
    p.add_argument("--identifiability", action="store_true",
                   help="demonstrate that lon0/lat1/lat2 are not recoverable, then exit")

    p = sub.add_parser("grid", help="build the tinshift grids")
    _add_common(p); _zone_arg(p)

    p = sub.add_parser("crs", help="write the QGIS Custom CRS definitions")
    _add_common(p); _zone_arg(p)

    p = sub.add_parser("figures", help="regenerate the paper's figures")
    _add_common(p)

    p = sub.add_parser("paper", help="render the paper to PDF")
    _add_common(p)
    p.add_argument("--lang", choices=("bg", "en"), help="only this edition")
    p.add_argument("--browser", help="path to a Chrome/Chromium binary")

    p = sub.add_parser("transform", help="transform points between CCS2005 and KS1970")
    _add_common(p)
    p.add_argument("--zone", choices=ZONES, required=True)
    p.add_argument("--method", choices=("grid", "analytic"), default="grid")
    p.add_argument("--direction", choices=("forward", "reverse"), default="forward",
                   help="forward: CCS2005 -> KS1970 K{n}. reverse: the other way.")
    p.add_argument("--no-fallback", action="store_true",
                   help="do not retry out-of-coverage points with the analytic model")
    p.add_argument("--xy", nargs=2, type=float, metavar=("X", "Y"))
    p.add_argument("--csv", type=Path, help="batch-transform a CSV")
    p.add_argument("--out", type=Path, help="output CSV (required with --csv)")
    p.add_argument("--x-field", default="x")
    p.add_argument("--y-field", default="y")

    p = sub.add_parser("reproject", help="reproject a whole vector file")
    _add_common(p)
    p.add_argument("--zone", choices=ZONES, required=True)
    p.add_argument("--direction", choices=("forward", "reverse"), default="forward")
    p.add_argument("--method", choices=("grid", "analytic"), default="grid")
    p.add_argument("--no-fallback", action="store_true")
    p.add_argument("--in", dest="in_path", type=Path, required=True)
    p.add_argument("--out", dest="out_path", type=Path, required=True)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    verbose = not getattr(args, "quiet", False)
    root = getattr(args, "root", None)
    zones = [args.zone] if getattr(args, "zone", None) else list(ZONES)
    cmd = args.command

    if cmd == "reconstruct":
        from .pipeline import reconstruct
        stages = tuple(s.strip() for s in args.stages.split(",") if s.strip())
        unknown = [s for s in stages if s not in STAGE_ORDER]
        if unknown:
            print(f"unknown stage(s): {unknown}; choose from {list(STAGE_ORDER)}",
                  file=sys.stderr)
            return 2
        results = reconstruct(zones, root, stages, verbose,
                              strict=not args.keep_going, browser=args.browser)
        return 0 if all(r.ok for r in results) else 1

    if cmd == "check":
        from .check import audit
        from .data import PRIMARY_DATASET
        problems = audit(args.dataset or PRIMARY_DATASET, root, verbose)
        return 1 if problems else 0

    if cmd == "fit":
        from .fit import fit_all, identifiability_report
        if args.identifiability:
            for z in zones:
                identifiability_report(z, root)
            return 0
        fit_all(zones, root, verbose)
        return 0

    if cmd == "grid":
        from .grid import build_all
        build_all(zones, root, verbose)
        return 0

    if cmd == "crs":
        from .crs import build_all
        build_all(zones, root, verbose)
        return 0

    if cmd == "figures":
        from .figures import build_figures
        build_figures(root, verbose)
        return 0

    if cmd == "paper":
        from .paper import LANGS, build_pdfs
        build_pdfs([args.lang] if args.lang else LANGS, root, args.browser, verbose)
        return 0

    if cmd == "transform":
        from .transform import transform
        fallback = not args.no_fallback
        if args.xy:
            X, Y = transform([args.xy[0]], [args.xy[1]], args.zone, args.direction,
                             args.method, fallback, root)
            print(f"{X[0]:.4f} {Y[0]:.4f}")
            return 0
        if args.csv:
            if not args.out:
                print("--csv requires --out", file=sys.stderr)
                return 2
            with args.csv.open(newline="") as fh:
                rows = list(csv.DictReader(fh))
            X, Y = transform([float(r[args.x_field]) for r in rows],
                             [float(r[args.y_field]) for r in rows],
                             args.zone, args.direction, args.method, fallback, root)
            fields = (list(rows[0].keys()) + ["x_out", "y_out"]) if rows else ["x_out", "y_out"]
            with args.out.open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields)
                w.writeheader()
                for r, xo, yo in zip(rows, X, Y):
                    r["x_out"], r["y_out"] = f"{xo:.4f}", f"{yo:.4f}"
                    w.writerow(r)
            if verbose:
                print(f"wrote {len(rows)} rows to {args.out}")
            return 0
        print("provide --xy X Y or --csv IN --out OUT", file=sys.stderr)
        return 2

    if cmd == "reproject":
        from .reproject import reproject_file
        reproject_file(args.in_path, args.out_path, args.zone, args.direction,
                       args.method, not args.no_fallback, root, verbose)
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
