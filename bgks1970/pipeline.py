"""
Reconstruct everything from the archived source data, in one run.

The stages are ordered by dependency, and each one only needs what the previous
ones wrote:

  1. audit    the input archives -- stops the run if an export is incomplete,
              because every number downstream would otherwise be quietly built on
              truncated data
  2. fit      the identifiable closed-form model per zone -> lcc_affine_fit_k*.json
  3. grid     the tinshift correction grid per zone      -> bg_k*_tinshift.json
              and the accuracy/seam summary              -> tinshift_accuracy.json
  4. crs      a QGIS-importable Custom CRS per zone      -> bg_k*_ks1970_fitted.wkt
              (reads the fit's JSON)
  5. figures  the paper's two SVG figures, injected into the paper
  6. paper    the paper rendered to PDF, one file per language

Nothing here needs the archives unpacked: GDAL reads the shapefiles straight out
of the zips.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable

from . import check, crs, fit, figures, grid, paper
from .data import ZONES

#: Stage name -> (description, callable). `paper` is last because it is the only
#: stage needing a browser, so a run without one still produces everything else.
STAGE_ORDER = ("audit", "fit", "grid", "crs", "figures", "paper")


@dataclass
class StageResult:
    name: str
    ok: bool
    seconds: float
    detail: str = ""
    result: object = field(default=None, repr=False)


def reconstruct(zones=ZONES, root=None, stages=STAGE_ORDER, verbose: bool = True,
                strict: bool = True, browser: str | None = None) -> list[StageResult]:
    """Run the whole reconstruction. Returns one StageResult per stage.

    `strict` stops at the first failing stage, which is what you want for the audit:
    fitting truncated data produces numbers that look fine and are not. Set it False
    to push through (for instance to regenerate figures when no browser is present
    for the PDF stage)."""
    say = print if verbose else (lambda *a, **k: None)
    results: list[StageResult] = []

    def run(name: str, fn: Callable[[], object], describe: Callable[[object], str]):
        say(f"\n{'='*70}\n== {name}\n{'='*70}")
        t0 = time.perf_counter()
        try:
            out = fn()
        except Exception as exc:  # noqa: BLE001 -- reported, then re-raised if strict
            dt = time.perf_counter() - t0
            results.append(StageResult(name, False, dt, f"{type(exc).__name__}: {exc}"))
            say(f"\n-- {name} FAILED after {dt:.1f}s: {exc}")
            if strict:
                raise
            return None
        dt = time.perf_counter() - t0
        results.append(StageResult(name, True, dt, describe(out), out))
        say(f"\n-- {name} ok in {dt:.1f}s")
        return out

    if "audit" in stages:
        def _audit():
            problems = check.audit(root=root, verbose=verbose)
            if problems:
                raise RuntimeError(
                    f"{len(problems)} problem(s) in the input data; fix them before "
                    f"fitting, or pass strict=False to continue anyway. First: {problems[0]}"
                )
            return problems
        run("audit", _audit, lambda p: "all layers complete and consistent")

    if "fit" in stages:
        run("fit", lambda: fit.fit_all(zones, root, verbose),
            lambda r: f"{len(r)} zones fitted")
    if "grid" in stages:
        run("grid", lambda: grid.build_all(zones, root, verbose),
            lambda r: f"{len(r)} grids built")
    if "crs" in stages:
        run("crs", lambda: crs.build_all(zones, root, verbose), lambda r: "WKT written")
    if "figures" in stages:
        run("figures", lambda: figures.build_figures(root, verbose),
            lambda r: "figures written and injected")
    if "paper" in stages:
        run("paper", lambda: paper.build_pdfs(root=root, browser=browser, verbose=verbose),
            lambda r: f"{len(r)} PDF(s)")

    if verbose:
        total = sum(r.seconds for r in results)
        say(f"\n{'='*70}\n== reconstruction summary\n{'='*70}")
        for r in results:
            mark = "ok  " if r.ok else "FAIL"
            say(f"  {mark}  {r.name:<9s} {r.seconds:6.1f}s  {r.detail}")
        say(f"  {'':4s}  {'total':<9s} {total:6.1f}s")
    return results
