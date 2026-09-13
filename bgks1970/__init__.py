"""
bgks1970 -- recovering Bulgaria's unpublished "1970 system" Lambert zones.

Bulgaria's cadastre agency (AGKK) uses a historical "Coordinate system 1970"
(KS1970) split into four Lambert Conformal Conic zones -- K3, K5, K7, K9 -- on the
Pulkovo 1942(58) datum. Its defining parameters were never published. This package
recovers the transformation from BGS2005/CCS2005 (EPSG:7801) empirically, by two
independent methods, from tens of thousands of exact point correspondences produced
by AGKK's own tool:

  * a closed-form model  -- an LCC parameterized by its cone constant, followed by a
    planar similarity; portable as a single PROJ pipeline string, ~0.14-0.21 m
  * a tinshift grid      -- piecewise-linear over the sampled lattice; sub-millimetre
    inside the sampled territory, and the recommended default

Quick start:

    import bgks1970 as bg
    bg.reconstruct()                      # rebuild every artifact, one call
    E, N = bg.transform([300000], [4800000], "k3")   # CCS2005 -> KS1970 K3

or from a shell:

    bgks1970 reconstruct

The sampled data is read straight out of the committed zip archives; nothing needs
unpacking.
"""
from __future__ import annotations

__version__ = "1.0.0"

from .data import QUADRANT, STEPS, ZONES, load_matched_points, resolve_split
from .datasets import Dataset
from .paths import output_dir, paper_dir, project_root
from .pipeline import reconstruct
from .transform import build_analytic_pipeline, build_grid_pipeline, transform

__all__ = [
    "__version__",
    "ZONES", "STEPS", "QUADRANT",
    "Dataset",
    "project_root", "output_dir", "paper_dir",
    "resolve_split", "load_matched_points",
    "transform", "build_grid_pipeline", "build_analytic_pipeline",
    "reconstruct",
]
