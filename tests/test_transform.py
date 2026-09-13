"""The published transforms, exercised through PROJ exactly as a user would."""
import json

import numpy as np

from bgks1970.data import ZONES, load_sources, resolve_split
from bgks1970.paths import output_dir
from bgks1970.transform import transform

# A point comfortably inside each zone's sampled territory.
INSIDE = {"k3": (300000, 4800000), "k5": (600000, 4650000),
          "k7": (600000, 4800000), "k9": (300000, 4650000)}


def _outputs_present() -> bool:
    return all((output_dir() / f"bg_{z}_tinshift.json").exists() for z in ZONES)


def test_grid_round_trips_to_nanometres():
    if not _outputs_present():
        return  # nothing built yet; `bgks1970 reconstruct` first
    for z in ZONES:
        x, y = INSIDE[z]
        E, N = transform([x], [y], z)
        xb, yb = transform(E, N, z, direction="reverse")
        assert np.hypot(xb[0] - x, yb[0] - y) < 1e-6, z


def test_analytic_round_trips_and_agrees_with_the_grid_to_its_stated_accuracy():
    if not _outputs_present():
        return
    for z in ZONES:
        x, y = INSIDE[z]
        Eg, Ng = transform([x], [y], z, method="grid")
        Ea, Na = transform([x], [y], z, method="analytic")
        xb, yb = transform(Ea, Na, z, direction="reverse", method="analytic")
        assert np.hypot(xb[0] - x, yb[0] - y) < 1e-6, z
        stated = json.loads((output_dir() / f"lcc_affine_fit_{z}.json").read_text())
        tol = 6 * stated["accuracy_m"]["out_of_sample_rmse"]
        assert np.hypot(Eg[0] - Ea[0], Ng[0] - Na[0]) < tol, z


def test_analytic_reproduces_its_recorded_accuracy():
    """The JSON's stated RMSE must be what the published pipeline actually delivers."""
    if not _outputs_present():
        return
    for z in ZONES:
        split = resolve_split(z)
        src, dst = load_sources(z, split.test)
        E, N = transform(src[:, 0], src[:, 1], z, method="analytic")
        rmse = float(np.sqrt(np.mean((E - dst[:, 0]) ** 2 + (N - dst[:, 1]) ** 2)))
        stated = json.loads((output_dir() / f"lcc_affine_fit_{z}.json").read_text())
        assert abs(rmse - stated["accuracy_m"]["out_of_sample_rmse"]) < 1e-6, z
