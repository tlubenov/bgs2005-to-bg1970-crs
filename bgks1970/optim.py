"""
Levenberg-Marquardt least squares, in plain numpy.

This package used to call scipy.optimize.least_squares(method="lm"). SciPy is a
heavy dependency to carry for one 7-parameter fit, and it is not always
installable next to a system GDAL/PROJ, so the ~60 lines below replace it.

The algorithm is the standard Marquardt formulation: at each iteration solve

    (JᵀJ + λ·diag(JᵀJ)) δ = -Jᵀr

for the step δ, accept it if it lowers the cost (and relax λ), otherwise reject
it and tighten λ. Scaling the damping by diag(JᵀJ) rather than by the identity
is what makes this work on parameter vectors as badly scaled as ours, where
lon0 ≈ 24 (degrees) sits next to tx ≈ 8.5·10⁶ (metres).

The Jacobian is estimated by forward differences with a per-parameter relative
step; for 7 parameters that is 7 extra residual evaluations per iteration, which
is nothing next to the cost of getting an analytic Jacobian right.
"""
from __future__ import annotations

from typing import Callable, NamedTuple

import numpy as np


class LMResult(NamedTuple):
    x: np.ndarray
    cost: float          # 0.5 * sum(r**2)
    n_iter: int
    converged: bool
    message: str


def _jacobian(fun, x, r0, args, rel_step):
    J = np.empty((r0.size, x.size), dtype=np.float64)
    for k in range(x.size):
        h = rel_step * max(abs(x[k]), 1.0)
        xp = x.copy()
        xp[k] += h
        J[:, k] = (fun(xp, *args) - r0) / h
    return J


def least_squares_lm(
    fun: Callable[..., np.ndarray],
    x0,
    args: tuple = (),
    max_iter: int = 300,
    ftol: float = 1e-14,
    xtol: float = 1e-14,
    rel_step: float = 1e-8,
    verbose: bool = False,
) -> LMResult:
    """Minimize 0.5·||fun(x, *args)||² over x. `fun` must return a 1-D residual array."""
    x = np.asarray(x0, dtype=np.float64).copy()
    r = np.asarray(fun(x, *args), dtype=np.float64)
    cost = float(r @ r)
    lam = 1e-3
    message = f"hit max_iter={max_iter} without meeting ftol/xtol"
    converged = False
    it = 0

    for it in range(1, max_iter + 1):
        J = _jacobian(fun, x, r, args, rel_step)
        JtJ = J.T @ J
        g = J.T @ r
        scale = np.maximum(np.diag(JtJ), 1e-30)

        accepted = False
        for _ in range(40):  # inner loop: grow λ until the step is downhill
            try:
                step = np.linalg.solve(JtJ + lam * np.diag(scale), -g)
            except np.linalg.LinAlgError:
                lam *= 10.0
                continue
            x_new = x + step
            r_new = np.asarray(fun(x_new, *args), dtype=np.float64)
            cost_new = float(r_new @ r_new)
            if np.isfinite(cost_new) and cost_new < cost:
                d_cost = cost - cost_new
                d_x = np.max(np.abs(step) / np.maximum(np.abs(x), 1.0))
                x, r, cost = x_new, r_new, cost_new
                lam = max(lam / 10.0, 1e-12)
                accepted = True
                if d_cost <= ftol * max(cost, 1e-300):
                    converged, message = True, f"ftol reached (Δcost={d_cost:.3e})"
                elif d_x <= xtol:
                    converged, message = True, f"xtol reached (Δx={d_x:.3e})"
                break
            lam *= 10.0
        else:
            message = "λ grew past 10⁴⁰ without finding a downhill step"

        if verbose:
            print(f"    iter {it:3d}  cost={cost:.6e}  lambda={lam:.2e}"
                  f"{'' if accepted else '  (no step accepted)'}")
        if converged or not accepted:
            if not accepted and not converged:
                converged, message = True, "no further downhill step exists (local minimum)"
            break

    return LMResult(x=x, cost=0.5 * cost, n_iter=it, converged=converged, message=message)


if __name__ == "__main__":
    # Self-test: recover known parameters of the exact model this project fits
    # (LCC 2SP + 2D similarity) from noise-free synthetic data.
    from .lcc import GRS80_A, GRS80_E2, lcc_2sp_forward, similarity_2d

    truth = np.array([24.5, 41.9, 43.7, 0.004, 1.00013, 8_500_000.0, 4_600_000.0])
    rng = np.random.default_rng(7)
    lat = np.radians(rng.uniform(41.3, 44.2, 4000))
    lon = np.radians(rng.uniform(22.4, 28.6, 4000))

    def model(p, lat_r, lon_r):
        lon0, lat1, lat2, th, s, tx, ty = p
        x, y = lcc_2sp_forward(lat_r, lon_r, np.radians(42.6678756833333),
                               np.radians(lat1), np.radians(lat2), np.radians(lon0),
                               a=GRS80_A, e2=GRS80_E2)
        return similarity_2d(x, y, th, s, tx, ty)

    E, N = model(truth, lat, lon)

    def resid(p, lat_r, lon_r, Eo, No):
        Em, Nm = model(p, lat_r, lon_r)
        return np.concatenate([Em - Eo, Nm - No])

    x0 = np.array([25.5, 42.0, 43.3, 0.0, 1.0, 8_400_000.0, 4_500_000.0])
    res = least_squares_lm(resid, x0, args=(lat, lon, E, N), verbose=True)
    Ef, Nf = model(res.x, lat, lon)
    d = np.hypot(Ef - E, Nf - N)
    print(f"\n  converged={res.converged} ({res.message}) in {res.n_iter} iters")
    print(f"  max positional error after recovery: {d.max():.6e} m")
    assert d.max() < 1e-6, "LM failed to recover a noise-free synthetic fit"
    print("  OK -- least_squares_lm recovers the exact model to sub-micrometre level.")
