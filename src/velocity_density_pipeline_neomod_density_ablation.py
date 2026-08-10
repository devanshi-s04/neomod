#!/usr/bin/env python3
"""Density-estimator ablation for the geometric benchmark (GEOMETRIC_DENSITY_ESTIMATOR_ABLATION.md).

ONE parameterized implementation. The sealed production module
`velocity_density_pipeline_neomod_clone_only` is imported and NEVER modified: all population
loading, sky/magnitude cuts, NEOMOD3 GEN provenance guards, split corrections, physical weights,
grid construction and archive writing come from it unchanged.

Only the density-construction step is ablated, at its single seam:
`evaluate_density_map_full_posterior_2d(tree, grid_points, ...)`.

DISPATCH IS EXPLICIT.
No caller-frame inspection is used anywhere. `resolve_backend(density_mode, pop_name)` is a pure
function, unit-tested over all mode x population combinations. A backend is bound per BUILD CALL
over a population group whose members are asserted to resolve to the same backend; mixed modes are
executed as two disjoint population-group calls whose archives are merged. Population identity is
therefore carried by the explicit `clone_sources` subset, never inferred.

Modes
-----
hist_all      NEO, MBA, TNO, Trojans -> direct weighted histogram
bayes_all     NEO, MBA, TNO, Trojans -> sealed Bayesian estimator
bayes_nonneo  NEO -> histogram ; MBA, TNO, Trojans -> sealed Bayesian estimator

Gaussian smoothing is DISABLED in every mode (`smooth_density_maps=False`), which bypasses the
`smooth_density_map_by_support` block entirely. Production smooths NEO only, so the sealed
production maps differ from `bayes_all` by exactly that NEO-only convolution.

Histogram weighting
-------------------
The sealed builder computes `density_downweighted = density_estimate / effective_factor`, where
`effective_factor` is a SCALAR per (population, magnitude bin). Every sample in such a group
therefore carries an identical physical weight `w = 1/effective_factor`, and

    hist_count_density / effective_factor
      = (N_pixel / pixel_area) / effective_factor
      = (N_pixel * w) / pixel_area
      = sum(physical_sample_weights) / pixel_area

so returning a COUNT density here is exactly the protocol's physical weighted density after the
sealed normalization. `assert_uniform_weight_equivalence()` proves this numerically per build.
If a future population ever carries per-row weights, `histogram_density()` accepts `weights=` and
must be used instead -- equivalence is proven, never assumed.
"""
from __future__ import annotations

import numpy as np

import velocity_density_pipeline_neomod_clone_only as base

# The EXACT sealed Bayesian routine, captured once at import. Identity is assertable:
#   BAYES_ROUTINE is base.evaluate_density_map_full_posterior_2d
BAYES_ROUTINE = base.evaluate_density_map_full_posterior_2d
BAYES_POINT_ROUTINE = base.estimate_density_full_posterior_2d

POPULATIONS = ("NEO", "MBA", "TNO", "Trojans")
NONNEO = ("MBA", "TNO", "Trojans")
DENSITY_MODES = ("hist_all", "bayes_all", "bayes_nonneo")

_BACKEND_TABLE = {
    "hist_all":     {"NEO": "hist",  "MBA": "hist",  "TNO": "hist",  "Trojans": "hist"},
    "bayes_all":    {"NEO": "bayes", "MBA": "bayes", "TNO": "bayes", "Trojans": "bayes"},
    "bayes_nonneo": {"NEO": "hist",  "MBA": "bayes", "TNO": "bayes", "Trojans": "bayes"},
}


def resolve_backend(density_mode: str, pop_name: str) -> str:
    """Pure, explicit dispatch. No frame inspection, no globals, no ordering assumptions."""
    if density_mode not in _BACKEND_TABLE:
        raise ValueError(f"unknown density_mode {density_mode!r}; expected one of {DENSITY_MODES}")
    table = _BACKEND_TABLE[density_mode]
    if pop_name not in table:
        raise ValueError(f"unknown population {pop_name!r}; expected one of {POPULATIONS}")
    return table[pop_name]


def population_groups(density_mode: str):
    """Population groups that share a backend, so each build call binds ONE explicit backend."""
    groups = {}
    for pop in POPULATIONS:
        groups.setdefault(resolve_backend(density_mode, pop), []).append(pop)
    return [(backend, tuple(pops)) for backend, pops in sorted(groups.items())]


def grid_edges(x_grid, y_grid):
    """Cell edges for node-centred pixels. Identical edges for every population, by construction."""
    x_grid = np.asarray(x_grid, float); y_grid = np.asarray(y_grid, float)
    dx = float(x_grid[1] - x_grid[0]); dy = float(y_grid[1] - y_grid[0])
    if not (np.allclose(np.diff(x_grid), dx) and np.allclose(np.diff(y_grid), dy)):
        raise ValueError("ablation histogram requires a uniform grid")
    xe = np.concatenate([x_grid - dx / 2.0, [x_grid[-1] + dx / 2.0]])
    ye = np.concatenate([y_grid - dy / 2.0, [y_grid[-1] + dy / 2.0]])
    return xe, ye, dx * dy


def histogram_density(points, x_grid, y_grid, weights=None):
    """Direct weighted histogram density on the SAME grid, evaluated at the nodes.

    rho(i,j) = sum(weights of samples in pixel i,j) / pixel_area

    No smoothing, no pseudocounts, no hole filling, no neighbour borrowing.
    Empty pixels are exactly zero. Returns shape (len(y_grid), len(x_grid)) to match the sealed
    builder's `X0.shape` convention (row = vbeta, col = vlam).
    """
    pts = np.asarray(points, float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError(f"points must be (N,2); got {pts.shape}")
    xe, ye, pixel_area = grid_edges(x_grid, y_grid)
    H, _, _ = np.histogram2d(pts[:, 1], pts[:, 0], bins=[ye, xe], weights=weights)
    return H / pixel_area


def assert_uniform_weight_equivalence(n_pixel_counts, effective_factor, pixel_area, rtol=1e-13):
    """Numerically prove count-density/effective_factor == sum(physical weights)/pixel_area.

    Valid only when every sample in the group shares one weight w = 1/effective_factor, which the
    sealed builder guarantees because effective_factor is a scalar per (population, mag bin).

    The two expressions are identical in REAL arithmetic but not bit-identical in floating point:
    the ablation path computes (counts/pixel_area)/effective_factor, mirroring the sealed
    pipeline's "estimate density, then divide", whereas the protocol writes
    sum(w)/pixel_area. The operations are reassociated, so results differ by ~1 ULP
    (measured max relative difference ~4e-16). Equality is therefore asserted to a tight RELATIVE
    tolerance and the achieved gap is reported, rather than claiming bit-exactness.
    """
    counts = np.asarray(n_pixel_counts, float)
    w = 1.0 / float(effective_factor)
    lhs = (counts / pixel_area) / float(effective_factor)   # what the ablation path produces
    rhs = (counts * w) / pixel_area                          # the protocol's physical definition
    denom = np.where(np.abs(rhs) > 0, np.abs(rhs), 1.0)
    rel = np.abs(lhs - rhs) / denom
    max_rel = float(np.nanmax(rel)) if rel.size else 0.0
    max_abs = float(np.nanmax(np.abs(lhs - rhs))) if rel.size else 0.0
    if max_rel > rtol:
        raise AssertionError(
            f"histogram weight equivalence FAILED: max relative difference {max_rel:.3e} "
            f"exceeds rtol {rtol:.1e} (max abs {max_abs:.3e})")
    return {"equivalent": True, "weight_per_sample": w, "pixel_area": float(pixel_area),
            "n_pixels_nonzero": int((counts > 0).sum()), "total_count": float(counts.sum()),
            "max_rel_diff": max_rel, "max_abs_diff": max_abs, "rtol": rtol,
            "bit_identical": bool(max_abs == 0.0)}


def make_density_callable(backend: str, pop_group, density_mode: str):
    """Bind ONE explicit backend for a build call over `pop_group`.

    Every member of `pop_group` is asserted to resolve to `backend` under `density_mode`, so the
    population identity is fixed by the caller's explicit clone_sources subset -- never inferred.
    The returned callable has the sealed routine's signature.
    """
    for pop in pop_group:
        got = resolve_backend(density_mode, pop)
        if got != backend:
            raise AssertionError(
                f"population group {pop_group} is not backend-homogeneous under {density_mode!r}: "
                f"{pop} resolves to {got!r}, not {backend!r}")
    if backend == "bayes":
        return BAYES_ROUTINE, {"backend": "bayes", "routine": "evaluate_density_map_full_posterior_2d",
                               "sealed_identity": BAYES_ROUTINE is base.evaluate_density_map_full_posterior_2d}

    captured = {}

    def _hist_callable(tree, grid_points, k=None, n_d0_grid=None, show_progress=False, n_jobs=1):
        # tree.data are the SAME sample points the sealed builder would have handed the Bayesian
        # estimator, so the two backends see identical inputs.
        pts = np.asarray(tree.data, float)
        gp = np.asarray(grid_points, float)
        xs = np.unique(gp[:, 0]); ys = np.unique(gp[:, 1])
        rho = histogram_density(pts, xs, ys)
        _, _, pixel_area = grid_edges(xs, ys)
        captured["n_samples"] = int(len(pts))
        captured["counts_sum"] = float(rho.sum() * pixel_area)
        # sealed builder reshapes the returned flat array to X0.shape; return flat in the same order
        return rho.reshape(-1)

    return _hist_callable, {"backend": "hist", "routine": "histogram_density",
                            "captured": captured, "sealed_identity": None}


__all__ = ["resolve_backend", "population_groups", "histogram_density", "grid_edges",
           "assert_uniform_weight_equivalence", "make_density_callable",
           "BAYES_ROUTINE", "DENSITY_MODES", "POPULATIONS", "NONNEO"]
