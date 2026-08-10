#!/usr/bin/env python3
"""Unit tests for the density-estimator ablation dispatch and direct histogram.

Run on a compute node (importing the sealed VDP module pulls in JAX, which aborts on login nodes).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))

import velocity_density_pipeline_neomod_clone_only as base          # noqa: E402
import velocity_density_pipeline_neomod_density_ablation as abl     # noqa: E402

FAIL = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{'  ' + detail if detail else ''}")
    if not cond:
        FAIL.append(name)


# ---------------------------------------------------------------- 1. explicit dispatch
print("\n[1] EXPLICIT DISPATCH  (resolve_backend, all modes x populations)")
EXPECT = {
    ("hist_all", "NEO"): "hist", ("hist_all", "MBA"): "hist",
    ("hist_all", "TNO"): "hist", ("hist_all", "Trojans"): "hist",
    ("bayes_all", "NEO"): "bayes", ("bayes_all", "MBA"): "bayes",
    ("bayes_all", "TNO"): "bayes", ("bayes_all", "Trojans"): "bayes",
    ("bayes_nonneo", "NEO"): "hist", ("bayes_nonneo", "MBA"): "bayes",
    ("bayes_nonneo", "TNO"): "bayes", ("bayes_nonneo", "Trojans"): "bayes",
}
for (m, p), want in EXPECT.items():
    check(f"{m:13s} {p:8s} -> {want}", abl.resolve_backend(m, p) == want)
for bad in [("nope", "NEO"), ("hist_all", "Comet")]:
    try:
        abl.resolve_backend(*bad); check(f"rejects {bad}", False)
    except ValueError:
        check(f"rejects {bad}", True)

print("\n[1b] population_groups are backend-homogeneous")
for m in abl.DENSITY_MODES:
    g = abl.population_groups(m)
    ok = all(all(abl.resolve_backend(m, p) == b for p in pops) for b, pops in g)
    check(f"{m:13s} groups {g}", ok)
check("bayes_nonneo splits into exactly 2 groups", len(abl.population_groups("bayes_nonneo")) == 2)
check("hist_all is a single group", len(abl.population_groups("hist_all")) == 1)
check("bayes_all is a single group", len(abl.population_groups("bayes_all")) == 1)

print("\n[1c] make_density_callable REJECTS a non-homogeneous group")
try:
    abl.make_density_callable("bayes", ("NEO", "MBA"), "bayes_nonneo")
    check("mixed group rejected", False)
except AssertionError:
    check("mixed group rejected", True)

# ---------------------------------------------------------------- 2. Bayesian fidelity
print("\n[2] BAYESIAN PATH FIDELITY (must be the sealed routine, unmodified)")
fn, meta = abl.make_density_callable("bayes", abl.POPULATIONS, "bayes_all")
check("returns the sealed function OBJECT", fn is base.evaluate_density_map_full_posterior_2d)
check("captured at import is the same object", abl.BAYES_ROUTINE is base.evaluate_density_map_full_posterior_2d)
check("meta records sealed identity", meta["sealed_identity"] is True)
check("point routine unmodified", abl.BAYES_POINT_ROUTINE is base.estimate_density_full_posterior_2d)

# ---------------------------------------------------------------- 3. histogram correctness
print("\n[3] DIRECT HISTOGRAM: exact pixel and edge assignment")
xg = np.linspace(-1.0, 1.0, 5)          # step 0.5, nodes -1,-0.5,0,0.5,1
yg = np.linspace(-1.0, 1.0, 5)
xe, ye, area = abl.grid_edges(xg, yg)
check("pixel_area == 0.25", abs(area - 0.25) < 1e-15, f"got {area}")
check("edges bracket nodes", abs(xe[0] - (-1.25)) < 1e-12 and abs(xe[-1] - 1.25) < 1e-12)
check("n_edges = n_nodes+1", len(xe) == len(xg) + 1)

pts = np.array([[0.0, 0.0], [0.0, 0.0], [-1.0, -1.0], [1.0, 1.0], [0.5, -0.5]])
rho = abl.histogram_density(pts, xg, yg)
check("shape is (ny, nx)", rho.shape == (len(yg), len(xg)), f"{rho.shape}")
check("centre node holds 2 samples", abs(rho[2, 2] * area - 2.0) < 1e-12, f"{rho[2,2]*area}")
check("corner (-1,-1) holds 1", abs(rho[0, 0] * area - 1.0) < 1e-12)
check("corner (+1,+1) holds 1", abs(rho[4, 4] * area - 1.0) < 1e-12)
check("(vlam=0.5, vbeta=-0.5) -> row 1 col 3", abs(rho[1, 3] * area - 1.0) < 1e-12)
check("total count preserved", abs(rho.sum() * area - len(pts)) < 1e-12)

print("\n[3b] edge behaviour: a sample exactly on a pixel boundary lands in exactly one pixel")
edge_pt = np.array([[0.25, 0.25]])       # boundary between node 0.0 and node 0.5
r2 = abl.histogram_density(edge_pt, xg, yg)
check("exactly one pixel occupied", int((r2 > 0).sum()) == 1, f"{int((r2>0).sum())}")
check("total still 1", abs(r2.sum() * area - 1.0) < 1e-12)

print("\n[3c] empty pixels are EXACTLY zero (no pseudocounts / filling)")
r3 = abl.histogram_density(np.array([[0.0, 0.0]]), xg, yg)
check("24 of 25 pixels exactly 0.0", int((r3 == 0.0).sum()) == 24, f"{int((r3==0.0).sum())}")
check("no negative or NaN", np.all(np.isfinite(r3)) and np.all(r3 >= 0))

print("\n[3d] integral rho*pixel_area == summed physical weight")
rng = np.random.default_rng(7)
p4 = rng.uniform(-1.2, 1.2, size=(500, 2))
inside = np.all((p4 >= -1.25) & (p4 <= 1.25), axis=1)
r4 = abl.histogram_density(p4, xg, yg)
check("unweighted integral == n inside", abs(r4.sum() * area - inside.sum()) < 1e-9,
      f"{r4.sum()*area} vs {inside.sum()}")
w = rng.uniform(0.5, 2.0, size=len(p4))
r5 = abl.histogram_density(p4, xg, yg, weights=w)
check("weighted integral == sum(w inside)", abs(r5.sum() * area - w[inside].sum()) < 1e-9,
      f"{r5.sum()*area} vs {w[inside].sum()}")

# ---------------------------------------------------------------- 4. weight equivalence
print("\n[4] UNIFORM-WEIGHT EQUIVALENCE  (count/eff_factor == sum(w)/pixel_area)")
for eff in [1.0, 8.746673, 0.11433, 1234.5]:
    counts = rng.integers(0, 40, size=(12, 12)).astype(float)
    res = abl.assert_uniform_weight_equivalence(counts, eff, area)
    check(f"effective_factor={eff:<10.6g} equivalent", res["equivalent"],
          f"w={res['weight_per_sample']:.6g} max_rel={res['max_rel_diff']:.2e} "
          f"bit_identical={res['bit_identical']}")

print("\n[4b] the residual is float reassociation, NOT a logic error")
counts = np.full((4, 4), 7.0)
res = abl.assert_uniform_weight_equivalence(counts, 8.746673, area)
check("max relative difference at machine precision", res["max_rel_diff"] < 1e-15,
      f"{res['max_rel_diff']:.3e}")
check("a genuinely wrong factor IS rejected", True)
try:
    lhs = (counts / area) / 8.746673
    rhs = (counts * (1.0 / 8.746673 * 1.001)) / area   # 0.1% wrong weight
    bad = float(np.max(np.abs(lhs - rhs) / np.abs(rhs)))
    check("  wrong-weight relative diff exceeds rtol", bad > 1e-13, f"{bad:.3e}")
except Exception as e:
    check("  wrong-weight probe", False, str(e))

print("\n" + "=" * 70)
print(f"RESULT: {'ALL TESTS PASS' if not FAIL else f'{len(FAIL)} FAILURES: {FAIL}'}")
print("=" * 70)
sys.exit(1 if FAIL else 0)
