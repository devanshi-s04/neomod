#!/usr/bin/env python3
"""Mechanism test for the |v| 0.25-0.5 false positives (docs §9.4).

Chain under test:
  removing MBA cloning (clone_factor 5 -> 1) thins MBA support ~5x
  -> many more cells fall below support_mask_min=1
  -> the support mask ZEROES rho_MBA there, but NEO is EXEMPT (_support_mask_skip)
  -> P(NEO) = rho_NEO / sum(rho) -> 1 in the velocity band where real MBAs live
  -> false positives, concentrated exactly at |v| 0.25-0.5.
"""
import os, sys, numpy as np
W = "/mmfs1/gscratch/dirac/ds2004/sorcha"
sys.path.insert(0, f"{W}/neomod"); sys.path.insert(0, f"{W}/neomod/src"); os.chdir(f"{W}/neomod")
import velocity_density_pipeline_neomod_clone_only as v
for name in ("_support_mask_skip", "DEFAULT_SUPPORT_MASK_SKIP", "SUPPORT_MASK_SKIP"):
    if hasattr(v, name): print(f"support-mask skip list ({name}):", getattr(v, name))
cen = "prob_maps_grid_dlon+020_lat-12.npz"
for tag, dirn in (("production (MBA cf=5, ±2)", f"{W}/prob_maps_grid_s3m_nbody"),
                  ("NEOMOD3    (MBA cf=1, ±5)", f"{W}/prob_maps_grid_neomod3_full")):
    z = np.load(f"{dirn}/{cen}", allow_pickle=True)
    x, y = z["x_grid"], z["y_grid"]; X, Y = np.meshgrid(x, y)
    R = np.maximum(np.abs(X), np.abs(Y)); ann = (R > 0.25) & (R <= 0.5)
    sM = np.asarray(z["support_count__MBA__mag23"], float)
    sN = np.asarray(z["support_count__NEO__mag23"], float)
    pm = v.ProbMapSet.from_npz(f"{dirn}/{cen}", support_mask_min=1, mask_radius_deg_per_day=np.inf)
    p = pm.score_visible(X[ann].ravel(), Y[ann].ravel(), np.full(int(ann.sum()), 23.5))["NEO"]
    print(f"\n  {tag}")
    print(f"     annulus cells                   : {int(ann.sum()):,}")
    print(f"     MBA support < 1 (=> MBA ZEROED) : {100*(sM[ann]<1).mean():5.1f}%")
    print(f"     NEO support < 1 (but EXEMPT)    : {100*(sN[ann]<1).mean():5.1f}%")
    print(f"     P(NEO): median {np.nanmedian(p):.3f}  frac>0.5 {100*np.nanmean(p>0.5):5.1f}%"
          f"  frac>0.99 {100*np.nanmean(p>0.99):5.1f}%")
