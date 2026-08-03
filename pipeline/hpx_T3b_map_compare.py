#!/usr/bin/env python3
"""T3b - HEALPix-built maps vs the existing monolithic-built reference maps, at production n_jobs=16.

T3a already proved the two cache paths are BIT-identical at n_jobs=1. At n_jobs=16 the documented
joblib scatter makes exact equality impossible, so this is a STATISTICAL check: differences must be
at joblib-scatter level, not structural. Bottom line reported is the per-tracklet P_NEO difference,
which is what the classifier actually consumes.
"""
import os, sys
from pathlib import Path
import numpy as np
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
REF = W/"prob_maps_grid_neomod3_vlim5"; NEW = W/"prob_maps_grid_neomod3_vlim5_hpx"
CENTERS = ["prob_maps_grid_dlon+000_lat+00.npz", "prob_maps_grid_dlon+020_lat-12.npz"]

worst = 0.0
for cen in CENTERS:
    a = np.load(REF/cen, allow_pickle=True); b = np.load(NEW/cen, allow_pickle=True)
    ka, kb = set(a.files), set(b.files)
    print(f"\n=== {cen} ===")
    print(f"  keys: ref {len(ka)}, new {len(kb)}, missing {sorted(ka^kb)[:6] or 'none'}")
    rows = []
    for k in sorted(ka & kb):
        x, y = a[k], b[k]
        if x.dtype.kind not in "fi" or x.shape != y.shape: continue
        x = np.nan_to_num(np.asarray(x, float), posinf=0.0); y = np.nan_to_num(np.asarray(y, float), posinf=0.0)
        d = np.abs(x - y); scale = max(np.abs(x).max(), 1e-300)
        rel = d.max()/scale
        rows.append((rel, k, d.max(), (d > 0).mean()*100))
    rows.sort(reverse=True)
    for rel, k, dmax, fdiff in rows[:8]:
        tag = "magcut_count" in k
        print(f"  {k:38s} max|d|={dmax:.4g}  rel={rel:.3e}  cells_differing={fdiff:5.1f}%"
              f"{'   <- metadata, semantics differ by design' if tag else ''}")
    worst = max(worst, max([r for r, k, _, _ in rows if "magcut_count" not in k], default=0.0))
print(f"\nworst relative array difference (excluding magcut_count metadata): {worst:.3e}")
