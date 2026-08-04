#!/usr/bin/env python3
"""AUDIT 1: verify the smoothing threshold actually propagated into all 64 E0 products.

The launch previously mutated vdp.DEFAULT_SMOOTH_SUPPORT_THRESHOLD after import, but
generate_probability_maps() binds `smooth_support_threshold=DEFAULT_SMOOTH_SUPPORT_THRESHOLD` as a
DEFAULT ARGUMENT -- evaluated once at definition time -- so the mutation was a no-op and every
candidate would have been built at the import-time default. Fixed by explicit keyword; this audit
proves it from the stored metadata rather than trusting the fix.

Reports a 16 x 4 table of requested vs stored threshold, and asserts every OTHER configuration
field is identical within each center.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
THRS = [2, 3, 5, 10]
CENTERS = [(0,0),(20,-12),(10,2),(-10,-2),(0,8),(50,-1),(50,1),(-50,1),
           (90,0),(140,0),(-140,0),(0,25),(0,-25),(0,50),(0,-50),(-50,-50)]
CFG_KEYS = ["smooth_density_maps", "smooth_support_scale_by_clone_factor", "smooth_sigma_pixels",
            "smooth_truncate_sigma", "max_sep_deg", "center_lon_deg", "center_lat_deg",
            "x_grid", "y_grid", "population_names", "mag_bin_mins", "mag_bin_maxs"]
rows, cfg_fail, missing = [], [], []
for dlon, lat in CENTERS:
    name = f"prob_maps_grid_dlon{dlon:+04d}_lat{lat:+03d}.npz"
    rec = {"center": f"{dlon:+04d}/{lat:+03d}"}
    ref = None
    for t in THRS:
        p = W/f"prob_maps_e0_thr{t}"/name
        if not p.exists():
            rec[f"thr{t}"] = "MISSING"; missing.append(str(p)); continue
        z = np.load(p, allow_pickle=True)
        stored = float(z["smooth_support_threshold"]) if "smooth_support_threshold" in z.files else None
        rec[f"thr{t}"] = stored
        cfg = {k: (np.asarray(z[k]).tobytes() if k in z.files else None) for k in CFG_KEYS}
        if ref is None:
            ref = cfg
        else:
            diff = [k for k in CFG_KEYS if cfg[k] != ref[k]]
            if diff:
                cfg_fail.append((rec["center"], t, diff))
    rows.append(rec)
d = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("=== requested vs STORED smoothing threshold (16 centers x 4 candidates) ===")
print(d.to_string(index=False))
ok_vals = True
for t in THRS:
    col = d[f"thr{t}"]
    good = col.apply(lambda v: isinstance(v, float) and abs(v - t) < 1e-9)
    print(f"  threshold {t:>2}: {int(good.sum())}/16 centers store exactly {t}")
    ok_vals &= bool(good.all())
distinct_ok = all(len({d.loc[i, f'thr{t}'] for t in THRS}) == 4 for i in d.index)
print(f"\n  every center has 4 DISTINCT thresholds: {distinct_ok}")
print(f"  other config fields identical within center: {not cfg_fail}")
for c, t, diff in cfg_fail[:5]:
    print(f"    MISMATCH {c} thr{t}: {diff}")
print(f"  missing products: {len(missing)}")
ok = ok_vals and distinct_ok and not cfg_fail and not missing
print(f"\n{'='*60}\nAUDIT 1: {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
