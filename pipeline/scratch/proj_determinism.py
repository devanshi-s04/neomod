#!/usr/bin/env python3
"""Is build_visible_subset_dataframe deterministic run-to-run IN THE SAME ENVIRONMENT?
If yes, the 1-ULP difference vs the stored shard is drift since the original build (ephemeris/IERS
or library version), not non-determinism. If no, the projection itself is non-reproducible.
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
import velocity_density_pipeline_gmm as vdp
import neomod3_sampler as nm3s
EPOCH = "2027-08-25T00:00:00"
rng = np.random.default_rng(42)
df = nm3s.sample_neomod3_orbits(50_000, EPOCH, rng=rng)
_, scorer = vdp.load_s3m_population("neo", verbose=False)
outs = []
for k in range(2):
    v = vdp.build_visible_subset_dataframe(df.copy(), obstime_str=EPOCH, scorer=scorer,
        max_sep_deg=180.0, chunk=200_000, show_progress=False,
        center_mode="custom_ecliptic", center_lon_deg=0.0, center_lat_deg=0.0)
    outs.append(v)
a, b = outs
print(f"rows: {len(a):,} vs {len(b):,}")
cols = [c for c in a.columns if a[c].dtype.kind == "f"]
worst = 0.0; worstc = ""
for c in cols:
    d = np.nanmax(np.abs(a[c].to_numpy(float) - b[c].to_numpy(float)))
    if d > worst: worst, worstc = d, c
print(f"same-process, same-env, two runs: worst |diff| = {worst:.6g}  ({worstc})")
print("VERDICT:", "DETERMINISTIC run-to-run" if worst == 0 else "NON-DETERMINISTIC run-to-run")
import astropy, erfa
from astropy.utils import iers
print(f"\nastropy {astropy.__version__}  erfa {erfa.__version__}")
try:
    print("IERS_Auto cache:", iers.IERS_Auto.iers_table is not None)
except Exception as e:
    print("IERS:", e)
