#!/usr/bin/env python3
"""Prove the frozen projection environment reproduces a fixed projection IDENTICALLY across
separate processes (E0 governance item 3). --write-ref then --check-ref in two distinct jobs.
"""
import argparse, hashlib, json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
import frozen_projection_env as fpe
info = fpe.activate()                      # MUST come before any coordinate work
print("frozen env active:", json.dumps({k: (v[:40] if isinstance(v, str) else v)
                                        for k, v in info.items()}, indent=2))
import velocity_density_pipeline_gmm as vdp
import neomod3_sampler as nm3s
EPOCH = "2027-08-25T00:00:00"
REF = W/"outputs/splits/frozen_projection_reference.parquet"
ap = argparse.ArgumentParser(); ap.add_argument("--write-ref", action="store_true")
a = ap.parse_args()
df = nm3s.sample_neomod3_orbits(20_000, EPOCH, rng=np.random.default_rng(777))
_, scorer = vdp.load_s3m_population("neo", verbose=False)
vis = vdp.build_visible_subset_dataframe(df, obstime_str=EPOCH, scorer=scorer, max_sep_deg=180.0,
    chunk=200_000, show_progress=False, center_mode="custom_ecliptic",
    center_lon_deg=0.0, center_lat_deg=0.0)
cols = sorted(c for c in vis.columns if vis[c].dtype.kind == "f")
digest = hashlib.sha256(np.ascontiguousarray(vis[cols].to_numpy(np.float64)).tobytes()).hexdigest()
print(f"\nrows {len(vis):,}   float64 digest {digest}")
if a.write_ref:
    vis.to_parquet(REF, index=False)
    (REF.with_suffix(".sha256")).write_text(digest)
    print(f"reference written -> {REF}")
else:
    want = REF.with_suffix(".sha256").read_text().strip()
    ref = pd.read_parquet(REF)
    exact = len(ref) == len(vis) and all(
        np.array_equal(ref[c].to_numpy(), vis[c].to_numpy(), equal_nan=True) for c in cols)
    print(f"reference digest {want}")
    print(f"digest match: {digest == want}   arrays bit-identical: {exact}")
    print("\nVERDICT:", "REPRODUCIBLE under the frozen environment"
          if (digest == want and exact) else "NOT REPRODUCIBLE")
    sys.exit(0 if (digest == want and exact) else 1)
