#!/usr/bin/env python3
"""Smoke test: build ONE mag bin of ONE center with the NEOMOD3 clone-only module."""
import os, sys, time
from pathlib import Path
import numpy as np
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
import copy, pandas as pd
import velocity_density_pipeline_neomod_clone_only as vnm
from astropy.time import Time
from astropy.coordinates import get_sun, GeocentricTrueEcliptic

EPOCH="2027-08-25T00:00:00"
t=Time(EPOCH,scale="utc"); sun=get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
anti=(sun.lon.deg+180.0)%360.0
clon, clat = (anti+20.0)%360.0, -12.0
print(f"center (dlon+020, lat-12) -> ({clon:.3f}, {clat})", flush=True)

cache=pd.read_parquet(W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet")
_, scorer = vnm.load_s3m_population("neo", verbose=False)
ps=copy.deepcopy(vnm.DEFAULT_POPULATION_SETTINGS)
clone_sources={}
for pop,cfg in ps.items():
    sub=cache[cache["population"]==pop].reset_index(drop=True)
    clone_sources[pop]=dict(df=sub, scorer=scorer, clone_factor=cfg["clone_factor"],
        use_conditional_cloner=True, scatter_size=cfg["scatter_size"],
        scatter_alpha=cfg["scatter_alpha"], _mag_app=sub["mag_app"].to_numpy(float))
    print(f"  {pop}: {len(sub):,}", flush=True)

x, y, X0, Y0, grid_points = vnm.make_default_grid(grid_lim=(-5.0, 5.0), grid_step=0.01)
print(f"grid {len(x)}x{len(y)} = {len(grid_points):,} eval points", flush=True)
t0=time.time()
res=vnm.build_cloned_maps_for_center_magbin(
    center_lon_deg=clon, center_lat_deg=clat, center_label="smoke",
    clone_sources=clone_sources, obstime_str=EPOCH, max_sep_deg=30.0,
    rng=np.random.default_rng(42), overlay_max=200_000,
    x_grid=x, y_grid=y, grid_points=grid_points, X0=X0,
    k_map=vnm.DEFAULT_K_MAP, n_d0_grid_map=vnm.DEFAULT_N_D0_GRID_MAP,
    mag_min=23.0, mag_max=24.0,
    smooth_density_maps=True, smooth_population_names=("NEO",),
    smooth_support_threshold=vnm.DEFAULT_SMOOTH_SUPPORT_THRESHOLD,
    smooth_sigma_pixels=vnm.DEFAULT_SMOOTH_SIGMA_PIXELS,
    smooth_truncate_sigma=vnm.DEFAULT_SMOOTH_TRUNCATE_SIGMA,
    smooth_support_scale_by_clone_factor=True, smooth_presmoothing_passes=None,
    n_jobs=16,
)
print(f"\nbuilt in {time.time()-t0:.0f}s", flush=True)
cell=(x[1]-x[0])*(y[1]-y[0])
print("\npop      magcut_count   int(rho)      support_sum   max|v| w/ support")
V=np.maximum(np.abs(X0),np.abs(Y0))
for p in ["MBA","NEO","TNO","Trojans"]:
    d=res["density_maps_downweighted_raw"][p]; sc=res["support_count_maps"][p]
    h=sc>=1
    print(f"{p:>8} {res['magcut_counts'][p]:>13,} {d.sum()*cell:>12,.1f} {sc.sum():>14,.0f} "
          f"{(V[h].max() if h.any() else 0):>18.2f}")
