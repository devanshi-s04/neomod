#!/usr/bin/env python3
"""T3 - map regression: monolithic vs HEALPix cache must give BIT-IDENTICAL maps (docs §8.6).

Builds the SAME (center, mag bin) twice with n_jobs=1 -- once reading the monolithic cache, once the
HEALPix-partitioned cache -- and requires every density/support/nearest array to match exactly.
n_jobs=1 removes the documented joblib parallel non-determinism so "bit-identical" is a legitimate
requirement; T2 already proved the input clone sets are element-wise identical, so ANY difference
here is a real wiring bug.
"""
import os, sys, time, copy
from pathlib import Path
import numpy as np, pandas as pd
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
from astropy.time import Time
from astropy.coordinates import get_sun, GeocentricTrueEcliptic

EPOCH="2027-08-25T00:00:00"; MAGLO,MAGHI=23.0,24.0
HPX=str(W/"outputs/neomod3_projection_cache/by_pixel")

def build(use_hpx):
    for m in list(sys.modules):
        if m.startswith("velocity_density_pipeline_neomod"): del sys.modules[m]
    os.environ["NEOMOD3_CACHE_DIR"] = HPX if use_hpx else ""
    import velocity_density_pipeline_neomod_clone_only as v
    t=Time(EPOCH,scale="utc"); sun=get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
    clon=((sun.lon.deg+180.0)%360.0 + 20.0)%360.0; clat=-12.0
    cache=pd.read_parquet(W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet")
    _, scorer = v.load_s3m_population("neo", verbose=False)
    ps=copy.deepcopy(v.DEFAULT_POPULATION_SETTINGS); cs={}
    for pop,cfg in ps.items():
        sub=cache[cache["population"]==pop].reset_index(drop=True)
        cs[pop]=dict(df=sub, scorer=scorer, clone_factor=cfg["clone_factor"],
            use_conditional_cloner=True, scatter_size=cfg["scatter_size"],
            scatter_alpha=cfg["scatter_alpha"], _mag_app=sub["mag_app"].to_numpy(float))
    x,y,X0,Y0,gp = v.make_default_grid(grid_lim=(-5.0,5.0), grid_step=0.01)
    t0=time.time()
    res=v.build_cloned_maps_for_center_magbin(
        center_lon_deg=clon, center_lat_deg=clat, center_label="T3",
        clone_sources=cs, obstime_str=EPOCH, max_sep_deg=30.0,
        rng=np.random.default_rng(42), overlay_max=200_000,
        x_grid=x, y_grid=y, grid_points=gp, X0=X0,
        k_map=v.DEFAULT_K_MAP, n_d0_grid_map=v.DEFAULT_N_D0_GRID_MAP,
        mag_min=MAGLO, mag_max=MAGHI,
        smooth_density_maps=True, smooth_population_names=("NEO",),
        smooth_support_threshold=v.DEFAULT_SMOOTH_SUPPORT_THRESHOLD,
        smooth_sigma_pixels=v.DEFAULT_SMOOTH_SIGMA_PIXELS,
        smooth_truncate_sigma=v.DEFAULT_SMOOTH_TRUNCATE_SIGMA,
        smooth_support_scale_by_clone_factor=True, smooth_presmoothing_passes=None,
        n_jobs=1)
    print(f"  [{'HEALPIX' if use_hpx else 'MONOLITHIC'}] built in {time.time()-t0:.0f}s", flush=True)
    return res

print("=== building MONOLITHIC (reference) ==="); A=build(False)
print("=== building HEALPIX ==="); B=build(True)
print("\n=== T3a bit-identical comparison (n_jobs=1) ===")
ok=True
for key in ["density_maps_downweighted_raw","support_count_maps","nearest_dist_maps"]:
    for pop in ["MBA","NEO","TNO","Trojans"]:
        a=np.asarray(A[key][pop]); b=np.asarray(B[key][pop])
        same=a.shape==b.shape and np.array_equal(np.nan_to_num(a,posinf=1e30),np.nan_to_num(b,posinf=1e30))
        ok&=same
        d=np.nanmax(np.abs(np.nan_to_num(a,posinf=0)-np.nan_to_num(b,posinf=0))) if a.shape==b.shape else float('nan')
        print(f"  {key:34s} {pop:>8}: identical={str(same):>5}  max|diff|={d:.3e}")
print(f"  magcut_counts: {A['magcut_counts']} vs {B['magcut_counts']}")
print(f"\n{'='*60}\nT3a OVERALL: {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
