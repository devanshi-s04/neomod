#!/usr/bin/env python3
"""BOUNDED TNO numerical-resolution diagnostic. NOT part of threshold selection.

Check A reported R = int(rho)/N_full of 3-20 for TNO while MBA sat at ~1.14. Hypothesis: TNO
velocities span ~0.01-0.05 deg/day, comparable to ONE 0.01 deg/day grid pixel, so summing rho*dA
over a hyper-concentrated distribution mis-integrates. If R converges toward the MBA value as the
grid refines, it is a coarse-grid integration artifact; if it stays high, the map is defective.

Rebuilds ONE representative (center, mag bin) over a SMALL velocity window at 0.01 / 0.005 /
0.0025 deg/day and reports integral, convergence, and the VDP-score change for CAL objects there.
"""
import copy, json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
os.environ.setdefault("NEOMOD3_CACHE_DIR", str(W/"outputs/neomod3_projection_cache/by_pixel"))
import velocity_density_pipeline_neomod_clone_only as v
from astropy.time import Time
from astropy.utils import iers; iers.conf.auto_max_age = None
from astropy.coordinates import get_sun, GeocentricTrueEcliptic

EPOCH = "2027-08-25T00:00:00"
prov = json.loads((W/"outputs/splits/split_provenance.json").read_text())
v.NONNEO_SPLIT_FRACTIONS = prov
BINS = {"14_16":(14.,16.),"16_18":(16.,18.),"18_20":(18.,20.),"mag20":(20.,21.),
        "mag21":(21.,22.),"mag22":(22.,23.),"mag23":(23.,24.),"mag24+":(24.,25.)}

# ---- 1. find the worst TNO (center, bin) from the existing threshold-10 maps ----
CENTERS = [(0,0),(20,-12),(10,2),(-10,-2),(0,8),(50,-1)]
worst = None
for dl, la in CENTERS:
    z = np.load(W/"prob_maps_e0_thr10"/f"prob_maps_grid_dlon{dl:+04d}_lat{la:+03d}.npz", allow_pickle=True)
    x, y = z["x_grid"], z["y_grid"]; dA = float((x[1]-x[0])*(y[1]-y[0]))
    for b in BINS:
        k = f"density_unsmoothed__TNO__{b}"; sk = f"support_count__TNO__{b}"
        if k not in z.files: continue
        sup = np.nan_to_num(np.asarray(z[sk], float)).sum()
        if sup < 1000: continue
        f = prov["f_by_population_magbin"]["TNO"].get(b) or prov["f_gen_by_population"]["TNO"]
        R = (np.nan_to_num(np.asarray(z[k], float)).sum()*dA)/(sup/f)
        if worst is None or R > worst["R"]:
            worst = dict(dlon=dl, lat=la, bin=b, R=R, support=sup, f=f)
print(f"worst TNO bin: center ({worst['dlon']:+d},{worst['lat']:+d}) bin {worst['bin']}  "
      f"R={worst['R']:.3f}  GEN support={worst['support']:,.0f}  f={worst['f']:.4f}")

# ---- 2. TNO velocity extent there -> choose a bounded window ----
t = Time(EPOCH, scale="utc")
sun = get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
clon = (sun.lon.deg + 180.0 + worst["dlon"]) % 360.0
man = pd.read_parquet(W/"outputs/splits/nonneo_split_manifest.parquet")
gen_ids = set(man.ObjID[man.split == "GEN"])
cache = pd.read_parquet(W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet")
cache = cache[cache.ObjID.isin(gen_ids) | (cache.population == "NEO")].reset_index(drop=True)
_, scorer = v.load_s3m_population("neo", verbose=False)
mlo, mhi = BINS[worst["bin"]]
tn = cache[(cache.population == "TNO") & (cache.mag_app >= mlo) & (cache.mag_app < mhi)]
vext = float(np.nanmax(np.maximum(tn.vlam.abs(), tn.vbeta.abs()))) if len(tn) else 0.1
LIM = max(0.15, min(1.0, 2.5*vext))
print(f"TNO |v| max in that bin (all-sky GEN): {vext:.4f} deg/day -> window +-{LIM:.3f}")

def build(step):
    cs = {}
    for pop, cfg in copy.deepcopy(v.DEFAULT_POPULATION_SETTINGS).items():
        sub = cache[cache.population == pop].reset_index(drop=True)
        cs[pop] = dict(df=sub, scorer=scorer, clone_factor=cfg["clone_factor"],
                       use_conditional_cloner=True, scatter_size=cfg["scatter_size"],
                       scatter_alpha=cfg["scatter_alpha"], _mag_app=sub["mag_app"].to_numpy(float))
    x, y, X0, Y0, gp = v.make_default_grid(grid_lim=(-LIM, LIM), grid_step=step)
    r = v.build_cloned_maps_for_center_magbin(
        center_lon_deg=clon, center_lat_deg=float(worst["lat"]), center_label=f"tno{step}",
        clone_sources=cs, obstime_str=EPOCH, max_sep_deg=30.0,
        rng=np.random.default_rng(1), overlay_max=1000,
        x_grid=x, y_grid=y, grid_points=gp, X0=X0,
        k_map=v.DEFAULT_K_MAP, n_d0_grid_map=v.DEFAULT_N_D0_GRID_MAP,
        mag_min=mlo, mag_max=mhi,
        smooth_density_maps=True, smooth_population_names=("NEO",),
        smooth_support_threshold=10.0, smooth_sigma_pixels=v.DEFAULT_SMOOTH_SIGMA_PIXELS,
        smooth_truncate_sigma=v.DEFAULT_SMOOTH_TRUNCATE_SIGMA,
        smooth_support_scale_by_clone_factor=False, smooth_presmoothing_passes=None, n_jobs=16)
    return r, x, y

print(f"\n{'step':>8} {'npix':>10} {'int(rho)_TNO':>14} {'GEN support':>12} {'N_full':>10} {'R':>8}")
out, dens = [], {}
for step in (0.01, 0.005, 0.0025):
    r, x, y = build(step)
    dA = float((x[1]-x[0])*(y[1]-y[0]))
    u = np.nan_to_num(np.asarray(r["density_maps_unsmoothed"]["TNO"], float))
    sup = np.nan_to_num(np.asarray(r["support_count_maps"]["TNO"], float)).sum()
    integ = u.sum()*dA
    Nfull = sup/worst["f"]
    R = integ/Nfull if Nfull else np.nan
    out.append(dict(step=step, npix=u.size, int_rho=integ, gen_support=sup, N_full=Nfull, R=R))
    dens[step] = {p: np.nan_to_num(np.asarray(r["density_maps_unsmoothed"][p], float))
                  for p in ["NEO","MBA","TNO","Trojans"]}
    dens[step]["_grid"] = (x, y)
    print(f"{step:>8.4f} {u.size:>10,} {integ:>14,.1f} {sup:>12,.0f} {Nfull:>10,.0f} {R:>8.3f}")

df = pd.DataFrame(out)
print(f"\nconvergence: R = {' -> '.join(f'{r:.3f}' for r in df.R)}")
print(f"  MBA reference (Check A median) = 1.1391 ; expected b = k/(k-1) = {10/9:.4f}")
print(f"  |R(0.0025) - b| = {abs(df.R.iloc[-1]-10/9):.3f}   |R(0.01) - b| = {abs(df.R.iloc[0]-10/9):.3f}")
converging = abs(df.R.iloc[-1]-10/9) < abs(df.R.iloc[0]-10/9)
print(f"  VERDICT: {'coarse-grid INTEGRATION ARTIFACT (R converges toward b)' if converging else 'NOT converging -- possible real map defect'}")

# ---- 3. VDP score change for CAL objects in that TNO region ----
cal = pd.read_parquet(W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet")
cen_name = f"prob_maps_grid_dlon{worst['dlon']:+04d}_lat{worst['lat']:+03d}.npz"
sel = cal[(cal.prob_map_file == cen_name) & (cal.mean_mag >= mlo) & (cal.mean_mag < mhi) &
          (cal.vlam.abs() <= LIM) & (cal.vbeta.abs() <= LIM)]
print(f"\nCAL rows in that TNO region: {len(sel):,}  (NEO {int((sel.population=='NEO').sum())})")
if len(sel):
    pn = {}
    for step in dens:
        x, y = dens[step]["_grid"]
        ix = np.clip(np.searchsorted(x, sel.vlam.to_numpy()) - 1, 0, len(x)-1)
        iy = np.clip(np.searchsorted(y, sel.vbeta.to_numpy()) - 1, 0, len(y)-1)
        tot = sum(dens[step][p][iy, ix] for p in ["NEO","MBA","TNO","Trojans"])
        pn[step] = np.where(tot > 0, dens[step]["NEO"][iy, ix]/np.where(tot > 0, tot, 1), 0.0)
    base = pn[0.01]
    for step in (0.005, 0.0025):
        dd = pn[step] - base
        print(f"  P(NEO) change 0.01 -> {step}:  mean|d|={np.abs(dd).mean():.3e}  "
              f"max|d|={np.abs(dd).max():.3e}  n_changed={int((dd!=0).sum())}/{len(dd)}")
df.to_csv(W/"outputs/e0_results/TNO_RESOLUTION_DIAGNOSTIC.csv", index=False)
json.dump(dict(worst=worst, window=LIM, rows=out, converging=bool(converging),
               b_expected=10/9, mba_reference=1.1391),
          open(W/"outputs/e0_results/TNO_RESOLUTION_DIAGNOSTIC.json","w"), indent=2, default=str)
print(f"\nwrote TNO_RESOLUTION_DIAGNOSTIC.csv/json")
