#!/usr/bin/env python3
"""E0 acceptance checks A-D on the GEN pilot maps.
Criteria FROZEN in docs/E0_PILOT_PREREGISTRATION.md (sha 2608c472c04f2d81) before any output.
A density normalisation | B support is statistical | C zero leakage | D abstention/coverage
"""
import json, sys, glob, os
from pathlib import Path
import numpy as np, pandas as pd
import astropy.units as u
from astropy.time import Time
from astropy.utils import iers; iers.conf.auto_max_age = None
from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic, get_sun
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
pd.set_option("display.width", 250)
GEN = W/"prob_maps_grid_neomod3_GEN"
EPOCH = "2027-08-25T00:00:00"; MAXSEP = 30.0
BINS = [(14.,16.,"14_16"),(16.,18.,"16_18"),(18.,20.,"18_20"),(20.,21.,"mag20"),
        (21.,22.,"mag21"),(22.,23.,"mag22"),(23.,24.,"mag23"),(24.,25.,"mag24+")]
prov = json.load(open(W/"outputs/splits/split_provenance.json"))
man = pd.read_parquet(W/"outputs/splits/nonneo_split_manifest.parquet")
fails = []

t = Time(EPOCH, scale="utc")
sun = get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
cache = pd.read_parquet(W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet",
                        columns=["ObjID","population","lam_deg","beta_deg","mag_app"])
cache = cache[cache.population != "NEO"].reset_index(drop=True)
obj = SkyCoord(lon=cache.lam_deg.to_numpy()*u.deg, lat=cache.beta_deg.to_numpy()*u.deg,
               frame=GeocentricTrueEcliptic(obstime=t))
splitmap = dict(zip(man.ObjID, man.split))
cache["split"] = cache.ObjID.map(splitmap)

CENTERS = [(0,0),(20,-12),(10,2),(-10,-2),(0,8),(50,-1),(50,1),(-50,1),
           (90,0),(140,0),(-140,0),(0,25),(0,-25),(0,50),(0,-50),(-50,-50)]

print("="*104); print("CHECK A -- density normalisation   PASS: 1.05 <= R <= 1.20, spread <= 0.10")
print("  R = int(rho)dA / TRUE objects in patch+bin (full population, not just GEN)")
print("  expect ~1.11 (kNN k/(k-1) bias, not yet applied). R~1.00 => split correction MISSING.")
rows = []
for dlon, lat in CENTERS:
    name = f"prob_maps_grid_dlon{dlon:+04d}_lat{lat:+03d}.npz"
    f = GEN/name
    if not f.exists(): continue
    z = np.load(f, allow_pickle=True)
    x, y = z["x_grid"], z["y_grid"]; dA = float((x[1]-x[0])*(y[1]-y[0]))
    clon = (sun.lon.deg + 180.0 + dlon) % 360.0
    ctr = SkyCoord(lon=clon*u.deg, lat=lat*u.deg, frame=GeocentricTrueEcliptic(obstime=t))
    inpatch = obj.separation(ctr).deg <= MAXSEP
    for pop in ["MBA","TNO","Trojans"]:
        for lo, hi, lab in BINS:
            k = f"density_raw__{pop}__{lab}"
            if k not in z.files: continue
            sel = inpatch & (cache.population.to_numpy()==pop) & \
                  (cache.mag_app.to_numpy()>=lo) & (cache.mag_app.to_numpy()<hi)
            n_true = int(sel.sum())
            if n_true < 1000: continue
            integ = np.nan_to_num(np.asarray(z[k], float)).sum()*dA
            rows.append(dict(center=f"{dlon:+04d}/{lat:+03d}", pop=pop, magbin=lab,
                             n_true=n_true, int_rho=integ, R=integ/n_true))
A = pd.DataFrame(rows)
summ = A.groupby("pop").R.agg(["count","min","median","max"])
summ["spread"] = summ["max"]-summ["min"]
print(summ.to_string(float_format=lambda v: f"{v:,.4f}"))
okA = bool(((A.R>=1.05)&(A.R<=1.20)).all() and (summ["spread"]<=0.10).all())
print(f"  gated bins (n_true>=1000): {len(A)}   in range: {int(((A.R>=1.05)&(A.R<=1.20)).sum())}")
if not okA:
    bad = A[(A.R<1.05)|(A.R>1.20)]
    print(f"  OUT OF RANGE ({len(bad)}):"); print(bad.head(12).to_string(index=False, float_format=lambda v:f"{v:,.4f}"))
    fails.append("A")
print(f"  CHECK A: {'PASS' if okA else 'FAIL'}")

print("\n"+"="*104); print("CHECK B -- support is statistical (integral, not inflated)")
bad_int = bad_scale = 0; ratios = []
for dlon, lat in CENTERS:
    name = f"prob_maps_grid_dlon{dlon:+04d}_lat{lat:+03d}.npz"
    if not (GEN/name).exists(): continue
    zg = np.load(GEN/name, allow_pickle=True)
    zf = np.load(W/"prob_maps_grid_neomod3_full"/name, allow_pickle=True) if (
        W/"prob_maps_grid_neomod3_full"/name).exists() else None
    for pop in ["MBA","TNO","Trojans"]:
        for _,_,lab in BINS:
            k = f"support_count__{pop}__{lab}"
            if k not in zg.files: continue
            sg = np.nan_to_num(np.asarray(zg[k], float))
            if not np.all(sg == np.floor(sg)): bad_int += 1
            if zf is not None and k in zf.files:
                sf = np.nan_to_num(np.asarray(zf[k], float))
                if sf.sum() > 1000:
                    ratios.append(sg.sum()/sf.sum())
r = np.array(ratios)
print(f"  all support values integral: {bad_int == 0}   ({bad_int} violations)")
print(f"  GEN/full support ratio over {len(r)} dense bins: median {np.median(r):.4f} "
      f"[{r.min():.4f}, {r.max():.4f}]   expect ~0.60 (support FALLS with the split)")
okB = (bad_int == 0) and (0.55 <= np.median(r) <= 0.65)
if not okB: fails.append("B")
print(f"  CHECK B: {'PASS' if okB else 'FAIL'}  (support must NOT rise -- that would mean the split "
      f"fraction leaked into the support path)")

print("\n"+"="*104); print("CHECK C -- zero leakage between GEN objects and CAL rows")
gen_ids = set(man.ObjID[man.split=="GEN"]); cal_ids = set(man.ObjID[man.split=="CAL"])
test_ids = set(man.ObjID[man.split=="TEST"])
ov_gc, ov_gt, ov_ct = len(gen_ids&cal_ids), len(gen_ids&test_ids), len(cal_ids&test_ids)
print(f"  GEN n CAL = {ov_gc}   GEN n TEST = {ov_gt}   CAL n TEST = {ov_ct}")
okC = (ov_gc == 0 and ov_gt == 0 and ov_ct == 0)
if not okC: fails.append("C")
print(f"  CHECK C: {'PASS' if okC else 'FAIL'}")

print("\n"+"="*104); print("CHECK D -- coverage / abstention (REPORTED, not gated)")
cov = []
for dlon, lat in CENTERS:
    name = f"prob_maps_grid_dlon{dlon:+04d}_lat{lat:+03d}.npz"
    if not (GEN/name).exists(): continue
    z = np.load(GEN/name, allow_pickle=True)
    for pop in ["NEO","MBA"]:
        for _,_,lab in BINS:
            k = f"support_count__{pop}__{lab}"
            if k not in z.files: continue
            s = np.nan_to_num(np.asarray(z[k], float))
            cov.append(dict(pop=pop, magbin=lab, pix_any=(s>0).sum(), pix_ge2=(s>=2).sum(),
                            frac_ge2=float((s>=2).sum()/max((s>0).sum(),1))))
D = pd.DataFrame(cov).groupby(["pop","magbin"]).agg(
    pix_any=("pix_any","mean"), pix_ge2=("pix_ge2","mean"), frac_ge2=("frac_ge2","mean"))
print(D.to_string(float_format=lambda v: f"{v:,.1f}"))
print("  frac_ge2 = fraction of occupied pixels meeting the smoothing threshold (2 raw clones)")

print("\n"+"="*104)
print(f"E0 A-D: {'ALL PASS' if not fails else 'FAILED: '+','.join(fails)}")
sys.exit(1 if fails else 0)
