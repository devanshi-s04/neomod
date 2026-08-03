#!/usr/bin/env python3
"""A1.2 / A1.3 proof on a real generated mini-map. A1.2 is NOT complete until this passes.

1. density_unsmoothed is captured BEFORE smoothing (differs from density_raw where smoothing acted,
   and equals it where it did not)
2. NPZ round-trip is EXACT (float64, bit-identical)
3. thresholds 2/3/5/10 produce bit-identical density_unsmoothed AND support_count arrays
   (only the post-smoothing product may differ)
4. each archive records its own requested threshold, and it matches
5. neo_provenance is persisted and reports source_role NEOMOD3_GEN
"""
import copy, json, os, sys, tempfile
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
man = pd.read_parquet(W/"outputs/splits/nonneo_split_manifest.parquet")
gen_ids = set(man.ObjID[man.split == "GEN"])
cache = pd.read_parquet(W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet")
cache = cache[cache.ObjID.isin(gen_ids) | (cache.population == "NEO")].reset_index(drop=True)
_, scorer = v.load_s3m_population("neo", verbose=False)
t = Time(EPOCH, scale="utc")
sun = get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
clon = (sun.lon.deg + 180.0) % 360.0

def build(thr):
    cs = {}
    for pop, cfg in copy.deepcopy(v.DEFAULT_POPULATION_SETTINGS).items():
        sub = cache[cache.population == pop].reset_index(drop=True)
        cs[pop] = dict(df=sub, scorer=scorer, clone_factor=cfg["clone_factor"],
                       use_conditional_cloner=True, scatter_size=cfg["scatter_size"],
                       scatter_alpha=cfg["scatter_alpha"],
                       _mag_app=sub["mag_app"].to_numpy(float))
    x, y, X0, Y0, gp = v.make_default_grid(grid_lim=(-1.0, 1.0), grid_step=0.05)  # small + fast
    return v.build_cloned_maps_for_center_magbin(
        center_lon_deg=clon, center_lat_deg=0.0, center_label=f"T{thr}",
        clone_sources=cs, obstime_str=EPOCH, max_sep_deg=30.0,
        rng=np.random.default_rng(1), overlay_max=1000,
        x_grid=x, y_grid=y, grid_points=gp, X0=X0,
        k_map=v.DEFAULT_K_MAP, n_d0_grid_map=100,
        mag_min=23.0, mag_max=24.0,
        smooth_density_maps=True, smooth_population_names=("NEO",),
        smooth_support_threshold=float(thr),
        smooth_sigma_pixels=v.DEFAULT_SMOOTH_SIGMA_PIXELS,
        smooth_truncate_sigma=v.DEFAULT_SMOOTH_TRUNCATE_SIGMA,
        smooth_support_scale_by_clone_factor=False, smooth_presmoothing_passes=None,
        n_jobs=1)

fails = []
res = {thr: build(thr) for thr in (2, 3, 5, 10)}
print("=== 1. density_unsmoothed captured BEFORE smoothing ===")
r = res[2]
for pop in ["NEO", "MBA"]:
    u = np.asarray(r["density_maps_unsmoothed"][pop], float)
    s = np.asarray(r["density_maps_downweighted_raw"][pop], float)
    same = np.array_equal(u, s)
    exp = "differ (NEO is smoothed)" if pop == "NEO" else "equal (MBA not smoothed)"
    ok = (not same) if pop == "NEO" else same
    print(f"  {pop:>4}: unsmoothed==raw? {same}   expected: {exp}   {'OK' if ok else 'FAIL'}")
    if not ok: fails.append(f"1/{pop}")

print("\n=== 3. thresholds share bit-identical unsmoothed density + support ===")
for pop in ["NEO", "MBA", "TNO", "Trojans"]:
    us = [np.asarray(res[k]["density_maps_unsmoothed"][pop], float) for k in (2, 3, 5, 10)]
    sp = [np.asarray(res[k]["support_count_maps"][pop], float) for k in (2, 3, 5, 10)]
    ok_u = all(np.array_equal(us[0], a) for a in us[1:])
    ok_s = all(np.array_equal(sp[0], a) for a in sp[1:])
    sm = [np.asarray(res[k]["density_maps_downweighted_raw"][pop], float) for k in (2, 3, 5, 10)]
    differs = not all(np.array_equal(sm[0], a) for a in sm[1:])
    print(f"  {pop:>8}: unsmoothed identical={ok_u}  support identical={ok_s}  "
          f"smoothed differs across thresholds={differs}")
    if not (ok_u and ok_s): fails.append(f"3/{pop}")

print("\n=== 2 + 4 + 5. NPZ round-trip, threshold metadata, provenance ===")
with tempfile.TemporaryDirectory() as td:
    p = Path(td)/"mini.npz"
    v.save_maps_to_npz(str(p), results={"mag23": res[5]},
        mag_bins=[{"label": "mag23", "mag_min": 23.0, "mag_max": 24.0}],
        x_grid=np.linspace(-1, 1, 41), y_grid=np.linspace(-1, 1, 41),
        obstime_str=EPOCH, center_lon_deg=clon, center_lat_deg=0.0,
        center_label="mini", max_sep_deg=30.0,
        population_names=["NEO", "MBA", "TNO", "Trojans"],
        smoothing={"enabled": True, "population_names": ("NEO",),
                   "support_scale_by_clone_factor": False, "support_threshold": 5.0,
                   "sigma_pixels": v.DEFAULT_SMOOTH_SIGMA_PIXELS,
                   "truncate_sigma": v.DEFAULT_SMOOTH_TRUNCATE_SIGMA, "passes": []},
        neo_provenance=res[5]["neo_provenance"], save_overlays=False)
    z = np.load(p, allow_pickle=True)
    for pop in ["NEO", "MBA"]:
        k = f"density_unsmoothed__{pop}__mag23"
        ok = k in z.files and np.array_equal(
            np.asarray(z[k], float), np.asarray(res[5]["density_maps_unsmoothed"][pop], float))
        print(f"  round-trip exact {pop}: {ok}   dtype={z[k].dtype if k in z.files else 'MISSING'}")
        if not ok: fails.append(f"2/{pop}")
    thr = float(z["smooth_support_threshold"])
    print(f"  recorded threshold = {thr}  (requested 5.0)  {'OK' if thr == 5.0 else 'FAIL'}")
    if thr != 5.0: fails.append("4")
    pj = json.loads(str(z["neo_provenance_json"])) if "neo_provenance_json" in z.files else {}
    print(f"  neo_provenance.source_role = {pj.get('source_role')!r}  seed={pj.get('seed')}")
    print(f"    n_s3m_neo_rows_discarded = {pj.get('n_s3m_neo_rows_discarded'):,}")
    if pj.get("source_role") != "NEOMOD3_GEN": fails.append("5")

print(f"\n{'='*60}\n{'ALL PASS' if not fails else 'FAILED: ' + ','.join(fails)}")
sys.exit(1 if fails else 0)
