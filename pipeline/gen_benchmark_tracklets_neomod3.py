#!/usr/bin/env python3
"""NEOMOD3-based benchmark tracklets -- the test set whose populations MATCH the maps (docs §11.1).

Why this exists: the S3M benchmark scores NEOMOD3-prior maps against S3M truth, and S3M contains no
H>25 NEOs at all, so the two describe different universes and no calibration/contamination result on
it can be trusted (§9.10).

Composition mirrors the maps exactly:
    NEO                -> fresh INDEPENDENT NEOMOD3 draw (H 15-28), propagated to the map epoch
    MBA/TNO/Trojans    -> Stage-0 n-body epoch cache (the same objects the maps are built from)

Independence matters for NEO: maps are evaluated on a fixed velocity grid, so a clone that IS in the
map raises the density of its own cell. Re-using cache clones as test objects would bias the test
optimistically for the one population under examination. A different seed removes that entirely.

Differences from gen_benchmark_tracklets_s3m.py, all deliberate:
  * magnitude cut 14 <= mag_app < 25 -- the maps' own span. v3 had none, so only 32.7% of it was
    scoreable (NEO just 9%).
  * population mix from ABSOLUTE expected counts (NEOMOD3's own normalisation for NEO, real S3M
    counts for the rest) rather than S3M catalogue proportions.
  * non-NEO taken from the Stage-0 cache instead of being re-propagated.

Schema is identical to v3 so every downstream scorer works unchanged.

    sbatch neomod/pipeline/slurm/benchmark_neomod3.sbatch          # NEO shards
    python neomod/pipeline/gen_benchmark_tracklets_neomod3.py build
"""
from __future__ import annotations
import argparse, json, os, sys, time, glob
from pathlib import Path
import numpy as np, pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"adam_core_stub")); sys.path.insert(0, str(W/"neomod"/"src"))
EPOCH = "2027-08-25T00:00:00"
DT_DAYS = 30.0/1440.0
MAG_MIN, MAG_MAX = 14.0, 25.0
LON_STEP, SUN_EXCLUSION = 10.0, 40.0
DLON_LIMIT = 180.0 - SUN_EXCLUSION
LAT_BASE = [0, 1, 2, 3, 4, 5, 8, 12, 18, 25, 35, 50]
# None => NO downscaling: the benchmark carries the TRUE absolute object counts for this epoch,
# so each sky direction holds the number really there. Set an int only if size forces it -- the
# factor is then applied uniformly to ALL populations (never per-population; that is the v1 bug).
TARGET_TOTAL = None
NEO_SEED = 20270825           # deliberately NOT the cache's seed 42
SPLIT_ROLE = os.environ.get("BM_SPLIT_ROLE", "TEST_UNSEALED")
OUT = Path(os.environ.get("BM_OUT_DIR", str(W/"outputs"/"benchmark_tracklets_neomod3")))
# EXPLICIT NEO source. Defaulting to OUT/"neo_shards" silently coupled the NEO realisation to the
# output directory: pointing BM_OUT_DIR at the CAL directory made the builder look for shards that
# were never there (it exited rather than mis-sourcing, but the coupling is the bug).
SHARDS = Path(os.environ.get("BM_NEO_SHARDS", str(OUT/"neo_shards")))
NEO_SEED_EXPECTED = int(os.environ.get("BM_NEO_SEED_EXPECTED", "0")) or None
EPOCH_CACHE = W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet"
META = json.load(open(W/"outputs/neomod3_projection_cache/cache_metadata.json"))
TOTAL_NEO_ABS = META["total_weight_absolute_NEO_count"]      # 11,432,918 NEOs, H 15-28


def build_grid():
    dl = [round(d, 6) for d in np.arange(-180.0, 180.0, LON_STEP) if abs(d) <= DLON_LIMIT + 1e-9]
    lats = sorted({float(v) for v in LAT_BASE} | {float(-v) for v in LAT_BASE})
    return np.array([(d, l) for l in lats for d in dl], dtype=float)


def map_filename(d, l): return f"prob_maps_grid_dlon{int(round(d)):+04d}_lat{int(round(l)):+03d}.npz"
def map_label(d, l):    return f"grid_dlon{int(round(d)):+04d}_lat{int(round(l)):+03d}"


def neo_shard(shard: int, nshards: int, n_orbits_total: int):
    """Sample + project one shard of an INDEPENDENT NEOMOD3 draw."""
    os.chdir(W/"neomod")
    import velocity_density_pipeline_gmm as vdp
    import neomod3_sampler as nm3s
    SHARDS.mkdir(parents=True, exist_ok=True)
    out = SHARDS/f"neo_shard_{shard:03d}.parquet"
    if out.exists():
        print(f"[skip] {out}"); return
    n = n_orbits_total//nshards
    rng = np.random.default_rng(NEO_SEED + shard)
    t0 = time.time()
    df = nm3s.sample_neomod3_orbits(n, EPOCH, rng=rng)
    print(f"shard {shard}: drew {n:,} -> {len(df):,} valid NEO orbits ({time.time()-t0:.0f}s)", flush=True)
    _, scorer = vdp.load_s3m_population("neo", verbose=False)
    vis = vdp.build_visible_subset_dataframe(df, obstime_str=EPOCH, scorer=scorer,
                                             max_sep_deg=180.0, chunk=200_000, show_progress=False)
    vis = vis[(vis.mag_app >= MAG_MIN) & (vis.mag_app < MAG_MAX)].reset_index(drop=True)
    vis["n_orbits_drawn"] = n
    vis.to_parquet(out, index=False)
    print(f"shard {shard}: {len(vis):,} rows in {MAG_MIN}<=mag<{MAG_MAX}  ({time.time()-t0:.0f}s)", flush=True)


def _tracklet_frame(sub, pop, lam, beta, dlon, grid_arr, mjd_ref, objid_prefix):
    d = dlon[:, None] - grid_arr[None, :, 0]
    b = beta[:, None] - grid_arr[None, :, 1]
    ci = np.argmin(d**2 + b**2, axis=1)
    ad, al = grid_arr[ci, 0], grid_arr[ci, 1]
    ra0 = sub.ra_deg.to_numpy(float); dec0 = sub.dec_deg.to_numpy(float)
    dra = sub.dra_deg_day.to_numpy(float); ddec = sub.ddec_deg_day.to_numpy(float)
    mag = sub.mag_app.to_numpy(float)
    o = pd.DataFrame()
    o["ObjID"] = [f"{objid_prefix}{i:08d}" for i in range(len(sub))]
    o["population"] = pop
    o["prob_map_file"] = [map_filename(x, y) for x, y in zip(ad, al)]
    o["prob_map"] = [map_label(x, y) for x, y in zip(ad, al)]
    o["n_det_per_night"] = 2
    o["mean_ra"], o["mean_dec"] = ra0, dec0
    o["mean_dra"], o["mean_ddec"], o["mean_mag"] = dra, ddec, mag
    o["ra0"], o["dec0"], o["mjd0_utc"], o["mag0"] = ra0, dec0, mjd_ref, mag
    o["ra1"] = (ra0 + dra*DT_DAYS) % 360.0
    o["dec1"] = np.clip(dec0 + ddec*DT_DAYS, -89.99, 89.99)
    o["mjd1_utc"], o["mag1"] = mjd_ref + DT_DAYS, mag
    o["lam_deg"], o["beta_deg"], o["dlon_from_antisun_deg"] = lam, beta, dlon
    for c in ["H", "vlam", "vbeta", "a_au", "e", "q_au"]:
        if c in sub.columns: o[c] = sub[c].to_numpy(float)
    return o


def build(_args):
    from astropy.time import Time
    from astropy.coordinates import GeocentricTrueEcliptic, SkyCoord, get_sun
    from astropy.utils import iers; iers.conf.auto_max_age = None
    import astropy.units as u
    grid_arr = build_grid()
    t = Time(EPOCH, scale="utc"); mjd_ref = t.mjd
    antisun = (get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t)).lon.deg + 180.0) % 360.0
    print(f"epoch {EPOCH}  antisun_lon={antisun:.4f}  grid cells={len(grid_arr)}")

    files = sorted(glob.glob(str(SHARDS/"neo_shard_*.parquet")))
    if not files:
        sys.exit(f"no NEO shards in {SHARDS} -- set BM_NEO_SHARDS explicitly")
    neo = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    n_drawn = int(pd.read_parquet(files[0], columns=["n_orbits_drawn"]).n_orbits_drawn.iloc[0])*len(files)
    w_new = TOTAL_NEO_ABS/n_drawn
    import hashlib as _hl
    _shard_sha = _hl.sha256(b"".join(
        _hl.sha256(Path(f).read_bytes()).digest() for f in files)).hexdigest()
    print(f"NEO: {len(neo):,} clones from {n_drawn:,} draws -> w_new = {w_new:.5f} objects/clone")
    print(f"NEO source dir   : {SHARDS}")
    print(f"NEO shards       : {len(files)}  combined sha256 = {_shard_sha}")
    print(f"NEO split role   : {SPLIT_ROLE}")

    cache = pd.read_parquet(EPOCH_CACHE, columns=["ObjID", "population", "ra_deg", "dec_deg", "dra_deg_day",
                                                  "ddec_deg_day", "mag_app", "lam_deg", "beta_deg",
                                                  "vlam", "vbeta", "H"])
    cache = cache[(cache.mag_app >= MAG_MIN) & (cache.mag_app < MAG_MAX)]
    print(f"epoch cache in mag range: {len(cache):,}")

    # --- role fraction, derived from the frozen split manifest counts ---
    _cnt = json.loads((W/"outputs/splits/split_provenance.json").read_text())["counts"]
    if SPLIT_ROLE in ("GEN", "CAL", "TEST"):
        _num = sum(v[SPLIT_ROLE] for v in _cnt.values())
        _den = sum(sum(v.values()) for v in _cnt.values())
        role_fraction = _num/_den
        print(f"role fraction ({SPLIT_ROLE}) derived from split manifest: "
              f"{_num:,}/{_den:,} = {role_fraction:.6f}")
    else:
        role_fraction = 1.0
    pops, expected = {}, {}
    for pop in ["NEO", "MBA", "TNO", "Trojans"]:
        if pop == "NEO":
            sub = neo.copy()
            c = SkyCoord(ra=sub.ra_deg.to_numpy()*u.deg, dec=sub.dec_deg.to_numpy()*u.deg, frame="icrs")
            e = c.transform_to(GeocentricTrueEcliptic(obstime=t))
            lam, beta = e.lon.deg, e.lat.deg
        else:
            sub = cache[cache.population == pop].reset_index(drop=True)
            if SPLIT_ROLE in ("GEN", "CAL", "TEST"):
                import pandas as _pd
                _man = _pd.read_parquet(W/"outputs/splits/nonneo_split_manifest.parquet")
                _keep = set(_man.ObjID[_man.split == SPLIT_ROLE])
                _n0 = len(sub)
                sub = sub[sub.ObjID.isin(_keep)].reset_index(drop=True)
                print(f"    [{SPLIT_ROLE}] {pop}: {len(sub):,} of {_n0:,} objects", flush=True)
            lam, beta = sub.lam_deg.to_numpy(float), sub.beta_deg.to_numpy(float)
        dlon = ((lam - antisun + 180.0) % 360.0) - 180.0
        keep = np.abs(dlon) <= DLON_LIMIT + 0.5
        sub, lam, beta, dlon = sub[keep].reset_index(drop=True), lam[keep], beta[keep], dlon[keep]
        pops[pop] = (sub, lam, beta, dlon)
        # NEO is an INDEPENDENT NEOMOD3 draw representing the FULL sky, while non-NEO is a
        # partition holding only fraction f of the objects. Without scaling NEO by the same f the
        # class ratio is wrong by 1/f (CAL v1: NEO 3.76% vs the physical 0.776%). The fraction is
        # DERIVED from split_provenance.json, never hardcoded.
        expected[pop] = len(sub)*w_new*role_fraction if pop == "NEO" else float(len(sub))
        print(f"  {pop:>8}: {len(sub):,} in-grid  -> expected TRUE objects {expected[pop]:,.0f}")

    tot = sum(expected.values())
    scale = 1.0 if TARGET_TOTAL is None else min(1.0, TARGET_TOTAL/tot)
    print(f"\ntotal expected true objects in-grid, {MAG_MIN}<=mag<{MAG_MAX}: {tot:,.0f}")
    print(f"scaling by {scale:.6f}" + ("  (none -- TRUE absolute counts)" if scale == 1.0
          else f" -> ~{TARGET_TOTAL:,} rows"))
    print(f"NEO share of the benchmark: {100*expected['NEO']/tot:.2f}%  "
          f"(S3M-proportional v3 was 1.76%)")

    rng = np.random.default_rng(7)
    frames = []
    OUT.mkdir(parents=True, exist_ok=True)
    for i, (pop, (sub, lam, beta, dlon)) in enumerate(pops.items()):
        n_rows = int(round(expected[pop]*scale))
        if n_rows > len(sub):
            # Clamping HERE would apply a per-population cap while the others keep the common
            # factor -- precisely the v1 benchmark bug (MBA capped 200k, TNO uncapped -> MBA
            # suppressed ~69x, TNO share inflated). Refuse rather than silently distort the mix.
            sys.exit(f"FATAL {pop}: need {n_rows:,} rows, only {len(sub):,} available.\n"
                     f"  Clamping would cap ONE population and break the ratios (the v1 bug).\n"
                     f"  Fix by drawing more NEOMOD3 orbits (--n-orbits-total) or lowering "
                     f"TARGET_TOTAL, then rerun.")
        idx = rng.choice(len(sub), size=n_rows, replace=False)
        f = _tracklet_frame(sub.iloc[idx].reset_index(drop=True), pop, lam[idx], beta[idx], dlon[idx],
                            grid_arr, mjd_ref, f"NM{i}")
        f.to_parquet(OUT/f"tracklets_{pop}.parquet", index=False)
        frames.append(f); print(f"  wrote {pop}: {len(f):,} rows")
    comb = pd.concat(frames, ignore_index=True)
    comb.to_parquet(OUT/"tracklets_benchmark_neomod3.parquet", index=False)
    print(f"\nTOTAL {len(comb):,} rows -> {OUT/'tracklets_benchmark_neomod3.parquet'}")
    print((comb.population.value_counts(normalize=True)*100).round(2).to_string())

    # --- ratio integrity: ONE common factor, so every share must survive the scaling ---
    print("\n=== ratio check (v1 bug = per-population caps distorting these) ===")
    chk = pd.DataFrame({
        "true_expected": pd.Series(expected),
        "true_share_%": pd.Series({k: 100*v/tot for k, v in expected.items()}),
        "benchmark_rows": comb.population.value_counts(),
    })
    chk["bench_share_%"] = 100*chk.benchmark_rows/len(comb)
    chk["share_error_pp"] = chk["bench_share_%"] - chk["true_share_%"]
    print(chk.to_string(float_format=lambda v: f"{v:,.3f}"))
    worst = chk.share_error_pp.abs().max()
    print(f"  worst share error: {worst:.4f} percentage points "
          f"({'OK' if worst < 0.05 else 'INVESTIGATE'})")

    # --- spatial realism: does each direction bin hold the right number for that sky position? ---
    print("\n=== per-direction realism (benchmark rows / scale vs true in-grid objects) ===")
    bins = [(0,20),(20,40),(40,70),(70,110),(110,141)]
    rows_dir = []
    for lo, hi in bins:
        m = (comb.dlon_from_antisun_deg.abs() >= lo) & (comb.dlon_from_antisun_deg.abs() < hi)
        true_here = 0.0
        for pop, (sub, lam, beta, dlon) in pops.items():
            k = (np.abs(dlon) >= lo) & (np.abs(dlon) < hi)
            true_here += k.sum()*(w_new if pop == "NEO" else 1.0)
        rows_dir.append(dict(direction=f"{lo}-{hi} deg", bench_rows=int(m.sum()),
                             implied_true=int(m.sum()/scale), true_objects=int(true_here),
                             neo_frac_pct=100*float(comb[m].population.eq("NEO").mean())))
    dd = pd.DataFrame(rows_dir)
    dd["ratio"] = dd.implied_true/dd.true_objects
    print(dd.to_string(index=False, float_format=lambda v: f"{v:,.3f}"))
    print("  ratio ~1 => each sky direction carries a realistic object count for this date.")
    json.dump(dict(epoch=EPOCH, split_role=SPLIT_ROLE, role_fraction=role_fraction,
                   neo_shards_dir=str(SHARDS), mag_min=MAG_MIN, mag_max=MAG_MAX, neo_seed=NEO_SEED,
                   n_neo_orbits_drawn=n_drawn, w_abs_new=w_new, expected_true=expected,
                   scale=scale, n_rows=len(comb)), open(OUT/"benchmark_metadata.json", "w"), indent=2)


if __name__ == "__main__":
    p = argparse.ArgumentParser(); s = p.add_subparsers(dest="cmd", required=True)
    a = s.add_parser("neo-shard"); a.add_argument("--shard", type=int, required=True)
    a.add_argument("--nshards", type=int, default=10)
    a.add_argument("--n-orbits-total", type=int, default=10_000_000)
    a.set_defaults(f=lambda x: neo_shard(x.shard, x.nshards, x.n_orbits_total))
    b = s.add_parser("build"); b.set_defaults(f=build)
    args = p.parse_args(); args.f(args)
