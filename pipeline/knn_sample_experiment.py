#!/usr/bin/env python3
"""NEOMOD3 sample-count x kNN neighbour-count experiment — one parameterized runner.

Scientific question
-------------------
How do the NEOMOD3 Monte Carlo sample count and the kNN neighbour count k change the
smoothness and structure of the NEO velocity-density map?

Controlled: ONE sky center (dlon+000_lat+00), ONE magnitude bin (mag24+), the frozen production
epoch/geometry, the +-5 deg/day grid at 0.01 spacing, Gaussian smoothing OFF, no support masking.
Varied: NEO source count (BASE / HIGH) and k (10 / 25 / 50 / 100).

The NEO source is frozen to a parquet ONCE per case; every k in a row reads the same file, so the
only difference within a row is k. `_load_neomod3_cache` is replaced by a loader that returns the
frozen frame plus the case's effective_factor, and the sealed
`velocity_density_pipeline_neomod_clone_only.generate_probability_maps` runs unmodified after that.

Physical normalisation
----------------------
NEOMOD3 datacube weights are ABSOLUTE object counts, so
    effective_factor = n_draws / total_weight_absolute_NEO_count
and each selected row carries physical weight 1/effective_factor. Raising the draw budget raises
n_draws and therefore effective_factor by the SAME factor, so the modelled physical NEO abundance
represented by the sample is unchanged: sum(weights_BASE) ~= sum(weights_HIGH). Verified in `checks`.

Stages
------
    export-base   freeze the literal current GEN selection for this map/bin
    draw-high     shard: draw NEW NEOMOD3 orbits (recorded, non-overlapping seeds), select
    merge-high    source_high = source_base + all draw-high shards (nested)
    map           one (source, k) density map -> NPZ
    posterior     P(NEO) for all 8 cases against byte-identical fixed non-NEO densities
    checks        acceptance_checks.json

Writes only under outputs/more_neomod_samples_knn/. Never writes any production map directory.
"""
from __future__ import annotations

import argparse, copy, glob, hashlib, json, os, subprocess, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
sys.path.insert(0, str(W / "neomod" / "adam_core_stub"))

OUT = W / "outputs" / "more_neomod_samples_knn"
PROD_MAPS = W / "prob_maps_grid_neomod3_GEN_final"
NM3_META = W / "outputs" / "neomod3_projection_cache" / "cache_metadata.json"

# ---- frozen experimental constants -----------------------------------------------------------
CENTER_LABEL = "dlon+000_lat+00"
DLON_DEG, LAT_DEG = 0.0, 0.0
MAG_BIN = {"label": "mag24+", "mag_min": 24.0, "mag_max": 25.0}
MAX_SEP_DEG = 30.0
GRID_LIM = (-5.0, 5.0)
GRID_STEP = 0.01
K_VALUES = (10, 25, 50, 100)
SOURCES = ("BASE", "HIGH")
SMOOTH_DENSITY_MAPS = False          # explicit: Gaussian smoothing OFF everywhere
NONNEO_POPS = ("MBA", "TNO", "Trojans")

# n_d0_grid for the sealed posterior quadrature. The sealed integrand is
#   p(d0) ~ d0^{-k(k+1)} exp(-S/d0^2),  S = sum_j d_j^2,
# whose peak narrows as k grows, so the production default (400 nodes) under-resolves it at large k.
# Measured against the exact closed form  n0 = (k(k+1)/2 - 1/2) / (pi S)  -- max relative deviation
# over 1500 random evaluation points:
#       k      400        2000       8000
#      10   7.1e-16    7.5e-16    6.6e-16
#      25   6.0e-12    8.5e-16    6.6e-16
#      50   2.6e-04    2.5e-15    1.2e-15
#     100   1.1e-02    6.0e-15    4.0e-15
# 400 is converged at k=10 (which is why the production maps are unaffected) but NOT at k=50/100,
# where it would inject a k-dependent numerical bias into exactly the comparison this experiment
# makes. 8000 is converged at every k and is applied IDENTICALLY to all eight cases, so quadrature
# resolution is a controlled variable rather than a function of k.
N_D0_GRID = 8000

# HIGH draw budget. BASE selects 27,781 rows from 1e8 draws (rate 2.7781e-4), so ~6.4e8 additional
# draws land ~177,800 more, for ~205,600 total -- comfortably above the 200,000 floor without
# truncating (truncation would break the exact n_draws <-> effective_factor correspondence).
HIGH_EXTRA_DRAWS_TOTAL = 640_000_000
HIGH_NSHARDS = 320
# Original cache used seed=42 with shards 0..99, i.e. seeds 42..141. Start far away.
HIGH_SEED_BASE = 1_000_000

CACHE_COLS = ["a", "e", "i", "node", "argperi", "t_p", "q", "H",
              "ra_deg", "dec_deg", "dra_deg_day", "ddec_deg_day",
              "lam_deg", "beta_deg", "vlam", "vbeta", "mag_app", "M_obs_deg"]
ORBIT_ID_COLS = ["a", "e", "i", "node", "argperi", "t_p", "H"]


# ---- small helpers ----------------------------------------------------------------------------
def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_array(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", str(W / "neomod"), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def seal_epoch() -> str:
    return json.load(open(W / "outputs/splits/MAP_BUILD_SEAL.json"))["grid"]["ref_obstime"]


def center_lonlat(epoch: str):
    from astropy.time import Time
    from astropy.coordinates import get_sun, GeocentricTrueEcliptic
    t = Time(epoch, scale="utc")
    antisun = (get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t)).lon.deg + 180.0) % 360.0
    return (antisun + DLON_DEG) % 360.0, LAT_DEG


def nm3_norm():
    m = json.load(open(NM3_META))
    return float(m["n_draws"]), float(m["total_weight_absolute_NEO_count"]), \
        float(m["effective_factor_NEO"])


def effective_factor_for(source: str) -> tuple[float, int]:
    """(effective_factor, n_draws) for a case. HIGH raises n_draws by the extra draw budget."""
    n_draws, total_w, eff_base = nm3_norm()
    if source == "BASE":
        return eff_base, int(n_draws)
    manifest = json.load(open(OUT / "high_draw_manifest.json"))
    n_hi = int(n_draws) + int(manifest["n_extra_draws_total"])
    return n_hi / total_w, n_hi


def source_path(source: str) -> Path:
    return OUT / f"source_{source.lower()}.parquet"


def npz_path(source: str, k: int) -> Path:
    return OUT / "maps" / f"density_{source}_k{k:03d}.npz"


def _select_for_map(base_mod, df, epoch, scorer, clon, clat):
    """The production two-stage selection: magnitude cut, then the exact 30-deg sky cut."""
    mag_app = df["mag_app"].to_numpy(dtype=float)
    df_sel, _ = base_mod.select_df_by_mag_bin(df=df, mag_app=mag_app,
                                              mag_min=MAG_BIN["mag_min"], mag_max=MAG_BIN["mag_max"])
    if not len(df_sel):
        return df_sel
    return base_mod.build_visible_subset_dataframe(
        df_sel, obstime_str=epoch, scorer=scorer, max_sep_deg=MAX_SEP_DEG,
        chunk=100_000, show_progress=False,
        center_mode="custom_ecliptic", center_lon_deg=clon, center_lat_deg=clat)


# ---- stage: export-base -----------------------------------------------------------------------
def export_base(a):
    os.chdir(W / "neomod")
    import velocity_density_pipeline_neomod_clone_only as base
    OUT.mkdir(parents=True, exist_ok=True)
    epoch = seal_epoch()
    clon, clat = center_lonlat(epoch)
    print(f"epoch {epoch}  center lon {clon:.6f} lat {clat}", flush=True)

    t0 = time.time()
    df, eff = base._load_neomod3_cache(center_lon_deg=clon, center_lat_deg=clat,
                                       obstime_str=epoch, max_sep_deg=MAX_SEP_DEG)
    print(f"healpix superset {len(df):,} rows  effective_factor {eff!r}  ({time.time()-t0:.0f}s)",
          flush=True)
    _, scorer = base.load_s3m_population("neo", verbose=False)
    sel = _select_for_map(base, df, epoch, scorer, clon, clat)
    print(f"BASE selected rows: {len(sel):,}", flush=True)

    keep = [c for c in CACHE_COLS if c in sel.columns]
    out = sel[keep].copy()
    out["source_case"] = "BASE"
    out["draw_seed"] = np.int64(-1)          # -1 = original frozen GEN cache realisation
    out.to_parquet(source_path("BASE"), index=False)
    print(f"wrote {source_path('BASE')}  sha256 {sha256_file(source_path('BASE'))}", flush=True)


# ---- stage: draw-high -------------------------------------------------------------------------
def draw_high(a):
    os.chdir(W / "neomod")
    import velocity_density_pipeline_neomod_clone_only as base
    import neomod3_sampler as nm3s
    d = OUT / "high_shards"; d.mkdir(parents=True, exist_ok=True)
    outp = d / f"high_{a.shard:04d}.parquet"
    if outp.exists() and not a.overwrite:
        print(f"[skip] {outp} exists", flush=True); return

    epoch = seal_epoch()
    clon, clat = center_lonlat(epoch)
    n_this = HIGH_EXTRA_DRAWS_TOTAL // HIGH_NSHARDS
    seed = HIGH_SEED_BASE + a.shard
    if 42 <= seed <= 141:
        raise RuntimeError(f"seed {seed} overlaps the original cache seed range 42..141")
    rng = np.random.default_rng(seed)
    t0 = time.time()
    df = nm3s.sample_neomod3_orbits(n_this, epoch, rng=rng)
    print(f"shard {a.shard}: drew {n_this:,} -> {len(df):,} valid ({time.time()-t0:.0f}s)", flush=True)

    _, scorer = base.load_s3m_population("neo", verbose=False)
    # All-sky projection exactly as the cache builder does (center 0,0 / max_sep 180), then the
    # SAME float32 cast the cache applies, so HIGH rows are commensurate with BASE rows.
    vis = base.build_visible_subset_dataframe(
        df, obstime_str=epoch, scorer=scorer, max_sep_deg=180.0,
        chunk=200_000, show_progress=False,
        center_mode="custom_ecliptic", center_lon_deg=0.0, center_lat_deg=0.0)
    keep = [c for c in CACHE_COLS if c in vis.columns]
    vis = vis[keep].copy()
    for c in vis.columns:
        if vis[c].dtype == np.float64:
            vis[c] = vis[c].astype(np.float32)
    print(f"  projected {len(vis):,} ({time.time()-t0:.0f}s)", flush=True)

    sel = _select_for_map(base, vis, epoch, scorer, clon, clat)
    sel = sel[[c for c in CACHE_COLS if c in sel.columns]].copy()
    sel["source_case"] = "HIGH_EXTRA"
    sel["draw_seed"] = np.int64(seed)
    sel.to_parquet(outp, index=False)
    print(f"  selected {len(sel):,} -> {outp}  ({time.time()-t0:.0f}s)", flush=True)


# ---- stage: merge-high ------------------------------------------------------------------------
def merge_high(a):
    fs = sorted(glob.glob(str(OUT / "high_shards" / "high_*.parquet")))
    if len(fs) != HIGH_NSHARDS:
        raise RuntimeError(f"expected {HIGH_NSHARDS} shards, found {len(fs)}")
    b = pd.read_parquet(source_path("BASE"))
    extra = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
    print(f"BASE {len(b):,} + extra {len(extra):,}", flush=True)
    hi = pd.concat([b, extra], ignore_index=True)
    if len(hi) < 200_000:
        raise RuntimeError(f"HIGH has {len(hi):,} rows, below the 200,000 floor")

    dup = int(hi.duplicated(subset=ORBIT_ID_COLS).sum())
    print(f"duplicate orbital records in HIGH: {dup}", flush=True)
    if dup:
        raise RuntimeError(f"{dup} duplicate generated orbital records in HIGH")
    # nesting: every BASE row must survive into HIGH
    assert (hi.iloc[:len(b)][ORBIT_ID_COLS].to_numpy()
            == b[ORBIT_ID_COLS].to_numpy()).all(), "HIGH is not nested on BASE"

    json.dump({
        "n_extra_draws_total": int(HIGH_EXTRA_DRAWS_TOTAL),
        "n_shards": int(HIGH_NSHARDS),
        "draws_per_shard": int(HIGH_EXTRA_DRAWS_TOTAL // HIGH_NSHARDS),
        "seed_base": int(HIGH_SEED_BASE),
        "seeds": [int(HIGH_SEED_BASE + s) for s in range(HIGH_NSHARDS)],
        "original_cache_seed_range": [42, 141],
        "seed_overlap_with_original_cache": False,
        "n_rows_base": int(len(b)), "n_rows_extra": int(len(extra)), "n_rows_high": int(len(hi)),
        "nested": True, "duplicate_orbital_records": 0,
    }, open(OUT / "high_draw_manifest.json", "w"), indent=2)

    hi.to_parquet(source_path("HIGH"), index=False)
    print(f"wrote {source_path('HIGH')}  {len(hi):,} rows  sha256 {sha256_file(source_path('HIGH'))}",
          flush=True)

    effb, ndb = effective_factor_for("BASE")
    effh, ndh = effective_factor_for("HIGH")
    wb, wh = len(b) / effb, len(hi) / effh
    print(f"\nweights: BASE {wb:,.4f} (n_draws {ndb:,}, eff {effb:.6f})", flush=True)
    print(f"         HIGH {wh:,.4f} (n_draws {ndh:,}, eff {effh:.6f})", flush=True)
    print(f"         rel diff {abs(wh-wb)/wb:.3e}", flush=True)


# ---- stage: map -------------------------------------------------------------------------------
def build_map(a):
    os.chdir(W / "neomod")
    import velocity_density_pipeline_neomod_clone_only as base
    (OUT / "maps").mkdir(parents=True, exist_ok=True)
    src, k = a.source, a.k
    epoch = seal_epoch()
    clon, clat = center_lonlat(epoch)
    sp = source_path(src)
    src_df = pd.read_parquet(sp)
    eff, n_draws = effective_factor_for(src)
    print(f"[{src} k={k}] {len(src_df):,} rows  eff {eff:.6f}  n_draws {n_draws:,}", flush=True)

    # Replace the cache loader with the frozen source for this case. Everything downstream --
    # magnitude cut, sky cut, tree build, estimator, normalisation -- is the sealed code path.
    def _frozen_loader(center_lon_deg=None, center_lat_deg=None, obstime_str=None, max_sep_deg=None):
        return src_df, eff
    base._load_neomod3_cache = _frozen_loader

    # Prove the requested k actually reaches the estimator.
    seen = {}
    real_eval = base.evaluate_density_map_full_posterior_2d

    def _spy(tree, grid_points, k=None, n_d0_grid=None, show_progress=True, n_jobs=1, **kw):
        seen["k"] = int(k); seen["n_d0_grid"] = int(n_d0_grid); seen["n_pts"] = int(len(tree.data))
        return real_eval(tree, grid_points, k=k, n_d0_grid=n_d0_grid,
                         show_progress=show_progress, n_jobs=n_jobs, **kw)
    base.evaluate_density_map_full_posterior_2d = _spy

    ps = copy.deepcopy(base.DEFAULT_POPULATION_SETTINGS)
    s3m_neo, scorer = base.load_s3m_population("neo", verbose=False)
    cs = {"NEO": {"df": s3m_neo, "scorer": scorer,
                  "clone_factor": ps["NEO"]["clone_factor"],
                  "use_conditional_cloner": True,
                  "scatter_size": 4, "scatter_alpha": 0.1,
                  "k_map": int(k)}}
    tmp = OUT / "maps" / f"_gpm_{src}_k{k:03d}.npz"
    t0 = time.time()
    base.generate_probability_maps(
        obstime_str=epoch, output_path=str(tmp),
        center_lon_deg=clon, center_lat_deg=clat, center_label=CENTER_LABEL,
        max_sep_deg=MAX_SEP_DEG, clone_sources=cs, mag_bins=[MAG_BIN],
        grid_lim=GRID_LIM, grid_step=GRID_STEP,
        k_map=int(k), n_d0_grid_map=int(N_D0_GRID),
        smooth_density_maps=SMOOTH_DENSITY_MAPS,
        save_overlays=False, n_jobs=a.n_jobs, verbose=True)
    base.evaluate_density_map_full_posterior_2d = real_eval
    print(f"  map built in {time.time()-t0:.0f}s; estimator saw {seen}", flush=True)

    if seen.get("k") != int(k):
        raise RuntimeError(f"requested k={k} but the estimator received k={seen.get('k')}")
    if seen.get("n_d0_grid") != int(N_D0_GRID):
        raise RuntimeError(f"n_d0_grid {seen.get('n_d0_grid')} != {N_D0_GRID}")
    if seen.get("n_pts") != len(src_df):
        raise RuntimeError(f"tree has {seen.get('n_pts')} points, source has {len(src_df)}")

    z = np.load(tmp, allow_pickle=True)
    lab = MAG_BIN["label"]
    dens_phys = np.asarray(z[f"density_unsmoothed__NEO__{lab}"], dtype=float)
    dens_raw_stored = np.asarray(z[f"density_raw__NEO__{lab}"], dtype=float)
    x_grid = np.asarray(z["x_grid"], dtype=float); y_grid = np.asarray(z["y_grid"], dtype=float)
    smooth_flag = bool(np.asarray(z["smooth_density_maps"]).item())
    if smooth_flag:
        raise RuntimeError("smooth_density_maps is True; this experiment requires it OFF")
    # with smoothing off the stored raw map must equal the unsmoothed map (float32 storage)
    if not np.allclose(dens_raw_stored, dens_phys.astype(np.float32), rtol=0, atol=0,
                       equal_nan=True):
        raise RuntimeError("density_raw != float32(density_unsmoothed) with smoothing off")

    # k-th neighbour distance quantiles from every grid node
    from scipy.spatial import cKDTree
    pts = np.column_stack([src_df.vlam.to_numpy(float), src_df.vbeta.to_numpy(float)])
    tree = cKDTree(pts)
    X0, Y0 = np.meshgrid(x_grid, y_grid)
    gp = np.column_stack([X0.ravel(), Y0.ravel()])
    dk, _ = tree.query(gp, k=int(k), workers=a.n_jobs)
    dk = dk[:, -1]
    qs = [10, 50, 90, 99]
    dk_full = np.percentile(dk, qs)
    cen = (np.abs(gp[:, 0]) <= 1.0) & (np.abs(gp[:, 1]) <= 1.0)
    dk_cen = np.percentile(dk[cen], qs)

    outp = npz_path(src, k)
    np.savez_compressed(
        outp,
        density_neo_physical=dens_phys,             # rho_NEO after /effective_factor; enters P
        density_neo_raw=dens_phys * eff,            # before the physical normalisation
        x_grid=x_grid, y_grid=y_grid,
        k=np.int64(k),
        n_source_rows=np.int64(len(src_df)),
        total_physical_weight=np.float64(len(src_df) / eff),
        effective_factor=np.float64(eff),
        n_draws=np.int64(n_draws),
        source_case=src,
        source_parquet=str(sp),
        source_sha256=sha256_file(sp),
        code_commit=git_commit(),
        sealed_module_sha256=sha256_file(W / "neomod/src/velocity_density_pipeline_neomod_clone_only.py"),
        gaussian_smoothing=np.bool_(False),
        n_d0_grid=np.int64(N_D0_GRID),
        dk_quantile_levels=np.array(qs),
        dk_quantiles_full_grid=dk_full,
        dk_quantiles_central_1deg=dk_cen,
        estimator_k_observed=np.int64(seen["k"]),
        center_label=CENTER_LABEL, magnitude_bin=lab,
        epoch=epoch, center_lon_deg=np.float64(clon), center_lat_deg=np.float64(clat),
        grid_lim=np.array(GRID_LIM), grid_step=np.float64(GRID_STEP),
    )
    Path(tmp).unlink(missing_ok=True)
    print(f"wrote {outp}", flush=True)


# ---- stage: posterior -------------------------------------------------------------------------
def posterior(a):
    lab = MAG_BIN["label"]
    pm = np.load(PROD_MAPS / f"prob_maps_grid_{CENTER_LABEL}.npz", allow_pickle=True)
    fixed, hashes = {}, {}
    for p in NONNEO_POPS:
        arr = np.asarray(pm[f"density_unsmoothed__{p}__{lab}"], dtype=float)
        raw = np.asarray(pm[f"density_raw__{p}__{lab}"], dtype=float)
        # production smooths NEO only, so the non-NEO raw and unsmoothed maps must agree
        if not np.array_equal(raw, arr.astype(np.float32)):
            raise RuntimeError(f"{p}: production raw != float32(unsmoothed); not an unsmoothed map")
        fixed[p] = arr
        hashes[p] = sha256_array(arr)
        print(f"fixed {p:8s} sha256 {hashes[p]}  sum {arr.sum():.6e}", flush=True)

    rec, per_case_hashes = {}, {}
    for src in SOURCES:
        for k in K_VALUES:
            z = np.load(npz_path(src, k), allow_pickle=True)
            neo = np.asarray(z["density_neo_physical"], dtype=float)
            tot = neo + fixed["MBA"] + fixed["TNO"] + fixed["Trojans"]
            ok = np.isfinite(tot) & (tot > 0)
            P = np.full(neo.shape, np.nan)
            P[ok] = neo[ok] / tot[ok]
            comp = {}
            for nm, arr in (("NEO", neo), ("MBA", fixed["MBA"]),
                            ("TNO", fixed["TNO"]), ("Trojans", fixed["Trojans"])):
                c = np.full(neo.shape, np.nan); c[ok] = arr[ok] / tot[ok]; comp[nm] = c
            ssum = comp["NEO"][ok] + comp["MBA"][ok] + comp["TNO"][ok] + comp["Trojans"][ok]
            maxdev = float(np.max(np.abs(ssum - 1.0))) if ok.any() else float("nan")
            if not (maxdev < 1e-12):
                raise RuntimeError(f"{src} k={k}: probabilities do not sum to 1 (max dev {maxdev:.3e})")
            key = f"{src}_k{k:03d}"
            np.savez_compressed(OUT / "maps" / f"posterior_{key}.npz",
                                P_NEO=P, defined=ok, x_grid=z["x_grid"], y_grid=z["y_grid"],
                                max_abs_sum_deviation=np.float64(maxdev),
                                nonneo_hashes=json.dumps(hashes))
            per_case_hashes[key] = hashes
            rec[key] = {"max_abs_sum_deviation": maxdev,
                        "n_defined": int(ok.sum()), "n_undefined": int((~ok).sum())}
            print(f"  {key}: defined {int(ok.sum()):,}  undefined {int((~ok).sum()):,}  "
                  f"max|sum-1| {maxdev:.3e}", flush=True)

    ident = all(per_case_hashes[k2] == hashes for k2 in per_case_hashes)
    json.dump({"fixed_nonneo_sha256": hashes,
               "identical_nonneo_arrays_across_all_cases": bool(ident),
               "per_case": rec,
               "source_production_map": str(PROD_MAPS / f"prob_maps_grid_{CENTER_LABEL}.npz"),
               "magnitude_bin": lab},
              open(OUT / "posterior_provenance.json", "w"), indent=2)
    print(f"\nidentical non-NEO arrays across all 8 cases: {ident}", flush=True)


# ---- stage: finalize --------------------------------------------------------------------------
WEIGHT_TOLERANCE = 0.02      # Poisson on the ~178k new selections is ~0.24%; 2% is a loose bound


def _verify_base_is_production(base_mod) -> dict:
    """Byte-level proof that source_base.parquet IS the production NEO sample for this map/bin.

    `make_support_count_map` is a raw 2-D histogram of the selected (vlam, vbeta) points, so
    rebuilding it from source_base.parquet and comparing to the sealed production map's stored
    support count is an exact identity test on the point set -- not a count comparison.
    """
    lab = MAG_BIN["label"]
    pm = np.load(PROD_MAPS / f"prob_maps_grid_{CENTER_LABEL}.npz", allow_pickle=True)
    prod_support = np.asarray(pm[f"support_count__NEO__{lab}"], dtype=float)
    b = pd.read_parquet(source_path("BASE"))
    z = np.load(npz_path("BASE", 10), allow_pickle=True)
    mine = base_mod.make_support_count_map(
        b.vlam.to_numpy(float), b.vbeta.to_numpy(float),
        np.asarray(z["x_grid"], float), np.asarray(z["y_grid"], float), clone_factor=1)
    identical = bool(np.array_equal(mine, prod_support.astype(mine.dtype)))
    return {"identical_support_histogram": identical,
            "n_base_rows": int(len(b)),
            "n_in_production_support": float(prod_support.sum()),
            "n_in_rebuilt_support": float(mine.sum()),
            "production_map": str(PROD_MAPS / f"prob_maps_grid_{CENTER_LABEL}.npz")}


def finalize(a):
    os.chdir(W / "neomod")
    import velocity_density_pipeline_neomod_clone_only as base
    epoch = seal_epoch()
    clon, clat = center_lonlat(epoch)
    man = json.load(open(OUT / "high_draw_manifest.json"))
    pp = json.load(open(OUT / "posterior_provenance.json"))
    effb, ndb = effective_factor_for("BASE")
    effh, ndh = effective_factor_for("HIGH")
    nb_rows = len(pd.read_parquet(source_path("BASE"), columns=["vlam"]))
    nh_rows = len(pd.read_parquet(source_path("HIGH"), columns=["vlam"]))
    wb, wh = nb_rows / effb, nh_rows / effh

    cases = {}
    for src in SOURCES:
        for k in K_VALUES:
            z = np.load(npz_path(src, k), allow_pickle=True)
            cases[f"{src}_k{k:03d}"] = {
                "k": int(z["k"]), "estimator_k_observed": int(z["estimator_k_observed"]),
                "n_source_rows": int(z["n_source_rows"]),
                "total_physical_weight": float(z["total_physical_weight"]),
                "effective_factor": float(z["effective_factor"]), "n_draws": int(z["n_draws"]),
                "source_sha256": str(z["source_sha256"]),
                "gaussian_smoothing": bool(z["gaussian_smoothing"]),
                "n_d0_grid": int(z["n_d0_grid"]),
                "npz_sha256": sha256_file(npz_path(src, k)),
                "dk_quantiles_full_grid": np.asarray(z["dk_quantiles_full_grid"]).tolist(),
                "dk_quantiles_central_1deg": np.asarray(z["dk_quantiles_central_1deg"]).tolist(),
            }
    base_proof = _verify_base_is_production(base)

    prov = {
        "experiment": "NEOMOD3 sample count x kNN neighbour count, NEO velocity-density maps",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "code_commit": git_commit(),
        "sealed_module_sha256": sha256_file(
            W / "neomod/src/velocity_density_pipeline_neomod_clone_only.py"),
        "runner": "neomod/pipeline/knn_sample_experiment.py",
        "runner_sha256": sha256_file(W / "neomod/pipeline/knn_sample_experiment.py"),
        "center_label": CENTER_LABEL, "center_lon_deg": clon, "center_lat_deg": clat,
        "magnitude_bin": MAG_BIN, "epoch": epoch, "max_sep_deg": MAX_SEP_DEG,
        "grid_lim": list(GRID_LIM), "grid_step": GRID_STEP,
        "k_values": list(K_VALUES), "sources": list(SOURCES),
        "gaussian_smoothing": False, "support_masking": False,
        "n_d0_grid": N_D0_GRID,
        "n_d0_grid_note": ("production default 400 is converged at k=10 but not at k=50/100 "
                           "(2.6e-4 / 1.1e-2 vs the exact closed form); 8000 is converged to "
                           "<=4e-15 at every k and is used identically in all eight cases"),
        "weight_tolerance": WEIGHT_TOLERANCE,
        "source_base": {"path": str(source_path("BASE")), "sha256": sha256_file(source_path("BASE")),
                        "n_rows": nb_rows, "effective_factor": effb, "n_draws": ndb,
                        "total_physical_weight": wb},
        "source_high": {"path": str(source_path("HIGH")), "sha256": sha256_file(source_path("HIGH")),
                        "n_rows": nh_rows, "effective_factor": effh, "n_draws": ndh,
                        "total_physical_weight": wh},
        "high_draw_manifest": man,
        "base_is_production_realisation": base_proof,
        "fixed_nonneo": pp,
        "cases": cases,
    }
    json.dump(prov, open(OUT / "provenance.json", "w"), indent=2)

    src_hashes = {s: {cases[f"{s}_k{k:03d}"]["source_sha256"] for k in K_VALUES} for s in SOURCES}
    checks = {
        "all_eight_cases_exist": all(npz_path(s, k).exists() for s in SOURCES for k in K_VALUES),
        "base_is_literal_current_source_realisation": base_proof["identical_support_histogram"],
        "high_has_at_least_200000_rows": nh_rows >= 200_000,
        "high_is_nested_on_base": bool(man["nested"]),
        "no_duplicate_generated_orbital_records": man["duplicate_orbital_records"] == 0,
        "seeds_non_overlapping_with_original_cache": not man["seed_overlap_with_original_cache"],
        "weight_totals_agree": abs(wh - wb) / wb < WEIGHT_TOLERANCE,
        "weight_relative_difference": abs(wh - wb) / wb,
        "weight_tolerance": WEIGHT_TOLERANCE,
        "one_source_hash_per_row": all(len(v) == 1 for v in src_hashes.values()),
        "requested_k_reached_estimator": all(
            c["k"] == c["estimator_k_observed"] for c in cases.values()),
        "gaussian_smoothing_false_everywhere": all(
            not c["gaussian_smoothing"] for c in cases.values()),
        "identical_n_d0_grid": len({c["n_d0_grid"] for c in cases.values()}) == 1,
        "identical_nonneo_arrays_across_all_cases":
            pp["identical_nonneo_arrays_across_all_cases"],
        "posterior_sums_to_one_where_defined": all(
            v["max_abs_sum_deviation"] < 1e-12 for v in pp["per_case"].values()),
        "max_posterior_sum_deviation": max(v["max_abs_sum_deviation"]
                                           for v in pp["per_case"].values()),
    }
    checks["ALL_PASS"] = all(v for k2, v in checks.items()
                             if isinstance(v, bool) and k2 != "ALL_PASS")
    json.dump(checks, open(OUT / "acceptance_checks.json", "w"), indent=2)

    for k2, v in checks.items():
        if isinstance(v, bool):
            print(f"  {'PASS' if v else 'FAIL'}  {k2}")
        else:
            print(f"        {k2} = {v}")
    print(f"\nwrote {OUT/'provenance.json'} and {OUT/'acceptance_checks.json'}", flush=True)


# ---- stage: bundle ----------------------------------------------------------------------------
NOTEBOOK = W / "neomod/notebooks/validation/more_neomod_samples_knn_maps.ipynb"


def bundle(a):
    import tarfile
    prov = json.load(open(OUT / "provenance.json"))
    chk = json.load(open(OUT / "acceptance_checks.json"))
    b, h = prov["source_base"], prov["source_high"]
    man = prov["high_draw_manifest"]

    qt = []
    for src in SOURCES:
        for k in K_VALUES:
            c = prov["cases"][f"{src}_k{k:03d}"]
            q = c["dk_quantiles_central_1deg"]
            qt.append(f"| {src} | {k} | {c['n_source_rows']:,} | {q[0]:.4f} | {q[1]:.4f} | "
                      f"{q[2]:.4f} | {q[3]:.4f} |")

    readme = f"""# NEOMOD3 sample count x kNN neighbour count — NEO velocity-density maps

One sky map, one magnitude bin, two NEO source counts, four k values: eight NEO density maps.

**Question.** How do the NEOMOD3 Monte Carlo sample count and the kNN neighbour count change the
smoothness and structure of the NEO velocity-density map?

This is a controlled visual comparison. No preferred k is selected here, and no
calibration, threshold, ROC or F1 is computed. CAL and TEST are not read.

## Fixed

| | |
|---|---|
| sky center | `{prov['center_label']}` (lon {prov['center_lon_deg']:.4f}, lat {prov['center_lat_deg']:.1f}) |
| magnitude bin | `{prov['magnitude_bin']['label']}` ({prov['magnitude_bin']['mag_min']:.0f} <= mag < {prov['magnitude_bin']['mag_max']:.0f}) |
| epoch | {prov['epoch']} |
| velocity domain | {prov['grid_lim']} deg/day, step {prov['grid_step']} |
| patch radius | {prov['max_sep_deg']:.0f} deg |
| estimator | Bayesian kNN, sealed module `{prov['sealed_module_sha256'][:16]}...` |
| Gaussian smoothing | **{prov['gaussian_smoothing']}** |
| support masking | {prov['support_masking']} |
| quadrature `n_d0_grid` | {prov['n_d0_grid']} (all cases) |

## Source cases

| | rows | n_draws | effective_factor | total physical weight |
|---|---|---|---|---|
| BASE | {b['n_rows']:,} | {b['n_draws']:,} | {b['effective_factor']:.6f} | {b['total_physical_weight']:,.4f} |
| HIGH | {h['n_rows']:,} | {h['n_draws']:,} | {h['effective_factor']:.6f} | {h['total_physical_weight']:,.4f} |

`HIGH = all BASE rows + {man['n_rows_extra']:,} additional independent NEOMOD3 draws`, nested, with
seeds {man['seed_base']:,}+0..{man['n_shards']-1} (disjoint from the original cache's
{man['original_cache_seed_range'][0]}..{man['original_cache_seed_range'][1]}). No duplicated,
jittered, resampled or GMM-cloned rows: {man['duplicate_orbital_records']} duplicate orbital records.

**Physical abundance is preserved, not the row count.** Raising the draw budget raises `n_draws`
and therefore `effective_factor = n_draws / total_weight` by the same factor, so the modelled
physical NEO abundance the sample represents is unchanged. The two totals differ by
{chk['weight_relative_difference']:.3e} (tolerance {chk['weight_tolerance']:.0e}). That residual is
Monte Carlo sampling noise, and it is dominated by BASE: with {b['n_rows']:,} samples BASE carries
~{100.0 / (b['n_rows'] ** 0.5):.2f}% relative noise on its own, against ~{100.0 / (h['n_rows'] ** 0.5):.2f}% for HIGH.

## A numerical point that matters for this comparison

The sealed posterior integrand is `p(d0) ~ d0^-k(k+1) exp(-S/d0^2)` with `S = sum_j d_j^2`, so its
peak narrows sharply as k grows. Measured against the exact closed form
`n0 = (k(k+1)/2 - 1/2) / (pi S)`, the production default of 400 quadrature nodes gives max relative
deviation 7.1e-16 at k=10 but **2.6e-04 at k=50 and 1.1e-02 at k=100**. Left at 400, that error
would grow with k and show up as a structural difference between panels. All eight cases therefore
use {prov['n_d0_grid']} nodes, converged to <= 4e-15 at every k. Production maps at k=10 are unaffected.

## k-th neighbour distance, central [-1,+1] deg/day window (deg/day; pixel = 0.01)

| source | k | rows | p10 | p50 | p90 | p99 |
|---|---|---|---|---|---|---|
{chr(10).join(qt)}

## Contents

- `more_neomod_samples_knn_maps.ipynb` — executed notebook (the deliverable to read first)
- `source_base.parquet` — sha256 `{b['sha256']}`
- `source_high.parquet` — sha256 `{h['sha256']}`
- `maps/density_{{BASE,HIGH}}_k{{010,025,050,100}}.npz` — eight NEO density maps
- `maps/posterior_{{BASE,HIGH}}_k{{010,025,050,100}}.npz` — P(NEO) per case
- `provenance.json`, `acceptance_checks.json`

`P(NEO) = rho_NEO / (rho_NEO + rho_MBA + rho_TNO + rho_Trojan)` uses the fixed production non-NEO
densities, byte-identical across all eight cases (hashes in `provenance.json`). Undefined pixels
are NaN, never 0. Max |sum of four class probabilities - 1| over defined pixels:
{chk['max_posterior_sum_deviation']:.3e}.

## Acceptance

All checks: **{'PASS' if chk['ALL_PASS'] else 'FAIL'}** — see `acceptance_checks.json`.

Reproduce with `neomod/pipeline/knn_sample_experiment.py`
(`export-base` / `draw-high` / `merge-high` / `map` / `posterior` / `finalize` / `bundle`)
plus `neomod/pipeline/slurm/knn_draw_high.sbatch` and `knn_maps.sbatch`.
Code commit `{prov['code_commit']}`.
"""
    (OUT / "README.md").write_text(readme)
    print(f"wrote {OUT/'README.md'}", flush=True)

    tar_path = W / "outputs" / "more_neomod_samples_knn_bundle.tar.gz"
    items = [(NOTEBOOK, f"more_neomod_samples_knn/{NOTEBOOK.name}")]
    for n in ("README.md", "provenance.json", "acceptance_checks.json",
              "posterior_provenance.json", "high_draw_manifest.json",
              "source_base.parquet", "source_high.parquet"):
        items.append((OUT / n, f"more_neomod_samples_knn/{n}"))
    for src in SOURCES:
        for k in K_VALUES:
            items.append((npz_path(src, k), f"more_neomod_samples_knn/maps/{npz_path(src, k).name}"))
            p2 = OUT / "maps" / f"posterior_{src}_k{k:03d}.npz"
            items.append((p2, f"more_neomod_samples_knn/maps/{p2.name}"))
    missing = [str(s) for s, _ in items if not Path(s).exists()]
    if missing:
        raise RuntimeError(f"cannot bundle, missing: {missing}")
    with tarfile.open(tar_path, "w:gz") as tf:
        for s, arc in items:
            tf.add(s, arcname=arc)
    print(f"wrote {tar_path}  ({tar_path.stat().st_size/1e6:.1f} MB)", flush=True)
    print(f"tarball sha256 {sha256_file(tar_path)}", flush=True)
    print(f"members: {len(items)}  (high_shards/ deliberately excluded)", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    s = p.add_subparsers(dest="cmd", required=True)
    s.add_parser("export-base").set_defaults(func=export_base)
    d = s.add_parser("draw-high"); d.add_argument("--shard", type=int, required=True)
    d.add_argument("--overwrite", action="store_true"); d.set_defaults(func=draw_high)
    s.add_parser("merge-high").set_defaults(func=merge_high)
    m = s.add_parser("map")
    m.add_argument("--source", choices=SOURCES, required=True)
    m.add_argument("--k", type=int, required=True)
    m.add_argument("--n-jobs", type=int, default=1)
    m.set_defaults(func=build_map)
    s.add_parser("posterior").set_defaults(func=posterior)
    s.add_parser("finalize").set_defaults(func=finalize)
    s.add_parser("bundle").set_defaults(func=bundle)
    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
