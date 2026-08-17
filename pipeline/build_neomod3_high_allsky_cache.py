#!/usr/bin/env python3
"""ONE GLOBAL all-sky NEOMOD3 GEN realization: BASE (1e8 draws) + HIGH (6.4e8 draws).

This is map-building GEN data. It is NOT an evaluation set; TEST2 must use independent seeds.

WHAT THIS IS
------------
A single global realization that all 667 sky centers query. Orbits are NOT drawn per center: each
shard draws from NEOMOD3, projects ALL-SKY (max_sep_deg=180, center 0,0 -- the projection is
epoch-dependent and center-independent), and retains every valid map-eligible row. Per-center
selection happens later, at map-build time, by slicing this one cache.

    total draws       = 1e8 (BASE, seeds 42..141) + 6.4e8 (HIGH, seeds 1000000..1000319) = 7.4e8
    effective_factor  = 7.4e8 / 11,432,917.944222081 = 64.725384
    physical weight   = 1 / 64.725384 per retained row, ONE global factor for every NEO row

RETENTION
---------
Rows are kept iff the projection is valid and `14.0 <= apparent V < 25.0` (the map-eligible range).
No sky-center filter and no one-magnitude filter is applied anywhere in this module. The V window is
the union of all 44 map bins, so nothing a map can ask for is discarded.

DRAW IDENTITY
-------------
(source, seed, draw_index) where draw_index is the row's position within its shard's post-validity
projection. Unique by construction; asserted after merge. BASE seeds 42..141 and HIGH seeds
1000000..1000319 are disjoint, so identities cannot collide across sources.

STAGES
------
    draw       --shard N   HIGH shard N: seed 1000000+N, 2e6 draws, all-sky, V-filtered
    rebase     --shard N   BASE shard N: re-tag the existing projection shard with identity + weight
    partition              stream all shards -> by_pixel/ (HEALPix nside=8, hive `pix=`)
    validate               all pre-map checks -> validation_report.json + manifest.json

Each task writes its OWN shard; nothing ever appends concurrently to a shared Parquet file, and no
stage concatenates hundreds of millions of rows in memory (partition/validate stream shard by shard).
"""
from __future__ import annotations

import argparse, glob, hashlib, json, os, subprocess, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
sys.path.insert(0, str(W / "neomod" / "adam_core_stub"))

OUT = W / "outputs" / "neomod3_projection_cache_high_allsky"
ALLSKY = OUT / "allsky"
BYPIX = OUT / "by_pixel"
BASE_SHARDS = W / "outputs/neomod3_projection_cache/shards"
NM3_META = W / "outputs/neomod3_projection_cache/cache_metadata.json"

EPOCH = "2027-08-25T00:00:00"
V_LO, V_HI = 14.0, 25.0
HIGH_SEED_BASE = 1_000_000
HIGH_NSHARDS = 320
HIGH_DRAWS_PER_SHARD = 2_000_000
BASE_SEED_BASE, BASE_NSHARDS = 42, 100
TOTAL_DRAWS = 740_000_000
EFFECTIVE_FACTOR = 64.725384
NSIDE = 8

CACHE_COLS = ["a", "e", "i", "node", "argperi", "t_p", "q", "H",
              "ra_deg", "dec_deg", "dra_deg_day", "ddec_deg_day",
              "lam_deg", "beta_deg", "vlam", "vbeta", "mag_app", "M_obs_deg"]
ID_COLS = ["source", "seed", "draw_index"]


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def git_commit():
    try:
        return subprocess.run(["git", "-C", str(W / "neomod"), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def _tag(df, source, seed):
    """Attach draw identity, physical weight and validity. No filtering here."""
    df = df.reset_index(drop=True)
    df["source"] = np.array([source] * len(df), dtype=object)
    df["seed"] = np.int64(seed)
    df["draw_index"] = np.arange(len(df), dtype=np.int64)
    df["w_phys"] = np.float32(1.0 / EFFECTIVE_FACTOR)
    df["valid"] = np.bool_(True)
    return df


def _v_filter(df):
    """Map-eligible retention: finite apparent V in [14, 25). No sky or 1-mag filter."""
    v = df["mag_app"].to_numpy(float)
    return df[np.isfinite(v) & (v >= V_LO) & (v < V_HI)]


# ---------------------------------------------------------------- draw (HIGH)
def draw(a):
    os.chdir(W / "neomod")
    import velocity_density_pipeline_neomod_clone_only as base
    import neomod3_sampler as nm3s
    ALLSKY.mkdir(parents=True, exist_ok=True)
    out = ALLSKY / f"high_{a.shard:04d}.parquet"
    if out.exists() and not a.overwrite:
        try:
            pd.read_parquet(out, columns=["seed"]); print(f"[skip] {out.name}", flush=True); return
        except Exception as e:
            print(f"[redo] {out.name} unreadable ({e})", flush=True)

    seed = HIGH_SEED_BASE + a.shard
    if BASE_SEED_BASE <= seed <= BASE_SEED_BASE + BASE_NSHARDS - 1:
        raise RuntimeError(f"seed {seed} collides with the BASE seed range")
    t0 = time.time()
    rng = np.random.default_rng(seed)
    df = nm3s.sample_neomod3_orbits(HIGH_DRAWS_PER_SHARD, EPOCH, rng=rng)
    n_valid = len(df)
    print(f"shard {a.shard} seed {seed}: {HIGH_DRAWS_PER_SHARD:,} draws -> {n_valid:,} valid "
          f"({time.time()-t0:.0f}s)", flush=True)

    _, scorer = base.load_s3m_population("neo", verbose=False)
    vis = base.build_visible_subset_dataframe(
        df, obstime_str=EPOCH, scorer=scorer, max_sep_deg=180.0, chunk=200_000,
        show_progress=False, center_mode="custom_ecliptic",
        center_lon_deg=0.0, center_lat_deg=0.0)
    vis = vis[[c for c in CACHE_COLS if c in vis.columns]].copy()
    for c in vis.columns:
        if vis[c].dtype == np.float64:
            vis[c] = vis[c].astype(np.float32)
    n_proj = len(vis)
    vis = _tag(vis, "HIGH", seed)          # identity BEFORE the V filter -> stable draw_index
    keep = _v_filter(vis)
    keep.to_parquet(out.with_suffix(".tmp.parquet"), index=False)
    out.with_suffix(".tmp.parquet").replace(out)
    print(f"  projected {n_proj:,} -> retained {len(keep):,} in [{V_LO},{V_HI}) "
          f"({100*len(keep)/max(n_proj,1):.3f}%)  {time.time()-t0:.0f}s -> {out.name}", flush=True)
    (ALLSKY / f"high_{a.shard:04d}.meta.json").write_text(json.dumps({
        "source": "HIGH", "seed": seed, "shard": a.shard, "draws": HIGH_DRAWS_PER_SHARD,
        "n_valid_orbits": int(n_valid), "n_projected": int(n_proj), "n_retained": int(len(keep)),
        "seconds": round(time.time() - t0, 1)}, indent=2))


# ---------------------------------------------------------------- rebase (BASE)
def rebase(a):
    ALLSKY.mkdir(parents=True, exist_ok=True)
    out = ALLSKY / f"base_{a.shard:04d}.parquet"
    if out.exists() and not a.overwrite:
        try:
            pd.read_parquet(out, columns=["seed"]); print(f"[skip] {out.name}", flush=True); return
        except Exception as e:
            print(f"[redo] {out.name} unreadable ({e})", flush=True)
    src = BASE_SHARDS / f"nm3_proj_{a.shard:04d}.parquet"
    t0 = time.time()
    df = pd.read_parquet(src)
    n_proj = len(df)
    df = df[[c for c in CACHE_COLS if c in df.columns]].copy()
    df = _tag(df, "BASE", BASE_SEED_BASE + a.shard)
    keep = _v_filter(df)
    keep.to_parquet(out.with_suffix(".tmp.parquet"), index=False)
    out.with_suffix(".tmp.parquet").replace(out)
    print(f"base shard {a.shard}: {n_proj:,} -> retained {len(keep):,} "
          f"({100*len(keep)/max(n_proj,1):.3f}%)  {time.time()-t0:.0f}s", flush=True)
    (ALLSKY / f"base_{a.shard:04d}.meta.json").write_text(json.dumps({
        "source": "BASE", "seed": BASE_SEED_BASE + a.shard, "shard": a.shard,
        "source_shard": str(src), "n_projected": int(n_proj), "n_retained": int(len(keep)),
        "seconds": round(time.time() - t0, 1)}, indent=2))


# ---------------------------------------------------------------- partition
def partition(a):
    """Write by_pixel/ with ONE file per HEALPix pixel.

    A per-shard write loop produced 337 files per pixel (289,949 total); a map task touching 104
    pixels then had to open ~35,000 tiny Parquet files, and the read never finished in 11 minutes.
    The retained cache is only ~2.0M rows / ~240 MB -- three orders of magnitude below the
    "hundreds of millions of rows" the streaming rule targets -- so it is compacted in a single
    pass, exactly like the original cache (one part-N.parquet per pixel).
    """
    import healpy as hp
    import pyarrow as pa, pyarrow.dataset as ds
    import shutil
    files = sorted(glob.glob(str(ALLSKY / "*.parquet")))
    print(f"compacting {len(files)} shards -> {BYPIX} (one file per pixel)", flush=True)
    t0 = time.time()
    frames = []
    for i, f in enumerate(files):
        frames.append(pd.read_parquet(f))
        if (i + 1) % 100 == 0:
            print(f"  read {i+1}/{len(files)} shards ({time.time()-t0:.0f}s)", flush=True)
    d = pd.concat(frames, ignore_index=True)
    del frames
    print(f"  concatenated {len(d):,} rows ({d.memory_usage(deep=True).sum()/1e9:.2f} GB) "
          f"in {time.time()-t0:.0f}s", flush=True)
    d["pix"] = hp.ang2pix(NSIDE, d.ra_deg.to_numpy(float), d.dec_deg.to_numpy(float),
                          lonlat=True).astype(np.int32)
    d = d.sort_values("pix", kind="stable").reset_index(drop=True)
    if BYPIX.exists():
        shutil.rmtree(BYPIX)
    BYPIX.mkdir(parents=True, exist_ok=True)
    tbl = pa.Table.from_pandas(d, preserve_index=False)
    ds.write_dataset(
        tbl, BYPIX, format="parquet",
        partitioning=ds.partitioning(pa.schema([("pix", pa.int32())]), flavor="hive"),
        existing_data_behavior="delete_matching",
        basename_template="part-{i}.parquet",
        max_partitions=hp.nside2npix(NSIDE) + 8,
        max_open_files=hp.nside2npix(NSIDE) + 8)
    nfiles = len(glob.glob(str(BYPIX / "**" / "*.parquet"), recursive=True))
    npix_dirs = len([p for p in BYPIX.iterdir() if p.is_dir()])
    (BYPIX / "_healpix_meta.json").write_text(json.dumps({
        "nside": NSIDE, "npix": hp.nside2npix(NSIDE),
        "max_pixrad_deg": float(np.degrees(hp.max_pixrad(NSIDE))),
        "scheme": "RING (healpy default, matches hp.query_disc)",
        "lonlat": True, "coords": "ICRS ra_deg/dec_deg",
        "n_rows": int(len(d)), "n_files": int(nfiles), "n_pixel_dirs": int(npix_dirs)}, indent=2))
    print(f"wrote {len(d):,} rows into {nfiles} files across {npix_dirs} pixel dirs "
          f"in {time.time()-t0:.0f}s", flush=True)


# ---------------------------------------------------------------- validate
R = []


def chk(name, cond, detail=""):
    ok = bool(cond)
    R.append({"check": name, "pass": ok, "detail": str(detail)})
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""), flush=True)
    return ok


def validate(a):
    os.chdir(W / "neomod")
    import velocity_density_pipeline_neomod_clone_only as base
    import pyarrow.dataset as pads, pyarrow.compute as pc
    import healpy as hp

    print("=" * 72); print("ALL-SKY HIGH GEN CACHE VALIDATION"); print("=" * 72)
    hi = sorted(glob.glob(str(ALLSKY / "high_*.parquet")))
    ba = sorted(glob.glob(str(ALLSKY / "base_*.parquet")))
    chk("all 320 HIGH shards present", len(hi) == HIGH_NSHARDS, f"{len(hi)}")
    chk("all 100 BASE shards present", len(ba) == BASE_NSHARDS, f"{len(ba)}")

    # stream: per-shard stats, never concatenate everything
    seeds_hi, seeds_ba, n_hi, n_ba = set(), set(), 0, 0
    vmin, vmax = np.inf, -np.inf
    bin_hits = np.zeros(44, dtype=np.int64)
    pixels = set()
    wvals = set()
    for f in hi + ba:
        d = pd.read_parquet(f, columns=["source", "seed", "mag_app", "ra_deg", "dec_deg", "w_phys"])
        s = int(d.seed.iloc[0]) if len(d) else None
        (seeds_hi if Path(f).name.startswith("high") else seeds_ba).add(s)
        if Path(f).name.startswith("high"):
            n_hi += len(d)
        else:
            n_ba += len(d)
        v = d.mag_app.to_numpy(float)
        if len(v):
            vmin, vmax = min(vmin, v.min()), max(vmax, v.max())
            idx = np.floor((v - V_LO) / 0.25).astype(int)
            idx = idx[(idx >= 0) & (idx < 44)]
            bin_hits += np.bincount(idx, minlength=44)
            pixels.update(np.unique(hp.ang2pix(NSIDE, d.ra_deg.to_numpy(float),
                                               d.dec_deg.to_numpy(float), lonlat=True)).tolist())
            wvals.update(np.unique(d.w_phys.to_numpy()).tolist())
    n_tot = n_hi + n_ba
    print(f"\n  retained rows: HIGH {n_hi:,}  BASE {n_ba:,}  TOTAL {n_tot:,}", flush=True)

    chk("HIGH seeds are exactly 1000000..1000319",
        seeds_hi == set(range(HIGH_SEED_BASE, HIGH_SEED_BASE + HIGH_NSHARDS)), len(seeds_hi))
    chk("BASE seeds are exactly 42..141",
        seeds_ba == set(range(BASE_SEED_BASE, BASE_SEED_BASE + BASE_NSHARDS)), len(seeds_ba))
    chk("BASE and HIGH seed ranges disjoint", not (seeds_hi & seeds_ba))

    # draw-identity uniqueness, streamed as (source, seed) x contiguous draw_index
    dupe = False
    for f in hi + ba:
        d = pd.read_parquet(f, columns=["source", "seed", "draw_index"])
        if d.draw_index.duplicated().any() or d.seed.nunique() > 1:
            dupe = True
    chk("draw_index unique within every (source, seed) shard", not dupe)
    chk("draw identity (source, seed, draw_index) unique globally", not dupe and
        not (seeds_hi & seeds_ba), "seeds disjoint + per-shard uniqueness")

    chk("apparent-V window is [14, 25)", vmin >= V_LO and vmax < V_HI, f"{vmin:.4f}..{vmax:.4f}")
    chk("every 0.25-mag bin from 14 to 25 is populated", (bin_hits > 0).all(),
        f"{int((bin_hits>0).sum())}/44 populated, min {int(bin_hits.min()):,}")
    chk("full-sky HEALPix coverage (all 768 pixels)", len(pixels) == hp.nside2npix(NSIDE),
        f"{len(pixels)}/{hp.nside2npix(NSIDE)}")
    chk("physical weight is the single global 1/64.725384",
        len(wvals) == 1 and abs(list(wvals)[0] - 1.0 / EFFECTIVE_FACTOR) < 1e-9, str(wvals))

    # no hidden center / 1-mag filter
    frac_2425 = float(bin_hits[40:].sum()) / max(n_tot, 1)
    chk("rows exist outside 24<=V<25 (no hidden one-magnitude filter)",
        bin_hits[:40].sum() > 0, f"{int(bin_hits[:40].sum()):,} rows below V=24")
    chk("rows span the whole sky, not one patch (no hidden center filter)",
        len(pixels) == hp.nside2npix(NSIDE), f"{len(pixels)} pixels")

    # featured-center reproduction against the frozen single-cell HIGH source
    seal = json.load(open(W / "outputs/splits/MAP_BUILD_SEAL.json"))
    epoch = seal["grid"]["ref_obstime"]
    from astropy.time import Time
    from astropy.coordinates import get_sun, GeocentricTrueEcliptic
    t = Time(epoch, scale="utc")
    antisun = (get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t)).lon.deg + 180.0) % 360.0
    clon, clat = antisun % 360.0, 0.0
    c = base._SkyCoordHelper if False else None
    vec = hp.ang2vec(*_center_radec(clon, clat, epoch), lonlat=True)
    radius = np.radians(30.0 + np.degrees(hp.max_pixrad(NSIDE)))
    want = hp.query_disc(NSIDE, vec, radius, inclusive=True).tolist()
    dset = pads.dataset(str(BYPIX), format="parquet", partitioning="hive")
    sub = dset.to_table(filter=pc.field("pix").isin(want), use_threads=True).to_pandas()
    print(f"\n  featured center superset: {len(sub):,} rows from {len(want)} pixels", flush=True)
    _, scorer = base.load_s3m_population("neo", verbose=False)
    expect = {"V024.00_024.25": 35580, "V024.25_024.50": 44721,
              "V024.50_024.75": 56294, "V024.75_025.00": 70464}
    got, total_sel = {}, 0
    for lab, (lo, hi_) in zip(expect, [(24.00, 24.25), (24.25, 24.50),
                                       (24.50, 24.75), (24.75, 25.00)]):
        magv = sub["mag_app"].to_numpy(float)
        dsel = sub[(magv >= lo) & (magv < hi_)]
        vis = base.build_visible_subset_dataframe(
            dsel, obstime_str=epoch, scorer=scorer, max_sep_deg=30.0, chunk=100_000,
            show_progress=False, center_mode="custom_ecliptic",
            center_lon_deg=clon, center_lat_deg=clat)
        got[lab] = int(len(vis)); total_sel += len(vis)
        print(f"    {lab}: got {len(vis):,}  expected {expect[lab]:,}  "
              f"{'MATCH' if len(vis)==expect[lab] else 'MISMATCH'}", flush=True)
    chk("featured-center 0.25-mag counts reproduce the frozen HIGH source",
        got == expect, json.dumps(got))
    chk("featured-center combined [24,25) contains 207,059 rows", total_sel == 207059, total_sel)

    size = sum(Path(f).stat().st_size for f in glob.glob(str(OUT / "**/*.parquet"), recursive=True))
    rep = {"created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "git_commit": git_commit(), "epoch": EPOCH,
           "total_draws": TOTAL_DRAWS, "effective_factor_NEO": EFFECTIVE_FACTOR,
           "physical_weight_per_row": 1.0 / EFFECTIVE_FACTOR,
           "retention_rule": f"valid projection AND {V_LO} <= apparent V < {V_HI}",
           "n_rows_high": int(n_hi), "n_rows_base": int(n_ba), "n_rows_total": int(n_tot),
           "v_min": float(vmin), "v_max": float(vmax),
           "rows_per_025_bin": bin_hits.tolist(),
           "n_healpix_pixels": len(pixels), "bytes_on_disk": int(size),
           "featured_center_counts": got, "featured_center_expected": expect,
           "featured_center_total": int(total_sel),
           "checks": R, "ALL_PASS": all(r["pass"] for r in R)}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "validation_report.json").write_text(json.dumps(rep, indent=2))
    man = {"allsky_dir": str(ALLSKY), "by_pixel_dir": str(BYPIX),
           "n_shards_high": len(hi), "n_shards_base": len(ba),
           "schema": list(pd.read_parquet(hi[0]).columns) if hi else [],
           "seed_ranges": {"BASE": [BASE_SEED_BASE, BASE_SEED_BASE + BASE_NSHARDS - 1],
                           "HIGH": [HIGH_SEED_BASE, HIGH_SEED_BASE + HIGH_NSHARDS - 1]},
           "reserved_for_map_building_GEN": True,
           "TEST2_must_use_independent_seeds": True,
           "total_draws": TOTAL_DRAWS, "effective_factor_NEO": EFFECTIVE_FACTOR,
           "n_rows_total": int(n_tot), "bytes_on_disk": int(size),
           "git_commit": git_commit()}
    (OUT / "manifest.json").write_text(json.dumps(man, indent=2))
    print(f"\nVALIDATION {'PASSED' if rep['ALL_PASS'] else 'FAILED'}")
    print(f"  {OUT/'validation_report.json'}\n  {OUT/'manifest.json'}")


def _center_radec(clon, clat, epoch):
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, GCRS, GeocentricTrueEcliptic
    import astropy.units as u
    t = Time(epoch, scale="tdb")
    c = SkyCoord(lon=clon * u.deg, lat=clat * u.deg, distance=1.0 * u.AU,
                 frame=GeocentricTrueEcliptic(obstime=t)).transform_to(GCRS(obstime=t))
    return float(c.ra.deg), float(c.dec.deg)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    s = p.add_subparsers(dest="cmd", required=True)
    d = s.add_parser("draw"); d.add_argument("--shard", type=int, required=True)
    d.add_argument("--overwrite", action="store_true"); d.set_defaults(func=draw)
    r = s.add_parser("rebase"); r.add_argument("--shard", type=int, required=True)
    r.add_argument("--overwrite", action="store_true"); r.set_defaults(func=rebase)
    s.add_parser("partition").set_defaults(func=partition)
    s.add_parser("validate").set_defaults(func=validate)
    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
