#!/usr/bin/env python3
"""Build the HEALPix-partitioned NEOMOD3 projection cache (docs/new_neomod_cloning.md §8).

WHY: every map task currently loads the whole 97M-row / 7.88 GB cache to extract the ~50k clones in
its own 30 deg patch -> 667 x 7.88 GB = 5.3 TB of I/O for a full-grid build, ~20 GB RAM per task, and
an expensive startup that makes ckpt preemption costly.

HOW: HEALPix tiles the sky into DISJOINT cells, so each clone is written EXACTLY ONCE (no duplication
-- unlike per-center slicing, where 667 overlapping 30 deg patches would duplicate every clone ~45x).
At read time a center pulls only the ~72 of 768 pixels its disc touches: ~0.74 GB instead of 7.88 GB.

SAFETY: this is a pure I/O reorganisation -- the clone values are byte-identical, only their file
location changes. Correctness is preserved because (1) the pixel query is a conservative SUPERSET
(inclusive=True plus a max_pixrad margin) and (2) the exact 30 deg angular cut still runs afterwards
in build_visible_subset_dataframe, so extra clones are harmless and under-reading is prevented.
The monolithic cache is KEPT; this is an addition, not a replacement.

Commands:
    build      partition the cache          (parallel; saturates the allocated cores)
    validate   T1 partition-integrity tests (row count, checksums, pixel assignment, file count)
"""
from __future__ import annotations
import argparse, json, os, shutil, sys, time
from pathlib import Path
import numpy as np

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
CACHE_DIR = W / "outputs" / "neomod3_projection_cache"
MONOLITHIC = CACHE_DIR / "neomod3_projection_20270825T000000.parquet"
SHARD_DIR = CACHE_DIR / "shards"
BY_PIXEL = CACHE_DIR / "by_pixel"
NSIDE = 8
CHECKSUM_COLS = ["vlam", "vbeta", "mag_app", "H", "ra_deg", "dec_deg"]


def _threads():
    n = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 8))
    import pyarrow as pa
    pa.set_cpu_count(n)
    pa.set_io_thread_count(n)
    print(f"[threads] using {n} cores (pyarrow cpu+io thread pools)", flush=True)
    return n


def build(args):
    import pyarrow as pa, pyarrow.parquet as pq, pyarrow.dataset as ds
    import healpy as hp
    n = _threads()
    t0 = time.time()
    if BY_PIXEL.exists():
        if not args.overwrite:
            sys.exit(f"{BY_PIXEL} exists; pass --overwrite")
        shutil.rmtree(BY_PIXEL)

    # Read the 100 existing shards as one dataset -- pyarrow parallelises the read across them.
    src = ds.dataset(str(SHARD_DIR), format="parquet")
    tbl = src.to_table(use_threads=True)
    print(f"[read] {tbl.num_rows:,} rows, {tbl.nbytes/1e9:.2f} GB in {time.time()-t0:.0f}s", flush=True)

    # HEALPix pixel from ICRS ra/dec in DEGREES (lonlat=True). Disjoint tiles -> no duplication.
    t1 = time.time()
    ra = tbl.column("ra_deg").to_numpy(zero_copy_only=False).astype(np.float64)
    dec = tbl.column("dec_deg").to_numpy(zero_copy_only=False).astype(np.float64)
    pix = hp.ang2pix(NSIDE, ra, dec, lonlat=True).astype(np.int32)
    tbl = tbl.append_column("pix", pa.array(pix))
    print(f"[ang2pix] nside={NSIDE} on {len(pix):,} points in {time.time()-t1:.0f}s "
          f"({pix.min()}..{pix.max()})", flush=True)

    # CRITICAL: sort by `pix` before writing. The cache is in projection order (random sky
    # positions), so without this every record batch scatters rows across all 768 open file
    # handles -- tiny simultaneous appends to 768 files, which on GPFS crawls at ~2 MB/s.
    # Sorted, each partition is one contiguous run: open, write, close. Measured 2 MB/s -> minutes.
    t15 = time.time()
    tbl = tbl.sort_by([("pix", "ascending")])
    print(f"[sort] by pix in {time.time()-t15:.0f}s (makes partition writes sequential)", flush=True)

    # One file per pixel: partition on `pix`, one output file per partition.
    t2 = time.time()
    ds.write_dataset(
        tbl, str(BY_PIXEL), format="parquet",
        partitioning=ds.partitioning(pa.schema([("pix", pa.int32())]), flavor="hive"),
        existing_data_behavior="overwrite_or_ignore",
        max_partitions=hp.nside2npix(NSIDE) + 8,
        max_open_files=hp.nside2npix(NSIDE) + 8,
        use_threads=True,
    )
    print(f"[write] {BY_PIXEL} in {time.time()-t2:.0f}s", flush=True)

    meta = dict(nside=NSIDE, npix=hp.nside2npix(NSIDE),
                max_pixrad_deg=float(np.degrees(hp.max_pixrad(NSIDE))),
                n_rows=int(tbl.num_rows), source=str(MONOLITHIC),
                coord="ICRS ra_deg/dec_deg, lonlat=True, degrees",
                note="Disjoint sky tiles; each clone written exactly once. Read with "
                     "hp.query_disc(nside, vec, radians(max_sep+max_pixrad), inclusive=True). "
                     "The exact angular cut still runs in build_visible_subset_dataframe.")
    (BY_PIXEL / "_healpix_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[done] total {time.time()-t0:.0f}s", flush=True)


def validate(args):
    """T1 — partition integrity (docs §8.6)."""
    import pyarrow.parquet as pq, pyarrow.dataset as ds
    import healpy as hp
    _threads()
    npix = hp.nside2npix(NSIDE)
    ok = True

    mono = pq.read_table(str(MONOLITHIC), columns=CHECKSUM_COLS, use_threads=True)
    part = ds.dataset(str(BY_PIXEL), format="parquet", partitioning="hive").to_table(
        columns=CHECKSUM_COLS + ["pix"], use_threads=True)

    print("\n=== T1.1 row count ===")
    print(f"  monolithic {mono.num_rows:,} | partitioned {part.num_rows:,}")
    r = mono.num_rows == part.num_rows
    ok &= r; print(f"  -> {'PASS' if r else 'FAIL'}")

    print("\n=== T1.2 column checksums (no rows lost, altered, or duplicated) ===")
    for c in CHECKSUM_COLS:
        a = mono.column(c).to_numpy(zero_copy_only=False).astype(np.float64)
        b = part.column(c).to_numpy(zero_copy_only=False).astype(np.float64)
        sa, sb = np.nansum(a), np.nansum(b)
        rel = abs(sa - sb) / max(abs(sa), 1e-12)
        good = rel < 1e-9 and np.isclose(np.nanmin(a), np.nanmin(b)) and np.isclose(np.nanmax(a), np.nanmax(b))
        ok &= good
        print(f"  {c:>8}: sum rel-diff {rel:.2e}  min/max match {np.isclose(np.nanmin(a),np.nanmin(b))}"
              f"/{np.isclose(np.nanmax(a),np.nanmax(b))}  -> {'PASS' if good else 'FAIL'}")

    print("\n=== T1.3 pixel assignment (catches ra/dec<->lonlat, deg/rad, frame errors) ===")
    rng = np.random.default_rng(0)
    n_s = min(1_000_000, part.num_rows)
    idx = rng.choice(part.num_rows, size=n_s, replace=False)
    ra = part.column("ra_deg").take(idx).to_numpy(zero_copy_only=False).astype(np.float64)
    dec = part.column("dec_deg").take(idx).to_numpy(zero_copy_only=False).astype(np.float64)
    stored = part.column("pix").take(idx).to_numpy(zero_copy_only=False)
    recomputed = hp.ang2pix(NSIDE, ra, dec, lonlat=True)
    nbad = int((stored != recomputed).sum())
    ok &= nbad == 0
    print(f"  sampled {n_s:,}: mismatches {nbad}  -> {'PASS' if nbad==0 else 'FAIL'}")

    print("\n=== T1.4 partition files ===")
    dirs = sorted(BY_PIXEL.glob("pix=*"))
    counts = np.bincount(part.column("pix").to_numpy(zero_copy_only=False), minlength=npix)
    n_pop = int((counts > 0).sum())
    print(f"  directories {len(dirs)} | pixels with rows {n_pop} of {npix}")
    print(f"  empty pixels: {npix-n_pop} (legitimate only if the sky there is genuinely unpopulated)")
    good = len(dirs) == n_pop
    ok &= good; print(f"  -> {'PASS' if good else 'FAIL'}")
    tot = sum(f.stat().st_size for f in BY_PIXEL.rglob('*.parquet'))
    print(f"  on disk: {tot/1e9:.2f} GB across {len(list(BY_PIXEL.rglob('*.parquet')))} files "
          f"(monolithic {MONOLITHIC.stat().st_size/1e9:.2f} GB)")

    print(f"\n{'='*60}\nT1 OVERALL: {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)



# ---------------------------------------------------------------------------- the read path
def read_patch(center_lon_ecl, center_lat_ecl, obstime_str, max_sep_deg=30.0, columns=None):
    """Read ONLY the HEALPix partitions whose pixels touch this center's disc.

    The query is a deliberate SUPERSET: inclusive=True returns every pixel that overlaps the disc
    even partially, and we add max_pixrad so no edge pixel can be missed. The EXACT angular cut is
    still applied afterwards by build_visible_subset_dataframe, so extra rows are harmless and
    under-reading is structurally prevented (docs §8.5).
    """
    import pyarrow.dataset as ds, pyarrow.compute as pc
    import healpy as hp
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, GCRS, GeocentricTrueEcliptic
    import astropy.units as u
    t_obs = Time(obstime_str, scale="tdb")
    # same center convention as _center_skycoord(center_mode="custom_ecliptic")
    c = SkyCoord(lon=center_lon_ecl*u.deg, lat=center_lat_ecl*u.deg, distance=1.0*u.AU,
                 frame=GeocentricTrueEcliptic(obstime=t_obs)).transform_to(GCRS(obstime=t_obs))
    vec = hp.ang2vec(float(c.ra.deg), float(c.dec.deg), lonlat=True)
    radius = np.radians(max_sep_deg + np.degrees(hp.max_pixrad(NSIDE)))
    pixels = np.sort(hp.query_disc(NSIDE, vec, radius, inclusive=True))
    dset = ds.dataset(str(BY_PIXEL), format="parquet", partitioning="hive")
    tbl = dset.to_table(filter=pc.field("pix").isin(pixels.tolist()),
                        columns=columns, use_threads=True)
    return tbl, pixels


def validate2(args):
    """T2 - clone-set equality per center vs the monolithic cache (docs §8.6)."""
    import pyarrow.parquet as pq
    import pandas as pd
    from astropy.time import Time
    from astropy.coordinates import get_sun, GeocentricTrueEcliptic
    _threads()
    EPOCH = "2027-08-25T00:00:00"
    t = Time(EPOCH, scale="tdb")
    sun = get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
    anti = (sun.lon.deg + 180.0) % 360.0
    # 12 centers spanning the extremes: antisun, sun-exclusion edges, ecliptic and both poles
    dlon_lat = [(0,0),(20,-12),(0,-1),(140,0),(-140,0),(0,50),(0,-50),
                (70,25),(-70,-25),(140,50),(-140,-50),(30,8)]
    cols = ["ra_deg","dec_deg","vlam","vbeta","mag_app"]
    print(f"loading monolithic cache once ...", flush=True)
    mono = pq.read_table(str(MONOLITHIC), columns=cols, use_threads=True).to_pandas()
    import healpy as hp
    from astropy.coordinates import SkyCoord, GCRS
    import astropy.units as u
    ok_all = True
    print(f"\n{'center':>16} {'pix':>5} {'partitioned':>12} {'monolithic':>11} {'equal?':>7}")
    for dl, la in dlon_lat:
        clon = (anti + dl) % 360.0
        tbl, pixels = read_patch(clon, la, EPOCH, 30.0, columns=cols)
        part = tbl.to_pandas()
        # exact 30 deg cut, identical convention, applied to BOTH
        c = SkyCoord(lon=clon*u.deg, lat=la*u.deg, distance=1.0*u.AU,
                     frame=GeocentricTrueEcliptic(obstime=t)).transform_to(GCRS(obstime=t))
        def cut(df):
            sc = SkyCoord(ra=df.ra_deg.to_numpy()*u.deg, dec=df.dec_deg.to_numpy()*u.deg,
                          frame=GCRS(obstime=t))
            return df[sc.separation(c).deg <= 30.0]
        A = cut(part); B = cut(mono)
        a = np.sort(A[["vlam","vbeta","mag_app"]].to_numpy(), axis=0)
        b = np.sort(B[["vlam","vbeta","mag_app"]].to_numpy(), axis=0)
        eq = (len(A) == len(B)) and a.shape == b.shape and np.array_equal(a, b)
        ok_all &= eq
        print(f"dlon{dl:+04d}_lat{la:+03d} {len(pixels):>5} {len(A):>12,} {len(B):>11,} {'YES' if eq else 'NO':>7}")
    print(f"\n{'='*60}\nT2 OVERALL: {'PASS' if ok_all else 'FAIL'}")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    s = p.add_subparsers(dest="cmd", required=True)
    b = s.add_parser("build"); b.add_argument("--overwrite", action="store_true"); b.set_defaults(f=build)
    v = s.add_parser("validate"); v.set_defaults(f=validate)
    v2 = s.add_parser("validate2"); v2.set_defaults(f=validate2)
    a = p.parse_args(); a.f(a)
