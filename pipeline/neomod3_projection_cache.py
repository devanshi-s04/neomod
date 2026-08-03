#!/usr/bin/env python3
"""NEOMOD3 projection cache — sample + project ONCE, slice per map center.

WHY THIS EXISTS
---------------
The NEO clone source is now NEOMOD3 (see docs/new_neomod_cloning.md). NEOMOD3 is a GLOBAL debiased
model, so clones must be drawn all-sky and then sky-cut to the map patch. Sampling per (center,
mag-bin) is hopeless: a 30 deg patch is ~6.7% of the sky, and the draw-to-in-map-range cost is
strongly H-dependent (measured, outputs/neomod3_H_cut_test):

    H band   draws per object landing in the maps' 14-25 apparent-mag range
    15-25          16
    25-26         189
    26-27         631
    27-28       2,208

BUT the projection is EPOCH-dependent and NOT CENTER-dependent. So we sample and project once,
globally, and every one of the 667 map centers just slices the same cache. This is exactly the
Stage 0 pattern from fixing_integrator.md §12.2 -- and it plugs into the SAME mechanism:
`build_visible_subset_dataframe` detects the cache columns (ra_deg, dec_deg, vlam, vbeta) and skips
straight to the center-dependent sky cut. No new slice code, no new sky-cut convention.

WHY NO H CUT
------------
H is ABSOLUTE magnitude; the map bins and LSST's limit are APPARENT. For NEOs these decouple
(mag_app ~ H + 5log10(r*Delta)), because NEOs come close. Measured on 2M projected samples: objects
landing in the map range span the FULL H range 15.23-28.00, and cutting at H<25 would delete 35.7%
of all detectable NEOs and **81.2% of the detectable |v|>2 NEOs** -- precisely the high-velocity
population whose absence caused the kNN-bleed / support-mask artifact in the +-5 grid pilot
(outputs/vdp_bleed_diagnostic_bundle). The faint population is twice as fast on average
(median |v| 1.55 vs 0.76 deg/day; 38.3% vs 4.9% above |v|=2). So the cache keeps 15 <= H <= 28.

USAGE
-----
    # 1. build shards (Slurm array; each shard draws n_total/nshards orbits)
    python neomod3_projection_cache.py build --shard $SLURM_ARRAY_TASK_ID --nshards 40 \
        --n-orbits-total 40000000 --epoch 2027-08-25T00:00:00
    # 2. merge
    python neomod3_projection_cache.py merge
    # 3. inspect
    python neomod3_projection_cache.py stats

Then map generation passes the merged parquet as the NEO clone source; the existing cache branch in
build_visible_subset_dataframe does the rest.

Read-only w.r.t. everything else. Writes only under outputs/neomod3_projection_cache/.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
os.chdir(W / "neomod")  # s3m_loader searches ['.', 'S3Mdata'] relative to CWD

DEFAULT_EPOCH = "2027-08-25T00:00:00"          # must match the map --ref-obstime
OUT_DIR = W / "outputs" / "neomod3_projection_cache"
H_RANGE = (15.0, 28.0)                          # NEOMOD3's calibrated range; deliberately NOT cut

# Columns the pipeline's cache branch needs, plus what the cloner/mag-binning needs downstream.
#   ra_deg, dec_deg, vlam, vbeta -> _CACHE_POSITION_COLS (triggers the skip-propagation path)
#   a, node, argperi, t_p        -> L_obs_deg is derived from these if absent
#   mag_app                      -> magnitude-bin assignment
CACHE_COLS = ["a", "e", "i", "node", "argperi", "t_p", "q", "H",
              "ra_deg", "dec_deg", "dra_deg_day", "ddec_deg_day",
              "lam_deg", "beta_deg", "vlam", "vbeta", "mag_app", "M_obs_deg"]


def _shard_path(shard: int) -> Path:
    return OUT_DIR / "shards" / f"nm3_proj_{shard:04d}.parquet"


def merged_path(epoch: str) -> Path:
    return OUT_DIR / f"neomod3_projection_{epoch.replace(':', '').replace('-', '')}.parquet"


# ---------------------------------------------------------------------------- build
def build(args) -> None:
    import velocity_density_pipeline_gmm as vdp
    import neomod3_sampler as nm3s

    t0 = time.time()
    (OUT_DIR / "shards").mkdir(parents=True, exist_ok=True)
    out = _shard_path(args.shard)
    if out.exists() and not args.overwrite:
        print(f"[skip] {out} exists", flush=True)
        return

    n_this = args.n_orbits_total // args.nshards
    # per-shard seed: reproducible, and independent across shards
    rng = np.random.default_rng(args.seed + args.shard)
    print(f"shard {args.shard}/{args.nshards}: drawing {n_this:,} NEOMOD3 orbits "
          f"(H in {H_RANGE}) at epoch {args.epoch}", flush=True)

    df = nm3s.sample_neomod3_orbits(n_this, args.epoch, rng=rng)
    print(f"  valid NEO orbits (q<1.3): {len(df):,}  ({time.time()-t0:.0f}s)", flush=True)
    if not len(df):
        raise RuntimeError("NEOMOD3 sampler returned no valid orbits")

    # observer geometry only -- the scorer is population-independent here (same trick the
    # n-body map driver uses to avoid loading the 13.9M-row MBA file per task)
    _, scorer = vdp.load_s3m_population("neo", verbose=False)

    # ALL-SKY projection: max_sep_deg=180 keeps everything, so the cache is center-agnostic.
    # Each map center applies its own 30 deg cut later via the same function.
    vis = vdp.build_visible_subset_dataframe(
        df, obstime_str=args.epoch, scorer=scorer, max_sep_deg=180.0,
        chunk=args.chunk, show_progress=False,
        center_mode="custom_ecliptic", center_lon_deg=0.0, center_lat_deg=0.0,
    )
    print(f"  projected: {len(vis):,}  ({time.time()-t0:.0f}s)", flush=True)

    keep = [c for c in CACHE_COLS if c in vis.columns]
    missing = [c for c in ("ra_deg", "dec_deg", "vlam", "vbeta", "mag_app") if c not in keep]
    if missing:
        raise RuntimeError(f"projection is missing required cache columns: {missing}")
    sub = vis[keep].copy()
    for c in sub.columns:                       # halve the on-disk size; well within precision
        if sub[c].dtype == np.float64:
            sub[c] = sub[c].astype(np.float32)
    tmp = out.with_suffix(".tmp.parquet")
    sub.to_parquet(tmp, index=False)
    tmp.replace(out)
    print(f"  wrote {len(sub):,} rows -> {out}  ({time.time()-t0:.0f}s)", flush=True)


# ---------------------------------------------------------------------------- merge
def merge(args) -> None:
    shards = sorted((OUT_DIR / "shards").glob("nm3_proj_*.parquet"))
    if not shards:
        raise FileNotFoundError(f"no shards in {OUT_DIR/'shards'}")
    print(f"merging {len(shards)} shards ...", flush=True)
    frames = []
    for p in shards:
        d = pd.read_parquet(p)
        frames.append(d)
        print(f"  {p.name}: {len(d):,}", flush=True)
    cache = pd.concat(frames, ignore_index=True)
    out = merged_path(args.epoch)
    tmp = out.with_suffix(".tmp.parquet")
    cache.to_parquet(tmp, index=False)
    tmp.replace(out)
    print(f"\nwrote {len(cache):,} rows -> {out}  ({out.stat().st_size/1e9:.2f} GB)", flush=True)
    _report(cache, args.epoch)


# ---------------------------------------------------------------------------- stats
def stats(args) -> None:
    p = merged_path(args.epoch)
    if not p.exists():
        raise FileNotFoundError(p)
    _report(pd.read_parquet(p), args.epoch)


def _report(cache: pd.DataFrame, epoch: str) -> None:
    v = np.maximum(cache.vlam.abs(), cache.vbeta.abs())
    inmap = cache.mag_app.between(14.0, 25.0)
    print(f"\n=== NEOMOD3 projection cache — epoch {epoch} ===", flush=True)
    print(f"  rows: {len(cache):,}", flush=True)
    print(f"  H      : {cache.H.min():.2f} .. {cache.H.max():.2f}   (no cut, by design)", flush=True)
    print(f"  mag_app: median {cache.mag_app.median():.2f}   in 14-25 map range: "
          f"{int(inmap.sum()):,} ({100*inmap.mean():.2f}%)", flush=True)
    print(f"  |v|    : median {np.median(v):.3f}  p99 {np.percentile(v, 99):.3f}  max {v.max():.1f} deg/day", flush=True)
    im = cache[inmap]
    if len(im):
        vi = np.maximum(im.vlam.abs(), im.vbeta.abs())
        print(f"  IN-MAP subset ({len(im):,}): |v|>2 -> {int((vi > 2).sum()):,} ({100*(vi > 2).mean():.1f}%)  "
              f"| from H>=25: {100*(im.H >= 25).mean():.1f}%", flush=True)
    # what one 30 deg patch gets (6.7% of sky)
    frac = 0.067
    print(f"\n  expected per 30-deg map center (~{100*frac:.1f}% of sky):", flush=True)
    print(f"    all clones ~{int(len(cache)*frac):,}   in-map-range ~{int(inmap.sum()*frac):,}", flush=True)
    print(f"    of which |v|>2 ~{int((np.maximum(im.vlam.abs(), im.vbeta.abs()) > 2).sum()*frac):,}"
          if len(im) else "", flush=True)
    print("\n  -> feed this parquet to the NEO clone source; build_visible_subset_dataframe detects\n"
          "     ra_deg/dec_deg/vlam/vbeta and applies ONLY the per-center sky cut (§12.2 path).", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build")
    b.add_argument("--shard", type=int, required=True)
    b.add_argument("--nshards", type=int, default=40)
    b.add_argument("--n-orbits-total", type=int, default=40_000_000)
    b.add_argument("--epoch", default=DEFAULT_EPOCH)
    b.add_argument("--chunk", type=int, default=200_000)
    b.add_argument("--seed", type=int, default=42)
    b.add_argument("--overwrite", action="store_true")
    b.set_defaults(func=build)

    m = sub.add_parser("merge")
    m.add_argument("--epoch", default=DEFAULT_EPOCH)
    m.set_defaults(func=merge)

    s = sub.add_parser("stats")
    s.add_argument("--epoch", default=DEFAULT_EPOCH)
    s.set_defaults(func=stats)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
