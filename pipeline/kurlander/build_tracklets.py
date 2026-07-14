#!/usr/bin/env python3
"""build_tracklets.py — Kurlander (2025) detections -> 2+-detection nightly tracklets.

Works identically on Arnor (Epyc NFS) and Hyak (rsync'd copy): pass --root.

Files are partitioned by OBJECT BATCH (not time), so any subset of chunks is a uniformly
random object subsample — and 69% of object-nights already have >=2 detections inside a
single chunk, so tracklets form without cross-file joins.

Tracklet = all detections of one ObjID on one night (floor(fieldMJD_TAI)), >=2 detections.
Rates are least-squares slopes in RA (RAW alpha_dot, wrap-safe) and Dec, matching the v5
parquet's `mean_dra` convention -> score with dra_cosdec=False.

Truth label: NEO iff q < 1.3 au (q from COM files directly; from a(1-e) for the KEP NEO file).
Uses NOISY astrometry (RA_deg/Dec_deg) and measured mag (trailedSourceMag) — same information
a real pipeline (and digest2) would see. Never reads Range/RangeRate (truth-only columns).
"""
import argparse, glob, os, sys, time
import numpy as np, pandas as pd

POPS = {  # subdir -> (glob, population label)
    "neomod":    ("outfiles/neo_output_1.h5",   "NEO"),
    "s3m":       ("outfiles/*_*.h5",            "MBA"),
    "trojanmod": ("outfiles/*.h5",              "Trojan"),
    "cfeps":     ("outfiles/*.h5",              "TNO"),
}
NEED = ["ObjID", "fieldMJD_TAI", "RA_deg", "Dec_deg", "trailedSourceMag", "optFilter", "Linked"]


def tracklets_from_df(df, pop):
    """Group by (ObjID, night); LSQ position+rate; one row per tracklet. Fully vectorised
    (a per-group Python loop costs ~40 s/file; this is ~2 s/file)."""
    df = df.dropna(subset=["RA_deg", "Dec_deg", "fieldMJD_TAI", "trailedSourceMag"])
    if len(df) == 0:
        return pd.DataFrame()
    df = df.assign(night=np.floor(df.fieldMJD_TAI.values).astype(int))
    g = df.groupby(["ObjID", "night"], sort=False)
    df = df[g["fieldMJD_TAI"].transform("size").values >= 2]
    if len(df) == 0:
        return pd.DataFrame()

    g = df.groupby(["ObjID", "night"], sort=False)
    t = df.fieldMJD_TAI.values
    # unwrap RA about each group's first detection to kill the 0/360 seam
    ra_first = g["RA_deg"].transform("first").values
    ra_u = ra_first + ((df.RA_deg.values - ra_first + 180.0) % 360.0 - 180.0)
    df = df.assign(_rau=ra_u)
    g = df.groupby(["ObjID", "night"], sort=False)

    tbar = g["fieldMJD_TAI"].transform("mean").values
    rbar = g["_rau"].transform("mean").values
    dbar = g["Dec_deg"].transform("mean").values
    dt = t - tbar
    df = df.assign(_num_ra=dt * (ra_u - rbar), _num_dec=dt * (df.Dec_deg.values - dbar), _den=dt**2)

    a = df.groupby(["ObjID", "night"], sort=False).agg(
        n_det=("fieldMJD_TAI", "size"), mjd_mean=("fieldMJD_TAI", "mean"),
        mean_ra=("_rau", "mean"), mean_dec=("Dec_deg", "mean"),
        mean_mag=("trailedSourceMag", "mean"), filter=("optFilter", "first"),
        linked=("Linked", "any"), tmin=("fieldMJD_TAI", "min"), tmax=("fieldMJD_TAI", "max"),
        num_ra=("_num_ra", "sum"), num_dec=("_num_dec", "sum"), den=("_den", "sum"),
    ).reset_index()
    a = a[a.den > 0]
    a["mean_dra"] = a.num_ra / a.den                 # RAW alpha_dot [deg/day]
    a["mean_ddec"] = a.num_dec / a.den
    a["mean_ra"] = a.mean_ra % 360.0
    a["span_min"] = (a.tmax - a.tmin) * 1440.0
    a["population"] = pop
    a["ObjID"] = a.ObjID.astype(str)
    return a[["ObjID", "population", "night", "n_det", "mjd_mean", "mean_ra", "mean_dec",
              "mean_dra", "mean_ddec", "mean_mag", "filter", "linked", "span_min"]]


def add_truth(tr, src):
    """Attach q (and a,e,inc) from the source frame; label NEO iff q<1.3."""
    cols = src.columns
    key = src.drop_duplicates("ObjID").set_index("ObjID")
    if "q" in cols:
        q = key["q"]
    else:                                     # KEP file (NEO): q = a(1-e)
        q = key["a"] * (1 - key["e"])
    tr["q_au"] = tr.ObjID.map(q)
    for c in ("e", "inc"):
        if c in cols:
            tr[c] = tr.ObjID.map(key[c])
    tr["is_neo"] = tr.q_au < 1.3
    return tr


def process_file(path, pop):
    df = pd.read_hdf(path)
    keep = [c for c in NEED if c in df.columns]
    tr = tracklets_from_df(df[keep], pop)
    if len(tr):
        tr = add_truth(tr, df)
    return tr


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/astro/users/jkurla/public_html/LSST_Sorcha_predictions")
    ap.add_argument("--pops", nargs="+", default=["neomod", "s3m", "trojanmod", "cfeps"])
    ap.add_argument("--max-files", type=int, default=0, help="per population; 0 = all")
    ap.add_argument("--file-stride", type=int, default=1, help="take every Nth file (random object subset)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--shard", type=int, default=-1, help="Slurm array id: process only this file index")
    ap.add_argument("--nshards", type=int, default=1)
    a = ap.parse_args()

    tasks = []
    for sub in a.pops:
        pat, pop = POPS[sub]
        files = sorted(glob.glob(os.path.join(a.root, sub, pat)))
        files = files[::a.file_stride]
        if a.max_files:
            files = files[:a.max_files]
        tasks += [(f, pop) for f in files]
    if a.shard >= 0:
        tasks = tasks[a.shard::a.nshards]
    print(f"{len(tasks)} files to process", flush=True)

    parts, t0 = [], time.time()
    for i, (f, pop) in enumerate(tasks):
        t = time.time()
        try:
            tr = process_file(f, pop)
        except Exception as ex:
            print(f"  FAIL {f}: {ex}", flush=True); continue
        if len(tr):
            parts.append(tr)
        print(f"  [{i+1}/{len(tasks)}] {pop:7s} {os.path.basename(f):20s} "
              f"-> {len(tr):>7,} tracklets ({time.time()-t:.1f}s)", flush=True)
    if not parts:
        sys.exit("no tracklets built")
    out = pd.concat(parts, ignore_index=True)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    out.to_parquet(a.out, index=False)
    print(f"\nwrote {a.out}: {len(out):,} tracklets in {time.time()-t0:.0f}s")
    print(out.population.value_counts().to_dict())
    print(f"NEO frac (q<1.3): {out.is_neo.mean():.4f}")
