#!/usr/bin/env python3
"""Targeted benchmark for the same-object velocity validation (Arnor request 2026-07-12).

Builds benchmark rates at the v3 epoch (MJD 61642 = 2027-08-25) for EXACTLY the S3M
objects Sorcha detected on night 61642 (union across case1/2/3, 11,612 ObjIDs), so
identity is SHARED BY CONSTRUCTION. Every row carries `s3m_objid` (the S* id), letting
Arnor join Sorcha<->benchmark on true identity instead of the failed (e,H) cKDTree match
(see docs/NOTE_FOR_HYAK_benchmark_objid_mapping.md).

No population caps -> 100% overlap with Sorcha's night-61642 set. Rates are computed with
the SAME `score_orbital_df` machinery as benchmark v3, so they are apples-to-apples with
the v3 parquet (only the object set and the identity column differ).

This is a velocity-layer validation set: it carries the rate/geometry columns but is NOT
VDP/digest2-scored (no P_NEO columns) — the same-object velocity comparison doesn't need them.

Run as a 4-task array (one per population), then --combine-only:
    sbatch neomod/pipeline/slurm/benchmark_night61642.sbatch
    python neomod/pipeline/gen_benchmark_night61642.py --combine-only
"""
from __future__ import annotations
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.coordinates import get_sun, GeocentricTrueEcliptic

WORKDIR = "/mmfs1/gscratch/dirac/ds2004/sorcha"
NEOMD = os.path.join(WORKDIR, "neomod")
sys.path.insert(0, os.path.join(NEOMD, "adam_core_stub"))
sys.path.insert(0, os.path.join(NEOMD, "src"))
sys.path.insert(0, os.path.join(NEOMD, "pipeline"))

import gen_benchmark_tracklets_s3m as bench  # reuse constants + helpers (build_grid, etc.)

OUT_DIR = os.path.join(WORKDIR, "outputs", "benchmark_night61642")
TARGET_IDS_FILE = os.path.join(OUT_DIR, "target_objids.txt")
COMBINED = os.path.join(OUT_DIR, "benchmark_night61642.parquet")


def load_targets() -> set[str]:
    with open(TARGET_IDS_FILE) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def generate(pop_id: int) -> None:
    pop_label, s3m_pop, _cap = bench.POPULATIONS[pop_id]
    out_path = os.path.join(OUT_DIR, f"tracklets_{pop_label}.parquet")
    targets = load_targets()

    t_ref = Time(bench.REF_OBSTIME, scale="utc")
    mjd_ref = t_ref.mjd
    sun = get_sun(t_ref).transform_to(GeocentricTrueEcliptic(obstime=t_ref))
    antisun_lon = (sun.lon.deg + 180.0) % 360.0
    print(f"[{pop_label}] epoch {bench.REF_OBSTIME} (MJD {mjd_ref:.3f})", flush=True)

    os.chdir(NEOMD)
    import velocity_density_pipeline_gmm as vdp
    seed = vdp.ProbMapSet.from_npz(
        os.path.join(bench.PROB_MAPS_DIR, "prob_maps_grid_dlon+000_lat+00.npz"))

    df, scorer = vdp.load_s3m_population(s3m_pop)
    df = df[df["OID"].astype(str).str.strip().isin(targets)].reset_index(drop=True)
    print(f"[{pop_label}] target objects present in this population: {len(df):,}", flush=True)
    os.makedirs(OUT_DIR, exist_ok=True)
    if len(df) == 0:
        pd.DataFrame().to_parquet(out_path, index=False)
        print(f"[{pop_label}] none — wrote empty {out_path}", flush=True)
        return

    _, pop_df = seed.score_orbital_df(
        df=df, scorer=scorer, obstime_str=bench.REF_OBSTIME,
        max_sep_deg=200.0, chunk=50_000, show_progress=False, return_visible=True,
    )
    pop_df = pop_df.reset_index(drop=True)
    print(f"[{pop_label}] {len(pop_df):,} propagated/visible", flush=True)

    if "lam_deg" in pop_df.columns and "beta_deg" in pop_df.columns:
        lam = pop_df["lam_deg"].to_numpy(float)
        beta = pop_df["beta_deg"].to_numpy(float)
    else:
        lam, beta = bench.equatorial_to_ecliptic(
            pop_df["ra_deg"].to_numpy(float), pop_df["dec_deg"].to_numpy(float), t_ref)

    dlon = ((lam - antisun_lon + 180.0) % 360.0) - 180.0
    in_grid = np.abs(dlon) <= bench.DLON_LIMIT + 0.5
    n_excl = int((~in_grid).sum())
    pop_df = pop_df[in_grid].reset_index(drop=True)
    dlon, beta, lam = dlon[in_grid], beta[in_grid], lam[in_grid]
    print(f"[{pop_label}] sun-exclusion removed {n_excl:,}, kept {len(pop_df):,}", flush=True)

    grid_arr = bench.build_grid()
    cell = bench.nearest_cell(dlon, beta, grid_arr)
    asgn_dlon, asgn_lat = grid_arr[cell, 0], grid_arr[cell, 1]

    ra0 = pop_df["ra_deg"].to_numpy(float)
    dec0 = pop_df["dec_deg"].to_numpy(float)
    dra = pop_df["dra_deg_day"].to_numpy(float)
    ddec = pop_df["ddec_deg_day"].to_numpy(float)
    mag = pop_df["mag_app"].to_numpy(float)
    ra1 = (ra0 + dra * bench.DT_DAYS) % 360.0
    dec1 = np.clip(dec0 + ddec * bench.DT_DAYS, -89.99, 89.99)

    n = len(pop_df)
    out = pd.DataFrame()
    out["s3m_objid"] = pop_df["OID"].astype(str).str.strip().to_numpy()   # <-- true identity
    out["ObjID"] = [f"BM{i:08d}" for i in range(n)]                        # BM id kept for parity
    out["population"] = pop_label
    out["prob_map_file"] = [bench.map_filename(d, l) for d, l in zip(asgn_dlon, asgn_lat)]
    out["prob_map"] = [bench.map_label(d, l) for d, l in zip(asgn_dlon, asgn_lat)]
    out["n_det_per_night"] = 2
    out["mean_ra"], out["mean_dec"] = ra0, dec0
    out["mean_dra"], out["mean_ddec"] = dra, ddec
    out["mean_mag"] = mag
    out["ra0"], out["dec0"], out["mjd0_utc"], out["mag0"] = ra0, dec0, mjd_ref, mag
    out["ra1"], out["dec1"], out["mjd1_utc"], out["mag1"] = ra1, dec1, mjd_ref + bench.DT_DAYS, mag
    out["lam_deg"], out["beta_deg"] = lam, beta
    out["dlon_from_antisun_deg"] = dlon
    for col in ["H", "vlam", "vbeta", "a_au", "e", "q_au"]:
        if col in pop_df.columns:
            out[col] = pop_df[col].to_numpy(float)

    out.to_parquet(out_path, index=False)
    print(f"[{pop_label}] wrote {len(out):,} rows to {out_path}", flush=True)


def combine() -> None:
    files = sorted(glob.glob(os.path.join(OUT_DIR, "tracklets_*.parquet")))
    frames = [pd.read_parquet(f) for f in files]
    frames = [f for f in frames if len(f)]
    out = pd.concat(frames, ignore_index=True)
    out.to_parquet(COMBINED, index=False)
    print(f"wrote {COMBINED}: {len(out):,} rows", flush=True)
    print("populations:", out["population"].value_counts().to_dict(), flush=True)
    print(f"unique s3m_objid: {out['s3m_objid'].nunique():,}", flush=True)
    # coverage vs the requested target set
    req = load_targets()
    got = set(out["s3m_objid"].unique())
    print(f"coverage: {len(got & req):,} / {len(req):,} requested objects "
          f"({100*len(got & req)/len(req):.1f}%); {len(req - got):,} not covered "
          f"(mostly outside the benchmark population's element box in load_s3m_population "
          f"-- Sorcha's population labels are broader than the strict S3M element cuts; "
          f"a smaller number are sun-excluded, |dlon|>140deg).", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pop-id", type=int, choices=[0, 1, 2, 3])
    ap.add_argument("--combine-only", action="store_true")
    a = ap.parse_args()
    if a.combine_only:
        combine()
    elif a.pop_id is not None:
        generate(a.pop_id)
    else:
        ap.error("provide --pop-id (0-3) or --combine-only")
