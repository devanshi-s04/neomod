#!/usr/bin/env python3
"""
Generate pure-S3M benchmark tracklets against the 667-map GMM grid.

Propagates one S3M population to the grid reference epoch (2026-01-01T00:00:00),
assigns each object to its nearest map cell in the 667-map antisun-relative grid,
and writes a per-population parquet in the same schema that sorcha_phase2.py expects.

Run as a 4-task SLURM array (one task per population):
    sbatch neomod/pipeline/slurm/benchmark_gen_tracklets_s3m.sh

After all 4 tasks complete, combine on the login node:
    python neomod/pipeline/gen_benchmark_tracklets_s3m.py --combine-only

Then score:
    sbatch neomod/pipeline/slurm/benchmark_score_vdp_s3m.sh
    sbatch neomod/pipeline/slurm/benchmark_run_digest2_s3m.sh
    python neomod/pipeline/sorcha_phase2.py combine \
        --work-dir outputs/phase2_benchmark_s3m \
        --outfile outputs/phase2_benchmark_s3m/benchmark_comparison_s3m.parquet \
        --overwrite

Grid geometry (must match sorcha_gen_maps_grid.py exactly):
  - ref epoch:     2026-01-01T00:00:00
  - lon step:      10 deg
  - sun exclusion: 40 deg  (|dlon| <= 140 deg)
  - lat base:      0,1,2,3,4,5,8,12,18,25,35,50  (symmetric, 23 latitudes)
  - total maps:    667

Population mapping (--pop-id):
  0 = NEO      (neo,    12 900)
  1 = MBA      (mba,    650 000)
  2 = TNO      (tno,    1 600)
  3 = Trojans  (trojan, 6 000)

Caps chosen 2026-07-06 (Option A, see WAGG_SORCHA_HYAK_CONTEXT.md "benchmark
population-cap discovery") to make the benchmark's population mix proportional
to the TRUE S3M catalog (NEO 1.87%, MBA 96.54%, TNO 0.34%, Trojan 1.25%), back-
solved through each population's own observed post-filter survival rate so the
FINAL scored counts land near those proportions at a similar total size
(~475,000 rows) to the original (arbitrarily-capped) benchmark. Original run
used NEO=no cap, MBA=200_000, TNO=no cap, Trojan=100_000, which suppressed MBA
~69x below its true share while leaving TNO uncapped -- inflating TNO's
apparent share of the benchmark relative to a real Sorcha survey.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from astropy.coordinates import GeocentricTrueEcliptic, SkyCoord, get_sun
from astropy.time import Time
import astropy.units as u

WORKDIR = "/mmfs1/gscratch/dirac/ds2004/sorcha"
NEOMD = os.path.join(WORKDIR, "neomod")
sys.path.insert(0, os.path.join(NEOMD, "adam_core_stub"))
sys.path.insert(0, os.path.join(NEOMD, "src"))

# v3 (2026-07-07): epoch matched to MJD 61642 (2027-08-25), the busiest single
# Sorcha night across all of case1/2/3 (4,316 / 4,490 / 4,353 tracklets resp.,
# highest of any night in every case) -- see sorcha_full_pipeline.md sec 9 and
# WAGG_SORCHA_HYAK_CONTEXT.md. Maps are antisun-relative / epoch-independent
# (sorcha_gen_maps_grid.py), so only REF_OBSTIME needs to change; grid geometry
# and cell assignment logic are unaffected. Combines with v2's proportional
# population caps (Option A) for a benchmark that is BOTH population- and
# epoch-representative of a real Sorcha survey night.
REF_OBSTIME = "2027-08-25T00:00:00"
DT_DAYS = 30.0 / 1440.0
PROB_MAPS_DIR = os.path.join(WORKDIR, "prob_maps_grid_s3m")
# v3: proportional-cap benchmark (same caps as v2) + busiest-night epoch match.
# Writes to a NEW directory so v1 and v2 are left untouched.
OUT_DIR = os.path.join(WORKDIR, "outputs", "benchmark_tracklets_s3m_v3")
COMBINED_PARQUET = os.path.join(OUT_DIR, "tracklets_benchmark_v3.parquet")

POPULATIONS = {
    0: ("NEO",     "neo",    12_900),
    1: ("MBA",     "mba",    650_000),
    2: ("TNO",     "tno",    1_600),
    3: ("Trojans", "trojan", 6_000),
}

LON_STEP = 10.0
LAT_BASE = [0, 1, 2, 3, 4, 5, 8, 12, 18, 25, 35, 50]
SUN_EXCLUSION = 40.0
DLON_LIMIT = 180.0 - SUN_EXCLUSION   # 140.0


def build_grid() -> np.ndarray:
    dlons = []
    d = -180.0
    while d < 180.0 - 1e-9:
        if abs(d) <= DLON_LIMIT + 1e-9:
            dlons.append(round(d, 6))
        d += LON_STEP
    lat_set = set()
    for v in LAT_BASE:
        lat_set.add(float(v))
        lat_set.add(float(-v))
    lats = sorted(lat_set)
    return np.array([(dlon, lat) for lat in lats for dlon in dlons], dtype=np.float64)


def map_filename(dlon: float, lat: float) -> str:
    return f"prob_maps_grid_dlon{int(round(dlon)):+04d}_lat{int(round(lat)):+03d}.npz"


def map_label(dlon: float, lat: float) -> str:
    return f"grid_dlon{int(round(dlon)):+04d}_lat{int(round(lat)):+03d}"


def nearest_cell(dlon_obj: np.ndarray, lat_obj: np.ndarray,
                 grid_arr: np.ndarray) -> np.ndarray:
    ddlon = dlon_obj[:, None] - grid_arr[None, :, 0]
    dlat  = lat_obj[:, None]  - grid_arr[None, :, 1]
    return np.argmin(ddlon**2 + dlat**2, axis=1)


def equatorial_to_ecliptic(ra_deg: np.ndarray, dec_deg: np.ndarray,
                            obstime: Time) -> tuple[np.ndarray, np.ndarray]:
    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    ecl = c.transform_to(GeocentricTrueEcliptic(obstime=obstime))
    return ecl.lon.deg, ecl.lat.deg


def generate_population(pop_id: int) -> pd.DataFrame:
    pop_label, s3m_pop, max_obj = POPULATIONS[pop_id]
    out_path = os.path.join(OUT_DIR, f"tracklets_{pop_label}.parquet")

    grid_arr = build_grid()

    t_ref = Time(REF_OBSTIME, scale="utc")
    mjd_ref = t_ref.mjd
    sun = get_sun(t_ref).transform_to(GeocentricTrueEcliptic(obstime=t_ref))
    antisun_lon = (sun.lon.deg + 180.0) % 360.0
    print(f"[{pop_label}] antisun_lon={antisun_lon:.4f} deg  MJD={mjd_ref:.3f}", flush=True)

    os.chdir(NEOMD)
    import velocity_density_pipeline_gmm as vdp  # noqa: E402

    seed_map_path = os.path.join(PROB_MAPS_DIR, "prob_maps_grid_dlon+000_lat+00.npz")
    seed_pms = vdp.ProbMapSet.from_npz(seed_map_path)

    rng = np.random.default_rng(42)
    df, scorer = vdp.load_s3m_population(s3m_pop)
    print(f"[{pop_label}] loaded {len(df):,} objects", flush=True)

    if max_obj is not None and len(df) > max_obj:
        idx = rng.choice(len(df), size=max_obj, replace=False)
        df = df.iloc[idx].reset_index(drop=True)
        print(f"[{pop_label}] subsampled to {len(df):,}", flush=True)

    print(f"[{pop_label}] propagating to {REF_OBSTIME} ...", flush=True)
    _, pop_df = seed_pms.score_orbital_df(
        df=df,
        scorer=scorer,
        obstime_str=REF_OBSTIME,
        max_sep_deg=200.0,
        chunk=50_000,
        show_progress=False,
        return_visible=True,
    )
    pop_df = pop_df.reset_index(drop=True)
    print(f"[{pop_label}] {len(pop_df):,} visible objects after propagation", flush=True)

    # Ecliptic coordinates
    if "lam_deg" in pop_df.columns and "beta_deg" in pop_df.columns:
        lam  = pop_df["lam_deg"].to_numpy(float)
        beta = pop_df["beta_deg"].to_numpy(float)
    else:
        print(f"[{pop_label}] computing ecliptic coords via astropy ...", flush=True)
        lam, beta = equatorial_to_ecliptic(
            pop_df["ra_deg"].to_numpy(float),
            pop_df["dec_deg"].to_numpy(float),
            t_ref,
        )

    dlon = ((lam - antisun_lon + 180.0) % 360.0) - 180.0

    in_grid = np.abs(dlon) <= DLON_LIMIT + 0.5
    n_excl = (~in_grid).sum()
    pop_df = pop_df[in_grid].reset_index(drop=True)
    dlon   = dlon[in_grid]
    beta   = beta[in_grid]
    lam    = lam[in_grid]
    print(f"[{pop_label}] sun exclusion removed {n_excl:,}, kept {len(pop_df):,}", flush=True)

    cell_idx   = nearest_cell(dlon, beta, grid_arr)
    asgn_dlon  = grid_arr[cell_idx, 0]
    asgn_lat   = grid_arr[cell_idx, 1]

    ra0  = pop_df["ra_deg"].to_numpy(float)
    dec0 = pop_df["dec_deg"].to_numpy(float)
    dra  = pop_df["dra_deg_day"].to_numpy(float)
    ddec = pop_df["ddec_deg_day"].to_numpy(float)
    mag  = pop_df["mag_app"].to_numpy(float)

    ra1  = (ra0 + dra * DT_DAYS) % 360.0
    dec1 = np.clip(dec0 + ddec * DT_DAYS, -89.99, 89.99)

    n = len(pop_df)
    out = pd.DataFrame()
    out["ObjID"]       = [f"BM{i:08d}" for i in range(n)]
    out["population"]  = pop_label
    out["prob_map_file"] = [map_filename(d, l) for d, l in zip(asgn_dlon, asgn_lat)]
    out["prob_map"]      = [map_label(d, l)    for d, l in zip(asgn_dlon, asgn_lat)]
    out["n_det_per_night"] = 2
    # VDP scoring columns
    out["mean_ra"]    = ra0
    out["mean_dec"]   = dec0
    out["mean_dra"]   = dra
    out["mean_ddec"]  = ddec
    out["mean_mag"]   = mag
    # digest2 columns
    out["ra0"]       = ra0
    out["dec0"]      = dec0
    out["mjd0_utc"]  = mjd_ref
    out["mag0"]      = mag
    out["ra1"]       = ra1
    out["dec1"]      = dec1
    out["mjd1_utc"]  = mjd_ref + DT_DAYS
    out["mag1"]      = mag
    # analysis columns
    out["lam_deg"]               = lam
    out["beta_deg"]              = beta
    out["dlon_from_antisun_deg"] = dlon
    for col in ["H", "vlam", "vbeta", "a_au", "e", "q_au"]:
        if col in pop_df.columns:
            out[col] = pop_df[col].to_numpy(float)

    os.makedirs(OUT_DIR, exist_ok=True)
    out.to_parquet(out_path, index=False)
    print(f"[{pop_label}] wrote {len(out):,} rows to {out_path}", flush=True)
    return out


def combine_populations() -> None:
    import glob
    pop_files = sorted(glob.glob(os.path.join(OUT_DIR, "tracklets_[MNTT]*.parquet")))
    if not pop_files:
        raise FileNotFoundError(f"No per-population tracklet files found in {OUT_DIR}")
    print(f"Combining {len(pop_files)} population files:", flush=True)
    frames = []
    for f in pop_files:
        df = pd.read_parquet(f)
        print(f"  {os.path.basename(f)}: {len(df):,} rows", flush=True)
        frames.append(df)
    combined = pd.concat(frames, ignore_index=True)
    print(f"Total: {len(combined):,} rows", flush=True)
    for pop, cnt in combined["population"].value_counts().sort_index().items():
        print(f"  {pop}: {cnt:,}", flush=True)
    abs_dlon = np.abs(combined["dlon_from_antisun_deg"].to_numpy(float))
    print("Direction bins:", flush=True)
    for label, lo, hi in [("0-20", 0, 20), ("20-40", 20, 40), ("40-70", 40, 70),
                           ("70-110", 70, 110), ("110-141", 110, 141)]:
        n = ((abs_dlon >= lo) & (abs_dlon < hi)).sum()
        print(f"  {label:>8}°: {n:,}", flush=True)
    combined.to_parquet(COMBINED_PARQUET, index=False)
    print(f"\nWrote combined file: {COMBINED_PARQUET}", flush=True)
    import math
    n_tasks = math.ceil(len(combined) / 5000)
    print(f"digest2 array size: 0-{n_tasks - 1}  (5000 rows/task)", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pop-id", type=int, choices=[0, 1, 2, 3],
                        help="Population index: 0=NEO 1=MBA 2=TNO 3=Trojans")
    parser.add_argument("--combine-only", action="store_true",
                        help="Skip generation; just combine existing per-pop files")
    args = parser.parse_args()

    if args.combine_only:
        combine_populations()
    elif args.pop_id is not None:
        generate_population(args.pop_id)
    else:
        parser.error("Provide --pop-id (0-3) or --combine-only")


if __name__ == "__main__":
    main()
