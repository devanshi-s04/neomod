#!/usr/bin/env python
"""
Extract a single-epoch antisun patch parquet for fair three-way ROC comparison.

Filters to: tracklets observed in May 2026 AND within 30° of the May 9 antisun
(ecliptic lon=228°, lat=0°). Joins S3M kNN scores with GMM scores so all three
classifiers (S3M kNN, GMM, digest2) appear on the same rows.

Output: outputs/phase2/wagg_sorcha_may2026_antisun_patch.parquet

Why NOT 100 sq deg (5.6° radius):
  A strict 5.6° radius around the antisun yields only 83 objects (6 NEOs) from the
  subsampled parquet — too few for a ROC curve. A 30° radius gives 14,873 objects
  (2,402 NEOs, 16.1% rate), which matches the full-comparison NEO rate and provides
  good ROC statistics. The "single epoch" property is preserved: all objects are
  observed in May 2026 and scored with the May 2026 antisun VDP maps.
"""

import sys
import numpy as np
import pandas as pd
from astropy.time import Time

WORKDIR = "/mmfs1/gscratch/astro/ds2004/sorcha"
S3M_PARQUET = f"{WORKDIR}/outputs/phase2/wagg_sorcha_comparison.parquet"
GMM_PARQUET = f"{WORKDIR}/outputs/phase2_gmm/wagg_sorcha_comparison_gmm.parquet"
OUT_PARQUET = f"{WORKDIR}/outputs/phase2/wagg_sorcha_may2026_antisun_patch.parquet"

ANTISUN_LON = 228.02   # May 9 2026 antisun ecliptic longitude (deg)
ANTISUN_LAT = 0.0
PATCH_RADIUS_DEG = 30.0

mjd_may1 = Time("2026-05-01T00:00:00", scale="utc").mjd
mjd_jun1 = Time("2026-06-01T00:00:00", scale="utc").mjd


def ang_sep_deg(lon1, lat1, lon2, lat2):
    l1, b1 = np.deg2rad(lon1), np.deg2rad(lat1)
    l2, b2 = np.deg2rad(lon2), np.deg2rad(lat2)
    cos_sep = np.sin(b1)*np.sin(b2) + np.cos(b1)*np.cos(b2)*np.cos(l1 - l2)
    return np.rad2deg(np.arccos(np.clip(cos_sep, -1.0, 1.0)))


print("Loading S3M parquet (needed columns only)...")
s3m_cols = [
    "tracklet_id", "population", "P_NEO_vdp", "P_NEO_d2",
    "vlam", "vbeta", "mean_mag", "ecl_lon", "ecl_lat", "mjd0_utc",
]
s3m = pd.read_parquet(S3M_PARQUET, columns=s3m_cols)
print(f"  S3M loaded: {len(s3m):,} rows")

print("Loading GMM parquet (P_NEO_vdp only)...")
gmm_cols = ["tracklet_id", "P_NEO_vdp"]
gmm = pd.read_parquet(GMM_PARQUET, columns=gmm_cols)
gmm = gmm.rename(columns={"P_NEO_vdp": "P_NEO_gmm"})
print(f"  GMM loaded: {len(gmm):,} rows")

# --- date filter: observed in May 2026 ---
date_mask = (s3m["mjd0_utc"] >= mjd_may1) & (s3m["mjd0_utc"] < mjd_jun1)
s3m_may = s3m[date_mask].copy()
print(f"\nMay 2026 observations: {len(s3m_may):,}")

# --- spatial filter: within 30° of May 9 antisun ---
sep = ang_sep_deg(
    s3m_may["ecl_lon"].values, s3m_may["ecl_lat"].values,
    ANTISUN_LON, ANTISUN_LAT,
)
patch_mask = sep <= PATCH_RADIUS_DEG
s3m_patch = s3m_may[patch_mask].copy()
print(f"Within {PATCH_RADIUS_DEG}° of antisun: {len(s3m_patch):,}")
print(f"  NEOs: {(s3m_patch.population == 'NEO').sum():,}")
print(f"  Population:\n{s3m_patch.population.value_counts().to_string()}")

# --- merge GMM scores ---
out = s3m_patch.merge(gmm, on="tracklet_id", how="left")
out = out.rename(columns={"P_NEO_vdp": "P_NEO_s3m"})

# --- report ---
print(f"\nFinal patch: {len(out):,} tracklets")
print(f"  NEOs: {(out.population == 'NEO').sum():,}  ({(out.population == 'NEO').mean():.1%} NEO rate)")
print(f"  P_NEO_s3m  range: [{out.P_NEO_s3m.min():.4f}, {out.P_NEO_s3m.max():.4f}]  NaN: {out.P_NEO_s3m.isna().sum()}")
print(f"  P_NEO_gmm  range: [{out.P_NEO_gmm.min():.4f}, {out.P_NEO_gmm.max():.4f}]  NaN: {out.P_NEO_gmm.isna().sum()}")
print(f"  P_NEO_d2   range: [{out.P_NEO_d2.min():.4f},  {out.P_NEO_d2.max():.4f}]   NaN: {out.P_NEO_d2.isna().sum()}")

out.to_parquet(OUT_PARQUET, index=False)
print(f"\nSaved: {OUT_PARQUET}")
