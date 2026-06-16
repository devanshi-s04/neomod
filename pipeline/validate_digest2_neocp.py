#!/usr/bin/env python3
"""
Validate our digest2 implementation against MPC's published NEOCP scores.

MPC publishes digest2 scores (0-100) for every object on the NEOCP.
We have those scores in neocp_data/tables/neocp_objects_history.csv.
We also have hourly ephemerides in neocp_data/ephemerides/.

Approach:
  For each NEOCP object, pick two consecutive ephemeris positions (1 hr apart)
  and construct a synthetic 2-detection MPC-format tracklet, then run our
  local digest2 binary on it.  Compare the resulting P(NEO) scores to MPC's.

Caveat: MPC scored each object on its ACTUAL submitted arc (nobs=4-8 detections,
arc_days=0.02-0.1).  We only supply a synthetic 2-detection tracklet from
ephemerides.  A perfect match is not expected, but high correlation validates
that our digest2 binary and MPC 80-col formatting are correct.
"""

import os
import sys
import subprocess
import tempfile
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
WORKDIR = "/mmfs1/gscratch/dirac/ds2004/sorcha"
NEOCP_DIR = f"{WORKDIR}/neomod/neocp_data"
HISTORY_CSV = f"{NEOCP_DIR}/tables/neocp_objects_history.csv"
EPHEM_DIR = f"{NEOCP_DIR}/ephemerides"

D2_EXEC = f"{WORKDIR}/digest2/digest2"
D2_DIR  = f"{WORKDIR}/digest2"

OUT_CSV = f"{WORKDIR}/outputs/validate_digest2_neocp.csv"
OUT_FIG = f"{WORKDIR}/outputs/validate_digest2_neocp.png"

OBSCODE = "T05"  # MPC code for temporary designation; use X05 (Rubin) or 500 (geocenter)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def utc_to_day_frac(dt_str: str):
    """Convert ISO UTC string to (year, month, day_decimal)."""
    dt = datetime.fromisoformat(dt_str.replace("Z", "+00:00"))
    day_frac = (dt.day
                + dt.hour / 24.0
                + dt.minute / 1440.0
                + dt.second / 86400.0)
    return dt.year, dt.month, day_frac


def format_mpc80(desig12, year, month, day_frac, ra_deg, dec_deg, mag,
                 obscode="X05"):
    """Return exactly 80-character MPC 80-column observation line."""
    ra_h  = ra_deg / 15.0
    rah   = int(ra_h);  ram_f = (ra_h - rah) * 60.0
    ram   = int(ram_f); ras   = min((ram_f - ram) * 60.0, 59.99)
    sign  = "+" if dec_deg >= 0 else "-"
    dec_a = abs(np.clip(dec_deg, -89.99, 89.99))
    decd  = int(dec_a);  decm_f = (dec_a - decd) * 60.0
    decm  = int(decm_f); decs   = min((decm_f - decm) * 60.0, 59.9)
    mag_v = max(0.0, min(99.9, float(mag) if np.isfinite(mag) else 21.0))
    line = (
        f"{desig12}"
        f"  C"
        f"{year:04d} {month:02d} {day_frac:08.5f}"
        f" {rah:02d} {ram:02d} {ras:05.2f} "
        f"{sign}{decd:02d} {decm:02d} {decs:04.1f}"
        f"          "
        f"{mag_v:4.1f} V      "
        f"{obscode:3s}"
    )
    assert len(line) == 80, f"MPC line length {len(line)} != 80: {line!r}"
    return line


def run_digest2(obs_lines: list[str]) -> dict[str, float]:
    """Run digest2 on a list of MPC 80-col obs lines and return {desig: score_0_1}."""
    cfg_text = "noheadings\nnorms\nNEO\n"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".config", delete=False) as f:
        f.write(cfg_text)
        cfg_path = f.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".obs", delete=False) as f:
        f.write("\n".join(obs_lines) + "\n")
        obs_path = f.name

    try:
        result = subprocess.run(
            [D2_EXEC, "-p", D2_DIR, "-c", cfg_path, obs_path],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            print(f"  digest2 stderr: {result.stderr[:500]}")
            return {}
        out = result.stdout.strip().splitlines()
    finally:
        os.unlink(cfg_path)
        os.unlink(obs_path)

    scores = {}
    for line in out:
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                scores[parts[0]] = int(parts[1]) / 100.0
            except ValueError:
                pass
    return scores


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("Loading NEOCP history …")
hist = pd.read_csv(HISTORY_CSV)
print(f"  {len(hist)} objects  columns: {list(hist.columns)}")

# Keep one row per designation (in case of duplicates), prefer highest score
hist = hist.sort_values("score", ascending=False).drop_duplicates("designation")
print(f"  {len(hist)} unique designations after dedup")

print("\nLoading ephemerides …")
ephem_files = sorted(os.listdir(EPHEM_DIR))
ephem_dfs = []
for fn in ephem_files:
    df = pd.read_csv(os.path.join(EPHEM_DIR, fn))
    if len(df) > 0:
        ephem_dfs.append(df)
ephem = pd.concat(ephem_dfs, ignore_index=True)
print(f"  {len(ephem)} total ephemeris rows from {len(ephem_dfs)} files")

# Dedup: multiple snapshot files may include the same (object, ephem_time) pair
ephem = ephem.drop_duplicates(subset=["designation", "ephem_time_utc"])
print(f"  {len(ephem)} rows after dedup on (designation, ephem_time_utc)")
print(f"  Unique designations in ephemerides: {sorted(ephem['designation'].unique())}")

# For each designation, take two consecutive observations (1 hr apart)
obs_pairs = []
missing = []
for desig in hist["designation"]:
    rows = ephem[ephem["designation"] == desig].sort_values("ephem_time_utc").reset_index(drop=True)
    if len(rows) < 2:
        missing.append(desig)
        continue
    r0 = rows.iloc[0]
    r1 = rows.iloc[1]
    obs_pairs.append({
        "designation": desig,
        "time0": r0["ephem_time_utc"],
        "time1": r1["ephem_time_utc"],
        "ra0": r0["ra_deg"],
        "dec0": r0["dec_deg"],
        "ra1": r1["ra_deg"],
        "dec1": r1["dec_deg"],
        "mag": r0["v_mag"],
    })

print(f"  {len(obs_pairs)} objects with ≥2 ephemeris entries")
if missing:
    print(f"  {len(missing)} objects missing ephemerides: {missing}")

# Build MPC 80-col lines
# digest2 truncates designations to 5 chars in output, causing collisions for
# long designations like ZTF10Df/ZTF10Di → both become ZTF10.
# Use synthetic index keys (T000000..T000016) to guarantee unique round-trip.
print("\nFormatting MPC 80-column tracklets …")
mpc_lines = []
desig_order = []
for idx, p in enumerate(obs_pairs):
    key12 = f"VA{idx:03d}".ljust(12)[:12]  # 5-char unique: VA000..VA016
    y0, m0, d0 = utc_to_day_frac(p["time0"])
    y1, m1, d1 = utc_to_day_frac(p["time1"])
    mpc_lines.append(format_mpc80(key12, y0, m0, d0, p["ra0"], p["dec0"], p["mag"]))
    mpc_lines.append(format_mpc80(key12, y1, m1, d1, p["ra1"], p["dec1"], p["mag"]))
    desig_order.append(p["designation"])

print(f"  {len(mpc_lines)//2} tracklets, {len(mpc_lines)} obs lines")
print(f"  Sample obs0: {mpc_lines[0]}")
print(f"  Sample obs1: {mpc_lines[1]}")

# Run digest2
print("\nRunning digest2 …")
scores = run_digest2(mpc_lines)
print(f"  Got {len(scores)} scores")
print(f"  All score keys: {sorted(scores.keys())}")

# Merge results: look up by synthetic key T{idx:06d}, possibly truncated by digest2
results = []
for idx, p in enumerate(obs_pairs):
    desig = p["designation"]
    # Try both full synthetic key and the 5-char truncation digest2 produces
    key_full = f"VA{idx:03d}"
    key_short = key_full[:5]
    our_score = scores.get(key_full, scores.get(key_short, np.nan))
    mpc_row = hist[hist["designation"] == desig].iloc[0]
    results.append({
        "designation": desig,
        "mpc_score_0_1": mpc_row["score_0_1"],
        "mpc_score_int": mpc_row["score"],
        "our_d2_score": our_score,
        "nobs": mpc_row.get("nobs", np.nan),
        "arc_days": mpc_row.get("arc_days", np.nan),
    })

df = pd.DataFrame(results)
print(f"\nResults summary ({len(df)} objects):")
print(f"  NaN our_d2_score: {df['our_d2_score'].isna().sum()}")
valid = df.dropna(subset=["our_d2_score", "mpc_score_0_1"])
print(f"  Valid pairs for comparison: {len(valid)}")

if len(valid) > 0:
    corr = valid[["mpc_score_0_1", "our_d2_score"]].corr().iloc[0, 1]
    spearman = valid[["mpc_score_0_1", "our_d2_score"]].corr(method="spearman").iloc[0, 1]
    mae = np.abs(valid["mpc_score_0_1"] - valid["our_d2_score"]).mean()
    high_mpc = valid[valid["mpc_score_0_1"] >= 0.90]
    high_agree = (high_mpc["our_d2_score"] >= 0.50).sum()
    print(f"  Pearson r:           {corr:.3f}")
    print(f"  Spearman rho:        {spearman:.3f}")
    print(f"  Mean absolute error: {mae:.3f}")
    print(f"  High-MPC (≥0.90) correctly scored ≥0.50 by ours: {high_agree}/{len(high_mpc)}")
    print(f"\n  Note: MPC scored these objects on their ACTUAL arc (nobs={int(valid['nobs'].min())}"
          f"–{int(valid['nobs'].max())}, arc_days={valid['arc_days'].min():.2f}"
          f"–{valid['arc_days'].max():.2f}).")
    print(f"  We used synthetic 2-detection tracklets (1-hr baseline) from ephemerides.")
    print(f"  Imperfect correlation is expected; it validates the binary is functioning.")
    print("\nPer-object comparison:")
    print(valid[["designation", "mpc_score_int", "mpc_score_0_1", "our_d2_score", "nobs", "arc_days"]]
          .sort_values("mpc_score_0_1", ascending=False).to_string(index=False))

# Save CSV
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
df.to_csv(OUT_CSV, index=False)
print(f"\nSaved: {OUT_CSV}")

# Plot
if len(valid) > 0:
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(valid["mpc_score_0_1"], valid["our_d2_score"],
               s=40, alpha=0.7, edgecolors="k", linewidths=0.4)
    ax.plot([0, 1], [0, 1], "r--", lw=1, label="1:1")
    ax.set_xlabel("MPC digest2 score (published)", fontsize=12)
    ax.set_ylabel("Our digest2 score (2-det synthetic tracklet)", fontsize=12)
    ax.set_title(f"NEOCP digest2 validation  |  n={len(valid)}  |  r={corr:.2f}  |  MAE={mae:.2f}", fontsize=11)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=150)
    print(f"Saved: {OUT_FIG}")
    plt.close(fig)
