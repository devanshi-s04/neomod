#!/usr/bin/env python3
"""Digest2 vs VDP (GMM pipeline) ROC comparison — pure-S3M single-epoch benchmark.

GMM-pipeline counterpart of run_digest2_comparison.py. It treats a pre-built GMM
probability-map file as a black box: it propagates the S3M populations to the map's
observation time, filters to the map's circular ecliptic sky patch, scores each object
with the GMM VDP maps via score_orbital_df, then builds X05 30-minute MPC tracklets and
scores the same objects with digest2. Saves a side-by-side ROC comparison.

This makes the single-epoch S3M benchmark methodologically identical to the Sorcha run
(same GMM cloner, same +/-2.0 deg/day velocity grid, same support-masked maps), so the
only difference between the two is Sorcha cadence vs. single-epoch census.

Run (on the Mac, where digest2 lives):
    python3 run_digest2_comparison_gmm.py \
        --prob-maps ../prob_maps_grid_s3m/prob_maps_grid_dlon+000_lat+00.npz \
        --d2-dir /path/to/mpcdev-digest2

By default the map's own obstime/center/grid are used, so you do not need to know them.
"""

import os
import argparse
import subprocess
import sys
import tempfile
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "adam_core_stub"))

# Score the pure-S3M census via the raw .s3m files (not the hybrid_elements.parquet).
# Must be set BEFORE importing the pipeline, which reads VDP_LOADER at import time.
os.environ.setdefault("VDP_LOADER", "s3m")
import velocity_density_pipeline_gmm as vdp     # GMM pipeline (the only change vs. the K|M script)

_SRC  = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_SRC)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Digest2 vs VDP (GMM) ROC comparison.")
parser.add_argument(
    "--prob-maps",
    default=os.path.join(_ROOT, "prob_maps_grid_s3m", "prob_maps_grid_dlon+000_lat+00.npz"),
    help="Path to the GMM probability-map .npz (the black box). "
         "Default: the pure-S3M antisun grid map.",
)
parser.add_argument(
    "--obstime-str",
    default=None,
    help="Observation time. Default: read from the map file (recommended).",
)
parser.add_argument(
    "--run-label",
    default="s3m_gmm",
    help="Label used in the output filenames.",
)
parser.add_argument(
    "--d2-dir",
    default="/Users/devanshisingh/Downloads/research/SolSys/digest3/mpcdev-digest2-278a31e734e4",
    help="Directory of the digest2 checkout (containing the 'digest2' binary + model files).",
)
args = parser.parse_args()

RUN_LABEL     = args.run_label
DT_DAYS       = 30.0 / 1440.0          # 30-min tracklet baseline
OBSCODE       = "X05"                   # Vera C. Rubin Observatory
DIGEST2_CHUNK_TRACKLETS = 5_000
DIGEST2_TIMEOUT_SEC     = 1_800

D2_DIR  = args.d2_dir
D2_EXEC = os.path.join(D2_DIR, "digest2")

PROB_MAPS_PATH = args.prob_maps
OUT_FIG        = os.path.join(_ROOT, f"roc_comparison_vdp_digest2_{RUN_LABEL}.png")
OUT_PARQUET    = os.path.join(_ROOT, f"s3m_digest2_comparison_{RUN_LABEL}.parquet")
OUT_VDP_INPUT  = os.path.join(_ROOT, f"s3m_digest2_comparison_vdp_input_{RUN_LABEL}.parquet")

# Populations: label -> (s3m_pop_name, max_objects)  -- same subsampling as the paper benchmark
POP_SETTINGS = {
    "NEO":     ("neo",    None),        # ~268 k
    "MBA":     ("mba",    200_000),     # subsample for speed
    "TNO":     ("tno",    None),        # ~48 k
    "Trojans": ("trojan", 100_000),     # subsample to avoid slow chunks
}


def format_mpc80(desig12, year, month, day_frac, ra_deg, dec_deg, mag):
    """Return exactly 80-character MPC 80-column observation line."""
    ra_h   = ra_deg / 15.0
    rah    = int(ra_h);  ram_f = (ra_h - rah) * 60.0
    ram    = int(ram_f); ras   = min((ram_f - ram) * 60.0, 59.99)
    sign   = "+" if dec_deg >= 0 else "-"
    dec_a  = abs(np.clip(dec_deg, -89.99, 89.99))
    decd   = int(dec_a);  decm_f = (dec_a - decd) * 60.0
    decm   = int(decm_f); decs   = min((decm_f - decm) * 60.0, 59.9)
    mag_v  = max(0.0, min(99.9, float(mag) if np.isfinite(mag) else 21.0))
    line = (
        f"{desig12}"                                     # [0-11]  12 chars
        f"  C"                                           # [12-14]  3 chars
        f"{year:04d} {month:02d} {day_frac:08.5f}"      # [15-30] 16 chars
        f" {rah:02d} {ram:02d} {ras:05.2f} "            # [31-43] 13 chars
        f"{sign}{decd:02d} {decm:02d} {decs:04.1f}"     # [44-54] 11 chars
        f"          "                                    # [55-64] 10 chars
        f"{mag_v:4.1f} V      "                         # [65-76] 12 chars
        f"{OBSCODE:3s}"                                  # [77-79]  3 chars
    )
    assert len(line) == 80
    return line


def run_digest2_chunk(obs_lines, cfg_text, chunk_tracklets):
    """Run digest2 on MPC 80-column lines in chunks and return stdout lines."""
    if len(obs_lines) % 2 != 0:
        raise ValueError("Expected exactly two observation lines per tracklet.")

    n_tracklets = len(obs_lines) // 2
    all_out_lines = []

    with tempfile.NamedTemporaryFile(mode="w", suffix=".config", delete=False) as cfg:
        cfg.write(cfg_text)
        cfg_path = cfg.name

    try:
        for start in range(0, n_tracklets, chunk_tracklets):
            stop = min(start + chunk_tracklets, n_tracklets)
            line_start = 2 * start
            line_stop = 2 * stop

            with tempfile.NamedTemporaryFile(mode="w", suffix=".obs", delete=False) as tmp:
                tmp.write("\n".join(obs_lines[line_start:line_stop]) + "\n")
                obs_path = tmp.name

            print(
                f"  digest2 chunk {start // chunk_tracklets + 1}/"
                f"{int(np.ceil(n_tracklets / chunk_tracklets))}: "
                f"{stop - start:,} tracklets",
                flush=True,
            )

            try:
                result = subprocess.run(
                    [D2_EXEC, "-p", D2_DIR, "-c", cfg_path, obs_path],
                    capture_output=True,
                    text=True,
                    timeout=DIGEST2_TIMEOUT_SEC,
                )
            finally:
                os.unlink(obs_path)

            if result.returncode != 0:
                raise RuntimeError(
                    f"digest2 failed on tracklets {start:,}:{stop:,} "
                    f"with exit code {result.returncode}\n"
                    f"stderr:\n{result.stderr[:1000]}"
                )

            if result.stderr.strip():
                print(f"    stderr: {result.stderr[:300]}", flush=True)

            out = result.stdout.splitlines()
            print(f"    output lines: {len(out):,}", flush=True)
            all_out_lines.extend(out)
    finally:
        os.unlink(cfg_path)

    return all_out_lines


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# Load GMM VDP probability maps (the black box)
print("Loading GMM VDP probability maps …")
print(f"  map file: {PROB_MAPS_PATH}")
prob_map_set = vdp.ProbMapSet.from_npz(PROB_MAPS_PATH, support_mask_min=1)
print(f"  populations: {prob_map_set.population_names}")
print(f"  map obstime: {prob_map_set.obstime_str}")
print(f"  map center:  ({prob_map_set.center_lon_deg:.1f}, {prob_map_set.center_lat_deg:.1f})")

# Use the map's own obstime unless explicitly overridden.
OBSTIME_STR = args.obstime_str or prob_map_set.obstime_str
if args.obstime_str and args.obstime_str != prob_map_set.obstime_str:
    raise ValueError(
        f"--obstime-str={args.obstime_str!r} does not match map obstime "
        f"{prob_map_set.obstime_str!r}. Omit --obstime-str to use the map's value."
    )
print(f"  using obstime: {OBSTIME_STR}")

_OBS_DT = datetime.fromisoformat(OBSTIME_STR)
YEAR0, MONTH0 = _OBS_DT.year, _OBS_DT.month
DAY0 = (
    _OBS_DT.day
    + _OBS_DT.hour / 24.0
    + _OBS_DT.minute / 1440.0
    + (_OBS_DT.second + _OBS_DT.microsecond / 1e6) / 86400.0
)
DAY1 = DAY0 + DT_DAYS

rng = np.random.default_rng(42)
all_frames = []

for pop_label, (s3m_pop, max_obj) in POP_SETTINGS.items():
    print(f"\n=== {pop_label} ({s3m_pop}) ===")
    df, scorer = vdp.load_s3m_population(s3m_pop)
    print(f"  Loaded {len(df):,} objects")

    if max_obj is not None and len(df) > max_obj:
        idx = rng.choice(len(df), size=max_obj, replace=False)
        df  = df.iloc[idx].reset_index(drop=True)
        print(f"  Subsampled to {len(df):,}")

    _, pop_df = prob_map_set.score_orbital_df(
        df=df,
        scorer=scorer,
        obstime_str=OBSTIME_STR,
        max_sep_deg=prob_map_set.max_sep_deg,
        chunk=50_000,
        show_progress=False,
        return_visible=True,
    )

    pop_df["true_population"] = pop_label
    if len(pop_df) > 0:
        print(f"  → {len(pop_df):,} total in map sky patch")
        all_frames.append(pop_df)

data = pd.concat(all_frames, ignore_index=True)
print(f"\nCombined dataset: {len(data):,} objects")
print(data["true_population"].value_counts().to_string())

# Filter to VDP map's operating region:
#   - within max_sep_deg of the map's ecliptic sky center
#   - within the map's velocity and magnitude grid bounds
x_min, x_max = float(prob_map_set.x_grid.min()), float(prob_map_set.x_grid.max())
y_min, y_max = float(prob_map_set.y_grid.min()), float(prob_map_set.y_grid.max())
mag_min = min(b["mag_min"] for b in prob_map_set.mag_bins)
mag_max = max(b["mag_max"] for b in prob_map_set.mag_bins)

in_region = (
    (data["sky_sep_deg"].to_numpy(float) <= prob_map_set.max_sep_deg)
    & (data["mag_app"].to_numpy(float) >= mag_min)
    & (data["mag_app"].to_numpy(float) <  mag_max)
    & (data["vlam"].to_numpy(float)  >= x_min) & (data["vlam"].to_numpy(float)  <= x_max)
    & (data["vbeta"].to_numpy(float) >= y_min) & (data["vbeta"].to_numpy(float) <= y_max)
)
data = data[in_region].reset_index(drop=True)
print(f"\nFiltered to VDP operating region (center=({prob_map_set.center_lon_deg:.0f}°, "
      f"{prob_map_set.center_lat_deg:.0f}°), radius={prob_map_set.max_sep_deg:.0f}°, "
      f"mag {mag_min:.0f}–{mag_max:.0f}): {len(data):,} objects")
print(data["true_population"].value_counts().to_string())

# VDP scoring was already done by ProbMapSet.score_orbital_df above.
data["P_NEO_vdp"] = data["P_NEO"]
print(f"  P_NEO range: [{data['P_NEO_vdp'].min():.4f}, {data['P_NEO_vdp'].max():.4f}]")
data.to_parquet(OUT_VDP_INPUT, index=False)
print(f"  Saved VDP-scored input checkpoint: {OUT_VDP_INPUT}")

# Build MPC 80-column tracklets
print("\nBuilding MPC tracklets …")
mpc_lines = []
for i in range(len(data)):
    row    = data.iloc[i]
    desig  = f"     D{i:06d}"           # 5 blanks + "D" + 6-digit int = 12 chars

    ra0, dec0 = float(row["ra_deg"]),  float(row["dec_deg"])
    dra, ddec = float(row["dra_deg_day"]), float(row["ddec_deg_day"])
    mag       = float(row["mag_app"])

    ra1  = (ra0 + dra * DT_DAYS) % 360.0
    dec1 = float(np.clip(dec0 + ddec * DT_DAYS, -89.99, 89.99))

    mpc_lines.append(format_mpc80(desig, YEAR0, MONTH0, DAY0, ra0, dec0, mag))
    mpc_lines.append(format_mpc80(desig, YEAR0, MONTH0, DAY1, ra1, dec1, mag))

print(f"  {len(mpc_lines)//2:,} tracklets, {len(mpc_lines):,} obs lines")
print(f"  Sample:\n    {mpc_lines[0]}\n    {mpc_lines[1]}")

print(f"\nRunning digest2 on {len(mpc_lines)//2:,} tracklets …", flush=True)
out_lines = run_digest2_chunk(
    obs_lines=mpc_lines,
    cfg_text="noheadings\nnorms\nNEO\n",
    chunk_tracklets=DIGEST2_CHUNK_TRACKLETS,
)
print(f"  Output lines: {len(out_lines)}")
if out_lines:
    print(f"  First 3: {out_lines[:3]}")

# Parse digest2 output
print("\nParsing digest2 output …")
d2_map = {}
duplicate_ids = []
unparsed_lines = []
for line in out_lines:
    parts = line.strip().split()
    if len(parts) >= 2:
        try:
            if parts[0] in d2_map:
                duplicate_ids.append(parts[0])
            d2_map[parts[0]] = int(parts[1]) / 100.0
        except ValueError:
            unparsed_lines.append(line)
    else:
        unparsed_lines.append(line)
print(f"  Parsed {len(d2_map):,} scores")
if duplicate_ids:
    print(f"  Warning: {len(duplicate_ids):,} duplicate digest2 IDs; first few: {duplicate_ids[:5]}")
if unparsed_lines:
    print(f"  Warning: {len(unparsed_lines):,} unparsed digest2 lines; first few: {unparsed_lines[:5]}")
missing_ids = [f"D{i:06d}" for i in range(len(data)) if f"D{i:06d}" not in d2_map]
if missing_ids:
    print(f"  Warning: {len(missing_ids):,} missing digest2 IDs; first few: {missing_ids[:5]}")
data["P_NEO_d2"] = [d2_map.get(f"D{i:06d}", 0.0) for i in range(len(data))]

# ROC curves
print("\nComputing ROC curves …")
is_neo = (data["true_population"] == "NEO").to_numpy()
N_neo  = is_neo.sum()
N_non  = (~is_neo).sum()
print(f"  N_NEO={N_neo:,}  N_non-NEO={N_non:,}")

thresholds = np.linspace(0.0, 1.0, 201)

def compute_roc(scores, is_neo, thresholds):
    N = is_neo.sum()
    comp = np.empty(len(thresholds))
    cont = np.empty(len(thresholds))
    for j, t in enumerate(thresholds):
        above = scores > t
        tp    = above[is_neo].sum()
        fp    = above[~is_neo].sum()
        comp[j] = tp / max(N, 1)
        cont[j] = fp / max(tp + fp, 1)
    return comp, cont

comp_vdp, cont_vdp = compute_roc(data["P_NEO_vdp"].to_numpy(), is_neo, thresholds)
comp_d2,  cont_d2  = compute_roc(data["P_NEO_d2"].to_numpy(),  is_neo, thresholds)

def best_f1(comp, cont, thresholds):
    purity = 1.0 - cont
    f1 = np.where(purity + comp > 0,
                  2 * purity * comp / (purity + comp + 1e-12), 0.0)
    i = np.argmax(f1)
    return thresholds[i], comp[i], cont[i], f1[i]

t_vdp, c_vdp, n_vdp, f_vdp = best_f1(comp_vdp, cont_vdp, thresholds)
t_d2,  c_d2,  n_d2,  f_d2  = best_f1(comp_d2,  cont_d2,  thresholds)

print(f"\n  VDP     optimal: t={t_vdp:.3f}  completeness={c_vdp:.1%}  contamination={n_vdp:.1%}  F1={f_vdp:.3f}")
print(f"  Digest2 optimal: t={t_d2:.3f}  completeness={c_d2:.1%}  contamination={n_d2:.1%}  F1={f_d2:.3f}")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

ax = axes[0]
ax.plot(comp_vdp * 100, cont_vdp * 100, lw=2.2, color="tab:blue",
        label="VDP (GMM, this work)")
ax.plot(comp_d2  * 100, cont_d2  * 100, lw=2.2, color="tab:orange",
        linestyle="--", label="digest2")
ax.scatter([c_vdp*100], [n_vdp*100], color="tab:blue",   s=100, zorder=5,
           label=f"VDP best (F1={f_vdp:.2f})")
ax.scatter([c_d2*100],  [n_d2*100],  color="tab:orange", s=100, zorder=5,
           marker="s", label=f"d2 best (F1={f_d2:.2f})")
ax.set_xlabel("Completeness (%)", fontsize=12)
ax.set_ylabel("Contamination (%)", fontsize=12)
ax.set_title(
    f"ROC Comparison: VDP (GMM) vs digest2\n"
    f"(pure S3M, ecliptic center=({prob_map_set.center_lon_deg:.0f}°, "
    f"{prob_map_set.center_lat_deg:.0f}°), radius={prob_map_set.max_sep_deg:.0f}°, "
    f"mag {mag_min:.0f}–{mag_max:.0f}, {OBSTIME_STR})\n"
    f"N_NEO={N_neo:,}, N_non={N_non:,}",
    fontsize=10,
)
ax.legend(fontsize=10); ax.grid(alpha=0.3)
ax.set_xlim(0, 100); ax.set_ylim(0, 100)

ax = axes[1]
ax.plot(thresholds, comp_vdp*100, lw=2, color="tab:blue",    label="VDP completeness")
ax.plot(thresholds, cont_vdp*100, lw=2, color="tab:blue",   linestyle=":", label="VDP contamination")
ax.plot(thresholds, comp_d2 *100, lw=2, color="tab:orange", linestyle="--", label="d2 completeness")
ax.plot(thresholds, cont_d2 *100, lw=2, color="tab:orange", linestyle="-.", label="d2 contamination")
ax.axvline(t_vdp, color="tab:blue",   alpha=0.4, lw=1.2)
ax.axvline(t_d2,  color="tab:orange", alpha=0.4, lw=1.2)
ax.set_xlabel("Score threshold", fontsize=12)
ax.set_ylabel("Completeness / Contamination (%)", fontsize=12)
ax.set_title("Score threshold vs. Performance", fontsize=11)
ax.legend(fontsize=9); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(OUT_FIG, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved ROC comparison: {OUT_FIG}")

# Summary
print("\n=== Summary ===")
print(f"  Dataset: {len(data):,} objects  "
      f"(map sky patch radius={prob_map_set.max_sep_deg:.1f}°, {OBSTIME_STR})")
print(f"  True NEO:     {N_neo:,} ({100*N_neo/len(data):.1f}%)")
print(f"  True non-NEO: {N_non:,} ({100*N_non/len(data):.1f}%)")
print()
print(f"  VDP    best F1={f_vdp:.3f}: completeness={c_vdp:.1%}, purity={1-n_vdp:.1%}")
print(f"  Digest2 best F1={f_d2:.3f}: completeness={c_d2:.1%}, purity={1-n_d2:.1%}")

data[["true_population", "lam_deg", "vlam", "vbeta", "mag_app", "H",
      "P_NEO_vdp", "P_NEO_d2"]].to_parquet(OUT_PARQUET, index=False)
print(f"\nSaved scored data: {OUT_PARQUET}")
