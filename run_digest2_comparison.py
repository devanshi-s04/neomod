#!/usr/bin/env python3
"""Digest2 vs VDP ROC curve comparison.

Propagates S3M objects to 2025-03-21, filters to the ecliptic strip
(|beta| < 3 deg), then scores each object with both the VDP probability
maps and the digest2 classifier.  Saves a side-by-side ROC comparison
figure.

Run:
    /path/to/venv/python3 run_digest2_comparison.py
"""

import os
import subprocess
import sys
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import neoscore as nsc
import velocity_density_pipeline as vdp

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OBSTIME_STR   = "2025-03-21T00:00:00"
BETA_MAX_DEG  = 3.0
DT_DAYS       = 30.0 / 1440.0          # 30-min tracklet baseline
OBSCODE       = "568"                   # Mauna Kea

D2_DIR  = "/Users/devanshisingh/Downloads/research/SolSys/digest3/mpcdev-digest2-278a31e734e4"
D2_EXEC = os.path.join(D2_DIR, "digest2")

_SELF = os.path.dirname(os.path.abspath(__file__))
PROB_MAPS_PATH = os.path.join(_SELF, "prob_maps_2025-03-21.npz")
OUT_FIG        = os.path.join(_SELF, "roc_comparison_vdp_digest2.png")
OUT_PARQUET    = os.path.join(_SELF, "s3m_digest2_comparison.parquet")

# Populations: label -> (s3m_pop_name, max_objects)
POP_SETTINGS = {
    "NEO":     ("neo",    None),        # ~268 k
    "MBA":     ("mba",    200_000),     # subsample for speed
    "TNO":     ("tno",    None),        # ~48 k
    "Trojans": ("trojan", 100_000),     # subsample to avoid slow chunks
}

YEAR0, MONTH0 = 2025, 3
DAY0 = 21.00000
DAY1 = DAY0 + DT_DAYS                  # 21.02083…

ECLIPTIC_OBL_DEG = 23.439291
AU_KM = 149_597_870.7

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def propagate_strip(df, scorer, obstime_str, beta_max_deg, chunk_size=50_000,
                    pop_label="?", verbose=True):
    """Propagate orbital elements, apply ecliptic-strip filter, return DataFrame.

    Avoids the SkyCoord.separation() call in build_visible_subset_dataframe
    (that call is only needed for the circular sky-patch cut, which we don't
    want here).  Instead we compute ecliptic beta directly and filter on that.

    Returns a DataFrame with columns:
        ra_deg, dec_deg, dra_deg_day, ddec_deg_day,
        lam_deg, beta_deg, vlam, vbeta, mag_app (all float64)
        + all original columns from df (a, e, i, node, argperi, t_p, H …)
    """
    eps   = np.deg2rad(ECLIPTIC_OBL_DEG)
    ceps  = np.cos(eps)
    seps  = np.sin(eps)

    _, rE, vE, r_tele = scorer._get_earth_and_observer(obstime_str)

    n_total  = len(df)
    n_chunks = int(np.ceil(n_total / chunk_size))
    chunks_out = []

    for ci, start in enumerate(range(0, n_total, chunk_size)):
        cdf = df.iloc[start : start + chunk_size].reset_index(drop=True)
        n_c = len(cdf)
        if verbose:
            print(f"  Chunk {ci+1}/{n_chunks} ({n_c:,}) … ", end="", flush=True)

        try:
            r_ecl, v_ecl = nsc.elements_to_helio_ecliptic_state(
                a_AU    = cdf["a"].to_numpy(float),
                e       = cdf["e"].to_numpy(float),
                inc_deg = cdf["i"].to_numpy(float),
                raan_deg= cdf["node"].to_numpy(float),
                argp_deg= cdf["argperi"].to_numpy(float),
                tp_mjd  = cdf["t_p"].to_numpy(float),
                obstime_str=obstime_str,
                method="newton", n_iter=10,
                chunk=min(chunk_size, n_c),
                show_progress=False,
            )
        except Exception as exc:
            print(f"SKIP (propagation error: {exc})")
            continue

        # r_ecl, v_ecl are in km, km/s — keep as-is; rE/vE/r_tele are also km/km/s

        # heliocentric equatorial (ecliptic → equatorial rotation)
        r_eq = np.column_stack([
            r_ecl[:, 0],
            ceps * r_ecl[:, 1] - seps * r_ecl[:, 2],
            seps * r_ecl[:, 1] + ceps * r_ecl[:, 2],
        ])
        v_eq = np.column_stack([
            v_ecl[:, 0],
            ceps * v_ecl[:, 1] - seps * v_ecl[:, 2],
            seps * v_ecl[:, 1] + ceps * v_ecl[:, 2],
        ])

        # geocentric relative position/velocity (km, km/s)
        r_rel = r_eq - (rE + r_tele)
        v_rel = v_eq - vE

        x, y, z   = r_rel[:, 0], r_rel[:, 1], r_rel[:, 2]
        vx, vy, vz = v_rel[:, 0], v_rel[:, 1], v_rel[:, 2]
        rho2 = x*x + y*y
        r2   = rho2 + z*z
        rho  = np.sqrt(rho2)

        good = (
            np.isfinite(rho) & (rho > 0) &
            np.isfinite(r2)  & (r2 > 0)
        )

        # geocentric ecliptic (for beta filter + vlam/vbeta)
        xe = r_rel[:, 0]
        ye =  ceps * r_rel[:, 1] + seps * r_rel[:, 2]
        ze = -seps * r_rel[:, 1] + ceps * r_rel[:, 2]
        vxe = v_rel[:, 0]
        vye =  ceps * v_rel[:, 1] + seps * v_rel[:, 2]
        vze = -seps * v_rel[:, 1] + ceps * v_rel[:, 2]

        rho2e = xe*xe + ye*ye
        r2e   = rho2e + ze*ze
        rhoe  = np.sqrt(rho2e)

        beta_rad = np.arctan2(ze, rhoe)
        lam_rad  = np.arctan2(ye, xe)
        beta_deg_arr = np.rad2deg(beta_rad)
        lam_deg_arr  = np.rad2deg(lam_rad) % 360.0

        dlam  = (xe*vye - ye*vxe) / np.maximum(rho2e, 1e-30)
        dbeta = (vze*rho2e - ze*(xe*vxe + ye*vye)) / np.maximum(r2e * rhoe, 1e-30)
        vlam  = np.rad2deg(dlam)  * 86400.0
        vbeta = np.rad2deg(dbeta) * 86400.0

        # equatorial RA/Dec and rates
        ra_rad  = np.arctan2(y, x)
        dec_rad = np.arctan2(z, rho)
        dra     = (x*vy - y*vx) / np.maximum(rho2, 1e-30)
        ddec    = (vz*rho2 - z*(x*vx + y*vy)) / np.maximum(r2 * rho, 1e-30)
        ra_deg  = np.rad2deg(ra_rad) % 360.0
        dec_deg = np.rad2deg(dec_rad)
        dra_deg_day  = np.rad2deg(dra)  * 86400.0
        ddec_deg_day = np.rad2deg(ddec) * 86400.0

        # ecliptic-strip filter
        m = good & (np.abs(beta_deg_arr) < beta_max_deg)
        n_kept = m.sum()
        if verbose:
            print(f"{n_kept:,} kept in |β|<{beta_max_deg}°")

        if n_kept == 0:
            continue

        out = cdf[m].copy().reset_index(drop=True)
        out["ra_deg"]       = ra_deg[m]
        out["dec_deg"]      = dec_deg[m]
        out["dra_deg_day"]  = dra_deg_day[m]
        out["ddec_deg_day"] = ddec_deg_day[m]
        out["lam_deg"]      = lam_deg_arr[m]
        out["beta_deg"]     = beta_deg_arr[m]
        out["vlam"]         = vlam[m]
        out["vbeta"]        = vbeta[m]

        # apparent magnitude (uses H from orbital elements)
        r_sun  = np.linalg.norm(r_eq[m], axis=1) / AU_KM   # km → AU
        Delta  = np.linalg.norm(r_rel[m], axis=1) / AU_KM  # km → AU
        u_sun  = -r_eq[m]  / np.linalg.norm(r_eq[m],  axis=1, keepdims=True)
        u_obs  = -r_rel[m] / np.linalg.norm(r_rel[m], axis=1, keepdims=True)
        cos_a  = np.clip(np.sum(u_sun * u_obs, axis=1), -1, 1)
        alpha  = np.arccos(cos_a)
        G = 0.15
        phi1 = np.exp(-3.33 * np.tan(0.5 * alpha)**0.63)
        phi2 = np.exp(-1.87 * np.tan(0.5 * alpha)**1.22)
        Phi  = np.maximum((1 - G) * phi1 + G * phi2, 1e-30)
        H_arr = out["H"].to_numpy(float)
        mag_app = H_arr + 5.0 * np.log10(r_sun * Delta) - 2.5 * np.log10(Phi)
        out["mag_app"]         = mag_app
        out["true_population"] = pop_label

        chunks_out.append(out)

    if not chunks_out:
        return pd.DataFrame()
    return pd.concat(chunks_out, ignore_index=True)


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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# Load VDP probability maps
print("Loading VDP probability maps …")
prob_map_set = vdp.ProbMapSet.from_npz(PROB_MAPS_PATH)
print(f"  populations: {prob_map_set.population_names}")

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

    pop_df = propagate_strip(df, scorer, OBSTIME_STR, BETA_MAX_DEG,
                             pop_label=pop_label, verbose=True)
    if len(pop_df) > 0:
        print(f"  → {len(pop_df):,} total in strip")
        all_frames.append(pop_df)

data = pd.concat(all_frames, ignore_index=True)
print(f"\nCombined dataset: {len(data):,} objects")
print(data["true_population"].value_counts().to_string())

# Filter to VDP map's operating region:
#   - within max_sep_deg (30°) of opposition (geocentric ecliptic lon=180°)
#   - within the map's velocity and magnitude grid bounds
#   - not masked as stationary (sky velocity >= mask_radius_deg_per_day)
lam    = data["lam_deg"].to_numpy(float)
v_sky  = np.sqrt(data["vlam"].to_numpy(float)**2 + data["vbeta"].to_numpy(float)**2)
x_min, x_max = float(prob_map_set.x_grid.min()), float(prob_map_set.x_grid.max())
y_min, y_max = float(prob_map_set.y_grid.min()), float(prob_map_set.y_grid.max())
mag_min = min(b["mag_min"] for b in prob_map_set.mag_bins)
mag_max = max(b["mag_max"] for b in prob_map_set.mag_bins)

in_region = (
    (np.abs(lam - 180.0) < prob_map_set.max_sep_deg)
    & (data["mag_app"].to_numpy(float) >= mag_min)
    & (data["mag_app"].to_numpy(float) <  mag_max)
    & (data["vlam"].to_numpy(float)  >= x_min) & (data["vlam"].to_numpy(float)  <= x_max)
    & (data["vbeta"].to_numpy(float) >= y_min) & (data["vbeta"].to_numpy(float) <= y_max)
    & (v_sky >= prob_map_set.mask_radius_deg_per_day)
)
data = data[in_region].reset_index(drop=True)
print(f"\nFiltered to VDP operating region (opposition ±{prob_map_set.max_sep_deg}°, "
      f"mag {mag_min:.0f}–{mag_max:.0f}): {len(data):,} objects")
print(data["true_population"].value_counts().to_string())

# VDP scoring
print("\nScoring with VDP …")
probs = prob_map_set.score_visible(
    vlam    = data["vlam"].to_numpy(float),
    vbeta   = data["vbeta"].to_numpy(float),
    mag_app = data["mag_app"].to_numpy(float),
)
data["P_NEO_vdp"] = probs["NEO"]
print(f"  P_NEO range: [{data['P_NEO_vdp'].min():.4f}, {data['P_NEO_vdp'].max():.4f}]")

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

# Write digest2 config
cfg_path = os.path.join(D2_DIR, "d2_roc_comparison.config")
with open(cfg_path, "w") as fh:
    fh.write("noheadings\nnorms\nNEO\n")

# Run digest2
with tempfile.NamedTemporaryFile(mode="w", suffix=".obs", delete=False) as tmp:
    tmp.write("\n".join(mpc_lines) + "\n")
    obs_path = tmp.name

print(f"\nRunning digest2 on {len(mpc_lines)//2:,} tracklets …")
result = subprocess.run(
    [D2_EXEC, "-p", D2_DIR, "-c", cfg_path, obs_path],
    capture_output=True, text=True, timeout=3600,
)
os.unlink(obs_path)
print(f"  Exit code: {result.returncode}")
if result.stderr.strip():
    print(f"  stderr: {result.stderr[:300]}")

out_lines = result.stdout.splitlines()
print(f"  Output lines: {len(out_lines)}")
if out_lines:
    print(f"  First 3: {out_lines[:3]}")

# Parse digest2 output
print("\nParsing digest2 output …")
d2_map = {}
for line in out_lines:
    parts = line.strip().split()
    if len(parts) >= 2:
        try:
            d2_map[parts[0]] = int(parts[1]) / 100.0
        except ValueError:
            pass
print(f"  Parsed {len(d2_map):,} scores")
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
        label="VDP (this work)")
ax.plot(comp_d2  * 100, cont_d2  * 100, lw=2.2, color="tab:orange",
        linestyle="--", label="digest2")
ax.scatter([c_vdp*100], [n_vdp*100], color="tab:blue",   s=100, zorder=5,
           label=f"VDP best (F1={f_vdp:.2f})")
ax.scatter([c_d2*100],  [n_d2*100],  color="tab:orange", s=100, zorder=5,
           marker="s", label=f"d2 best (F1={f_d2:.2f})")
ax.set_xlabel("Completeness (%)", fontsize=12)
ax.set_ylabel("Contamination (%)", fontsize=12)
ax.set_title(f"ROC Comparison: VDP vs digest2\n(S3M, opposition ±30°, mag 14–26, 2025-03-21)\nN_NEO={N_neo:,}, N_non={N_non:,}", fontsize=10)
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
print(f"  Dataset: {len(data):,} objects  (|beta|<{BETA_MAX_DEG}°, 2025-03-21)")
print(f"  True NEO:     {N_neo:,} ({100*N_neo/len(data):.1f}%)")
print(f"  True non-NEO: {N_non:,} ({100*N_non/len(data):.1f}%)")
print()
print(f"  VDP    best F1={f_vdp:.3f}: completeness={c_vdp:.1%}, purity={1-n_vdp:.1%}")
print(f"  Digest2 best F1={f_d2:.3f}: completeness={c_d2:.1%}, purity={1-n_d2:.1%}")

data[["true_population", "lam_deg", "vlam", "vbeta", "mag_app", "H",
      "P_NEO_vdp", "P_NEO_d2"]].to_parquet(OUT_PARQUET, index=False)
print(f"\nSaved scored data: {OUT_PARQUET}")
