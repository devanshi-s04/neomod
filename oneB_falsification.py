#!/usr/bin/env python3
"""
1B falsification test — does VDP's coarse sky-conditioning (667-map grid, 10 deg
lon spacing, nearest-centre snapping) cost F1 at intermediate elongation?

Re-scores the 40-70 deg antisun band of the eval subsample three ways:
  (1) assigned  = nearest grid map                       [baseline]
  (2) mis-1step = one lon step wrong (dlon +/- 10 deg)   [conditioning stress]
  (3) interp3   = 3-nearest inverse-distance weighted    [does interp help?]

Verdict:
  Cond2 ~ Cond1  -> maps locally flat over 10 deg; conditioning NOT the loss -> go to 1A.
  Cond3 >  Cond1 -> interpolation recovers signal      -> 1B-interp worth building.

Pure numpy + the .npz maps. No sklearn (broken here), no propagation, no Hyak.
Run: /astro/users/ds2004/.conda/envs/neofast_py310/bin/python oneB_falsification.py
"""
import os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, "src")
import velocity_density_pipeline_gmm as vdp

PARQUET   = os.environ.get("ONEB_PARQUET", "outputs/phase2/sorcha_comparison_v5_masked.parquet")
MAPDIR    = os.environ.get("ONEB_MAPDIR", "prob_maps_grid")
SUPP_MIN  = 1
BAND      = (40.0, 70.0)          # antisun elongation band to interrogate
GLON      = np.arange(-140, 141, 10)                                   # 29 lon centres
GLAT      = np.array([0,1,-1,2,-2,3,-3,4,-4,5,-5,8,-8,12,-12,18,-18,   # 23 lat centres
                      25,-25,35,-35,50,-50])
REPORT_ALL_BANDS = True           # also print baseline vs interp for every band
OUTFILE   = os.environ.get("ONEB_OUTFILE", "outputs/oneB_falsification_scores.parquet")


# ----------------------------------------------------------------------------- utils
def antisun_dlon(ecl_lon_deg, mjd_utc):
    T = (mjd_utc + 2400000.5) - 2451545.0
    lam_sun = (280.46 + 0.9856474 * T) % 360.0
    anti = (lam_sun + 180.0) % 360.0
    return ((ecl_lon_deg - anti + 180.0) % 360.0) - 180.0


def map_name(dlon, lat):
    return f"{MAPDIR}/prob_maps_grid_dlon{int(dlon):+04d}_lat{int(lat):+03d}.npz"


def nearest_centre(dlon, lat):
    """Vectorised nearest (GLON, GLAT) grid centre for each row (haversine on dlon/lat)."""
    dl = np.deg2rad(dlon[:, None] - GLON[None, :])                     # (N, 29)
    # lat term evaluated after we pick lon; do full 2D nearest instead:
    LON2, LAT2 = np.meshgrid(GLON, GLAT, indexing="ij")               # (29, 23)
    cl = LON2.ravel(); ca = LAT2.ravel()                              # (667,)
    dphi = np.deg2rad(lat[:, None] - ca[None, :])
    dlam = np.deg2rad(dlon[:, None] - cl[None, :])
    a = np.sin(dphi/2)**2 + np.cos(np.deg2rad(lat[:, None]))*np.cos(np.deg2rad(ca[None, :]))*np.sin(dlam/2)**2
    d = 2*np.arcsin(np.sqrt(np.clip(a, 0, 1)))                        # (N, 667) angular dist
    order = np.argsort(d, axis=1)
    return cl, ca, d, order


def score_by_map(df_idx, vlam, vbeta, mag, map_lon, map_lat):
    """Score rows given a per-row (map_lon, map_lat) assignment. Group by map,
    load-score-evict. Returns P(NEO) array aligned to df_idx order."""
    P = np.full(len(df_idx), np.nan)
    key = np.char.add(map_lon.astype(int).astype(str), np.char.add("_", map_lat.astype(int).astype(str)))
    uniq = pd.unique(key)
    for k in uniq:
        m = key == k
        dl, la = (int(x) for x in k.split("_"))
        path = map_name(dl, la)
        if not os.path.exists(path):
            continue
        pm = vdp.ProbMapSet.from_npz(path, support_mask_min=SUPP_MIN,
                                     mask_radius_deg_per_day=np.inf)
        out = pm.score_visible(vlam[m], vbeta[m], mag[m])
        P[m] = out["NEO"]
        del pm
    return P


def best_f1(y, s):
    """Hand-rolled best-F1 (sklearn broken here). Returns (F1, compl, contam, thr)."""
    ok = np.isfinite(s)
    y = y[ok]; s = s[ok]
    order = np.argsort(-s, kind="mergesort")
    ys = y[order].astype(float); ss = s[order]
    tp = np.cumsum(ys); fp = np.cumsum(1.0 - ys)
    P = tp / np.maximum(tp + fp, 1e-12)                # precision
    R = tp / max(y.sum(), 1)                           # recall / completeness
    f1 = 2 * P * R / np.maximum(P + R, 1e-12)
    bi = int(np.argmax(f1))
    return f1[bi], R[bi], 1.0 - P[bi], ss[bi]


# ----------------------------------------------------------------------------- load
import re
def parse_centre(pmf):
    m = re.search(r"dlon([+-]\d+)_lat([+-]\d+)", str(pmf))
    return (int(m.group(1)), int(m.group(2))) if m else (np.nan, np.nan)

t0 = time.time()
df = pd.read_parquet(PARQUET, columns=["ecl_lon", "ecl_lat", "mjd0_utc",
                                       "vlam", "vbeta", "mean_mag", "q_au",
                                       "prob_map_file", "P_NEO_vdp", "P_NEO_d2"])
df["dlon"] = antisun_dlon(df.ecl_lon.values, df.mjd0_utc.values)
df["absdlon"] = np.abs(df.dlon.values)
df["is_neo"] = (df.q_au.values < 1.3)
cen = np.array([parse_centre(p) for p in df.prob_map_file.values], dtype=float)
df["cen_lon"] = cen[:, 0]; df["cen_lat"] = cen[:, 1]
df = df[np.isfinite(df.cen_lon.values)].reset_index(drop=True)
print(f"loaded {len(df):,} tracklets in {time.time()-t0:.1f}s | "
      f"NEO={int(df.is_neo.sum()):,}")

# restrict to the 40-70 gap band (use canonical assigned centre lon, matches parquet)
band = df[(np.abs(df.cen_lon.values) >= BAND[0]) & (np.abs(df.cen_lon.values) < BAND[1])].reset_index(drop=True)
vl = band.vlam.values; vb = band.vbeta.values; mg = band.mean_mag.values
y = band.is_neo.values
print(f"\n=== {BAND[0]:.0f}-{BAND[1]:.0f} deg band: {len(band):,} tracklets, "
      f"{int(y.sum()):,} NEO ({100*y.mean():.1f}%) ===")

# 3-nearest grid centres to each tracklet (for interp + neighbour lookup)
cl, ca, dmat, order = nearest_centre(band.dlon.values, band.ecl_lat.values)

# ---------------------------------------------------------------- condition 1: CANONICAL
# score with the STORED prob_map_file assignment (reproduces parquet P_NEO_vdp)
assigned_lon = band.cen_lon.values; assigned_lat = band.cen_lat.values
t = time.time()
P1 = score_by_map(band.index, vl, vb, mg, assigned_lon, assigned_lat)
f1a = best_f1(y, P1); f1d = best_f1(y, band.P_NEO_d2.values)
f1_stored = best_f1(y, band.P_NEO_vdp.values)
print(f"\n[cond1 canonical]  scored in {time.time()-t:.1f}s | "
      f"check: my F1={f1a[0]:.3f} vs stored P_NEO_vdp F1={f1_stored[0]:.3f} "
      f"corr={np.corrcoef(np.nan_to_num(P1), band.P_NEO_vdp.values)[0,1]:.3f}")

# ------------------------------------------------------- condition 2: one lon step wrong
# shift dlon by +/-10 deg (toward Sun / toward antisun), re-snap lat unchanged
def step_lon(lon, sign):
    out = lon + sign * 10.0
    out = np.clip(out, GLON.min(), GLON.max())
    # snap to nearest existing lon centre
    return GLON[np.argmin(np.abs(out[:, None] - GLON[None, :]), axis=1)]

t = time.time()
P2p = score_by_map(band.index, vl, vb, mg, step_lon(assigned_lon, +1), assigned_lat)
P2m = score_by_map(band.index, vl, vb, mg, step_lon(assigned_lon, -1), assigned_lat)
f2p = best_f1(y, P2p); f2m = best_f1(y, P2m)
print(f"[cond2 mis-1step]  scored in {time.time()-t:.1f}s")

# --------------------------------------------------- condition 3: canonical + 2 neighbours
# centre-0 = STORED map (P1); add the 2 nearest grid centres that differ from stored.
t = time.time()
flat_index = {(int(cl[j]), int(ca[j])): j for j in range(len(cl))}
N = len(band)
j_store = np.array([flat_index[(int(assigned_lon[k]), int(assigned_lat[k]))] for k in range(N)])
d0 = dmat[np.arange(N), j_store]
# 2 nearest neighbours excluding the stored centre (vectorised; order is dist-sorted)
c = order[:, :4]
keep = c != j_store[:, None]
sel = np.argsort(~keep, axis=1, kind="stable")   # keep=True cols first, dist order preserved
pick = np.take_along_axis(c, sel, axis=1)
nb1 = pick[:, 0]; nb2 = pick[:, 1]
d1 = dmat[np.arange(N), nb1]; d2 = dmat[np.arange(N), nb2]
w0 = 1.0/np.maximum(d0, 1e-6); w1 = 1.0/np.maximum(d1, 1e-6); w2 = 1.0/np.maximum(d2, 1e-6)
wsum = w0 + w1 + w2
Pb = score_by_map(band.index, vl, vb, mg, cl[nb1], ca[nb1])
Pc = score_by_map(band.index, vl, vb, mg, cl[nb2], ca[nb2])
P3 = (w0*np.nan_to_num(P1) + w1*np.nan_to_num(Pb) + w2*np.nan_to_num(Pc)) / wsum
f3 = best_f1(y, P3)
print(f"[cond3 interp3]    scored in {time.time()-t:.1f}s (canonical + 2 nearest neighbours)")

# ----------------------------------------------------------------------------- report
def line(name, r):
    print(f"  {name:22s} F1={r[0]:.3f}  compl={100*r[1]:.1f}%  contam={100*r[2]:.1f}%  thr={r[3]:.3f}")

print(f"\n===== 1B FALSIFICATION RESULT ({BAND[0]:.0f}-{BAND[1]:.0f} deg band) =====")
line("cond1 assigned (VDP)", f1a)
line("cond2 mis +1 lon step", f2p)
line("cond2 mis -1 lon step", f2m)
line("cond3 interp 3-nearest", f3)
line("digest2 (reference)", f1d)

dF_mis = min(f2p[0], f2m[0]) - f1a[0]
dF_int = f3[0] - f1a[0]
print(f"\n  mis-1step   dF1 = {dF_mis:+.3f}  (large negative -> 10 deg resolution matters)")
print(f"  interp3     dF1 = {dF_int:+.3f}  (positive -> 1B interpolation worth building)")
print(f"  gap to d2   dF1 = {f1a[0]-f1d[0]:+.3f}")
if abs(dF_mis) < 0.005 and dF_int < 0.005:
    print("\n  VERDICT: conditioning is NOT the bottleneck. Skip 5deg regen -> go to 1A.")
elif dF_int >= 0.005:
    print("\n  VERDICT: interpolation helps. Build 1B-interp (Arnor); consider 5deg regen (Hyak).")
else:
    print("\n  VERDICT: mis-step hurts but interp doesn't recover it -> denser grid (Hyak).")

# save per-condition scores for later inspection
out = band[["dlon", "ecl_lat", "mean_mag", "q_au", "is_neo", "P_NEO_d2"]].copy()
out["P_assigned"] = P1; out["P_mis_plus"] = P2p; out["P_mis_minus"] = P2m; out["P_interp3"] = P3
out_dir = os.path.dirname(OUTFILE)
if out_dir:
    os.makedirs(out_dir, exist_ok=True)
out.to_parquet(OUTFILE, index=False)
print(f"\nsaved {OUTFILE} | total {time.time()-t0:.0f}s")
