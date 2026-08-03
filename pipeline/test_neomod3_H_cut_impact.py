#!/usr/bin/env python3
"""Does cutting NEOMOD3 at H<25 remove DETECTABLE objects? And does it remove the high-|v| tail?

H is ABSOLUTE magnitude; the VDP map bins (14-25) and LSST's ~24.5 limit are APPARENT magnitude.
For NEOs these decouple, because NEOs come close:  mag_app ~ H + 5log10(r*Delta) + phase.
At Delta ~ 0.05 AU an H=28 object appears at mag ~21.5 -- easily detected.

Angular velocity also scales ~ v_tan/Delta, so the SMALL-Delta (close) objects are BOTH the
detectable faint-H ones AND the high-|v| ones. This script measures both effects directly by
projecting NEOMOD3 samples through the same machinery the map builder uses.

Read-only. Output: outputs/neomod3_H_cut_test/
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
os.chdir(W / "neomod")
OUT = W / "outputs" / "neomod3_H_cut_test"
OBSTIME = "2027-08-25T00:00:00"
N_DRAW = int(os.environ.get("N_DRAW", 2_000_000))
MAP_MAG_LO, MAP_MAG_HI = 14.0, 25.0     # the VDP maps' apparent-mag range
LSST_LIM = 24.5
H_CUT = 25.0
V_LIM = 2.0                              # the production velocity-grid boundary

import velocity_density_pipeline_gmm as vdp
import neomod3_sampler as nm3s


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)

    print(f"sampling {N_DRAW:,} NEOMOD3 orbits (full H range) ...", flush=True)
    df = nm3s.sample_neomod3_orbits(N_DRAW, OBSTIME, rng=rng)
    print(f"  valid NEO orbits: {len(df):,}", flush=True)

    # project through the SAME machinery the map builder uses (all-sky, so we measure the
    # intrinsic apparent-mag / velocity distribution rather than one patch's sky cut)
    _, scorer = vdp.load_s3m_population("neo", verbose=False)   # observer geometry only
    print("projecting to apparent mag + sky velocity ...", flush=True)
    vis = vdp.build_visible_subset_dataframe(
        df, obstime_str=OBSTIME, scorer=scorer, max_sep_deg=180.0,
        chunk=200_000, show_progress=False,
        center_mode="custom_ecliptic", center_lon_deg=0.0, center_lat_deg=0.0,
    )
    print(f"  projected: {len(vis):,}", flush=True)

    vis = vis.copy()
    vis["maxabs_v"] = np.maximum(vis.vlam.abs(), vis.vbeta.abs())
    vis["faint_H"] = vis.H >= H_CUT
    in_map = vis.mag_app.between(MAP_MAG_LO, MAP_MAG_HI)
    vis["in_map_range"] = in_map
    vis.to_parquet(OUT / "neomod3_projected.parquet", index=False)

    n = len(vis); nf = int(vis.faint_H.sum()); nb = n - nf
    print("\n================ RESULT ================", flush=True)
    print(f"projected NEOMOD3 objects: {n:,}   H<{H_CUT:.0f}: {nb:,}   H>={H_CUT:.0f}: {nf:,}", flush=True)

    # --- Q1: do faint-H objects reach detectable apparent magnitudes? ---
    faint = vis[vis.faint_H]
    bright = vis[~vis.faint_H]
    print(f"\n[Q1] Do H>={H_CUT:.0f} objects reach DETECTABLE apparent magnitudes?")
    print(f"  H>={H_CUT:.0f} with mag_app in map range [{MAP_MAG_LO},{MAP_MAG_HI}]: "
          f"{int(faint.in_map_range.sum()):,} / {len(faint):,} ({100*faint.in_map_range.mean():.2f}%)")
    print(f"  H>={H_CUT:.0f} with mag_app < {LSST_LIM} (LSST limit):        "
          f"{int((faint.mag_app < LSST_LIM).sum()):,} ({100*(faint.mag_app < LSST_LIM).mean():.2f}%)")
    print(f"  their mag_app: min {faint.mag_app.min():.2f}  p1 {faint.mag_app.quantile(.01):.2f}  "
          f"median {faint.mag_app.median():.2f}")

    # --- Q2: of everything IN the map range, how much would the H cut delete? ---
    inmap = vis[in_map]
    print(f"\n[Q2] Of objects that LAND IN the map range [{MAP_MAG_LO},{MAP_MAG_HI}]:")
    print(f"  total: {len(inmap):,}")
    if len(inmap):
        print(f"  from H>={H_CUT:.0f} (would be DELETED by the cut): {int(inmap.faint_H.sum()):,} "
              f"({100*inmap.faint_H.mean():.1f}%)")

    # --- Q3: the high-velocity tail (the kNN-bleed / +-5 grid problem) ---
    fast = vis[vis.maxabs_v > V_LIM]
    fast_inmap = vis[(vis.maxabs_v > V_LIM) & in_map]
    print(f"\n[Q3] The HIGH-VELOCITY tail (|v| > {V_LIM} deg/day) — the +-5-grid / kNN-bleed region:")
    print(f"  all |v|>{V_LIM}: {len(fast):,}   of which H>={H_CUT:.0f}: {int(fast.faint_H.sum()):,} "
          f"({100*fast.faint_H.mean():.1f}%)" if len(fast) else "  none")
    if len(fast_inmap):
        print(f"  |v|>{V_LIM} AND in map range: {len(fast_inmap):,}   of which H>={H_CUT:.0f}: "
              f"{int(fast_inmap.faint_H.sum()):,} ({100*fast_inmap.faint_H.mean():.1f}%)")
        print(f"  --> cutting H<{H_CUT:.0f} would remove {100*fast_inmap.faint_H.mean():.1f}% of the "
              f"detectable high-velocity NEO clones")
    # velocity reach
    for lab, sub in [("H<25", bright), ("H>=25", faint)]:
        s = sub[sub.mag_app.between(MAP_MAG_LO, MAP_MAG_HI)]
        if len(s):
            print(f"  {lab:6s} in-map |v|: median {s.maxabs_v.median():.3f}  p99 {s.maxabs_v.quantile(.99):.3f}  "
                  f"max {s.maxabs_v.max():.3f}   frac|v|>{V_LIM}: {100*(s.maxabs_v>V_LIM).mean():.2f}%")

    # --- Q4: geometry — why (Delta) ---
    print(f"\n[Q4] Geocentric distance (why faint objects are visible at all):")
    for lab, sub in [("H<25 in-map", bright[bright.mag_app.between(MAP_MAG_LO, MAP_MAG_HI)]),
                     ("H>=25 in-map", faint[faint.mag_app.between(MAP_MAG_LO, MAP_MAG_HI)])]:
        if len(sub) and "delta_au" in sub.columns:
            print(f"  {lab:14s} Delta: median {sub.delta_au.median():.4f} AU  p10 {sub.delta_au.quantile(.1):.4f}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    ax = axes[0]
    ax.hist(bright.mag_app, bins=80, alpha=.55, label=f"H<{H_CUT:.0f}", color="tab:blue")
    ax.hist(faint.mag_app, bins=80, alpha=.55, label=f"H>={H_CUT:.0f} (cut would delete)", color="tab:orange")
    for x, c, l in [(MAP_MAG_HI, "k", "map limit 25"), (LSST_LIM, "r", "LSST 24.5")]:
        ax.axvline(x, color=c, ls="--", lw=1.2, label=l)
    ax.set_yscale("log"); ax.set_xlabel("apparent magnitude"); ax.set_ylabel("count")
    ax.legend(fontsize=8); ax.set_title("Faint-H objects DO reach detectable apparent mags", fontsize=10)

    ax = axes[1]
    b = bright[bright.mag_app.between(MAP_MAG_LO, MAP_MAG_HI)]
    f_ = faint[faint.mag_app.between(MAP_MAG_LO, MAP_MAG_HI)]
    bins = np.linspace(0, 6, 80)
    ax.hist(b.maxabs_v, bins=bins, alpha=.55, label=f"H<{H_CUT:.0f}", color="tab:blue", density=True)
    ax.hist(f_.maxabs_v, bins=bins, alpha=.55, label=f"H>={H_CUT:.0f}", color="tab:orange", density=True)
    ax.axvline(V_LIM, color="k", ls="--", lw=1.2, label=f"map grid |v|={V_LIM}")
    ax.set_yscale("log"); ax.set_xlabel("max(|vlam|,|vbeta|) deg/day"); ax.set_ylabel("density")
    ax.legend(fontsize=8); ax.set_title("in-map objects: velocity reach by H", fontsize=10)

    ax = axes[2]
    sub = vis[in_map]
    sc = ax.scatter(sub.delta_au if "delta_au" in sub.columns else sub.mag_app,
                    sub.maxabs_v, c=sub.H, s=2, alpha=.35, cmap="viridis")
    ax.axhline(V_LIM, color="k", ls="--", lw=1.2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("geocentric distance Delta [AU]" if "delta_au" in sub.columns else "apparent mag")
    ax.set_ylabel("max(|vlam|,|vbeta|) deg/day")
    plt.colorbar(sc, ax=ax, label="H")
    ax.set_title("close objects are fast AND faint-H-detectable", fontsize=10)
    fig.suptitle(f"Impact of cutting NEOMOD3 at H<{H_CUT:.0f}  (N={len(vis):,} projected, epoch {OBSTIME})")
    fig.savefig(OUT / "neomod3_H_cut_impact.png", dpi=150, bbox_inches="tight")
    print(f"\nwrote -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
