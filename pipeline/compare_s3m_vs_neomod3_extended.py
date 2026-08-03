#!/usr/bin/env python3
"""Extended S3M-vs-NEOMOD3 orbital-element comparison (reads the saved parquets — no reload/resample).

Adds to the base 1D comparison:
  * H<25 MATCHED cut — S3M's synthetic catalog stops at H=25 while NEOMOD3 runs to H=28, so the
    full-sample KS on q/e/i mixes "orbital-shape disagreement" with "different H regimes". Cutting
    both at H<25 isolates the pure orbital-shape difference.
  * 2D distributions (a-e, a-i, q-e) as side-by-side density panels + overlaid contours.
  * What NEOMOD3 ADDS: the H 25-28 faint population S3M has no objects for at all.

Red = S3M source NEO (real catalog objects). Blue = NEOMOD3 (debiased model, direct sample).
Read-only. Output: outputs/neomod3_vs_s3m_comparison/
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
OUT = W / "outputs" / "neomod3_vs_s3m_comparison"
OBSTIME_STR = "2027-08-25T00:00:00"
H_CUT = 25.0
RED, BLUE = "tab:red", "tab:blue"

PLOT_COLS = [("a", "semi-major axis a [AU]"), ("q", "perihelion q [AU]"), ("e", "eccentricity e"),
             ("i", "inclination i [deg]"), ("H", "absolute magnitude H"),
             ("M_obs_deg", f"M at {OBSTIME_STR}")]


def ks_table(s3m, nm3, tag):
    rows = []
    for col, label in PLOT_COLS:
        s = s3m[col].to_numpy(float); n = nm3[col].to_numpy(float)
        s = s[np.isfinite(s)]; n = n[np.isfinite(n)]
        if len(s) < 2 or len(n) < 2:
            continue
        k = ks_2samp(s, n)
        rows.append(dict(sample=tag, column=col, label=label, n_s3m=len(s), n_neomod3=len(n),
                         s3m_median=np.median(s), neomod3_median=np.median(n),
                         s3m_mean=s.mean(), neomod3_mean=n.mean(),
                         s3m_std=s.std(), neomod3_std=n.std(),
                         ks_stat=k.statistic, ks_pvalue=k.pvalue))
    return pd.DataFrame(rows)


def plot_1d(s3m, nm3, title, fname):
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for ax, (col, label) in zip(axes.ravel(), PLOT_COLS):
        ax.hist(s3m[col], bins=80, density=True, alpha=0.45, label="source NEO (S3M)", color=RED)
        ax.hist(nm3[col], bins=80, density=True, alpha=0.45, label="NEOMOD3 (direct sample)", color=BLUE)
        ax.set_xlabel(label); ax.set_ylabel("density"); ax.legend(fontsize=8)
    fig.suptitle(title)
    fig.savefig(OUT / fname, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_2d(s3m, nm3, pairs, title, fname):
    """Row 1: S3M 2D density. Row 2: NEOMOD3. Row 3: overlaid contours (red=S3M, blue=NEOMOD3)."""
    fig, axes = plt.subplots(3, len(pairs), figsize=(5.2 * len(pairs), 13), constrained_layout=True)
    for j, (cx, cy, lx, ly, rng) in enumerate(pairs):
        sx, sy = s3m[cx].to_numpy(float), s3m[cy].to_numpy(float)
        nx, ny = nm3[cx].to_numpy(float), nm3[cy].to_numpy(float)
        for row, (xx, yy, nm, cmap) in enumerate([(sx, sy, "S3M source NEO", "Reds"),
                                                  (nx, ny, "NEOMOD3", "Blues")]):
            ax = axes[row, j]
            ax.hist2d(xx, yy, bins=70, range=rng, cmap=cmap, density=True)
            ax.set_xlabel(lx); ax.set_ylabel(ly); ax.set_title(f"{nm}: {cx}–{cy}", fontsize=10)
        # contour overlay
        ax = axes[2, j]
        for xx, yy, color, nm in [(sx, sy, RED, "S3M"), (nx, ny, BLUE, "NEOMOD3")]:
            Hh, xe, ye = np.histogram2d(xx, yy, bins=60, range=rng, density=True)
            Xc = 0.5 * (xe[:-1] + xe[1:]); Yc = 0.5 * (ye[:-1] + ye[1:])
            lv = np.percentile(Hh[Hh > 0], [50, 75, 90, 97]) if (Hh > 0).any() else None
            if lv is not None:
                ax.contour(Xc, Yc, Hh.T, levels=lv, colors=color, linewidths=1.3, alpha=0.85)
            ax.plot([], [], color=color, label=nm)
        ax.set_xlim(rng[0]); ax.set_ylim(rng[1])
        ax.set_xlabel(lx); ax.set_ylabel(ly); ax.legend(fontsize=9)
        ax.set_title(f"overlay: {cx}–{cy}", fontsize=10)
    fig.suptitle(title)
    fig.savefig(OUT / fname, dpi=150, bbox_inches="tight"); plt.close(fig)


def main():
    s3m = pd.read_parquet(OUT / "s3m_source_neo.parquet")
    nm3 = pd.read_parquet(OUT / "neomod3_samples.parquet")
    print(f"loaded S3M {len(s3m):,} | NEOMOD3 {len(nm3):,}", flush=True)

    # ---------- H<25 matched cut ----------
    s3m_c = s3m[s3m.H < H_CUT].copy()
    nm3_c = nm3[nm3.H < H_CUT].copy()
    print(f"H<{H_CUT} cut: S3M {len(s3m_c):,} ({100*len(s3m_c)/len(s3m):.1f}%) | "
          f"NEOMOD3 {len(nm3_c):,} ({100*len(nm3_c)/len(nm3):.1f}%)", flush=True)

    stats_full = ks_table(s3m, nm3, "full")
    stats_cut = ks_table(s3m_c, nm3_c, f"H<{H_CUT:.0f}")
    stats = pd.concat([stats_full, stats_cut], ignore_index=True)
    stats.to_csv(OUT / "s3m_vs_neomod3_stats_full_and_Hcut.csv", index=False)
    pd.set_option("display.width", 220)
    print("\n=== KS: full sample ===\n" + stats_full[
        ["column", "n_s3m", "n_neomod3", "s3m_median", "neomod3_median", "ks_stat"]].to_string(index=False))
    print(f"\n=== KS: H<{H_CUT:.0f} matched ===\n" + stats_cut[
        ["column", "n_s3m", "n_neomod3", "s3m_median", "neomod3_median", "ks_stat"]].to_string(index=False))
    print("\n=== KS change (full -> H<25), orbital elements only ===")
    for col in ["a", "q", "e", "i"]:
        f = stats_full[stats_full.column == col].ks_stat.iloc[0]
        c = stats_cut[stats_cut.column == col].ks_stat.iloc[0]
        print(f"  {col}: {f:.4f} -> {c:.4f}  ({c - f:+.4f})")

    # ---------- plots ----------
    plot_1d(s3m_c, nm3_c,
            f"NEOMOD3 vs S3M source NEO — 1D distributions, H<{H_CUT:.0f} MATCHED\n"
            f"N_S3M={len(s3m_c):,}  N_NEOMOD3={len(nm3_c):,}  epoch={OBSTIME_STR}",
            f"s3m_vs_neomod3_1D_Hcut{int(H_CUT)}.png")

    pairs = [("a", "e", "semi-major axis a [AU]", "eccentricity e", [[0, 4.2], [0, 1]]),
             ("a", "i", "semi-major axis a [AU]", "inclination i [deg]", [[0, 4.2], [0, 88]]),
             ("q", "e", "perihelion q [AU]", "eccentricity e", [[0, 1.3], [0, 1]])]
    plot_2d(s3m, nm3, pairs,
            "S3M vs NEOMOD3 — 2D orbital-element distributions (full sample)",
            "s3m_vs_neomod3_2D_full.png")
    plot_2d(s3m_c, nm3_c, pairs,
            f"S3M vs NEOMOD3 — 2D orbital-element distributions (H<{H_CUT:.0f} matched)",
            f"s3m_vs_neomod3_2D_Hcut{int(H_CUT)}.png")

    # ---------- what NEOMOD3 adds: the faint tail S3M lacks ----------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    ax = axes[0]
    ax.hist(s3m.H, bins=90, density=False, alpha=0.5, label=f"S3M (max H={s3m.H.max():.2f})", color=RED)
    ax.hist(nm3.H, bins=90, density=False, alpha=0.5, label=f"NEOMOD3 (max H={nm3.H.max():.2f})", color=BLUE)
    ax.axvline(H_CUT, color="k", ls="--", lw=1.2, label=f"S3M cutoff H={H_CUT:.0f}")
    ax.set_xlabel("absolute magnitude H"); ax.set_ylabel("count"); ax.set_yscale("log")
    ax.legend(fontsize=8); ax.set_title("H coverage: S3M truncates at 25, NEOMOD3 reaches 28", fontsize=10)
    ax = axes[1]
    faint = nm3[nm3.H >= H_CUT]
    ax.hist(faint.q, bins=70, density=True, color=BLUE, alpha=0.65,
            label=f"NEOMOD3 H>{H_CUT:.0f} (n={len(faint):,})")
    ax.hist(nm3_c.q, bins=70, density=True, histtype="step", lw=1.6, color="k",
            label=f"NEOMOD3 H<{H_CUT:.0f}")
    ax.set_xlabel("perihelion q [AU]"); ax.set_ylabel("density"); ax.legend(fontsize=8)
    ax.set_title("the faint population S3M has NO objects for", fontsize=10)
    fig.suptitle("What NEOMOD3 adds beyond S3M")
    fig.savefig(OUT / "neomod3_adds_faint_tail.png", dpi=150, bbox_inches="tight"); plt.close(fig)

    frac_faint = 100.0 * len(nm3[nm3.H >= H_CUT]) / len(nm3)
    print(f"\nNEOMOD3 samples with H>={H_CUT:.0f} (absent from S3M entirely): "
          f"{len(nm3[nm3.H >= H_CUT]):,} ({frac_faint:.1f}%)")
    print(f"S3M H range: [{s3m.H.min():.2f}, {s3m.H.max():.2f}] | "
          f"NEOMOD3 H range: [{nm3.H.min():.2f}, {nm3.H.max():.2f}]")
    print(f"\nwrote plots + s3m_vs_neomod3_stats_full_and_Hcut.csv -> {OUT}")


if __name__ == "__main__":
    main()
