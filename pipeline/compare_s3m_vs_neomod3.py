#!/usr/bin/env python3
"""S3M (source NEO) vs NEOMOD3 (direct sample) orbital-element comparison.

Recreates the "NEO GMM clone sanity checks: 1D distributions" plot from
notebooks/dev/cloning_test_ZI_v2.ipynb (cell 14), but with NEOMOD3 samples (drawn
DIRECTLY from the debiased datacube, no GMM fit/smoothing) in place of GMM clones.
This is the first check before building velocity_density_pipeline_neomod_clone_only.py:
does NEOMOD3, sampled straight (no GMM), give a physically sane match/complement to the
real S3M NEO orbital distribution?

Red = S3M source NEO (real objects). Blue = NEOMOD3 (debiased model, direct sample).
Panels: a, q, e, i, H, M_obs_deg (M at the current n-body epoch, not the notebook's old
two-body epoch).

Read-only: does not touch production maps or parquets.
Output: outputs/neomod3_vs_s3m_comparison/
"""
from __future__ import annotations
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.time import Time

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
os.chdir(W / "neomod")  # s3m_loader searches ['.', 'S3Mdata'] relative to CWD
OUT = W / "outputs" / "neomod3_vs_s3m_comparison"
# current production n-body epoch (Stage-0 cache), NOT the old notebook's stale 2026-05-09
OBSTIME_STR = "2027-08-25T00:00:00"
RNG_SEED = 42
MAX_SOURCE_OBJECTS = None  # None = use the full S3M NEO catalog (268,511 objects; cheap to load)
NEOMOD3_OVERSAMPLE = 10  # NEOMOD3 samples = oversample x len(source), like the existing
                         # _clone_neo_gmm augmentation default (max(500, len(train)*10))

import velocity_density_pipeline_gmm as vdp
import neomod3_sampler as nm3


def as_q(df):
    if "q" in df.columns:
        return df["q"].to_numpy(float)
    return df["a"].to_numpy(float) * (1.0 - df["e"].to_numpy(float))


def mean_anomaly_deg(df, obstime_str):
    t_obs = Time(obstime_str, scale="tdb")
    a = df["a"].to_numpy(float); tp = df["t_p"].to_numpy(float)
    n_rad_day = np.sqrt(vdp.MU_SUN / ((a * vdp.AU_KM) ** 3)) * 86400.0
    M_rad = np.mod((t_obs.mjd - tp) * n_rad_day, 2.0 * np.pi)
    return np.degrees(M_rad)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    print("Loading S3M NEO source population ...", flush=True)
    df_neo, _scorer = vdp.load_s3m_population("neo", verbose=True)
    src = df_neo.copy().reset_index(drop=True)
    if MAX_SOURCE_OBJECTS is not None and len(src) > MAX_SOURCE_OBJECTS:
        idx = rng.choice(len(src), size=MAX_SOURCE_OBJECTS, replace=False)
        src = src.iloc[idx].reset_index(drop=True)
    src["q"] = as_q(src)
    src["M_obs_deg"] = mean_anomaly_deg(src, OBSTIME_STR)
    print(f"  S3M source NEOs: {len(src):,}", flush=True)

    n_nm3 = min(500_000, max(20_000, len(src) * NEOMOD3_OVERSAMPLE))
    print(f"Sampling {n_nm3:,} orbits directly from NEOMOD3 (no GMM) ...", flush=True)
    array4D, centers, widths = nm3.load_neomod3_array()
    nm3_df = nm3.sample_neomod3_orbits(n_nm3, OBSTIME_STR, rng=rng,
                                       array4D=array4D, centers=centers, widths=widths)
    nm3_df["M_obs_deg"] = mean_anomaly_deg(nm3_df, OBSTIME_STR)  # recompute for a consistent panel
    print(f"  NEOMOD3 samples (valid NEO orbits, q<1.3): {len(nm3_df):,}", flush=True)

    # ---- quantitative comparison: KS statistic + summary stats per column ----
    from scipy.stats import ks_2samp
    plot_cols = [("a", "semi-major axis a"), ("q", "perihelion q"), ("e", "eccentricity e"),
                 ("i", "inclination i"), ("H", "absolute magnitude H"),
                 ("M_obs_deg", f"M at {OBSTIME_STR}")]
    rows = []
    for col, label in plot_cols:
        s = src[col].to_numpy(float); n = nm3_df[col].to_numpy(float)
        s = s[np.isfinite(s)]; n = n[np.isfinite(n)]
        ks = ks_2samp(s, n)
        rows.append(dict(column=col, label=label,
                         s3m_median=np.median(s), s3m_mean=s.mean(), s3m_std=s.std(),
                         neomod3_median=np.median(n), neomod3_mean=n.mean(), neomod3_std=n.std(),
                         ks_stat=ks.statistic, ks_pvalue=ks.pvalue))
    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(OUT / "s3m_vs_neomod3_stats.csv", index=False)
    pd.set_option("display.width", 200)
    print("\n=== S3M vs NEOMOD3 per-column comparison ===")
    print(stats_df.to_string(index=False))

    # ---- the 6-panel plot (same layout/style as cloning_test_ZI_v2.ipynb cell 14) ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    axes = axes.ravel()
    for ax, (col, label) in zip(axes, plot_cols):
        ax.hist(src[col], bins=80, density=True, alpha=0.45, label="source NEO (S3M)", color="tab:red")
        ax.hist(nm3_df[col], bins=80, density=True, alpha=0.45, label="NEOMOD3 (direct sample)", color="tab:blue")
        ax.set_xlabel(label); ax.set_ylabel("density"); ax.legend(fontsize=8)
    fig.suptitle(f"NEOMOD3 (direct, no GMM) vs S3M source NEO — 1D orbital-element distributions\n"
                f"N_S3M={len(src):,}  N_NEOMOD3={len(nm3_df):,}  epoch={OBSTIME_STR}")
    fig.savefig(OUT / "s3m_vs_neomod3_1D.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    src[["a", "q", "e", "i", "H", "M_obs_deg"]].to_parquet(OUT / "s3m_source_neo.parquet", index=False)
    nm3_df[["a", "q", "e", "i", "H", "M_obs_deg", "node", "argperi", "t_p"]].to_parquet(
        OUT / "neomod3_samples.parquet", index=False)
    print(f"\nwrote s3m_vs_neomod3_1D.png, s3m_vs_neomod3_stats.csv, "
          f"s3m_source_neo.parquet, neomod3_samples.parquet -> {OUT}")


if __name__ == "__main__":
    main()
