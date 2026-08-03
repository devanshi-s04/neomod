#!/usr/bin/env python3
"""What does the pipeline's GMM cloner DO to a NEOMOD3 input distribution?

Feeds NEOMOD3 orbits into the production cloner `_clone_neo_gmm` (the exact function
velocity_density_pipeline_gmm.py uses for NEO) and compares input vs output orbital elements.

This isolates the GMM SMOOTHING itself: `_clone_neo_gmm` already augments its own training set with
NEOMOD3 (n_neomod3 = max(500, 10*len(train))), so feeding it NEOMOD3 makes the training set
essentially pure NEOMOD3 -- any difference between input and clones is what the mixture model does to
the distribution, not a source mismatch.

Red  = NEOMOD3 source (what we want the clones to reproduce)
Blue = GMM clones (what the current pipeline would actually put in the maps)

Read-only. Output: outputs/neomod3_vs_s3m_comparison/
"""
from __future__ import annotations
import os, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import ks_2samp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
os.chdir(W / "neomod")
OUT = W / "outputs" / "neomod3_vs_s3m_comparison"
EPOCH = "2027-08-25T00:00:00"
N_SOURCE = int(os.environ.get("N_SOURCE", 20_000))    # ~ a real map center's visible NEO count
CLONE_FACTOR = int(os.environ.get("CLONE_FACTOR", 10))
RED, BLUE = "tab:red", "tab:blue"

PLOT_COLS = [("a", "semi-major axis a [AU]"), ("q", "perihelion q [AU]"), ("e", "eccentricity e"),
             ("i", "inclination i [deg]"), ("H", "absolute magnitude H"),
             ("M_obs_deg", f"M at {EPOCH}")]

import velocity_density_pipeline_gmm as vdp


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    nm3 = pd.read_parquet(OUT / "neomod3_samples.parquet")
    src = nm3.sample(n=min(N_SOURCE, len(nm3)), random_state=42).reset_index(drop=True)
    n_clones = len(src) * CLONE_FACTOR
    print(f"source NEOMOD3 orbits: {len(src):,}   -> requesting {n_clones:,} GMM clones", flush=True)

    t0 = time.time()
    clones, gmm, scaler, diag = vdp._clone_neo_gmm(
        src, n_clones=n_clones, obstime_str=EPOCH, n_components=80, random_state=42,
    )
    print(f"GMM clones produced: {len(clones):,}  in {time.time()-t0:.0f}s", flush=True)
    print(f"  converged={diag.get('converged')}  acceptance={diag.get('acceptance_fraction'):.3f}", flush=True)
    src.to_parquet(OUT / "gmm_test_neomod3_source.parquet", index=False)
    clones.to_parquet(OUT / "gmm_test_neomod3_clones.parquet", index=False)

    # ---- KS: does the GMM reproduce its own input? ----
    rows = []
    for col, label in PLOT_COLS:
        if col not in src.columns or col not in clones.columns:
            continue
        a = src[col].to_numpy(float); b = clones[col].to_numpy(float)
        a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
        k = ks_2samp(a, b)
        rows.append(dict(column=col, label=label, n_source=len(a), n_clones=len(b),
                         source_median=np.median(a), clone_median=np.median(b),
                         source_std=a.std(), clone_std=b.std(), ks_stat=k.statistic, ks_pvalue=k.pvalue))
    stats = pd.DataFrame(rows)
    stats.to_csv(OUT / "gmm_test_neomod3_stats.csv", index=False)
    pd.set_option("display.width", 220)
    print("\n=== GMM clones vs their own NEOMOD3 input ===")
    print(stats[["column", "source_median", "clone_median", "source_std", "clone_std", "ks_stat"]].to_string(index=False))

    # ---- plot ----
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for ax, (col, label) in zip(axes.ravel(), PLOT_COLS):
        if col not in src.columns or col not in clones.columns:
            ax.set_visible(False); continue
        ax.hist(src[col], bins=80, density=True, alpha=0.45, color=RED, label="NEOMOD3 source (ours)")
        ax.hist(clones[col], bins=80, density=True, alpha=0.45, color=BLUE, label="GMM clones")
        ax.set_xlabel(label); ax.set_ylabel("density"); ax.legend(fontsize=8)
    fig.suptitle("GMM cloner applied to NEOMOD3 — does the mixture reproduce its own input?\n"
                 f"source {len(src):,}  ->  {len(clones):,} clones   (production `_clone_neo_gmm`, "
                 f"80 components, epoch {EPOCH})")
    fig.savefig(OUT / "gmm_clones_vs_neomod3_1D.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote gmm_clones_vs_neomod3_1D.png + parquets/stats -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
