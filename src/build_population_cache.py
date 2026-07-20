#!/usr/bin/env python3
"""
build_population_cache.py — one-time build of the population density tables the
1A ranging engine looks up. Loading the ~11M-row S3M MBA census every run is
wasteful, so we histogram all four S3M classes onto a common wide grid once and
cache to outputs/pop_cache_wide.npz.

Grid G_wide = (H, log10 a, e, i), chosen wide enough to hold every non-NEO class
(MBA a~1.8-3.3, Trojan a~5.2, TNO a~30-70) as well as NEOs. Axes are independent
of the NEOMOD3 native grid (H,a,e,i; a<=4.2) — the engine looks up the NEOMOD3
numerator and the S3M denominator on their own grids, so they need not share axes.

Outputs (outputs/pop_cache_wide.npz):
  H_edges, loga_edges, e_edges, i_edges   — grid edges
  hist_neo, hist_mba, hist_trojan, hist_tno   — raw COUNT histograms per class
  (densities = counts / bin_volume are formed at load time; see ranging_engine)

Design notes / knobs (documented for D2 review):
  - COUNTS are stored, not densities: the class ratio needs consistent absolute
    units, and counts on one census are automatically consistent. Density (per
    unit element volume) is formed in the engine by dividing by bin volume so the
    ratio is grid-resolution-independent.
  - log10 a axis: MBA fine structure (the dominant contaminant, 96%) is well
    resolved; Trojan/TNO tails are coarse but they are slow-moving and rarely
    velocity-confusable with NEOs anyway.
  - bin counts (NH, Nloga, Ne, Ni) are CLI-overridable knobs.

Usage:
  python src/build_population_cache.py            # default grid
  python src/build_population_cache.py --NH 112   # finer H
"""
import argparse, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import s3m_loader as s3

# ---- default G_wide grid (all knobs) ----
DEF = dict(
    NH=56,     H_min=0.0,  H_max=28.0,     # covers TNO(0.8-10), MBA(12-24), NEO(10-25)
    Nloga=50,  loga_min=0.0, loga_max=2.0,  # a in [1, 100] au
    Ne=25,     e_min=0.0,  e_max=1.0,
    Ni=30,     i_min=0.0,  i_max=90.0,
)


def build_grid(p):
    return dict(
        H_edges=np.linspace(p["H_min"], p["H_max"], p["NH"] + 1),
        loga_edges=np.linspace(p["loga_min"], p["loga_max"], p["Nloga"] + 1),
        e_edges=np.linspace(p["e_min"], p["e_max"], p["Ne"] + 1),
        i_edges=np.linspace(p["i_min"], p["i_max"], p["Ni"] + 1),
    )


def hist_pop(df, edges):
    """4D count histogram of a population on G_wide = (H, log10 a, e, i)."""
    a = df["a"].values
    ok = np.isfinite(a) & (a > 0)
    pts = np.vstack([
        df["H"].values[ok],
        np.log10(a[ok]),
        df["e"].values[ok],
        df["i"].values[ok],
    ]).T
    H, _ = np.histogramdd(
        pts, bins=[edges["H_edges"], edges["loga_edges"], edges["e_edges"], edges["i_edges"]]
    )
    return H


def main():
    ap = argparse.ArgumentParser()
    for k, v in DEF.items():
        ap.add_argument(f"--{k}", type=type(v), default=v)
    ap.add_argument("--out", default="outputs/pop_cache_wide.npz")
    args = vars(ap.parse_args())
    outpath = args.pop("out")

    edges = build_grid(args)
    Path(outpath).parent.mkdir(parents=True, exist_ok=True)

    hists = {}
    t0 = time.time()
    for pop in ["neo", "trojan", "tno", "mba"]:  # mba last (slowest, ~11M rows)
        t = time.time()
        df = s3.define_s3m(pop=pop, verbose=False)
        hists[f"hist_{pop}"] = hist_pop(df, edges)
        print(f"  {pop:6s}: {len(df):>9,} rows -> hist sum {hists[f'hist_{pop}'].sum():.3e} "
              f"({time.time()-t:.0f}s)", flush=True)
        del df

    np.savez_compressed(outpath, **edges, **hists,
                        grid_params=np.array(list(args.items()), dtype=object))
    print(f"saved {outpath} in {time.time()-t0:.0f}s | "
          f"shape {hists['hist_mba'].shape}")


if __name__ == "__main__":
    main()
