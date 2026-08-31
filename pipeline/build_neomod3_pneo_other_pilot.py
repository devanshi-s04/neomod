#!/usr/bin/env python3
"""Selected-centre NEO-vs-OTHER pilot map builder.

Builds rho_NEO and a DIRECTLY POOLED, physically weighted rho_OTHER for a handful of centres.
Every geometric and statistical setting is taken from the frozen builder by importing it, so the
only scientific change is the non-NEO aggregation.

WRITES ONLY under outputs/neomod3_pneo_other_pilot/pilot_maps/. The frozen map root is opened
read-only (for rho_NEO) and never modified.

THE OTHER ESTIMATOR
-------------------
The frozen builder estimates each non-NEO population separately and applies one scalar physical
normalisation:

    rho_pop(x) = n0_pop(x) / eff_pop ,      eff_pop = f_split(pop, magnitude_bin)

Pooling MBA + TNO + Trojans is NOT a matter of concatenating points, because f_split differs by
population and by magnitude bin (measured: MBA ~0.60, TNO 0.52-0.74, Trojans 0.49-0.64; all 42
populated bins differ). Every pooled sample therefore carries its own physical weight
w_i = 1 / f_split(pop_i, bin), and an unweighted kNN density divided by a single scalar would be
wrong.

This pilot keeps the sealed closed-form kNN estimator untouched and applies the local mean physical
weight of the same k neighbours it already queries:

    n0(x)        = ( k(k+1)/2 - 1/2 ) / ( pi * sum_j d_j^2 )      <- sealed, unchanged
    wbar(x)      = (1/k) * sum_{j in kNN(x)} w_j                  <- local mean physical weight
    rho_OTHER(x) = n0(x) * wbar(x)

When every pooled sample shares one weight w this is exactly n0(x) * w = n0(x) / eff, i.e. it
reduces identically to the frozen single-population line `dens = dens / eff`. The pilot estimator is
a strict generalisation of the frozen one, not a different estimator.

VALIDITY -- identical rule to the frozen builder, applied to the pooled sample:
    n_pooled >= k + 1  -> VALID
    2 <= n <= k        -> INVALID "insufficient_support"
    n < 2              -> INVALID "no_samples"/"below_min"
A population with no split fraction in a bin cannot be weighted, so its samples are EXCLUDED from
the pool and the exclusion is recorded. Nothing is zero-filled and no epsilon floor is applied.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"pipeline"))
sys.path.insert(0, str(W/"neomod"/"src"))
sys.path.insert(0, str(W/"neomod"/"adam_core_stub"))

import build_neomod3_mag025_k150_maps as FROZEN     # geometry/settings come from here

FROZEN_MAPS = W/"outputs"/"neomod3_mag025_k150_maps_v2"
OUT_ROOT    = W/"outputs"/"neomod3_pneo_other_pilot"/"pilot_maps"
K_OTHER     = 10          # default: same k as every frozen non-NEO population
NONNEO      = ("MBA", "TNO", "Trojans")


def pooled_other_density(pts, wts, grid_points, k, workers):
    """Sealed closed-form kNN density times the local mean physical weight of the same k
    neighbours. Reduces exactly to n0/eff when all weights are equal."""
    from scipy.spatial import cKDTree
    tree = cKDTree(pts)
    d, idx = tree.query(grid_points, k=int(k), workers=int(workers))
    S  = np.einsum("ij,ij->i", d, d)
    n0 = (k*(k+1)/2.0 - 0.5) / (np.pi * S)          # sealed closed form, verbatim
    wbar = wts[idx].mean(axis=1)                     # local mean physical weight
    return n0 * wbar


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--center", required=True)
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--out-dir", default=str(OUT_ROOT))
    ap.add_argument("--k-other", type=int, default=K_OTHER,
                    help="k for the pooled OTHER density. Default 10, matching every frozen "
                         "non-NEO population. THE ONLY QUANTITY a sensitivity study may vary.")
    a = ap.parse_args()

    k_other = int(a.k_other)
    # main() chdirs to neomod/ below, so a relative --out-dir would resolve against the wrong
    # directory. Make it absolute FIRST.
    outdir_abs = Path(a.out_dir).resolve()
    import os; os.chdir(W/"neomod")
    import velocity_density_pipeline_neomod_clone_only as base

    seal  = json.load(open(W/"outputs/splits/MAP_BUILD_SEAL.json"))
    epoch = seal["grid"]["ref_obstime"]                       # frozen map reference epoch
    label = a.center
    dlon, lat  = FROZEN.parse_center(label)
    clon, clat = FROZEN.center_lonlat(epoch, dlon, lat)
    bins   = FROZEN.mag025_bins()
    workers= int(a.n_jobs)

    n_side = int(round((FROZEN.GRID_LIM[1]-FROZEN.GRID_LIM[0])/FROZEN.GRID_STEP)) + 1
    x_grid = np.linspace(FROZEN.GRID_LIM[0],
                         FROZEN.GRID_LIM[0]+(n_side-1)*FROZEN.GRID_STEP, n_side)
    y_grid = x_grid.copy()
    X0, Y0 = np.meshgrid(x_grid, y_grid)
    grid_points = np.column_stack([X0.ravel(), Y0.ravel()])
    pixel_area  = FROZEN.GRID_STEP**2

    f_by = json.load(open(FROZEN.SPLIT_MAG025))["f_by_population_magbin"]
    _, scorer = base.load_s3m_population("neo", verbose=False)

    # ---- non-NEO GEN parents, exactly as the frozen builder selects them -----------
    cache = pd.read_parquet(FROZEN.EPOCH_CACHE)
    man   = pd.read_parquet(FROZEN.SPLIT_MANIFEST, columns=["ObjID","split"])
    gen   = set(man.ObjID[man.split == "GEN"])
    src   = {p: cache[(cache.population==p) & (cache.ObjID.isin(gen))].reset_index(drop=True)
             for p in NONNEO}
    del cache
    for p in NONNEO:
        print(f"[{label}] {p}: {len(src[p]):,} GEN parents", flush=True)

    # ---- rho_NEO copied verbatim from the frozen maps -----------------------------
    frozen_npz = FROZEN_MAPS/f"mag025_k150_{label}.npz"
    fz = np.load(frozen_npz, allow_pickle=True)
    frozen_keys = set(fz.keys())
    print(f"[{label}] frozen map {frozen_npz.name}: "
          f"{sum(1 for k in frozen_keys if k.startswith('density__NEO__'))} NEO bins", flush=True)

    arrays, coverage = {}, []
    t0 = time.time()
    for b in bins:
        lab = b["label"]; lo, hi = b["mag_min"], b["mag_max"]

        neo_key = f"density__NEO__{lab}"
        has_neo = neo_key in frozen_keys
        if has_neo:
            arrays[f"rho_NEO__{lab}"] = np.asarray(fz[neo_key], dtype=np.float32)

        # pooled OTHER
        pts_l, wts_l, per_pop, excluded = [], [], {}, {}
        for pop in NONNEO:
            f = f_by.get(pop, {}).get(lab)
            if f is None:
                excluded[pop] = "no_split_fraction"; per_pop[pop] = 0; continue
            s = src[pop]
            magv = s["mag_app"].to_numpy(float)
            sel = np.isfinite(magv) & (magv >= lo) & (magv < hi)
            ds = s[sel]
            if not len(ds):
                per_pop[pop] = 0; continue
            vis = base.build_visible_subset_dataframe(
                ds, obstime_str=epoch, scorer=scorer, max_sep_deg=FROZEN.MAX_SEP_DEG,
                chunk=100_000, show_progress=False, center_mode="custom_ecliptic",
                center_lon_deg=clon, center_lat_deg=clat)
            per_pop[pop] = int(len(vis))
            if len(vis):
                pts_l.append(np.column_stack([vis["vlam"].to_numpy(float),
                                              vis["vbeta"].to_numpy(float)]))
                wts_l.append(np.full(len(vis), 1.0/float(f)))

        n_pool = int(sum(per_pop.values()))
        reason, valid = "", False
        if n_pool == 0:      reason = "no_samples"
        elif n_pool < 2:     reason = "below_min"
        elif n_pool <= k_other: reason = "insufficient_support"
        else:                valid = True

        rec = dict(center=label, population="OTHER", magnitude_bin=lab, mag_min=lo, mag_max=hi,
                   n_pooled_visible=n_pool, k_requested=k_other,
                   k_effective=(k_other if valid else None), valid=bool(valid),
                   reason=(reason or "ok"),
                   n_MBA=per_pop.get("MBA",0), n_TNO=per_pop.get("TNO",0),
                   n_Trojans=per_pop.get("Trojans",0),
                   excluded_populations=(",".join(sorted(excluded)) or ""),
                   f_MBA=f_by.get("MBA",{}).get(lab), f_TNO=f_by.get("TNO",{}).get(lab),
                   f_Trojans=f_by.get("Trojans",{}).get(lab),
                   neo_density_available=bool(has_neo))
        if valid:
            pts = np.vstack(pts_l); wts = np.concatenate(wts_l)
            dens = pooled_other_density(pts, wts, grid_points, k_other, workers).reshape(X0.shape)
            arrays[f"rho_OTHER__{lab}"] = dens.astype(np.float32)
            rec["density_integral"] = float(np.nansum(dens)*pixel_area)
            rec["total_physical_weight"] = float(wts.sum())
            rec["mean_physical_weight"] = float(wts.mean())
            print(f"[{label}] OTHER {lab} n={n_pool:>7,} "
                  f"(MBA {per_pop.get('MBA',0)} TNO {per_pop.get('TNO',0)} "
                  f"Tro {per_pop.get('Trojans',0)}) {time.time()-t0:6.1f}s", flush=True)
        coverage.append(rec)

    outdir = outdir_abs; outdir.mkdir(parents=True, exist_ok=True)
    meta = dict(center=label, epoch=epoch, center_lon_deg=clon, center_lat_deg=clat,
                grid_lim=list(FROZEN.GRID_LIM), grid_step=FROZEN.GRID_STEP,
                k_OTHER=k_other, k_NEO_from_frozen_maps=150,
                gaussian_smoothing=False, support_masking=False,
                density_engine="closed_form + local mean physical weight",
                max_sep_deg=FROZEN.MAX_SEP_DEG,
                bin_scheme=dict(lo=FROZEN.MAG_LO, hi=FROZEN.MAG_HI, step=FROZEN.MAG_STEP,
                                semantics="half-open [lo, lo+step)"),
                frozen_map_source=str(frozen_npz),
                rho_NEO_origin="copied verbatim from the frozen map set",
                rho_OTHER_origin="directly pooled MBA+TNO+Trojans, per-sample physical weight 1/f_split",
                built_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    # np.savez_compressed appends ".npz" when the name does not already end in it, so the
    # temporary name must itself end in .npz or the rename target will not exist.
    tmp = outdir/f"pilot_{label}.tmp.npz"
    np.savez_compressed(tmp, x_grid=x_grid, y_grid=y_grid,
                        meta_json=json.dumps(meta),
                        coverage_json=json.dumps(coverage), **arrays)
    assert tmp.exists(), f"savez did not produce {tmp}"
    tmp.rename(outdir/f"pilot_{label}.npz")
    pd.DataFrame(coverage).to_csv(outdir/f"coverage_{label}.csv", index=False)
    nv = sum(1 for c in coverage if c["valid"])
    print(f"\n[{label}] OTHER valid in {nv}/{len(bins)} bins; "
          f"NEO present in {sum(1 for c in coverage if c['neo_density_available'])}/{len(bins)}; "
          f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
