#!/usr/bin/env python3
"""NEOMOD3 0.25-mag / k=150 full-grid map builder — ONE parameterized path.

See neomod/docs/NEOMOD3_MAG025_K150_FULLGRID_RUNBOOK.md.

Frozen: 667 centers, epoch 2027-08-25T00:00:00, 30 deg patches, +-5 deg/day at 0.01,
apparent V, 0.25-mag bins over 14 <= V < 25 (44 bins, half-open [lo, lo+0.25)),
NEO k=150, MBA/TNO/Trojans k=10, NO Gaussian smoothing, no support masking.

DENSITY ENGINE (--density-engine)
---------------------------------
`closed_form` (default) and `sealed_quadrature` are two evaluations of the SAME estimator, not two
estimators. The sealed posterior for the local spacing scale d0 is

    p(d0) ∝ prod_j (2/(d0 Γ(j))) exp[-(d_j/d0)^2] (d_j/d0)^(2j-1)
          = C(d) · d0^(-k(k+1)) · exp(-S/d0^2),        S = sum_j d_j^2

(k factors of 1/d0, and sum_j (2j-1) = k^2). Under the flat prior on d0 that the sealed code uses --
it normalizes with np.trapezoid(p, d0_grid), i.e. measure dd0 -- the posterior mean of
n0 = 1/(pi d0^2) is exactly

    n0 = ( k(k+1)/2 - 1/2 ) / ( pi * sum_j d_j^2 )

so `closed_form` is the analytic value of the integral `sealed_quadrature` computes numerically. It
is faster AND has no discretisation error: the sealed default n_d0_grid=400 is converged at k=10 but
reaches 2.6e-04 at k=50 and 1.1e-02 at k=100, which at NEO k=150 would be a k-dependent bias.
`--density-engine sealed_quadrature` is retained so the two can be compared (gate check 8).

INSUFFICIENT-SUPPORT RULE (single, explicit, never widens a slice)
-----------------------------------------------------------------
    n >= k_req + 1   -> VALID,   effective_k = k_req
    2 <= n <= k_req  -> INVALID, effective_k = None, reason "insufficient_support"
    n < 2            -> INVALID, effective_k = None, reason "no_samples"/"below_min"
    split fraction f missing for a non-NEO cell -> INVALID, reason "no_split_fraction"

The requested k is NEVER reduced and the magnitude slice is NEVER widened. INVALID densities are
NaN, never 0, and INVALID cells are not stored as all-NaN arrays -- they are recorded in the
coverage table with n, k_requested, k_effective, valid, reason.

Posteriors are formed only from VALID populations, and every bin records `populations_used` and
`populations_invalid` plus a `partial_denominator` flag, so a reduced denominator is explicit rather
than silent.

Writes ONLY under outputs/neomod3_mag025_k150_maps_v2/. Atomic: writes .tmp then renames, and the
.ok success marker is created only after the per-center validation passes.
"""
from __future__ import annotations

import argparse, copy, hashlib, json, os, subprocess, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
sys.path.insert(0, str(W / "neomod" / "adam_core_stub"))

OUT_ROOT = W / "outputs" / "neomod3_mag025_k150_maps_v2"
CENTER_LIST = W / "outputs/geometric_density_estimator_ablation/frozen_center_list.txt"
SPLIT_MAG025 = W / "outputs/splits/split_provenance_mag025.json"
EPOCH_CACHE = W / "outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet"
SPLIT_MANIFEST = W / "outputs/splits/nonneo_split_manifest.parquet"

MAX_SEP_DEG = 30.0
GRID_LIM = (-5.0, 5.0)
GRID_STEP = 0.01
K_BY_POP = {"NEO": 150, "MBA": 10, "TNO": 10, "Trojans": 10}
POPS = ("NEO", "MBA", "TNO", "Trojans")
NONNEO_POPS = ("MBA", "TNO", "Trojans")
SMOOTHING = False
MAG_LO, MAG_HI, MAG_STEP = 14.0, 25.0, 0.25
N_D0_GRID_REF = 8000          # only for --density-engine sealed_quadrature


def mag025_bins():
    n = int(round((MAG_HI - MAG_LO) / MAG_STEP))
    return [{"label": f"V{MAG_LO+i*MAG_STEP:06.2f}_{MAG_LO+(i+1)*MAG_STEP:06.2f}",
             "mag_min": round(MAG_LO + i * MAG_STEP, 2),
             "mag_max": round(MAG_LO + (i + 1) * MAG_STEP, 2), "index": i} for i in range(n)]


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def git_commit():
    try:
        return subprocess.run(["git", "-C", str(W / "neomod"), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


# ---------------------------------------------------------------- density engines
def density_closed_form(tree, grid_points, k, workers=1):
    """Exact posterior mean of n0 = 1/(pi d0^2): (k(k+1)/2 - 1/2) / (pi * sum_j d_j^2)."""
    d, _ = tree.query(grid_points, k=int(k), workers=int(workers))
    S = np.einsum("ij,ij->i", d, d)          # sum_j d_j^2 without materialising d**2
    return (k * (k + 1) / 2.0 - 0.5) / (np.pi * S)


def density_sealed_quadrature(tree, grid_points, k, workers=1, n_d0_grid=N_D0_GRID_REF):
    import velocity_density_pipeline_neomod_clone_only as base
    return base.evaluate_density_map_full_posterior_2d(
        tree, grid_points, k=int(k), n_d0_grid=int(n_d0_grid),
        show_progress=False, n_jobs=int(workers))


def evaluate_density(engine, tree, grid_points, k, workers):
    if engine == "closed_form":
        return density_closed_form(tree, grid_points, k, workers)
    if engine == "sealed_quadrature":
        return density_sealed_quadrature(tree, grid_points, k, workers)
    raise ValueError(f"unknown density engine {engine!r}")


# ---------------------------------------------------------------- sources
def load_center_list():
    return [c.strip() for c in CENTER_LIST.read_text().split() if c.strip()]


def center_lonlat(epoch, dlon, lat):
    from astropy.time import Time
    from astropy.coordinates import get_sun, GeocentricTrueEcliptic
    t = Time(epoch, scale="utc")
    antisun = (get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t)).lon.deg + 180.0) % 360.0
    return (antisun + dlon) % 360.0, lat


def parse_center(label):
    import re
    m = re.match(r"dlon([+-]\d+)_lat([+-]\d+)", label)
    if not m:
        raise ValueError(f"unparseable center label {label!r}")
    return float(m.group(1)), float(m.group(2))


def build_one_center(a):
    import re
    os.chdir(W / "neomod")
    import velocity_density_pipeline_neomod_clone_only as base
    from scipy.spatial import cKDTree

    seal = json.load(open(W / "outputs/splits/MAP_BUILD_SEAL.json"))
    epoch = seal["grid"]["ref_obstime"]
    centers = load_center_list()
    label = a.center or centers[a.task_id]
    dlon, lat = parse_center(label)
    clon, clat = center_lonlat(epoch, dlon, lat)
    workers = int(a.n_jobs)
    bins = mag025_bins()
    if a.only_bins:
        keep = set(a.only_bins.split(","))
        bins = [b for b in bins if b["label"] in keep]

    outdir = Path(a.out_dir) if a.out_dir else OUT_ROOT
    outdir.mkdir(parents=True, exist_ok=True)
    final = outdir / f"mag025_k150_{label}.npz"
    marker = outdir / f"mag025_k150_{label}.ok"
    if marker.exists() and not a.overwrite:
        try:
            json.load(open(marker))
            np.load(final, allow_pickle=True)["x_grid"]     # open-and-read, not just exists
            print(f"[skip] {label} already complete and readable", flush=True)
            return
        except Exception as e:
            print(f"[redo] {label} marker/file unreadable ({e}); rebuilding", flush=True)

    t_center0 = time.time()
    x_grid = np.arange(GRID_LIM[0], GRID_LIM[1] + 1e-9, GRID_STEP)
    n_side = int(round((GRID_LIM[1] - GRID_LIM[0]) / GRID_STEP)) + 1
    x_grid = np.linspace(GRID_LIM[0], GRID_LIM[0] + (n_side - 1) * GRID_STEP, n_side)
    y_grid = x_grid.copy()
    X0, Y0 = np.meshgrid(x_grid, y_grid)
    grid_points = np.column_stack([X0.ravel(), Y0.ravel()])
    pixel_area = GRID_STEP ** 2

    split025 = json.load(open(SPLIT_MAG025))
    f_by = split025["f_by_population_magbin"]

    # ---- sources -------------------------------------------------------------
    _, scorer = base.load_s3m_population("neo", verbose=False)
    neo_meta = json.load(open(W / "outputs/neomod3_projection_cache/cache_metadata.json"))

    if a.neo_source == "high_cell":
        neo_df = pd.read_parquet(a.neo_source_path)
        eff_neo = float(a.neo_effective_factor)
        neo_src_hash = sha256_file(a.neo_source_path)
        neo_src_desc = f"HIGH single-cell parquet {a.neo_source_path}"
    else:
        neo_df, eff_neo = base._load_neomod3_cache(
            center_lon_deg=clon, center_lat_deg=clat, obstime_str=epoch, max_sep_deg=MAX_SEP_DEG)
        neo_src_hash = "healpix_by_pixel_cache"
        neo_src_desc = f"NEOMOD3 by_pixel GEN cache (n_draws {neo_meta['n_draws']:,})"
    print(f"[{label}] NEO source: {neo_src_desc}  rows {len(neo_df):,}  eff {eff_neo:.6f}",
          flush=True)

    cache = pd.read_parquet(EPOCH_CACHE)
    man = pd.read_parquet(SPLIT_MANIFEST, columns=["ObjID", "split"])
    gen_ids = set(man.ObjID[man.split == "GEN"])
    nonneo = {}
    for pop in NONNEO_POPS:
        sub = cache[(cache.population == pop) & (cache.ObjID.isin(gen_ids))].reset_index(drop=True)
        nonneo[pop] = sub
        print(f"[{label}] {pop}: {len(sub):,} GEN parents", flush=True)
    del cache

    # ---- per (population, bin) ----------------------------------------------
    arrays, coverage, integrals = {}, [], []
    for pop in POPS:
        k_req = K_BY_POP[pop]
        src = neo_df if pop == "NEO" else nonneo[pop]
        magv = src["mag_app"].to_numpy(float)
        for b in bins:
            t0 = time.time()
            lab, lo, hi = b["label"], b["mag_min"], b["mag_max"]
            sel = np.isfinite(magv) & (magv >= lo) & (magv < hi)
            df_sel = src[sel]
            reason, valid, k_eff, n_vis = "", False, None, 0
            f_split = 1.0
            if pop != "NEO":
                f = f_by.get(pop, {}).get(lab)
                if f is None:
                    reason = "no_split_fraction"
                else:
                    f_split = float(f)
            if not reason and len(df_sel):
                vis = base.build_visible_subset_dataframe(
                    df_sel, obstime_str=epoch, scorer=scorer, max_sep_deg=MAX_SEP_DEG,
                    chunk=100_000, show_progress=False, center_mode="custom_ecliptic",
                    center_lon_deg=clon, center_lat_deg=clat)
                n_vis = len(vis)
            elif not reason:
                vis = df_sel
                n_vis = 0
            if not reason:
                if n_vis < 2:
                    reason = "no_samples" if n_vis == 0 else "below_min"
                elif n_vis <= k_req:
                    reason = "insufficient_support"
                else:
                    valid, k_eff = True, k_req
            eff = eff_neo if pop == "NEO" else (1.0 * f_split)
            rec = dict(center=label, population=pop, magnitude_bin=lab, mag_min=lo, mag_max=hi,
                       n_selected=int(len(df_sel)), n_visible=int(n_vis),
                       k_requested=int(k_req), k_effective=(int(k_eff) if k_eff else None),
                       valid=bool(valid), reason=(reason or "ok"),
                       effective_factor=float(eff), split_fraction=float(f_split),
                       physical_weight_per_sample=float(1.0 / eff) if eff else None,
                       total_physical_weight=(float(n_vis / eff) if (valid and eff) else None))
            if not valid:
                coverage.append(rec)
                continue
            pts = np.column_stack([vis["vlam"].to_numpy(float), vis["vbeta"].to_numpy(float)])
            tree = cKDTree(pts)
            dens = evaluate_density(a.density_engine, tree, grid_points, k_req, workers)
            dens = (dens / eff).reshape(X0.shape)          # physical normalisation
            arrays[f"density__{pop}__{lab}"] = dens.astype(np.float32)
            integrals.append(dict(center=label, population=pop, magnitude_bin=lab,
                                  density_integral=float(np.nansum(dens) * pixel_area),
                                  n_visible=int(n_vis),
                                  total_physical_weight=float(n_vis / eff)))
            rec["seconds"] = round(time.time() - t0, 3)
            coverage.append(rec)
            print(f"[{label}] {pop:8s} {lab} n={n_vis:>7,} k={k_req:<4d} "
                  f"{time.time()-t0:6.1f}s", flush=True)

    # ---- posteriors ----------------------------------------------------------
    post_meta = []
    for b in bins:
        lab = b["label"]
        present = [p for p in POPS if f"density__{p}__{lab}" in arrays]
        missing = [p for p in POPS if p not in present]
        if "NEO" not in present:
            post_meta.append(dict(center=label, magnitude_bin=lab, defined=False,
                                  populations_used=present, populations_invalid=missing,
                                  partial_denominator=bool(missing),
                                  n_defined_pixels=0, max_abs_sum_deviation=None,
                                  reason="NEO invalid"))
            continue
        tot = np.zeros(X0.shape, dtype=np.float64)
        for p in present:
            tot += arrays[f"density__{p}__{lab}"].astype(np.float64)
        ok = np.isfinite(tot) & (tot > 0)
        P = np.full(X0.shape, np.nan)
        P[ok] = arrays[f"density__NEO__{lab}"].astype(np.float64)[ok] / tot[ok]
        ssum = np.zeros(int(ok.sum()))
        for p in present:
            ssum += arrays[f"density__{p}__{lab}"].astype(np.float64)[ok] / tot[ok]
        maxdev = float(np.max(np.abs(ssum - 1.0))) if ok.any() else float("nan")
        if ok.any() and not (maxdev < 1e-12):
            raise RuntimeError(f"{label}/{lab}: posterior does not sum to 1 (max dev {maxdev:.3e})")
        arrays[f"P_NEO__{lab}"] = P.astype(np.float32)
        post_meta.append(dict(center=label, magnitude_bin=lab, defined=True,
                              populations_used=present, populations_invalid=missing,
                              partial_denominator=bool(missing),
                              n_defined_pixels=int(ok.sum()),
                              max_abs_sum_deviation=maxdev, reason="ok"))

    # ---- write atomically ----------------------------------------------------
    meta = dict(
        center_label=label, center_lon_deg=clon, center_lat_deg=clat, epoch=epoch,
        max_sep_deg=MAX_SEP_DEG, grid_lim=list(GRID_LIM), grid_step=GRID_STEP,
        k_by_population=K_BY_POP, gaussian_smoothing=SMOOTHING, support_masking=False,
        density_engine=a.density_engine, magnitude_quantity="apparent V (HG, G=0.15)",
        bin_scheme=dict(lo=MAG_LO, hi=MAG_HI, step=MAG_STEP, n_bins=len(mag025_bins()),
                        semantics="half-open [lo, lo+step)"),
        neo_source=neo_src_desc, neo_source_sha256=neo_src_hash,
        neo_effective_factor=eff_neo,
        split_provenance_mag025_sha256=sha256_file(SPLIT_MAG025),
        sealed_module_sha256=sha256_file(W / "neomod/src/velocity_density_pipeline_neomod_clone_only.py"),
        code_commit=git_commit(), build_seconds=round(time.time() - t_center0, 1),
    )
    tmp = final.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, x_grid=x_grid, y_grid=y_grid,
                        meta_json=json.dumps(meta),
                        coverage_json=json.dumps(coverage),
                        integrals_json=json.dumps(integrals),
                        posterior_json=json.dumps(post_meta), **arrays)
    z = np.load(tmp, allow_pickle=True)           # open-and-read before promoting
    _ = z["x_grid"], z["meta_json"]
    for kk in list(arrays)[:3]:
        _ = z[kk]
    tmp.replace(final)
    pd.DataFrame(coverage).to_parquet(outdir / f"coverage_{label}.parquet", index=False)
    pd.DataFrame(integrals).to_parquet(outdir / f"integrals_{label}.parquet", index=False)
    marker.write_text(json.dumps({
        "center": label, "sha256": sha256_file(final),
        "bytes": final.stat().st_size, "n_arrays": len(arrays),
        "n_valid_cells": int(sum(1 for c in coverage if c["valid"])),
        "n_cells": len(coverage), "build_seconds": meta["build_seconds"],
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, indent=2))
    print(f"[{label}] DONE {final.stat().st_size/1e6:.1f} MB in {meta['build_seconds']:.1f}s "
          f"({sum(1 for c in coverage if c['valid'])}/{len(coverage)} cells valid)", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--task-id", type=int, default=0)
    p.add_argument("--center", default=None)
    p.add_argument("--n-jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    p.add_argument("--density-engine", choices=("closed_form", "sealed_quadrature"),
                   default="closed_form")
    p.add_argument("--neo-source", choices=("by_pixel", "high_cell"), default="by_pixel")
    p.add_argument("--neo-source-path", default=str(
        W / "outputs/more_neomod_samples_knn/source_high.parquet"))
    p.add_argument("--neo-effective-factor", type=float, default=64.725384)
    p.add_argument("--only-bins", default=None, help="comma-separated bin labels")
    p.add_argument("--out-dir", default=None)
    p.add_argument("--overwrite", action="store_true")
    a = p.parse_args()
    build_one_center(a)


if __name__ == "__main__":
    main()
