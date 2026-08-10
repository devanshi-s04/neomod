#!/usr/bin/env python3
"""Build ONE ablation map (all magnitude bins) for one sky center and one density mode.

Mirrors the sealed driver `sorcha_gen_maps_grid.py` for loading, cuts, split corrections, weights,
grid and provenance. The ONLY difference is the density backend, dispatched explicitly.

Per-(population, magnitude-bin) build calls
-------------------------------------------
Every mode runs ONE `generate_probability_maps` call per (population, magnitude bin), passing a
single-entry `clone_sources` AND a single-entry `mag_bins`. Both identities are therefore explicit
call arguments; nothing is inferred from callback ordering. This is what makes the physical
normalization audit per (population, map, bin) well defined.

Merge rule (assertion 1)
------------------------
Only RAW per-population density arrays and their diagnostics are merged. Partial posteriors are
never averaged or concatenated. After all four populations are assembled the posterior is
recomputed from scratch and verified against an independent second computation.

Invalid densities are PRESERVED (assertion 2)
---------------------------------------------
`nan_to_num` is never applied to a population density. A missing or non-finite population density
makes the posterior NaN with a reason code; it is not silently treated as zero.

Shared metadata (assertion 3)
-----------------------------
Every common metadata key is asserted identical across all four partial archives before being
copied. `density_raw == density_unsmoothed` is asserted for every population and bin, since
smooth_density_maps=False must make them identical.

Integral ratios are reported WITHOUT attributing a cause (assertion 4).
"""
from __future__ import annotations
import argparse, hashlib, json, os, re, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
sys.path.insert(0, str(W / "neomod" / "adam_core_stub"))

POPS = ("NEO", "MBA", "TNO", "Trojans")
MAX_SEP = 30.0
# reason codes for an invalid posterior pixel
R_OK, R_NONFINITE, R_NEGATIVE, R_ZEROTOTAL = 0, 1, 2, 3


def sha256(p, ch=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(ch), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["hist_all", "bayes_all", "bayes_nonneo"])
    ap.add_argument("--center", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n-jobs", type=int, default=8)
    a = ap.parse_args()
    # Resolve BEFORE chdir: the sealed module requires cwd=neomod/, so a relative --out-dir
    # would silently land under neomod/ instead of the repo root.
    a.out_dir = str(Path(a.out_dir).resolve())
    os.chdir(W / "neomod")
    import velocity_density_pipeline_neomod_clone_only as base
    import velocity_density_pipeline_neomod_density_ablation as abl

    OUT = Path(a.out_dir); OUT.mkdir(parents=True, exist_ok=True)
    final = OUT / f"prob_maps_grid_{a.center}.npz"
    if final.exists():
        print(f"[skip] {final} exists"); return
    m = re.match(r"dlon([+-]\d+)_lat([+-]\d+)", a.center)
    dlon, lat = float(m.group(1)), float(m.group(2))
    seal = json.load(open(W / "outputs/splits/MAP_BUILD_SEAL.json"))
    ref_obstime = seal["grid"]["ref_obstime"]
    MAGBINS = [{"label": l, "mag_min": float(lo), "mag_max": float(hi)}
               for l, lo, hi in seal["magnitude_bins"]]
    # The velocity domain is NOT recorded in MAP_BUILD_SEAL; the production driver passed it via
    # --velocity-grid-limit. Read it from a PRODUCTION map so the ablation cannot drift from the
    # frozen domain. Module defaults are (-2,2), which would silently build a different experiment.
    _pm = W / seal["output_dir"] / f"prob_maps_grid_{a.center}.npz"
    with np.load(_pm, allow_pickle=True) as _z:
        _xg = np.asarray(_z["x_grid"], float)
    GRID_LIM = (float(_xg.min()), float(_xg.max()))
    GRID_STEP = float(round(_xg[1] - _xg[0], 12))
    if len(_xg) != 1001 or abs(GRID_STEP - 0.01) > 1e-12:
        raise AssertionError(f"unexpected production grid: n={len(_xg)} step={GRID_STEP}")
    print(f"[grid] from production map: lim={GRID_LIM} step={GRID_STEP} n={len(_xg)}", flush=True)

    from astropy.time import Time
    from astropy.coordinates import get_sun, GeocentricTrueEcliptic
    t = Time(ref_obstime, scale="utc")
    antisun = (get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t)).lon.deg + 180.0) % 360.0
    center_lon = (antisun + dlon) % 360.0

    base.NONNEO_SPLIT_FRACTIONS = json.load(open(W / "outputs/splits/split_provenance.json"))
    neo_meta = json.load(open(W / "outputs/neomod3_projection_cache/cache_metadata.json"))

    cache = pd.read_parquet(W / "outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet")
    man = pd.read_parquet(W / "outputs/splits/nonneo_split_manifest.parquet")
    # EXACTLY the sealed driver's split rule: GEN non-NEO rows, NEO passes through because it is
    # replaced by the NEOMOD3 GEN cache inside the builder.
    keep = set(man.ObjID[man.split == "GEN"])
    n0 = len(cache)
    cache = cache[cache.ObjID.isin(keep) | (cache.population == "NEO")].reset_index(drop=True)
    print(f"[split] GEN: kept {len(cache):,} of {n0:,} cache rows", flush=True)
    # EXACTLY the sealed driver's population settings and shared observer-only scorer.
    import copy as _copy
    pop_settings = _copy.deepcopy(base.DEFAULT_POPULATION_SETTINGS)
    _, shared_scorer = base.load_s3m_population("neo", verbose=False)
    clone_sources = {}
    for pop, cfg in pop_settings.items():
        sub = cache[cache["population"] == pop].reset_index(drop=True)
        clone_sources[pop] = {
            "df": sub, "scorer": shared_scorer,
            "clone_factor": cfg["clone_factor"],
            "use_conditional_cloner": cfg.get("use_conditional_cloner", True),
            "scatter_size": cfg.get("scatter_size", 4),
            "scatter_alpha": cfg.get("scatter_alpha", 0.1),
            "_mag_app": sub["mag_app"].to_numpy(dtype=float),
        }
    print(f"[cache] " + "  ".join(f"{p}={len(clone_sources[p]['df']):,}" for p in POPS), flush=True)

    def effective_factor_for(pop, mb):
        """Explicit reconstruction of the sealed builder's scalar effective_factor."""
        if pop == "NEO":
            for k in ("effective_factor_NEO", "effective_factor"):
                if k in neo_meta:
                    return float(neo_meta[k])
            return float(neo_meta["n_draws"]) / float(neo_meta["total_weight_absolute_NEO_count"])
        return 1.0 * float(base._nonneo_split_fraction(pop, mb["mag_min"], mb["mag_max"]))

    accounting = []
    skipped = []          # (population, bin) combinations with no samples -- builder wrote zeros
    orig = base.evaluate_density_map_full_posterior_2d
    partials = {}
    t0 = time.time()

    # ---------------- one call per (population, magnitude bin): both identities EXPLICIT --------
    for pop in POPS:
        backend = abl.resolve_backend(a.mode, pop)
        for mb in MAGBINS:
            # effective_factor is resolved LAZILY inside the wrapper. Some (population, bin)
            # combinations have no objects at all -- e.g. TNO/14_16 has no split-provenance
            # fraction because no TNO is that bright -- and the sealed builder returns zeros
            # early without ever computing a factor. Resolving eagerly would raise on bins the
            # production build legitimately skips.
            ctx = {"center": a.center, "mode": a.mode, "population": pop,
                   "magnitude_bin": mb["label"], "backend": backend}

            def _wrapper(tree, grid_points, k=None, n_d0_grid=None, show_progress=False,
                         n_jobs=1, _ctx=ctx, _backend=backend, _mb=mb, _pop=pop):
                _eff = effective_factor_for(_pop, _mb)
                _ctx = dict(_ctx)
                _ctx["effective_factor"] = _eff
                _ctx["physical_weight_per_sample"] = 1.0 / _eff
                pts = np.asarray(tree.data, float)
                gp = np.asarray(grid_points, float)
                xs = np.unique(gp[:, 0]); ys = np.unique(gp[:, 1])
                xe, ye, pixel_area = abl.grid_edges(xs, ys)
                inside = ((pts[:, 0] >= xe[0]) & (pts[:, 0] <= xe[-1]) &
                          (pts[:, 1] >= ye[0]) & (pts[:, 1] <= ye[-1]))
                w = 1.0 / _eff
                rec = dict(_ctx)
                rec.update({
                    "n_input_samples": int(len(pts)),
                    "n_inside_edges": int(inside.sum()),
                    "n_outside_edges": int((~inside).sum()),
                    "physical_weight_inside": float(inside.sum() * w),
                    "physical_weight_outside": float((~inside).sum() * w),
                    "expected_in_domain_physical_count": float(inside.sum() * w),
                    "pixel_area": float(pixel_area),
                })
                if _backend == "hist":
                    rho = abl.histogram_density(pts, xs, ys)
                    ib = float(rho.sum() * pixel_area)
                    out = rho.reshape(-1)
                else:
                    out = abl.BAYES_ROUTINE(tree, grid_points, k=k, n_d0_grid=n_d0_grid,
                                            show_progress=show_progress, n_jobs=n_jobs)
                    ib = float(np.nansum(out) * pixel_area)
                    rec["delegated_to_sealed"] = bool(
                        abl.BAYES_ROUTINE is base.evaluate_density_map_full_posterior_2d)
                rec["integral_before_physical_normalization"] = ib
                rec["integral_after_physical_normalization"] = ib / _eff
                exp = rec["expected_in_domain_physical_count"]
                rec["ratio_after_norm_to_expected"] = (ib / _eff / exp) if exp > 0 else np.nan
                accounting.append(rec)
                return out

            tmp = OUT / f"_partial_{a.center}_{pop}_{mb['label']}.npz"
            base.evaluate_density_map_full_posterior_2d = _wrapper
            try:
                base.generate_probability_maps(
                    obstime_str=ref_obstime, output_path=str(tmp),
                    center_lon_deg=center_lon, center_lat_deg=lat, center_label=a.center,
                    max_sep_deg=MAX_SEP, n_jobs=a.n_jobs, save_overlays=False,
                    grid_lim=GRID_LIM, grid_step=GRID_STEP,
                    smooth_density_maps=False,
                    clone_sources={pop: clone_sources[pop]},
                    mag_bins=[mb],
                )
            finally:
                base.evaluate_density_map_full_posterior_2d = orig
            partials[(pop, mb["label"])] = tmp
            if not any(r["population"] == pop and r["magnitude_bin"] == mb["label"]
                       for r in accounting):
                skipped.append({"center": a.center, "mode": a.mode, "population": pop,
                                "magnitude_bin": mb["label"], "backend": backend,
                                "reason": "no samples in bin; sealed builder wrote zeros"})
        print(f"[{a.mode}] {pop} via {backend} done ({time.time()-t0:.0f}s)", flush=True)

    # ---------------- ASSERTION 3: inventory keys, verify shared metadata agreement -------------
    DENS_PFX = ("density_raw__", "density_unsmoothed__", "support_count__", "nearest_dist__")
    inventories, metas, merged = {}, {}, {}
    for key, p in partials.items():
        z = np.load(p, allow_pickle=True)
        inventories[key] = list(z.files)
        meta_keys = [k for k in z.files if not any(k.startswith(x) for x in DENS_PFX)]
        metas[key] = {k: z[k] for k in meta_keys}
        z.close()

    common = set.intersection(*(set(v.keys()) for v in metas.values()))
    disagree = []
    for k in sorted(common):
        vals = [metas[key][k] for key in metas]
        ref = vals[0]
        for v in vals[1:]:
            same = (np.array_equal(ref, v) if isinstance(ref, np.ndarray)
                    else bool(np.all(ref == v)))
            if not same:
                disagree.append(k); break
    # EXACT key names. Substring matching is unsafe: "lat" is a substring of "popuLATion_names",
    # which legitimately differs between single-population partials.
    STRICT_KEYS = {"x_grid", "y_grid", "center_lon_deg", "center_lat_deg", "center_label",
                   "obstime_str", "grid_lim", "grid_step", "max_sep_deg",
                   "delta_lon_from_antisun_deg", "grid_lat_deg", "ref_obstime"}
    STRICT = sorted(STRICT_KEYS & common)
    EXPECT_DIFFER = {"population_names"}   # one population per partial, by construction
    strict_bad = [k for k in STRICT if k in disagree]
    if strict_bad:
        raise AssertionError(f"shared metadata DISAGREES across partial archives: {strict_bad}")
    print(f"[assert] {len(common)} common metadata keys; {len(STRICT)} strict keys "
          f"({', '.join(STRICT)}) agree across all {len(partials)} partials  OK", flush=True)
    unexpected = sorted(set(disagree) - EXPECT_DIFFER - set(STRICT))
    if unexpected:
        print(f"[note] metadata keys differing across partials (per-population, not copied): "
              f"{unexpected}", flush=True)
    for k in sorted(common):
        if k not in disagree:
            merged[k] = metas[list(metas)[0]][k]

    for (pop, lab), p in partials.items():
        z = np.load(p, allow_pickle=True)
        for k in z.files:
            if any(k.startswith(x) for x in DENS_PFX) and f"__{pop}__{lab}" in k:
                if k in merged:
                    raise AssertionError(f"duplicate density array {k}")
                merged[k] = z[k]
        z.close()

    bins = [b["label"] for b in MAGBINS]
    for b in bins:
        present = [p for p in POPS if f"density_raw__{p}__{b}" in merged]
        if len(present) != 4:
            raise AssertionError(f"bin {b}: expected 4 population densities, got {present}")
    print(f"[assert] exactly 4 density arrays per bin for all {len(bins)} bins  OK", flush=True)

    # smoothing disabled => raw and unsmoothed must be identical
    for b in bins:
        for p in POPS:
            r = merged.get(f"density_raw__{p}__{b}"); u = merged.get(f"density_unsmoothed__{p}__{b}")
            if u is None:
                continue
            # The sealed archive stores density_raw as float32 and density_unsmoothed as float64,
            # so bit-equality is impossible by construction. With smoothing disabled they must
            # agree to float32 precision: raw == float32(unsmoothed), NaN positions identical.
            ra = np.asarray(r); ua = np.asarray(u)
            if not np.array_equal(np.isnan(ra), np.isnan(ua)):
                raise AssertionError(f"NaN pattern differs for {p}/{b}")
            if not np.array_equal(ra, ua.astype(np.float32), equal_nan=True):
                fin = np.isfinite(ra) & np.isfinite(ua)
                mx = float(np.max(np.abs(ra[fin].astype(np.float64) - ua[fin]))) if fin.any() else 0.0
                raise AssertionError(f"density_raw != float32(density_unsmoothed) for {p}/{b} "
                                     f"(max|diff| {mx:.3e}): smoothing or rebinding occurred")
    print("[assert] density_raw == density_unsmoothed for every population and bin  OK", flush=True)

    # ---------------- ASSERTION 2: preserve invalid densities, explicit reason codes ------------
    for b in bins:
        arrs = {p: np.asarray(merged[f"density_raw__{p}__{b}"], float) for p in POPS}
        finite = np.ones_like(arrs["NEO"], dtype=bool)
        nonneg = np.ones_like(arrs["NEO"], dtype=bool)
        for p in POPS:
            finite &= np.isfinite(arrs[p])
            nonneg &= (arrs[p] >= 0) | ~np.isfinite(arrs[p])
        rho_total = arrs["NEO"] + arrs["MBA"] + arrs["TNO"] + arrs["Trojans"]   # NO nan_to_num
        valid = finite & nonneg & np.isfinite(rho_total) & (rho_total > 0)
        P = np.full(arrs["NEO"].shape, np.nan)
        np.divide(arrs["NEO"], rho_total, out=P, where=valid)
        reason = np.full(arrs["NEO"].shape, R_OK, dtype=np.int8)
        reason[~finite] = R_NONFINITE
        reason[finite & ~nonneg] = R_NEGATIVE
        reason[finite & nonneg & ~(np.isfinite(rho_total) & (rho_total > 0))] = R_ZEROTOTAL
        # independent recomputation must agree exactly on the valid set
        P2 = np.full(arrs["NEO"].shape, np.nan)
        tot2 = sum(np.asarray(merged[f"density_raw__{p}__{b}"], float) for p in POPS)
        np.divide(np.asarray(merged[f"density_raw__NEO__{b}"], float), tot2, out=P2, where=valid)
        if not np.array_equal(P[valid], P2[valid]):
            raise AssertionError(f"posterior for bin {b} is not a fresh normalization")
        if np.any(np.isfinite(P) & ~valid):
            raise AssertionError(f"bin {b}: finite posterior outside the valid mask")
        merged[f"P_NEO__{b}"] = P
        merged[f"rho_total__{b}"] = rho_total
        merged[f"posterior_valid__{b}"] = valid
        merged[f"posterior_reason__{b}"] = reason
    print("[assert] posterior recomputed fresh; invalid densities preserved as NaN + reason codes  OK",
          flush=True)

    # Explicit grid provenance in the ARCHIVE itself: downstream stages must never have to
    # reconstruct the velocity domain from a production map.
    merged["ablation_grid_lim"] = np.asarray(GRID_LIM, dtype=np.float64)
    merged["ablation_grid_step"] = np.float64(GRID_STEP)
    merged["ablation_n_grid"] = np.int64(len(merged["x_grid"]))
    merged["ablation_density_mode"] = np.asarray(a.mode)
    merged["ablation_backend_by_population"] = np.asarray(
        [f"{p_}={abl.resolve_backend(a.mode, p_)}" for p_ in POPS])
    merged["ablation_gaussian_smoothing"] = np.asarray("disabled")
    np.savez_compressed(final, **merged)
    for p in partials.values():
        Path(p).unlink(missing_ok=True)

    acc = pd.DataFrame(accounting)
    if len(acc):
        hist = acc[acc.backend == "hist"]
        bad = hist[(hist.expected_in_domain_physical_count > 0) &
                   (~np.isclose(hist.ratio_after_norm_to_expected, 1.0, rtol=1e-12, atol=0))]
        if len(bad):
            raise AssertionError(f"histogram ratio_after_norm_to_expected != 1 for {len(bad)} rows:\n"
                                 f"{bad[['population','magnitude_bin','ratio_after_norm_to_expected']]}")
        print(f"[assert] histogram ratio == 1 within 1e-12 for all {len(hist)} histogram rows  OK",
              flush=True)
    acc.to_parquet(OUT / f"domain_accounting_{a.center}.parquet", index=False)
    if skipped:
        pd.DataFrame(skipped).to_parquet(OUT / f"skipped_bins_{a.center}.parquet", index=False)
        print(f"[note] {len(skipped)} (population, bin) combinations had no samples: "
              + ", ".join(f"{r['population']}/{r['magnitude_bin']}" for r in skipped), flush=True)

    (OUT / f"manifest_{a.center}.json").write_text(json.dumps({
        "mode": a.mode, "center": a.center, "center_lon_deg": center_lon, "center_lat_deg": lat,
        "grid": {"grid_lim": list(GRID_LIM), "grid_step": GRID_STEP,
                 "n_grid": int(len(merged["x_grid"])),
                 "source": "read from the production map for this center and asserted "
                           "n=1001, step=0.01; NOT the module default (-2,2)"},
        "dlon_from_antisun_deg": dlon, "ref_obstime": ref_obstime,
        "gaussian_smoothing": {"enabled": False, "argument": "smooth_density_maps=False",
                               "verified": "density_raw == density_unsmoothed asserted per pop/bin"},
        "dispatch": {p: abl.resolve_backend(a.mode, p) for p in POPS},
        "bayes_routine_is_sealed": bool(abl.BAYES_ROUTINE is base.evaluate_density_map_full_posterior_2d),
        "density_units": "objects per (deg/day)^2 after division by effective_factor",
        "normalization": "density_estimate / effective_factor (scalar per population, mag bin)",
        "interpolation": "not applied at build time; density-first interpolation happens at scoring",
        "posterior_rule": "rho_NEO/(sum of four raw densities) where all finite, >=0 and total>0; "
                          "else NaN with reason code",
        "reason_codes": {"0": "valid", "1": "non-finite density", "2": "negative density",
                         "3": "zero or non-finite total"},
        "merge_rule": "raw per-population densities only; partial posteriors never merged",
        "n_partial_builds": len(partials),
        "source_hashes": {
            "ablation_module": sha256(W/"neomod/src/velocity_density_pipeline_neomod_density_ablation.py"),
            "sealed_module": sha256(W/"neomod/src/velocity_density_pipeline_neomod_clone_only.py"),
            "runner": sha256(W/"neomod/pipeline/run_density_ablation.py")},
        "input_hashes": {
            "epoch_cache": sha256(W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet"),
            "split_manifest": sha256(W/"outputs/splits/nonneo_split_manifest.parquet"),
            "split_provenance": sha256(W/"outputs/splits/split_provenance.json"),
            "map_build_seal": sha256(W/"outputs/splits/MAP_BUILD_SEAL.json")},
        "output_sha256": sha256(final), "elapsed_s": time.time()-t0,
    }, indent=2))
    print(f"[done] {final}  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
