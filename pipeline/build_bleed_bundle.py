#!/usr/bin/env python3
"""Build the VDP high-velocity kNN-bleed / asymmetric-support-mask diagnostic bundle.

Explains why the 28 out-of-±2 night-61642 NEOs recover to P_NEO=1.0 under the ±5 pilot maps:
NEO clones do not reach |v|>2 (support ends ~|v|<=2), the raw density there is kNN bleed, and the
support mask zeros MBA/TNO/Trojan (support<1) but SKIPS the smoothed NEO population, so NEO wins the
ratio by default. All quantities are derived from the stored map arrays (density_raw, support_count,
nearest_dist) + the Stage-1 deterministic parquets — no --save-overlays clone coords needed, and NO
further map regeneration.

Output: outputs/vdp_bleed_diagnostic_bundle/  (+ .tar.gz, excluding the ±5 map directory).
"""
from __future__ import annotations
import os, sys, glob, json, hashlib, subprocess, tarfile, time
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
PILOT = W / "prob_maps_grid_s3m_nbody_vlim5_pilot"
PROD = W / "prob_maps_grid_s3m_nbody"
S1 = W / "outputs" / "mag245_nbody_deterministic_rescore"
OUT = W / "outputs" / "vdp_bleed_diagnostic_bundle"
NIGHT = 61642
GRID_LIM_PROD = 2.0
STEP = 0.01
RADII = [0.05, 0.10, 0.15, 0.25]
K_TH = 10  # kth-nearest (DEFAULT_K_MAP)
NJOBS = int(os.environ.get("BUNDLE_JOBS", "8"))
POPS = ["MBA", "NEO", "TNO", "Trojans"]
FIVE = ["S0000WD6a", "S0000WxWa", "S00006EPa", "S00002Qma", "S0000NWGa"]
import velocity_density_pipeline_gmm as vdp


# ---------------------------------------------------------------- load the pilot-center objects
def load_objects():
    b = pd.read_parquet(S1 / "bench_deterministic.parquet", columns=[
        "s3m_objid", "population", "vlam", "vbeta", "prob_map_file", "mean_mag", "P_NEO_vdp", "P_NEO_d2_det"])
    b["s3m_objid"] = b.s3m_objid.astype(str)
    b = b.rename(columns={"s3m_objid": "ObjID", "mean_mag": "mag", "P_NEO_vdp": "vdp_prod"}); b["pipeline"] = "benchmark"
    s = pd.read_parquet(S1 / "sorcha_vcorr_deterministic.parquet", columns=[
        "ObjID", "night", "population", "vlam", "vbeta", "prob_map_file", "mean_mag_V", "P_NEO_vdp_Vband", "P_NEO_d2_det"])
    s["ObjID"] = s.ObjID.astype(str); s = s[s.night == NIGHT].copy()
    s = s.rename(columns={"mean_mag_V": "mag", "P_NEO_vdp_Vband": "vdp_prod"}); s["pipeline"] = "sorcha_vcorr"
    cols = ["ObjID", "pipeline", "population", "vlam", "vbeta", "prob_map_file", "mag", "vdp_prod", "P_NEO_d2_det"]
    return pd.concat([b[cols], s[cols]], ignore_index=True)


# ---------------------------------------------------------------- per-center scoring (parallel)
def score_center(cen, sub):
    """Return sub with ±5 masked and unmasked P_NEO added."""
    pm5m = vdp.ProbMapSet.from_npz(str(PILOT / cen), support_mask_min=1, mask_radius_deg_per_day=np.inf)
    pm5u = vdp.ProbMapSet.from_npz(str(PILOT / cen), support_mask_min=None, mask_radius_deg_per_day=np.inf)
    vl = sub.vlam.to_numpy(float); vb = sub.vbeta.to_numpy(float); mg = sub.mag.to_numpy(float)
    sub = sub.copy()
    sub["vdp_pilot5_masked"] = pm5m.score_visible(vl, vb, mg)["NEO"]
    sub["vdp_pilot5_unmasked"] = pm5u.score_visible(vl, vb, mg)["NEO"]
    return sub


# ---------------------------------------------------------------- per-cell / neighbourhood diagnostics
def magbin_of(z, mag):
    for k in range(len(z["mag_bin_labels"])):
        if z["mag_bin_mins"][k] <= mag < z["mag_bin_maxs"][k]:
            return str(z["mag_bin_labels"][k])
    return None


def cell_index(z, vlam, vbeta):
    xg, yg = z["x_grid"], z["y_grid"]
    ix = int(np.clip(round((vlam - xg[0]) / (xg[1] - xg[0])), 0, len(xg) - 1))
    iy = int(np.clip(round((vbeta - yg[0]) / (yg[1] - yg[0])), 0, len(yg) - 1))
    return ix, iy


def target_diagnostics(z, mb, vlam, vbeta):
    """All per-cell + neighbourhood quantities for one target, from the pilot NPZ arrays."""
    xg, yg = z["x_grid"], z["y_grid"]
    ix, iy = cell_index(z, vlam, vbeta)
    d = dict(cell_ix=ix, cell_iy=iy)
    dens = {p: z[f"density_raw__{p}__{mb}"] for p in POPS}
    sup = {p: z[f"support_count__{p}__{mb}"] for p in POPS}
    near = {p: z[f"nearest_dist__{p}__{mb}"] for p in POPS}
    # masked density: skip NEO (smoothed), others zeroed where support<1
    masked = {p: (dens[p] if p == "NEO" else np.where(sup[p] >= 1, dens[p], 0.0)) for p in POPS}
    for p in POPS:
        d[f"raw_density_{p}"] = float(dens[p][iy, ix])
        d[f"support_{p}"] = float(sup[p][iy, ix])
        d[f"masked_density_{p}"] = float(masked[p][iy, ix])
        d[f"nearest_clone_dist_{p}"] = float(near[p][iy, ix])
    # local clone counts (sum support within radius) + kth-nearest, per pop, from the support grid
    xx, yy = np.meshgrid(xg, yg)
    dist = np.hypot(xx - vlam, yy - vbeta)
    for p in POPS:
        s = sup[p]
        for R in RADII:
            d[f"clones_r{int(R*100):02d}_{p}"] = float(s[dist <= R].sum())
        hit = s >= 1
        if hit.any():
            dd = np.sort(dist[hit]); ss = s[hit][np.argsort(dist[hit])]
            cum = np.cumsum(ss)
            d[f"kth_nearest_dist_{p}"] = float(dd[np.searchsorted(cum, K_TH)]) if cum[-1] >= K_TH else float("nan")
        else:
            d[f"kth_nearest_dist_{p}"] = float("nan")
    # max |v| reached by genuinely supported NEO clones
    hitN = sup["NEO"] >= 1
    vmag = np.maximum(np.abs(xx), np.abs(yy))
    d["neo_supported_vmax"] = float(vmag[hitN].max()) if hitN.any() else 0.0
    d["genuine_support_at_query"] = bool(sup["NEO"][iy, ix] >= 1)
    return d


def main():
    t0 = time.time()
    (OUT / "arrays").mkdir(parents=True, exist_ok=True)
    obj = load_objects()
    obj["maxabs"] = np.maximum(obj.vlam.abs(), obj.vbeta.abs())
    obj["total_speed"] = np.hypot(obj.vlam, obj.vbeta)
    centers = [os.path.basename(f) for f in sorted(glob.glob(str(PILOT / "*.npz")))]
    inp = obj[obj.prob_map_file.isin(centers)].reset_index(drop=True)
    # targets = benchmark-NEO matched, night61642, |v|>2 in either pipeline
    b = obj[obj.pipeline == "benchmark"]; bneo = set(b[b.population == "NEO"].ObjID)
    s = obj[obj.pipeline == "sorcha_vcorr"]; matched = s[s.ObjID.isin(bneo)]
    tgt_ids = set(matched[matched.maxabs > GRID_LIM_PROD].ObjID) | set(
        b[b.ObjID.isin(set(matched.ObjID)) & (b.maxabs > GRID_LIM_PROD)].ObjID)
    print(f"objects at {len(centers)} centers: {len(inp):,} | targets: {len(tgt_ids)}", flush=True)

    # ---- parallel ±5 masked+unmasked scoring for every object ----
    groups = [(cen, inp[inp.prob_map_file == cen]) for cen in centers]
    scored = pd.concat(Parallel(n_jobs=NJOBS, verbose=5)(delayed(score_center)(c, g) for c, g in groups), ignore_index=True)
    scored["out_of_grid_production"] = scored.maxabs > GRID_LIM_PROD
    scored["is_target"] = scored.ObjID.isin(tgt_ids)

    # ---- ITEM 1: targets_28.parquet ----
    tgt = scored[scored.is_target].reset_index(drop=True)
    zc = {}
    def zof(c):
        if c not in zc: zc[c] = np.load(PILOT / c, allow_pickle=True)
        return zc[c]
    diag_rows = []
    for _, r in tgt.iterrows():
        z = zof(r.prob_map_file); mb = magbin_of(z, r.mag)
        base = dict(ObjID=r.ObjID, pipeline=r.pipeline, map_file=r.prob_map_file, mag_bin=mb,
                    vlam=r.vlam, vbeta=r.vbeta, total_speed=r.total_speed, magnitude=r.mag,
                    vdp_prod_pm2=r.vdp_prod, vdp_pilot5_masked=r.vdp_pilot5_masked,
                    vdp_pilot5_unmasked=r.vdp_pilot5_unmasked, digest2_det=r.P_NEO_d2_det,
                    out_of_grid_production=bool(r.out_of_grid_production))
        base.update(target_diagnostics(z, mb, r.vlam, r.vbeta) if mb else {})
        diag_rows.append(base)
    targets28 = pd.DataFrame(diag_rows)
    targets28.to_parquet(OUT / "targets_28.parquet", index=False)

    # ---- ITEM 2 + 3: clone support points + diagnostic arrays for the 5 cases ----
    csp_rows = []
    arrays = {}
    for i, oid in enumerate(FIVE):
        row = tgt[(tgt.ObjID == oid)].sort_values("pipeline").iloc[0]  # prefer benchmark (sorts first)
        cen = row.prob_map_file; z = zof(cen); mb = magbin_of(z, row.mag)
        xg, yg = z["x_grid"], z["y_grid"]; ix, iy = cell_index(z, row.vlam, row.vbeta)
        pm_m = vdp.ProbMapSet.from_npz(str(PILOT / cen), support_mask_min=1, mask_radius_deg_per_day=np.inf)
        pm_u = vdp.ProbMapSet.from_npz(str(PILOT / cen), support_mask_min=None, mask_radius_deg_per_day=np.inf)
        pre = f"case{i}_{oid}"
        arrays[f"{pre}__vlam_grid"] = xg.astype(np.float32)
        arrays[f"{pre}__vbeta_grid"] = yg.astype(np.float32)
        arrays[f"{pre}__masked_P_NEO"] = pm_m.get_probability_map(mb, "NEO").astype(np.float32)
        arrays[f"{pre}__unmasked_P_NEO"] = pm_u.get_probability_map(mb, "NEO").astype(np.float32)
        arrays[f"{pre}__query"] = np.array([row.vlam, row.vbeta], np.float32)
        arrays[f"{pre}__query_cell_ixiy"] = np.array([ix, iy], np.int32)
        arrays[f"{pre}__meta"] = np.array([oid, row.pipeline, cen, mb], dtype="<U40")
        for p in POPS:
            dens = z[f"density_raw__{p}__{mb}"].astype(np.float32); sc = z[f"support_count__{p}__{mb}"].astype(np.float32)
            arrays[f"{pre}__density_raw_{p}"] = dens
            arrays[f"{pre}__support_{p}"] = sc
            arrays[f"{pre}__masked_density_{p}"] = (dens if p == "NEO" else np.where(sc >= 1, dens, 0.0)).astype(np.float32)
            # clone support points (cell centres with support>=1)
            yy, xx = np.where(sc >= 1)
            for jy, jx in zip(yy, xx):
                csp_rows.append(dict(case=pre, ObjID=oid, population=p, clone_vlam=float(xg[jx]),
                                     clone_vbeta=float(yg[jy]), clone_weight_support=float(sc[jy, jx]),
                                     map_file=cen, mag_bin=mb))
        # production ±2 P_NEO for the same center+bin (overlapping domain)
        pcen = PROD / cen
        if pcen.exists():
            pm2 = vdp.ProbMapSet.from_npz(str(pcen), support_mask_min=1, mask_radius_deg_per_day=np.inf)
            z2 = np.load(pcen, allow_pickle=True)
            arrays[f"{pre}__prod_vlam_grid"] = z2["x_grid"].astype(np.float32)
            arrays[f"{pre}__prod_vbeta_grid"] = z2["y_grid"].astype(np.float32)
            arrays[f"{pre}__prod_P_NEO_pm2"] = pm2.get_probability_map(mb, "NEO").astype(np.float32)
        del pm_m, pm_u
    pd.DataFrame(csp_rows).to_parquet(OUT / "clone_support_points.parquet", index=False)
    np.savez_compressed(OUT / "arrays" / "diagnostic_map_arrays.npz", **arrays)

    # ---- ITEM 4: overlap_controls.parquet (|v|<1.9) ----
    ov = scored[scored.maxabs < 1.9].copy()
    ov_mb = ov.apply(lambda r: magbin_of(zof(r.prob_map_file), r.mag), axis=1)
    overlap = pd.DataFrame(dict(ObjID=ov.ObjID, population=ov.population, vlam=ov.vlam, vbeta=ov.vbeta,
                                prod_pm2=ov.vdp_prod, pilot5_masked=ov.vdp_pilot5_masked,
                                dP=ov.vdp_pilot5_masked - ov.vdp_prod, map_file=ov.prob_map_file, mag_bin=ov_mb.values))
    overlap.to_parquet(OUT / "overlap_controls.parquet", index=False)
    ov_median = float(overlap.dP.abs().median())

    # ---- ITEM 5: nonneo_controls.parquet ----
    nn = scored[scored.population != "NEO"].copy()
    nn_mb = nn.apply(lambda r: magbin_of(zof(r.prob_map_file), r.mag), axis=1)
    THR = 0.5
    nonneo = pd.DataFrame(dict(ObjID=nn.ObjID, population=nn.population, vlam=nn.vlam, vbeta=nn.vbeta,
                               prod_pm2=nn.vdp_prod, pilot5_masked=nn.vdp_pilot5_masked, pilot5_unmasked=nn.vdp_pilot5_unmasked,
                               class_prod=(nn.vdp_prod >= THR).astype(int), class_pilot=(nn.vdp_pilot5_masked >= THR).astype(int),
                               map_file=nn.prob_map_file, mag_bin=nn_mb.values))
    nonneo.to_parquet(OUT / "nonneo_controls.parquet", index=False)
    fp_prod = int(nonneo.class_prod.sum()); fp_pilot = int(nonneo.class_pilot.sum())

    # ---- ITEM 6: configuration_and_provenance.json ----
    def sha(p):
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for c in iter(lambda: fh.read(1 << 20), b""):
                h.update(c)
        return h.hexdigest()
    try:
        commit = subprocess.check_output(["git", "-C", str(W / "neomod"), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        commit = "unknown"
    probe_map = tgt.iloc[0].prob_map_file
    z0 = zof(probe_map)
    pm_probe = vdp.ProbMapSet.from_npz(str(PILOT / probe_map), support_mask_min=1, mask_radius_deg_per_day=np.inf)
    src = {p: sha(p) for p in [W/"neomod/pipeline/build_bleed_bundle.py", W/"neomod/pipeline/sorcha_gen_maps_grid.py",
           W/"neomod/src/velocity_density_pipeline_gmm.py", S1/"bench_deterministic.parquet",
           S1/"sorcha_vcorr_deterministic.parquet", PILOT/probe_map] if Path(p).exists()}
    cfg = dict(
        git_commit=commit, generated=time.strftime("%Y-%m-%dT%H:%M:%S"),
        map_generator_command=("python neomod/pipeline/sorcha_gen_maps_grid.py --delta-lon D --lat L "
                               "--ref-obstime 2027-08-25T00:00:00 --cache outputs/epoch_state_cache/"
                               "epoch_state_2027-08-25T000000.parquet --prob-maps-dir prob_maps_grid_s3m_nbody_vlim5_pilot "
                               "--mba-clone-factor 5 --velocity-grid-limit 5.0 --velocity-grid-step 0.01 --n-jobs 8 --overwrite"),
        pilot_grid=dict(limit_deg_day=5.0, step_deg_day=STEP, shape=[len(z0["x_grid"]), len(z0["y_grid"])]),
        production_grid=dict(limit_deg_day=GRID_LIM_PROD, step_deg_day=STEP, shape=[401, 401]),
        knn_params=dict(k_map=int(getattr(vdp, "DEFAULT_K_MAP", 10)), n_d0_grid=int(getattr(vdp, "DEFAULT_N_D0_GRID_MAP", 400)),
                        estimator="bayesian_knn_full_posterior_2d (Ivezic 2001 Appendix B)"),
        support_mask_min=1,
        support_mask_skip_populations=list(map(str, getattr(pm_probe, "_support_mask_skip", ["NEO"]))),
        smoothing=dict(smooth_density_maps=bool(z0["smooth_density_maps"]),
                       smooth_population_names=list(map(str, z0["smooth_population_names"])),
                       smooth_support_threshold=float(z0["smooth_support_threshold"]),
                       smooth_sigma_pixels=float(z0["smooth_sigma_pixels"])),
        clone_method=dict(NEO="GMM (log a,q,sin/cos angles,sin/cos M_obs; H sampled empirically)",
                          others="conditional K|M"),
        clone_factors=dict(NEO=80, MBA=5, TNO=10, Trojans=5),
        probability_formula="P(pop) = density_pop / sum_pops density_pop, per (vlam,vbeta) cell, bilinear-interpolated; out-of-grid -> 0",
        magnitude_bins=list(map(str, z0["mag_bin_labels"])),
        ref_obstime=str(z0["ref_obstime_str"]), populations=POPS,
        note="NEO is in support_mask_skip (it is the smoothed population), so support_mask_min zeros MBA/TNO/Trojan "
             "density where support<1 but leaves NEO's kNN-bleed density -> P_NEO defaults to 1.0 at |v|>2 where no "
             "population has genuine clone support. NEO clone support ends near |v|~2 in these maps.",
        source_file_sha256={str(Path(k).relative_to(W)): v for k, v in src.items()},
    )
    (OUT / "configuration_and_provenance.json").write_text(json.dumps(cfg, indent=2))

    # ---- ITEM 7: README.md + copy the build script into the bundle ----
    (OUT / "README.md").write_text(README.format(
        n_targets=targets28.ObjID.nunique(), n_rows=len(targets28), n_overlap=len(overlap),
        ov_median=ov_median, n_nonneo=len(nonneo), fp_prod=fp_prod, fp_pilot=fp_pilot, five=", ".join(FIVE)))
    import shutil
    shutil.copy(W / "neomod/pipeline/build_bleed_bundle.py", OUT / "build_bleed_bundle.py")

    # ---- validation ----
    n_rec = int(((tgt.vdp_prod < 1e-6) & (tgt.vdp_pilot5_masked > 0.5)).sum())
    print("\n=== VALIDATION ===", flush=True)
    print(f"  targets present: {targets28.ObjID.nunique()} unique ({len(targets28)} rows)  [expect 28 / 56]", flush=True)
    print(f"  targets recovered (prod~0 -> pilot>0.5): {n_rec}/{len(tgt)}", flush=True)
    print(f"  5 detailed cases with full arrays: {sum(1 for i in range(len(FIVE)) if any(k.startswith(f'case{i}_') for k in arrays))}/5", flush=True)
    print(f"  overlap median |dP| = {ov_median:.4f}  [expect ~0.0003]", flush=True)
    print(f"  non-NEO false-positives @0.5: prod {fp_prod} -> pilot {fp_pilot}  [expect 131 -> 140]", flush=True)

    # ---- tar.gz (exclude the ±5 map dir; bundle contains no maps) ----
    tarpath = W / "outputs" / "vdp_bleed_diagnostic_bundle.tar.gz"
    with tarfile.open(tarpath, "w:gz") as t:
        t.add(OUT, arcname="vdp_bleed_diagnostic_bundle")
    print(f"\nwrote bundle -> {OUT}\ntarball -> {tarpath} ({tarpath.stat().st_size/1e6:.1f} MB)  ({time.time()-t0:.0f}s)", flush=True)


README = """# VDP high-velocity kNN-bleed / support-mask diagnostic bundle

Explains why the 28 out-of-±2 night-61642 NEOs recover to P_NEO=1.0 under the ±5 velocity-grid pilot
maps. Short version: **NEO clones do not reach |v|>2** (support ends near |v|~2); the raw density out
there is **kNN bleed** from clones at lower |v|; the **support mask zeros MBA/TNO/Trojan** (support<1)
but **skips the smoothed NEO population**, so NEO wins the ratio by default → P_NEO=1.0 without any
genuine NEO clone support. All quantities derive from stored map arrays + Stage-1 parquets; no map
regeneration.

## Files

- **targets_28.parquet** — one row per (target object × pipeline; {n_rows} rows, {n_targets} objects).
  Columns: ObjID, pipeline; map_file, mag_bin; vlam, vbeta, total_speed (=hypot), magnitude;
  vdp_prod_pm2 (production ±2 score), vdp_pilot5_masked (±5, support_mask_min=1),
  vdp_pilot5_unmasked (±5, no mask); raw_density_{{MBA,NEO,TNO,Trojans}} at the query cell;
  support_{{pop}} (clone count at the cell); masked_density_{{pop}} (raw × support≥1, except NEO which
  is skipped); nearest_clone_dist_{{pop}} (from the map's nearest_dist array);
  clones_r{{05,10,15,25}}_{{pop}} (Σ support within 0.05/0.10/0.15/0.25 deg/day); kth_nearest_dist_{{pop}}
  (k=10, from cumulative support vs distance); neo_supported_vmax (max |v| with NEO support≥1);
  out_of_grid_production (|v|>2), genuine_support_at_query (NEO support≥1 at the cell).

- **clone_support_points.parquet** — for the 5 representative cases ({five}), the support-derived clone
  points: every grid cell with support≥1, per population: case, ObjID, population, clone_vlam,
  clone_vbeta, clone_weight_support (= support count in the cell), map_file, mag_bin. (Cell centres
  weighted by support; the maps were not built with --save-overlays so exact clone coords are not
  stored — the support grid is the complete record of where clones landed.)

- **arrays/diagnostic_map_arrays.npz** — for the 5 cases (keys `case{{i}}_{{ObjID}}__*`):
  vlam_grid, vbeta_grid; density_raw_{{pop}}, support_{{pop}}, masked_density_{{pop}}; masked_P_NEO,
  unmasked_P_NEO (full ±5, 1001²); prod_vlam_grid, prod_vbeta_grid, prod_P_NEO_pm2 (production ±2,
  401²); query (vlam,vbeta), query_cell_ixiy, meta [ObjID,pipeline,map,bin]. Sufficient to plot raw
  kNN density, clone support, the support-mask boundary, and masked vs unmasked probability without the
  original maps.

- **overlap_controls.parquet** — the |v|<1.9 overlap sample: ObjID, population, vlam, vbeta, prod_pm2,
  pilot5_masked, dP, map_file, mag_bin. Median |dP| = {ov_median:.4f}.

- **nonneo_controls.parquet** — the {n_nonneo} non-NEO contamination-control objects: ObjID,
  population, vlam, vbeta, prod_pm2, pilot5_masked, pilot5_unmasked, class_prod / class_pilot
  (threshold 0.5), map_file, mag_bin. False positives @0.5: prod {fp_prod} → pilot {fp_pilot}.

- **configuration_and_provenance.json** — git commit, generator command, grid limits/spacing, kNN
  params, support_mask_min, support-mask skip populations, smoothing, clone method/factors, probability
  formula, and sha256 of every source file used.

- **build_bleed_bundle.py** — the script that produced this bundle (copied in).

## Reproduces
28 targets recovered; overlap median |ΔP| ≈ 0.0003; non-NEO false positives 131 → 140.
"""

if __name__ == "__main__":
    main()
