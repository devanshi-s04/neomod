#!/usr/bin/env python3
"""STAGE B — explain every digest2=1 / VDP=0 zero for the night-61642 out-of-±2 NEOs.

For the 794 identity-matched night-61642 NEOs (benchmark pop==NEO ⋈ Sorcha), find every object whose
max(|vlam|,|vbeta|) exceeds the ACTUAL map boundary (±2.0 deg/day, verified in Stage A), and for
benchmark and Sorcha separately write the full diagnostic: map file, mag bin, vlam/vbeta, in/out-of-grid,
raw NEO density + total density + support + nearest NEO-clone distance at the nearest in-grid cell,
unmasked P_NEO, production-setting P_NEO, and an explicit zero reason.

Production settings (Stage A): mask_radius_deg_per_day=inf (nearest-dist mask OFF), support_mask_min=1,
maps prob_maps_grid_s3m_nbody, V-corrected mags on the Sorcha side.

Read-only. Output: outputs/mag245_nbody_per_object_digest2_audit/digest2_one_vdp_zero_summary.csv
"""
from __future__ import annotations
import os, sys, glob
from pathlib import Path
import numpy as np, pandas as pd

WORK = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(WORK / "neomod" / "src"))
MAPDIR = WORK / "prob_maps_grid_s3m_nbody"
S1 = WORK / "outputs" / "mag245_nbody_deterministic_rescore"
OUT = WORK / "outputs" / "mag245_nbody_per_object_digest2_audit"
GRID_LIM = 2.0          # verified actual boundary
NIGHT = 61642

import velocity_density_pipeline_gmm as vdp  # for the authoritative production/unmasked P_NEO


def load_matched():
    b = pd.read_parquet(S1 / "bench_deterministic.parquet", columns=[
        "s3m_objid", "population", "vlam", "vbeta", "prob_map_file", "mag_bin_label",
        "P_NEO_d2_det", "P_NEO_vdp", "mean_mag"])
    b["s3m_objid"] = b.s3m_objid.astype(str)
    b = b[b.population == "NEO"].drop_duplicates("s3m_objid").set_index("s3m_objid")
    s = pd.read_parquet(S1 / "sorcha_vcorr_deterministic.parquet", columns=[
        "ObjID", "night", "population", "vlam", "vbeta", "prob_map_file", "mag_bin_label_Vband",
        "P_NEO_d2_det", "P_NEO_vdp_Vband", "mean_mag_V"])
    s["ObjID"] = s.ObjID.astype(str)
    s = s[s.night == NIGHT].copy()
    m = s[s.ObjID.isin(b.index)].reset_index(drop=True)  # benchmark-NEO matched, 794
    return m, b


def edge_cell_lookup(z, pop, magbin, vlam, vbeta):
    """Nearest in-grid cell values (clamped) for density/support/nearest, + the raw ratio there.
    Returns (neo_dens, total_dens, support, nearest, unmasked_ratio_at_cell, ix, iy, in_grid)."""
    xg = z["x_grid"]; yg = z["y_grid"]
    nx, ny = len(xg), len(yg)
    ix = int(np.clip(np.round((vlam - xg[0]) / (xg[1] - xg[0])), 0, nx - 1))
    iy = int(np.clip(np.round((vbeta - yg[0]) / (yg[1] - yg[0])), 0, ny - 1))
    in_grid = (abs(vlam) <= GRID_LIM) and (abs(vbeta) <= GRID_LIM)
    pops = [p for p in z["population_names"]]
    dens = {p: z[f"density_raw__{p}__{magbin}"] for p in pops}
    total = sum(dens.values())
    neo = dens["NEO"][iy, ix]
    tot = total[iy, ix]
    sup = z[f"support_count__NEO__{magbin}"][iy, ix]
    near = z[f"nearest_dist__NEO__{magbin}"][iy, ix]
    ratio = (neo / tot) if tot > 0 else 0.0
    return float(neo), float(tot), float(sup), float(near), float(ratio), ix, iy, in_grid


def zero_reason(row, magbin_valid, map_present, in_grid, neo_dens, tot_dens, support, p_prod):
    if not map_present:
        return "missing_map"
    if not magbin_valid:
        return "magnitude_out_of_range"
    if not in_grid:
        return "out_of_grid"
    # in-grid but zero:
    if tot_dens <= 0:
        return "zero_total_density"
    if neo_dens <= 0:
        return "zero_NEO_density"
    if support < 1:
        return "support_mask"
    if p_prod == 0:
        return "other"
    return "nonzero"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    m, b = load_matched()
    bm = b.reindex(m.ObjID.to_numpy() if False else m.ObjID.astype(str).to_numpy()) \
        if False else b.reindex(m.ObjID.astype(str).to_numpy())
    # benchmark side arrays (aligned to m)
    m["b_vlam"] = bm.vlam.to_numpy(); m["b_vbeta"] = bm.vbeta.to_numpy()
    m["b_prob_map_file"] = bm.prob_map_file.to_numpy(); m["b_mag_bin"] = bm.mag_bin_label.to_numpy()
    m["b_d2"] = bm.P_NEO_d2_det.to_numpy(); m["b_vdp"] = bm.P_NEO_vdp.to_numpy()
    m["b_maxabs"] = np.maximum(np.abs(m.b_vlam), np.abs(m.b_vbeta))
    m["s_maxabs"] = np.maximum(np.abs(m.vlam), np.abs(m.vbeta))

    sel = (m.b_maxabs > GRID_LIM) | (m.s_maxabs > GRID_LIM)
    tgt = m[sel].copy()
    print(f"matched NEOs: {len(m)} | >±{GRID_LIM} in benchmark: {int((m.b_maxabs>GRID_LIM).sum())} "
          f"| in Sorcha: {int((m.s_maxabs>GRID_LIM).sum())} | both: {int(((m.b_maxabs>GRID_LIM)&(m.s_maxabs>GRID_LIM)).sum())} "
          f"| union (selected): {len(tgt)}", flush=True)

    # unique maps needed
    maps_needed = set(tgt.b_prob_map_file) | set(tgt.prob_map_file)
    print(f"unique map centers referenced: {len(maps_needed)}", flush=True)
    cache = {}
    def getz(name):
        if name not in cache:
            p = MAPDIR / name
            cache[name] = (np.load(p, allow_pickle=True) if p.exists() else None)
        return cache[name]

    rows = []
    for _, r in tgt.iterrows():
        for pipe, vlam, vbeta, pmf, magbin, d2, vdps, maxabs in [
            ("benchmark", r.b_vlam, r.b_vbeta, r.b_prob_map_file, r.b_mag_bin, r.b_d2, r.b_vdp, r.b_maxabs),
            ("sorcha_vcorr", r.vlam, r.vbeta, r.prob_map_file, r.mag_bin_label_Vband, r.P_NEO_d2_det, r.P_NEO_vdp_Vband, r.s_maxabs)]:
            z = getz(str(pmf))
            map_present = z is not None
            magbin_valid = isinstance(magbin, str) and magbin in (list(z["mag_bin_labels"]) if z is not None else [])
            if map_present and magbin_valid:
                neo, tot, sup, near, ratio, ix, iy, in_grid = edge_cell_lookup(z, "NEO", magbin, vlam, vbeta)
            else:
                neo = tot = sup = near = ratio = np.nan; ix = iy = -1
                in_grid = (abs(vlam) <= GRID_LIM) and (abs(vbeta) <= GRID_LIM)
            # authoritative production + unmasked P_NEO via ProbMapSet (out-of-grid -> 0 by design)
            p_prod = p_unmasked = np.nan
            if map_present:
                pm_prod = vdp.ProbMapSet.from_npz(str(MAPDIR / pmf), support_mask_min=1,
                                                  mask_radius_deg_per_day=np.inf)
                p_prod = float(pm_prod.score_visible(np.array([vlam]), np.array([vbeta]),
                                                     np.array([r.mean_mag if pipe=="benchmark" else r.mean_mag_V]))["NEO"][0])
                pm_raw = vdp.ProbMapSet.from_npz(str(MAPDIR / pmf), support_mask_min=None,
                                                 mask_radius_deg_per_day=np.inf)
                p_unmasked = float(pm_raw.score_visible(np.array([vlam]), np.array([vbeta]),
                                                        np.array([r.mean_mag if pipe=="benchmark" else r.mean_mag_V]))["NEO"][0])
            reason = zero_reason(r, magbin_valid, map_present, in_grid, neo, tot, sup, p_prod)
            rows.append(dict(
                ObjID=r.ObjID, pipeline=pipe, map_file=str(pmf), mag_bin=magbin,
                vlam=round(float(vlam),4), vbeta=round(float(vbeta),4), max_abs_component=round(float(maxabs),4),
                in_grid=bool(in_grid), edge_ix=ix, edge_iy=iy,
                raw_NEO_density=neo, total_density=tot, support_count=sup, nearest_NEO_clone_dist=near,
                P_NEO_unmasked=round(p_unmasked,6) if np.isfinite(p_unmasked) else np.nan,
                P_NEO_production=round(p_prod,6) if np.isfinite(p_prod) else np.nan,
                digest2_det=round(float(d2),4), vdp_stored=round(float(vdps),6),
                zero_reason=reason))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "digest2_one_vdp_zero_summary.csv", index=False)
    pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
    print(f"\n=== zero-reason breakdown (per pipeline row) ===", flush=True)
    print(df.groupby(["pipeline", "zero_reason"]).size().to_string(), flush=True)
    print(f"\n  objects geometrically OUT of ±{GRID_LIM}: "
          f"benchmark {int((df[df.pipeline=='benchmark'].in_grid==False).sum())}, "
          f"sorcha {int((df[df.pipeline=='sorcha_vcorr'].in_grid==False).sum())}", flush=True)
    print(f"  max |component| among selected: benchmark {df[df.pipeline=='benchmark'].max_abs_component.max():.4f}, "
          f"sorcha {df[df.pipeline=='sorcha_vcorr'].max_abs_component.max():.4f}", flush=True)
    print(f"  mag_bin valid (all in 14-25 range)? {df.mag_bin.notna().all()}  | all maps present? {df.map_file.map(lambda f:(MAPDIR/f).exists()).all()}", flush=True)
    print(f"  production P_NEO all zero for these? {np.allclose(df.P_NEO_production.fillna(-1), 0)}", flush=True)
    print(f"\nwrote {OUT/'digest2_one_vdp_zero_summary.csv'}  ({len(df)} rows, {df.ObjID.nunique()} objects, "
          f"{df.map_file.nunique()} unique maps)", flush=True)


if __name__ == "__main__":
    main()
