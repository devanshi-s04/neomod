#!/usr/bin/env python3
"""STAGE E — validate the ±5 velocity-domain pilot (29 centers) vs the ±2 production maps.

Parallel over the 29 centers (joblib, one worker per center). For every benchmark and Sorcha object
assigned to a pilot center, re-score VDP against the pilot ±5 map (production settings:
support_mask_min=1, nearest-dist mask OFF), compare to the stored production (±2) score, and check:

  1. all 28 out-of-±2 target NEOs are now geometrically in-grid (|v|<5) and score nonzero;
  2. recovered P_NEO comes from GENUINE populated NEO density (support>0), not extrapolation;
  3. existing in-domain (|v|<2) scores stay stable (overlap continuity: pilot±5 vs prod±2);
  4. no meaningful contamination increase in the pilot-center population.

Read-only. Outputs -> outputs/mag245_nbody_vlim5_pilot_validation/
"""
from __future__ import annotations
import os, sys, glob
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed

WORK = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(WORK / "neomod" / "src"))
PILOT = WORK / "prob_maps_grid_s3m_nbody_vlim5_pilot"
PROD = WORK / "prob_maps_grid_s3m_nbody"
S1 = WORK / "outputs" / "mag245_nbody_deterministic_rescore"
OUT = WORK / "outputs" / "mag245_nbody_vlim5_pilot_validation"
NIGHT = 61642
NJOBS = int(os.environ.get("STAGE_E_JOBS", "8"))


def load():
    b = pd.read_parquet(S1 / "bench_deterministic.parquet", columns=[
        "s3m_objid", "population", "vlam", "vbeta", "prob_map_file", "mean_mag", "P_NEO_vdp", "P_NEO_d2_det"])
    b["s3m_objid"] = b.s3m_objid.astype(str)
    b = b.rename(columns={"s3m_objid": "ObjID", "mean_mag": "mag", "P_NEO_vdp": "vdp_prod"}); b["pipeline"] = "benchmark"
    s = pd.read_parquet(S1 / "sorcha_vcorr_deterministic.parquet", columns=[
        "ObjID", "night", "population", "vlam", "vbeta", "prob_map_file", "mean_mag_V", "P_NEO_vdp_Vband", "P_NEO_d2_det"])
    s["ObjID"] = s.ObjID.astype(str); s = s[s.night == NIGHT].copy()
    s = s.rename(columns={"mean_mag_V": "mag", "P_NEO_vdp_Vband": "vdp_prod"}); s["pipeline"] = "sorcha_vcorr"
    cols = ["ObjID", "pipeline", "population", "vlam", "vbeta", "prob_map_file", "mag", "vdp_prod", "P_NEO_d2_det"]
    return b[cols], s[cols]


def process_center(cen, sub, tgt_ids):
    """Score one center's objects against pilot±5 (and prod±2 in the overlap). Returns (rows, overlap)."""
    import velocity_density_pipeline_gmm as vdp
    pm5 = vdp.ProbMapSet.from_npz(str(PILOT / cen), support_mask_min=1, mask_radius_deg_per_day=np.inf)
    vl = sub.vlam.to_numpy(float); vb = sub.vbeta.to_numpy(float); mg = sub.mag.to_numpy(float)
    new5 = pm5.score_visible(vl, vb, mg)["NEO"]
    ov = np.maximum(np.abs(vl), np.abs(vb)) < 1.9
    overlap = None
    if ov.sum() >= 5:
        pm2 = vdp.ProbMapSet.from_npz(str(PROD / cen), support_mask_min=1, mask_radius_deg_per_day=np.inf)
        p2 = pm2.score_visible(vl[ov], vb[ov], mg[ov])["NEO"]
        d = np.abs(new5[ov] - p2)
        overlap = dict(center=cen, n_overlap=int(ov.sum()), median_abs_dP=float(np.median(d)),
                       p95_abs_dP=float(np.percentile(d, 95)), max_abs_dP=float(d.max()))
    # detailed density lookup ONLY for the target objects at this center
    is_t = sub.ObjID.isin(tgt_ids).to_numpy()
    zz = np.load(str(PILOT / cen), allow_pickle=True) if is_t.any() else None
    rows = []
    for j in range(len(sub)):
        r = sub.iloc[j]
        neo_d = tot_d = supp = np.nan; mb = None
        if is_t[j] and zz is not None:
            for k in range(len(zz["mag_bin_labels"])):
                if zz["mag_bin_mins"][k] <= r.mag < zz["mag_bin_maxs"][k]:
                    mb = str(zz["mag_bin_labels"][k]); break
            if mb is not None:
                xg = zz["x_grid"]; yg = zz["y_grid"]
                ix = int(np.clip(round((r.vlam - xg[0]) / (xg[1] - xg[0])), 0, len(xg) - 1))
                iy = int(np.clip(round((r.vbeta - yg[0]) / (yg[1] - yg[0])), 0, len(yg) - 1))
                neo_d = float(zz[f"density_raw__NEO__{mb}"][iy, ix])
                tot_d = float(sum(zz[f"density_raw__{p}__{mb}"][iy, ix] for p in zz["population_names"]))
                supp = float(zz[f"support_count__NEO__{mb}"][iy, ix])
        rows.append(dict(ObjID=r.ObjID, pipeline=r.pipeline, population=r.population, center=cen,
                         vlam=round(r.vlam, 4), vbeta=round(r.vbeta, 4),
                         maxabs=round(max(abs(r.vlam), abs(r.vbeta)), 4), in_grid5=bool(max(abs(r.vlam), abs(r.vbeta)) < 5.0),
                         vdp_prod=round(float(r.vdp_prod), 6), vdp_pilot5=round(float(new5[j]), 6),
                         digest2_det=round(float(r.P_NEO_d2_det), 4), is_target=bool(is_t[j]),
                         mag_bin=mb, NEO_density=neo_d, total_density=tot_d, support=supp))
    return rows, overlap


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    centers = [os.path.basename(f) for f in sorted(glob.glob(str(PILOT / "*.npz")))]
    b, s = load()
    b["maxabs"] = np.maximum(b.vlam.abs(), b.vbeta.abs()); s["maxabs"] = np.maximum(s.vlam.abs(), s.vbeta.abs())
    both = pd.concat([b, s], ignore_index=True)
    inp = both[both.prob_map_file.isin(centers)].copy()
    bneo = set(b[b.population == "NEO"].ObjID); matched = s[s.ObjID.isin(bneo)]
    tgt_ids = set(matched[matched.maxabs > 2.0].ObjID) | set(b[b.ObjID.isin(set(matched.ObjID)) & (b.maxabs > 2.0)].ObjID)
    print(f"objects at 29 centers: {len(inp):,} ({inp.pipeline.value_counts().to_dict()}) | targets: {len(tgt_ids)}", flush=True)

    groups = [(cen, inp[inp.prob_map_file == cen].reset_index(drop=True)) for cen in centers]
    out = Parallel(n_jobs=NJOBS, verbose=5)(delayed(process_center)(cen, sub, tgt_ids) for cen, sub in groups)
    rows = [r for rr, _ in out for r in rr]
    overlap = [o for _, o in out if o is not None]
    res = pd.DataFrame(rows)
    res.to_parquet(OUT / "vlim5_pilot_rescore.parquet", index=False)
    ov = pd.DataFrame(overlap); ov.to_csv(OUT / "vlim5_overlap_continuity.csv", index=False)
    res[res.is_target].to_csv(OUT / "vlim5_recovery_targets.csv", index=False)

    print("\n============ ACCEPTANCE GATE ============", flush=True)
    tt = res[res.is_target]
    print(f"[1] targets: {tt.ObjID.nunique()} unique | all in_grid(<5): {bool(tt.in_grid5.all())} | max maxabs: {tt.maxabs.max():.3f}", flush=True)
    rec_ok = (tt.vdp_prod < 1e-6) & (tt.vdp_pilot5 > 1e-6)
    print(f"    recovered (prod≈0 -> pilot>0): {int(rec_ok.sum())}/{len(tt)} rows | pilot P_NEO median {tt.vdp_pilot5.median():.3f}, >0.5: {int((tt.vdp_pilot5>0.5).sum())}", flush=True)
    genuine = (tt.support >= 1) & (tt.NEO_density > 0)
    print(f"[2] recovered from GENUINE NEO density (support≥1 & NEO_dens>0): {int(genuine.sum())}/{len(tt)}", flush=True)
    print(f"[3] overlap continuity (|v|<1.9, pilot±5 vs prod±2): median|dP| {ov.median_abs_dP.median():.4f}, p95 {ov.p95_abs_dP.median():.4f}, worst-center max {ov.max_abs_dP.max():.4f}", flush=True)
    nn = res[res.population != "NEO"]
    fp_p = int((nn.vdp_prod >= 0.5).sum()); fp_n = int((nn.vdp_pilot5 >= 0.5).sum())
    print(f"[4] non-NEO false-positives @P≥0.5: prod {fp_p} -> pilot {fp_n} (Δ {fp_n-fp_p}) of {len(nn)} non-NEO", flush=True)
    g = res.groupby("population").apply(lambda x: pd.Series({
        "n": len(x), "median_prod": x.vdp_prod.median(), "median_pilot": x.vdp_pilot5.median(),
        "mean_dP": (x.vdp_pilot5 - x.vdp_prod).mean(), "n_dP>0.1": int(((x.vdp_pilot5 - x.vdp_prod).abs() > 0.1).sum())}))
    print("\nper-population P_NEO change (pilot − prod) at the 29 centers:\n" + g.to_string(), flush=True)
    print(f"\nwrote -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
