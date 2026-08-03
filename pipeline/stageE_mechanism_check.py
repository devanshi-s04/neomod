#!/usr/bin/env python3
"""Definitive check: WHY do the recovered targets score P_NEO=1.0 when NEO clone support≈0 there?
Compare raw vs support-masked density at the query cell, the total (all-pop) density, and how far
NEO clone support actually extends in |v|."""
import numpy as np, pandas as pd, sys
from pathlib import Path
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
OUT = W / "outputs/mag245_nbody_vlim5_pilot_validation"
PILOT = W / "prob_maps_grid_s3m_nbody_vlim5_pilot"
import velocity_density_pipeline_gmm as vdp

rec = pd.read_csv(OUT / "vlim5_recovery_targets.csv")
sample = rec[rec.vdp_pilot5 > 0.5].drop_duplicates("ObjID").head(5)
for _, r in sample.iterrows():
    cen = r.center; mb = r.mag_bin
    z = np.load(PILOT / cen, allow_pickle=True)
    xg = z["x_grid"]; yg = z["y_grid"]
    ix = int(np.clip(round((r.vlam - xg[0]) / (xg[1] - xg[0])), 0, len(xg) - 1))
    iy = int(np.clip(round((r.vbeta - yg[0]) / (yg[1] - yg[0])), 0, len(yg) - 1))
    pops = list(z["population_names"])
    neo_raw = float(z[f"density_raw__NEO__{mb}"][iy, ix])
    tot_raw = float(sum(z[f"density_raw__{p}__{mb}"][iy, ix] for p in pops))
    sup = float(z[f"support_count__NEO__{mb}"][iy, ix])
    # per-pop density at the cell
    perpop = {p: float(z[f"density_raw__{p}__{mb}"][iy, ix]) for p in pops}
    # masked vs unmasked prob at the cell
    pmM = vdp.ProbMapSet.from_npz(str(PILOT / cen), support_mask_min=1, mask_radius_deg_per_day=np.inf)
    pmR = vdp.ProbMapSet.from_npz(str(PILOT / cen), support_mask_min=None, mask_radius_deg_per_day=np.inf)
    pM = float(pmM.get_probability_map(mb, "NEO")[iy, ix])
    pR = float(pmR.get_probability_map(mb, "NEO")[iy, ix])
    # how far does NEO clone support reach in |v|?
    sc = z[f"support_count__NEO__{mb}"]
    xx, yy = np.meshgrid(xg, yg)
    vmag = np.maximum(np.abs(xx), np.abs(yy))
    has = sc >= 1
    max_v_support = float(vmag[has].max()) if has.any() else 0.0
    n_sup_gt2 = int((has & (vmag > 2.0)).sum())
    print(f"\n{r.ObjID} [{r.pipeline}] center={cen} mag_bin={mb}  |v|={r.maxabs:.2f}")
    print(f"  query cell: NEO_support={sup:.0f}  NEO_dens_raw={neo_raw:.4f}  total_dens_raw={tot_raw:.4f}")
    print(f"  per-pop density at cell: " + "  ".join(f"{p}={perpop[p]:.4f}" for p in pops))
    print(f"  P_NEO at cell: masked(support_min=1)={pM:.4f}   unmasked={pR:.4f}   scored={r.vdp_pilot5:.3f}")
    print(f"  NEO clone support in this map reaches |v|max={max_v_support:.2f} deg/day; cells with support≥1 & |v|>2: {n_sup_gt2}")
    del pmM, pmR
