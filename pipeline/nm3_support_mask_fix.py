#!/usr/bin/env python3
"""Test candidate fixes for the |v| 0.25-0.5 false positives (docs §9.4).

The support mask is applied at SCORE time (ProbMapSet.__init__), not baked into the npz -- so every
variant here is a pure re-score of the SAME 667 maps. No rebuild.

Variants:
  baseline    support_mask_min=1, NEO EXEMPT (_support_mask_skip={'NEO'})   <- current production
  symmetric   support_mask_min=1 applied to ALL populations incl. NEO       <- "stop exempting NEO"
  nomask      support_mask_min=None                                          <- control
  mba_relax   mask everyone except NEO, but MBA needs only support>0        <- "restore MBA in tail"

Reported per variant: pooled AUC/F1/completeness/contamination, the |v| 0.25-0.5 FP count (77% of
the problem), and what happens to the 7,727 |v|>2 NEOs that motivated the +-5 grid -- a fix that
recovers contamination by destroying the fast-NEO win is NOT a fix.
"""
from __future__ import annotations
import os, sys, glob
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score, precision_recall_curve
W = "/mmfs1/gscratch/dirac/ds2004/sorcha"
sys.path.insert(0, f"{W}/neomod/src"); os.chdir(f"{W}/neomod")
import velocity_density_pipeline_neomod_clone_only as vnm
MAPS = f"{W}/prob_maps_grid_neomod3_full"
OUT = f"{W}/outputs/neomod3_fullgrid"
VARIANTS = ["baseline", "symmetric", "nomask", "mba_relax"]


def rebuild_prob_maps(pm, variant):
    """Recompute pm._prob_maps under a different masking rule (mirrors __init__ lines 2544-2569)."""
    if variant == "baseline":
        return pm
    for mb in pm.mag_bins:
        label = mb["label"]
        density = pm.results[label]["density_maps_downweighted_raw"]
        nearest = pm.results[label]["nearest_dist_maps"]
        support = pm.results[label].get("support_count_maps", {})
        density = {p: np.asarray(a).copy() for p, a in density.items()}
        if variant != "nomask":
            for pop in pm.population_names:
                sc = support.get(pop)
                if sc is None:
                    continue
                if variant == "symmetric":
                    thr = 1.0                                   # NEO no longer exempt
                elif variant == "mba_relax":
                    if pop == "NEO":
                        continue                                # NEO still exempt
                    thr = 1e-12 if pop == "MBA" else 1.0        # MBA: any support at all counts
                density[pop] = np.where(np.asarray(sc) >= thr, density[pop], 0.0)
        total = sum(density.values())
        pf = {}
        for pop in pm.population_names:
            frac = np.zeros_like(density[pop], dtype=np.float64)
            good = total > 0
            frac[good] = density[pop][good] / total[good]
            frac[np.asarray(nearest[pop]) > pm.mask_radius_deg_per_day] = 0.0
            pf[pop] = frac
        pm._prob_maps[label] = pf
    return pm


def score_center(cen, sub):
    out = {}
    for var in VARIANTS:
        pm = vnm.ProbMapSet.from_npz(f"{MAPS}/{cen}",
                                     support_mask_min=(None if var == "nomask" else 1),
                                     mask_radius_deg_per_day=np.inf)
        pm = rebuild_prob_maps(pm, var)
        out[var] = np.asarray(pm.score_visible(sub.vlam.to_numpy(float), sub.vbeta.to_numpy(float),
                                               sub.mean_mag_V.to_numpy(float))["NEO"], float)
        del pm
    return sub.index.to_numpy(), out


d = pd.read_parquet(f"{OUT}/fullgrid_scored_tracklets.parquet")
have = {os.path.basename(f) for f in glob.glob(f"{MAPS}/*.npz")}
d = d[d.prob_map_file.isin(have)].reset_index(drop=True)
print(f"re-scoring {len(d):,} tracklets x {len(VARIANTS)} variants over {d.prob_map_file.nunique()} centers")
res = Parallel(n_jobs=20, verbose=5)(delayed(score_center)(c, g) for c, g in d.groupby("prob_map_file"))
for var in VARIANTS:
    d[var] = np.nan
for idx, out in res:
    for var in VARIANTS:
        d.loc[idx, var] = out[var]
d["d2"] = d.P_NEO_d2.fillna(0.0)
y = d.is_neo.to_numpy()

def report(name, s):
    s = np.nan_to_num(np.asarray(s, float))
    p, r, t = precision_recall_curve(y, s)
    f1 = np.divide(2*p*r, p+r, out=np.zeros_like(p), where=(p+r) > 0)
    i = int(np.argmax(f1[:-1])); T = float(t[i]) if i < len(t) else 1.0
    band = (d.maxabs_v > 0.25) & (d.maxabs_v <= 0.5)
    fp_band = int(((s >= T) & (y == 0) & band).sum())
    fast = d.maxabs_v > 2.0
    return dict(variant=name, AUC=roc_auc_score(y, s), bestF1=f1[i], completeness=r[i]*100,
                contamination=(1-p[i])*100, FP_total=int(((s >= T) & (y == 0)).sum()),
                FP_v025_05=fp_band,
                fastNEO_zero_pct=100*float((s[fast.to_numpy()] == 0).mean()),
                fastNEO_median=float(np.median(s[fast.to_numpy()])))

rows = [report(v, d[v].to_numpy()) for v in VARIANTS] + [report("digest2", d.d2.to_numpy())]
summ = pd.DataFrame(rows)
summ.to_csv(f"{OUT}/support_mask_variants.csv", index=False)
pd.set_option("display.width", 250)
print("\n" + "="*100)
print(summ.to_string(index=False, float_format=lambda v: f"{v:,.4g}"))
print("\nfastNEO_* = the 7,727 |v|>2 NEOs the ±5 grid exists to recover "
      "(zero_pct should stay LOW, median should stay HIGH)")
