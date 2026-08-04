#!/usr/bin/env python3
"""Direct apples-to-apples map-resolution experiment: step 0.01 vs 0.005.

Both map sets are scored EXCLUSIVELY through the unchanged production path
(ProbMapSet.score_visible). No scorer modification, no density substitution. The classifier's
frozen semantics -- pixelwise probability normalisation then bilinear interpolation -- are
identical for both; only the grid step differs.

Gate: rescoring the 0.01 maps must reproduce the archived E0 metrics and coverage exactly.
"""
import json, os, sys, time, resource
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, precision_recall_curve
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
import velocity_density_pipeline_neomod_clone_only as v
OUT = W/"outputs/e0_results"
A_DIR, B_DIR = W/"prob_maps_e0_thr10", W/"prob_maps_e0_thr10_step005"
CENTERS = [(0,0),(20,-12),(10,2),(-10,-2),(0,8),(50,-1),(50,1),(-50,1),
           (90,0),(140,0),(-140,0),(0,25),(0,-25),(0,50),(0,-50),(-50,-50)]
POPS = ["NEO","MBA","TNO","Trojans"]; THRS=[2,3,5,10]
t_start = time.time()

cal = pd.read_parquet(W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet")
names = {f"prob_maps_grid_dlon{a:+04d}_lat{b:+03d}.npz" for a,b in CENTERS}
cal = cal[cal.prob_map_file.isin(names)].reset_index(drop=True)

def score_dir(dirp):
    s = np.full(len(cal), np.nan); c = np.zeros(len(cal), bool)
    for cen, g in cal.groupby("prob_map_file"):
        pm = v.ProbMapSet.from_npz(str(dirp/cen), support_mask_min=1,
                                   mask_radius_deg_per_day=np.inf)
        r = pm.score_visible(g.vlam.to_numpy(float), g.vbeta.to_numpy(float),
                             g.mean_mag.to_numpy(float))
        idx = g.index.to_numpy()
        s[idx] = np.asarray(r["NEO"], float)
        c[idx] = np.sum([np.nan_to_num(np.asarray(r[p], float)) for p in POPS], axis=0) > 0
        del pm
    return s, c

# frozen common set (all four thresholds at 0.01) -- unchanged definition
cov = {}
for t in THRS:
    _, cov[t] = score_dir(W/f"prob_maps_e0_thr{t}")
common = np.logical_and.reduce([cov[t] for t in THRS])
y = (cal.population == "NEO").astype(int).to_numpy()
print(f"frozen common rows {int(common.sum()):,}  NEO {int(y[common].sum())}  "
      f"abstentions {int((~common).sum())} (NEO {int(y[~common].sum())})")

sA, cA = score_dir(A_DIR)
def mets(s, m):
    yy, ss = y[m], np.nan_to_num(s[m])
    p_, r_, t_ = precision_recall_curve(yy, ss)
    f1 = np.divide(2*p_*r_, p_+r_, out=np.zeros_like(p_), where=(p_+r_)>0); i = int(np.argmax(f1[:-1]))
    return dict(ROC_AUC=roc_auc_score(yy,ss), pAUC=roc_auc_score(yy,ss,max_fpr=0.01),
                F1=f1[i], thr=float(t_[i]) if i < len(t_) else 1.0)
MA = mets(sA, common)
arch = pd.read_csv(OUT/"E0_SIMPLE_RESULTS.csv"); a10 = arch[arch.threshold == 10].iloc[0]
print("\n=== GATE: rescoring 0.01 must reproduce archived E0 ===")
for k, got, want in [("ROC_AUC",MA["ROC_AUC"],a10.ROC_AUC),("pAUC",MA["pAUC"],a10.pAUC_fpr01_norm),
                     ("best_F1",MA["F1"],a10.best_F1),
                     ("coverage_pct",100*cA.mean(),a10.coverage_pct)]:
    print(f"  {k:12s} {got:.9f} vs archived {want:.9f}   d={got-want:+.3e}")
gate = (abs(MA["ROC_AUC"]-a10.ROC_AUC)<1e-9 and abs(MA["pAUC"]-a10.pAUC_fpr01_norm)<1e-9
        and abs(MA["F1"]-a10.best_F1)<1e-9 and abs(100*cA.mean()-a10.coverage_pct)<1e-6)
print(f"  GATE: {'PASS' if gate else 'FAIL'}")
if not gate:
    sys.exit("baseline reproduction failed -- refusing to interpret the comparison")

sB, cB = score_dir(B_DIR)
MB = mets(sB, common)
m = common
a, b = np.nan_to_num(sA[m]), np.nan_to_num(sB[m])
d = np.abs(b-a); thrA = MA["thr"]
flips = (a >= thrA) != (b >= thrA)
print(f"\n=== 0.01 vs 0.005, frozen {int(m.sum()):,} rows / {int(y[m].sum())} NEOs ===")
print(f"  changed {int((b!=a).sum()):,} ({100*(b!=a).mean():.2f}%)")
print(f"  |dP| median {np.median(d):.3e}  p95 {np.percentile(d,95):.3e}  "
      f"p99 {np.percentile(d,99):.3e}  max {d.max():.3e}")
print(f"  coverage {100*cA.mean():.4f}% -> {100*cB.mean():.4f}%")
print(f"  ROC AUC {MA['ROC_AUC']:.6f} -> {MB['ROC_AUC']:.6f}  (d={MB['ROC_AUC']-MA['ROC_AUC']:+.2e})")
print(f"  pAUC    {MA['pAUC']:.6f} -> {MB['pAUC']:.6f}  (d={MB['pAUC']-MA['pAUC']:+.2e})")
print(f"  best F1 {MA['F1']:.6f} -> {MB['F1']:.6f}  (d={MB['F1']-MA['F1']:+.2e})")
print(f"  flips at frozen 0.01 best-F1 threshold ({thrA:.6g}): {int(flips.sum())}")
S = pd.DataFrame({"truth":cal.population.to_numpy()[m],"mag":cal.mean_mag.to_numpy()[m],
                  "vmax":np.maximum(cal.vlam.abs(),cal.vbeta.abs()).to_numpy()[m],
                  "absd":d,"flip":flips})
for lab, key in [("truth population","truth"),
                 ("magnitude bin", pd.cut(S.mag,[14,20,22,23,24,25])),
                 ("velocity band", pd.cut(S.vmax,[0,0.25,0.5,1.0,2.0,5.1]))]:
    print(f"\n  by {lab}:")
    print(S.groupby(key, observed=True).agg(n=("absd","size"),median=("absd","median"),
          p99=("absd",lambda z:np.percentile(z,99)),max=("absd","max"),
          flips=("flip","sum")).to_string(float_format=lambda z:f"{z:,.3e}"))
top = S.assign(P_001=a, P_0005=b, mapf=cal.prob_map_file.to_numpy()[m],
               vlam=cal.vlam.to_numpy()[m], vbeta=cal.vbeta.to_numpy()[m]
               ).nlargest(10,"absd")
print("\n  ten largest score changes:")
print(top[["mapf","mag","vlam","vbeta","truth","P_001","P_0005","absd"]].to_string(
    index=False, float_format=lambda z:f"{z:,.5g}"))

szA = sum(f.stat().st_size for f in A_DIR.glob("*.npz"))
szB = sum(f.stat().st_size for f in B_DIR.glob("*.npz"))
peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20
print(f"\n=== cost ===")
print(f"  archive 16 maps: 0.01 = {szA/2**30:.2f} GiB   0.005 = {szB/2**30:.2f} GiB  "
      f"({szB/szA:.2f}x)")
print(f"  projected 667 maps: 0.01 = {szA/16*667/2**30:.0f} GiB   0.005 = {szB/16*667/2**30:.0f} GiB")
print(f"  comparison runtime {time.time()-t_start:.0f}s   peak RSS {peak:.1f} GiB")
json.dump(dict(gate=bool(gate), rows=int(m.sum()), neos=int(y[m].sum()),
   changed=int((b!=a).sum()), median=float(np.median(d)), p95=float(np.percentile(d,95)),
   p99=float(np.percentile(d,99)), max=float(d.max()), flips=int(flips.sum()),
   coverage=[100*float(cA.mean()),100*float(cB.mean())],
   roc=[MA["ROC_AUC"],MB["ROC_AUC"]], pauc=[MA["pAUC"],MB["pAUC"]], f1=[MA["F1"],MB["F1"]],
   archive_gib=[szA/2**30,szB/2**30], projected_667_gib=[szA/16*667/2**30,szB/16*667/2**30]),
   open(OUT/"MAP_RESOLUTION_COMPARISON.json","w"), indent=2)
top.to_csv(OUT/"MAP_RESOLUTION_TOP10.csv", index=False)
print(f"\nwrote MAP_RESOLUTION_COMPARISON.json and MAP_RESOLUTION_TOP10.csv")
