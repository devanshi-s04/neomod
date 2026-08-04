#!/usr/bin/env python3
"""Mechanism confirmation: is the 0.01-vs-0.005 disagreement caused by per-pixel occupancy masking?

DIAGNOSTIC ONLY -- no resolution is selected here and the production default is not changed.
Hypothesis: support_mask_min=1 is a PER-PIXEL occupancy test, so halving the pixel size creates MBA
support holes; the MBA denominator term is zeroed there and true MBAs are promoted to P(NEO)=1.
Test: rescore both map sets with support_mask_min=0 (mask effectively disabled) and see whether the
disagreement collapses.
"""
import json, os, sys
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
POPS = ["NEO","MBA","TNO","Trojans"]; THRS = [2,3,5,10]
BINS = [("14_16",14.,16.),("16_18",16.,18.),("18_20",18.,20.),("mag20",20.,21.),
        ("mag21",21.,22.),("mag22",22.,23.),("mag23",23.,24.),("mag24+",24.,25.)]

cal = pd.read_parquet(W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet")
names = {f"prob_maps_grid_dlon{a:+04d}_lat{b:+03d}.npz" for a,b in CENTERS}
cal = cal[cal.prob_map_file.isin(names)].reset_index(drop=True)
y = (cal.population == "NEO").astype(int).to_numpy()

def score_dir(dirp, smin):
    s = np.full(len(cal), np.nan); c = np.zeros(len(cal), bool)
    for cen, g in cal.groupby("prob_map_file"):
        pm = v.ProbMapSet.from_npz(str(dirp/cen), support_mask_min=smin,
                                   mask_radius_deg_per_day=np.inf)
        r = pm.score_visible(g.vlam.to_numpy(float), g.vbeta.to_numpy(float),
                             g.mean_mag.to_numpy(float))
        i = g.index.to_numpy(); s[i] = np.asarray(r["NEO"], float)
        c[i] = np.sum([np.nan_to_num(np.asarray(r[p], float)) for p in POPS], axis=0) > 0
        del pm
    return s, c

# frozen common set: unchanged definition (coverage across the four 0.01 thresholds, mask=1)
CACHE = OUT/"_maskmech_cache.npz"
cov = {}
for t in THRS:
    _, cov[t] = score_dir(W/f"prob_maps_e0_thr{t}", 1)
common = np.logical_and.reduce([cov[t] for t in THRS])
print(f"frozen rows {int(common.sum()):,}  NEO {int(y[common].sum())}  "
      f"abstentions {int((~common).sum())} (all NEO: {int(y[~common].sum())==int((~common).sum())})")

def mets(s, m):
    yy, ss = y[m], np.nan_to_num(s[m])
    p_, r_, t_ = precision_recall_curve(yy, ss)
    f1 = np.divide(2*p_*r_, p_+r_, out=np.zeros_like(p_), where=(p_+r_)>0); i = int(np.argmax(f1[:-1]))
    return dict(ROC=roc_auc_score(yy,ss), pAUC=roc_auc_score(yy,ss,max_fpr=0.01), F1=f1[i],
                thr=float(t_[i]) if i < len(t_) else 1.0)

out = {}
for smin, tag in [(1, "mask_on"), (0, "mask_off")]:
    sA, cA = score_dir(A_DIR, smin); sB, cB = score_dir(B_DIR, smin)
    m = common
    a, b = np.nan_to_num(sA[m]), np.nan_to_num(sB[m])
    d = np.abs(b-a); MA, MB = mets(sA, m), mets(sB, m)
    fl = (a >= MA["thr"]) != (b >= MA["thr"])
    out[tag] = dict(median=float(np.median(d)), p95=float(np.percentile(d,95)),
                    p99=float(np.percentile(d,99)), max=float(d.max()),
                    changed=int((b!=a).sum()), flips=int(fl.sum()),
                    covA=100*float(cA.mean()), covB=100*float(cB.mean()),
                    ROC=[MA["ROC"],MB["ROC"]], pAUC=[MA["pAUC"],MB["pAUC"]], F1=[MA["F1"],MB["F1"]])
    if smin == 1:
        flips_idx = np.where(m)[0][fl]; sA1, sB1 = sA, sB
    else:
        sA0, sB0 = sA, sB
    print(f"\n=== support_mask_min={smin} ({tag}) : 0.01 vs 0.005 ===")
    print(f"  changed {out[tag]['changed']:,}   flips {out[tag]['flips']}")
    print(f"  |dP| median {out[tag]['median']:.3e}  p95 {out[tag]['p95']:.3e}  "
          f"p99 {out[tag]['p99']:.3e}  max {out[tag]['max']:.3e}")
    print(f"  coverage {out[tag]['covA']:.4f}% -> {out[tag]['covB']:.4f}%")
    print(f"  ROC  {MA['ROC']:.6f} -> {MB['ROC']:.6f}  (d={MB['ROC']-MA['ROC']:+.2e})")
    print(f"  pAUC {MA['pAUC']:.6f} -> {MB['pAUC']:.6f}  (d={MB['pAUC']-MA['pAUC']:+.2e})")
    print(f"  F1   {MA['F1']:.6f} -> {MB['F1']:.6f}  (d={MB['F1']-MA['F1']:+.2e})")

print(f"\n>>> flip count {out['mask_on']['flips']} (mask on) vs {out['mask_off']['flips']} (mask off)")
print(f">>> pAUC disagreement {out['mask_on']['pAUC'][1]-out['mask_on']['pAUC'][0]:+.3e} (on) "
      f"vs {out['mask_off']['pAUC'][1]-out['mask_off']['pAUC'][0]:+.3e} (off)")

# ---- corner-level export for the original flips ----
def corners(z, lab, xs, ys):
    xg, yg = z["x_grid"], z["y_grid"]
    ix = int(np.clip(np.searchsorted(xg, xs)-1, 0, len(xg)-2)[0])
    iy = int(np.clip(np.searchsorted(yg, ys)-1, 0, len(yg)-2)[0])
    rec = {}
    for p in POPS:
        sc = np.nan_to_num(np.asarray(z[f"support_count__{p}__{lab}"], float))
        dr = np.nan_to_num(np.asarray(z[f"density_raw__{p}__{lab}"], float))
        nd = np.nan_to_num(np.asarray(z[f"nearest_dist__{p}__{lab}"], float), posinf=9e9)
        c4 = [(iy,ix),(iy,ix+1),(iy+1,ix),(iy+1,ix+1)]
        rec[f"{p}_sup_corners"] = [float(sc[a1,b1]) for a1,b1 in c4]
        rec[f"{p}_sup_min"] = float(min(sc[a1,b1] for a1,b1 in c4))
        rec[f"{p}_dens_raw_mean"] = float(np.mean([dr[a1,b1] for a1,b1 in c4]))
        rec[f"{p}_dens_masked_mean"] = float(np.mean(
            [dr[a1,b1] if sc[a1,b1] >= 1 or p == "NEO" else 0.0 for a1,b1 in c4]))
        rec[f"{p}_nearest_min"] = float(min(nd[a1,b1] for a1,b1 in c4))
    return rec

rows = []
for gi in flips_idx:
    r = cal.iloc[gi]
    lab = next(l for l,lo,hi in BINS if lo <= r.mean_mag < hi)
    e = dict(idx=int(gi), truth=r.population, map=r.prob_map_file, magbin=lab,
             vlam=float(r.vlam), vbeta=float(r.vbeta),
             P_001_maskon=float(sA1[gi]), P_0005_maskon=float(sB1[gi]),
             P_001_maskoff=float(sA0[gi]), P_0005_maskoff=float(sB0[gi]))
    for dirp, tag in [(A_DIR,"c001"), (B_DIR,"c0005")]:
        z = np.load(dirp/r.prob_map_file, allow_pickle=True)
        for k, vv in corners(z, lab, np.array([r.vlam]), np.array([r.vbeta])).items():
            e[f"{tag}_{k}"] = vv
    rows.append(e)
F = pd.DataFrame(rows)
F.to_csv(OUT/"SUPPORT_MASK_FLIP_CORNERS.csv", index=False)
has_coarse = F.c001_MBA_sup_min > 0
zero_fine = F.c0005_MBA_sup_min == 0
print(f"\n=== flip anatomy ({len(F)} flips) ===")
print(f"  MBA support > 0 at ALL coarse corners : {int(has_coarse.sum())}")
print(f"  MBA support == 0 at >=1 fine corner   : {int((F.c0005_MBA_sup_min==0).sum())}")
print(f"  BOTH (coarse supported -> fine hole)  : {int((has_coarse & zero_fine).sum())}  "
      f"({100*(has_coarse & zero_fine).mean():.1f}%)")
print(f"  truth MBA among flips: {int((F.truth=='MBA').sum())}/{len(F)}")

print("\n=== ten largest MBA changes, with and without the mask ===")
mb = F[F.truth == "MBA"].assign(d_on=lambda t: (t.P_0005_maskon-t.P_001_maskon).abs(),
                                d_off=lambda t: (t.P_0005_maskoff-t.P_001_maskoff).abs())
print(mb.nlargest(10, "d_on")[["map","magbin","vlam","vbeta","P_001_maskon","P_0005_maskon",
      "d_on","P_001_maskoff","P_0005_maskoff","d_off","c001_MBA_sup_min","c0005_MBA_sup_min"]
      ].to_string(index=False, float_format=lambda z: f"{z:,.5g}"))
json.dump(dict(mask_on=out["mask_on"], mask_off=out["mask_off"],
               flips_coarse_supported_fine_hole=int((has_coarse & zero_fine).sum()),
               n_flips=int(len(F))), open(OUT/"SUPPORT_MASK_MECHANISM.json","w"), indent=2)
print(f"\nwrote SUPPORT_MASK_MECHANISM.json and SUPPORT_MASK_FLIP_CORNERS.csv")
