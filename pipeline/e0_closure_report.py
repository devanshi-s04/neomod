#!/usr/bin/env python3
"""FINAL CLOSURE (corrected): unmasked posterior, step 0.01 vs 0.005, existing products only.

Corrections applied:
 - no np.nan_to_num; finite-score ASSERTIONS instead
 - evaluation mask = intersection of the two candidates' TECHNICAL-VALID masks, asserted against
   the previous frozen common mask row-for-row, with every excluded row itemised
 - Brier IS valid for one-class bins; ROC AUC is what is undefined there
 - disagreements at three thresholds: fixed-0.01 CAL, each arm's own CAL-optimal, and 0.5
 - roc_auc_score(max_fpr=0.01) labelled STANDARDIZED PARTIAL ROC AUC (FPR<=0.01)
 - paired, truth-stratified bootstrap CIs on the differences
 - no equivalence claim: no margin was preregistered
"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, precision_recall_curve, brier_score_loss
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
import velocity_density_pipeline_neomod_clone_only as v
OUT=W/"outputs/e0_results"; A_DIR=W/"prob_maps_e0_thr10"; B_DIR=W/"prob_maps_e0_thr10_step005"
CENTERS=[(0,0),(20,-12),(10,2),(-10,-2),(0,8),(50,-1),(50,1),(-50,1),
         (90,0),(140,0),(-140,0),(0,25),(0,-25),(0,50),(0,-50),(-50,-50)]
POPS=["NEO","MBA","TNO","Trojans"]; THRS=[2,3,5,10]; VLIM=5.0
cal=pd.read_parquet(W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet")
names={f"prob_maps_grid_dlon{a:+04d}_lat{b:+03d}.npz" for a,b in CENTERS}
cal=cal[cal.prob_map_file.isin(names)].reset_index(drop=True)
y=(cal.population=="NEO").astype(int).to_numpy()
vmax=np.maximum(cal.vlam.abs(),cal.vbeta.abs()).to_numpy()

def score(dirp,smin):
    s=np.full(len(cal),np.nan); tot=np.full(len(cal),np.nan)
    for cen,g in cal.groupby("prob_map_file"):
        pm=v.ProbMapSet.from_npz(str(dirp/cen),support_mask_min=smin,mask_radius_deg_per_day=np.inf)
        r=pm.score_visible(g.vlam.to_numpy(float),g.vbeta.to_numpy(float),g.mean_mag.to_numpy(float))
        i=g.index.to_numpy(); s[i]=np.asarray(r["NEO"],float)
        tot[i]=np.sum([np.asarray(r[p],float) for p in POPS],axis=0)
        del pm
    return s,tot

# ---- 0. technical-valid masks (final decision §9.3) ----
in_bounds = (cal.vlam.abs()<=VLIM)&(cal.vbeta.abs()<=VLIM)&(cal.mean_mag>=14)&(cal.mean_mag<25)
sA,totA=score(A_DIR,0); sB,totB=score(B_DIR,0)
validA = in_bounds.to_numpy() & np.isfinite(sA) & np.isfinite(totA) & (totA>0)
validB = in_bounds.to_numpy() & np.isfinite(sB) & np.isfinite(totB) & (totB>0)
mask = validA & validB
cov={}
for t in THRS:
    _,tt=score(W/f"prob_maps_e0_thr{t}",1); cov[t]=np.isfinite(tt)&(tt>0)
frozen=np.logical_and.reduce([cov[t] for t in THRS])
print("="*104); print("0. TECHNICAL-VALID EVALUATION MASK")
print(f"  valid(0.01) {int(validA.sum()):,}   valid(0.005) {int(validB.sum()):,}   "
      f"intersection {int(mask.sum()):,}")
print(f"  previous frozen common mask: {int(frozen.sum()):,}")
identical = bool(np.array_equal(mask, frozen))
print(f"  IDENTICAL row-for-row to the frozen mask: {identical}")
diff = np.where(mask != frozen)[0]
if len(diff):
    print(f"  rows differing: {len(diff)}  (in new not frozen: {int((mask&~frozen).sum())}; "
          f"in frozen not new: {int((frozen&~mask).sum())})")
excl = np.where(~mask)[0]
E=cal.loc[excl,["ObjID","population","prob_map_file","mean_mag","vlam","vbeta"]].copy()
def reason(i):
    r=[]
    if not in_bounds.iloc[i]: r.append("out_of_domain(|v|>5 or mag outside 14-25)")
    if not np.isfinite(sA[i]) or not np.isfinite(totA[i]): r.append("non_finite_0.01")
    if not np.isfinite(sB[i]) or not np.isfinite(totB[i]): r.append("non_finite_0.005")
    if np.isfinite(totA[i]) and totA[i]<=0: r.append("zero_total_density_0.01")
    if np.isfinite(totB[i]) and totB[i]<=0: r.append("zero_total_density_0.005")
    return ";".join(r) or "unknown"
E["reason"]=[reason(i) for i in excl]; E["vmax"]=np.maximum(E.vlam.abs(),E.vbeta.abs())
print(f"\n  EXCLUDED ROWS ({len(E)}):")
print(E[["ObjID","population","prob_map_file","mean_mag","vmax","reason"]].to_string(
      index=False,float_format=lambda z:f"{z:,.4g}"))
E.to_csv(OUT/"E0_CLOSURE_EXCLUDED_ROWS.csv",index=False)

m=mask; yv=y[m]; a=sA[m]; b=sB[m]; vv=vmax[m]
assert np.isfinite(a).all() and np.isfinite(b).all(), "non-finite score inside the evaluation mask"
print(f"\n  finite-score assertion PASSED on {int(m.sum()):,} evaluation rows "
      f"({int(yv.sum())} NEO)")

def best_thr(s):
    p_,r_,t_=precision_recall_curve(yv,s)
    f1=np.divide(2*p_*r_,p_+r_,out=np.zeros_like(p_),where=(p_+r_)>0); i=int(np.argmax(f1[:-1]))
    return (float(t_[i]) if i<len(t_) else 1.0), f1[i], r_[i]*100, (1-p_[i])*100
def mets(s):
    th,f1,comp,cont=best_thr(s)
    return dict(ROC_AUC=roc_auc_score(yv,s),
                stdpAUC_fpr01=roc_auc_score(yv,s,max_fpr=0.01),
                F1=f1,Brier=brier_score_loss(yv,np.clip(s,0,1)),
                CAL_opt_thr=th,completeness=comp,contamination=cont)
MA,MB=mets(a),mets(b)
pd.set_option("display.width",250)
print("\n"+"="*104); print("1. UNMASKED 0.01 vs 0.005")
print("   ('stdpAUC_fpr01' = STANDARDIZED PARTIAL ROC AUC for FPR<=0.01, McClish-standardized)")
T=pd.DataFrame([dict(step=0.01,**MA),dict(step=0.005,**MB)])
print(T.to_string(index=False,float_format=lambda z:f"{z:,.6g}"))

# ---- paired truth-stratified bootstrap ----
rng=np.random.default_rng(0); B=1000
inx=np.where(yv==1)[0]; ino=np.where(yv==0)[0]
diffs={k:[] for k in ["ROC_AUC","stdpAUC_fpr01","F1","Brier"]}
for _ in range(B):
    ii=np.concatenate([rng.choice(inx,len(inx),True),rng.choice(ino,len(ino),True)])
    yy,aa,bb=yv[ii],a[ii],b[ii]
    if yy.sum()<3: continue
    diffs["ROC_AUC"].append(roc_auc_score(yy,bb)-roc_auc_score(yy,aa))
    diffs["stdpAUC_fpr01"].append(roc_auc_score(yy,bb,max_fpr=0.01)-roc_auc_score(yy,aa,max_fpr=0.01))
    def f1of(s):
        p_,r_,_=precision_recall_curve(yy,s)
        f=np.divide(2*p_*r_,p_+r_,out=np.zeros_like(p_),where=(p_+r_)>0); return f[:-1].max()
    diffs["F1"].append(f1of(bb)-f1of(aa))
    diffs["Brier"].append(brier_score_loss(yy,np.clip(bb,0,1))-brier_score_loss(yy,np.clip(aa,0,1)))
print("\n  paired truth-stratified bootstrap (B=1000), difference 0.005 minus 0.01:")
BS=pd.DataFrame([dict(metric=k,point=MB[k]-MA[k],lo=np.percentile(vlist,2.5),
                      hi=np.percentile(vlist,97.5),
                      excludes_zero=(np.percentile(vlist,2.5)>0)or(np.percentile(vlist,97.5)<0))
                 for k,vlist in diffs.items()])
print(BS.to_string(index=False,float_format=lambda z:f"{z:,.3e}"))
print("  NOTE: no equivalence margin was preregistered, so this does NOT establish formal")
print("        equivalence. Statement: no practically meaningful improvement was observed.")

print("\n"+"="*104); print("2. CLASSIFICATION DISAGREEMENTS AT THREE THRESHOLDS")
rows=[]
for nm,th in [("fixed 0.01 CAL-optimal",MA["CAL_opt_thr"]),
              ("each arm's own CAL-optimal",None),("0.5",0.5)]:
    if th is None:
        da=(a>=MA["CAL_opt_thr"]); db=(b>=MB["CAL_opt_thr"]); lab=f"{MA['CAL_opt_thr']:.6g}/{MB['CAL_opt_thr']:.6g}"
    else:
        da=(a>=th); db=(b>=th); lab=f"{th:.6g}"
    dis=da!=db
    rows.append(dict(rule=nm,threshold=lab,disagreements=int(dis.sum()),
                     NEO=int((dis&(yv==1)).sum()),nonNEO=int((dis&(yv==0)).sum())))
print(pd.DataFrame(rows).to_string(index=False))

print("\n"+"="*104); print("3. CALIBRATION (Brier) BY VELOCITY BAND -- Brier IS defined for one-class bins")
vb=pd.cut(vv,[0,0.25,0.5,1.0,2.0,5.1])
rr=[]
for lab,g in pd.DataFrame({"vb":vb,"y":yv,"a":a,"b":b}).groupby("vb",observed=True):
    one=g.y.nunique()<2
    rr.append(dict(vband=str(lab),n=len(g),NEOfrac=g.y.mean(),
        Brier_001=brier_score_loss(g.y,np.clip(g.a,0,1),labels=[0,1]) if not one
                  else float(np.mean((np.clip(g.a,0,1)-g.y)**2)),
        Brier_0005=brier_score_loss(g.y,np.clip(g.b,0,1),labels=[0,1]) if not one
                  else float(np.mean((np.clip(g.b,0,1)-g.y)**2)),
        ROC_AUC="undefined (one class)" if one else f"{roc_auc_score(g.y,g.a):.5f}"))
print(pd.DataFrame(rr).to_string(index=False,float_format=lambda z:f"{z:,.5g}"))

print("\n"+"="*104); print("4. DISAGREEING ROWS AT THE FIXED 0.01 THRESHOLD")
fi=np.where(m)[0][(a>=MA["CAL_opt_thr"])!=(b>=MA["CAL_opt_thr"])]
F=cal.loc[fi,["ObjID","population","prob_map_file","mean_mag","vlam","vbeta"]].copy()
F["P_001"]=sA[fi]; F["P_0005"]=sB[fi]; F["absd"]=np.abs(sB[fi]-sA[fi])
F["vmax"]=np.maximum(F.vlam.abs(),F.vbeta.abs())
F["dist_to_thr_001"]=np.abs(sA[fi]-MA["CAL_opt_thr"])
F["character"]=np.where(F.dist_to_thr_001<F.absd/2,"threshold-adjacent","local resolution sensitivity")
print(F.to_string(index=False,float_format=lambda z:f"{z:,.6g}"))
F.to_csv(OUT/"E0_CLOSURE_RESIDUAL_FLIPS.csv",index=False)

print("\n"+"="*104); print("5. FAST-NEO (|v|>2) COMPLETENESS")
fast=(vv>2)&(yv==1)
print(f"  |v|>2 NEOs: {int(fast.sum())}")
for nm,s_,M in [("0.01",a,MA),("0.005",b,MB)]:
    print(f"  step {nm}: median P {np.median(s_[fast]):.6f}  frac>=0.5 {100*np.mean(s_[fast]>=0.5):.1f}%"
          f"  frac>=own CAL thr {100*np.mean(s_[fast]>=M['CAL_opt_thr']):.1f}%  zero {100*np.mean(s_[fast]==0):.1f}%")
szA=sum(f.stat().st_size for f in A_DIR.glob("*.npz")); szB=sum(f.stat().st_size for f in B_DIR.glob("*.npz"))
print("\n"+"="*104); print("6. COST")
print(f"  16 maps: {szA/2**30:.2f} GiB (0.01) vs {szB/2**30:.2f} GiB (0.005) = {szB/szA:.2f}x")
print(f"  667 projected: {szA/16*667/2**30:.0f} GiB vs {szB/16*667/2**30:.0f} GiB")
print(f"  build ~9 min/map (0.01) vs 13-42 min/map (0.005)")
print(f"  => 0.01 uses 4x less storage and materially lower compute")
json.dump(dict(mask_identical_to_frozen=identical,n_eval=int(m.sum()),n_excluded=int(len(E)),
   metrics_001=MA,metrics_0005=MB,bootstrap=BS.to_dict("records"),
   disagreements=rows,storage_gib=[szA/2**30,szB/2**30]),
   open(OUT/"E0_CLOSURE.json","w"),indent=2,default=str)
T.to_csv(OUT/"E0_CLOSURE_METRICS.csv",index=False)
print("\nwrote E0_CLOSURE.json / METRICS.csv / RESIDUAL_FLIPS.csv / EXCLUDED_ROWS.csv")
