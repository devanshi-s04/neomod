#!/usr/bin/env python3
"""FINAL CLOSURE: unmasked posterior, 0.01 vs 0.005, from existing products only.

Decision recorded upstream: occupancy masking REJECTED, standardized-Z gate REJECTED, no further
CAL-derived support gate to be sought. Posteriors are built from all unmasked population densities
symmetrically; abstention only on explicit technical domain failure.
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
POPS=["NEO","MBA","TNO","Trojans"]; THRS=[2,3,5,10]
cal=pd.read_parquet(W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet")
names={f"prob_maps_grid_dlon{a:+04d}_lat{b:+03d}.npz" for a,b in CENTERS}
cal=cal[cal.prob_map_file.isin(names)].reset_index(drop=True)
y=(cal.population=="NEO").astype(int).to_numpy()
vmax=np.maximum(cal.vlam.abs(),cal.vbeta.abs()).to_numpy()

def score(dirp,smin):
    s=np.full(len(cal),np.nan); c=np.zeros(len(cal),bool)
    for cen,g in cal.groupby("prob_map_file"):
        pm=v.ProbMapSet.from_npz(str(dirp/cen),support_mask_min=smin,mask_radius_deg_per_day=np.inf)
        r=pm.score_visible(g.vlam.to_numpy(float),g.vbeta.to_numpy(float),g.mean_mag.to_numpy(float))
        i=g.index.to_numpy(); s[i]=np.asarray(r["NEO"],float)
        c[i]=np.sum([np.nan_to_num(np.asarray(r[p],float)) for p in POPS],axis=0)>0
        del pm
    return s,c
cov={}
for t in THRS: _,cov[t]=score(W/f"prob_maps_e0_thr{t}",1)
common=np.logical_and.reduce([cov[t] for t in THRS])
sA,cA=score(A_DIR,0); sB,cB=score(B_DIR,0)
m=common; yv=y[m]; a=np.nan_to_num(sA[m]); b=np.nan_to_num(sB[m]); vv=vmax[m]
def mets(s):
    p_,r_,t_=precision_recall_curve(yv,s)
    f1=np.divide(2*p_*r_,p_+r_,out=np.zeros_like(p_),where=(p_+r_)>0); i=int(np.argmax(f1[:-1]))
    return dict(ROC=roc_auc_score(yv,s),pAUC=roc_auc_score(yv,s,max_fpr=0.01),F1=f1[i],
                Brier=brier_score_loss(yv,np.clip(s,0,1)),thr=float(t_[i]) if i<len(t_) else 1.0,
                completeness=r_[i]*100,contamination=(1-p_[i])*100)
MA,MB=mets(a),mets(b)
pd.set_option("display.width",250)
print("="*104); print("1. UNMASKED 0.01 vs 0.005 -- frozen rows", f"{int(m.sum()):,} / {int(yv.sum())} NEOs")
T=pd.DataFrame([dict(step=0.01,**MA),dict(step=0.005,**MB)])
print(T.to_string(index=False,float_format=lambda z:f"{z:,.6g}"))
print(f"  coverage {100*cA.mean():.4f}% vs {100*cB.mean():.4f}%")

d=np.abs(b-a); thr=MA["thr"]; fl=(a>=thr)!=(b>=thr)
print("\n"+"="*104); print("2. SCORE DIFFERENCES")
S=pd.DataFrame({"truth":cal.population.to_numpy()[m],"mag":cal.mean_mag.to_numpy()[m],
                "vmax":vv,"absd":d,"flip":fl})
for lab,key in [("population","truth"),("magnitude",pd.cut(S.mag,[14,20,22,23,24,25])),
                ("velocity",pd.cut(S.vmax,[0,0.25,0.5,1.0,2.0,5.1]))]:
    print(f"\n  by {lab}:")
    print(S.groupby(key,observed=True).agg(n=("absd","size"),median=("absd","median"),
        p99=("absd",lambda z:np.percentile(z,99)),max=("absd","max"),flips=("flip","sum")
        ).to_string(float_format=lambda z:f"{z:,.3e}"))
print("\n"+"="*104); print("3. CALIBRATION (Brier) BY VELOCITY BAND")
vb=pd.cut(S.vmax,[0,0.25,0.5,1.0,2.0,5.1])
rows=[]
for lab,g in pd.DataFrame({"vb":vb,"y":yv,"a":a,"b":b}).groupby("vb",observed=True):
    if g.y.nunique()<2: 
        rows.append(dict(vband=str(lab),n=len(g),NEOfrac=g.y.mean(),
                         Brier_001=brier_score_loss(g.y,np.clip(g.a,0,1)) if g.y.nunique()>1 else np.nan,
                         Brier_0005=np.nan)); continue
    rows.append(dict(vband=str(lab),n=len(g),NEOfrac=g.y.mean(),
        Brier_001=brier_score_loss(g.y,np.clip(g.a,0,1)),
        Brier_0005=brier_score_loss(g.y,np.clip(g.b,0,1))))
print(pd.DataFrame(rows).to_string(index=False,float_format=lambda z:f"{z:,.5g}"))
print("\n"+"="*104); print("4. THE TWO RESIDUAL FLIPS (unmasked)")
fi=np.where(m)[0][fl]
F=cal.loc[fi,["ObjID","population","prob_map_file","mean_mag","vlam","vbeta"]].copy()
F["P_001"]=sA[fi]; F["P_0005"]=sB[fi]; F["absd"]=np.abs(sB[fi]-sA[fi])
F["vmax"]=np.maximum(F.vlam.abs(),F.vbeta.abs()); F["threshold"]=thr
print(F.to_string(index=False,float_format=lambda z:f"{z:,.6g}"))
print(f"\n  both sides of threshold {thr:.6g}: "
      f"{[(round(float(x),6),round(float(z),6)) for x,z in zip(F.P_001,F.P_0005)]}")
print("\n"+"="*104); print("5. FAST-NEO (|v|>2) COMPLETENESS")
fast=(vv>2)&(yv==1)
print(f"  |v|>2 NEOs in frozen set: {int(fast.sum())}")
for nm,s_,M in [("0.01",a,MA),("0.005",b,MB)]:
    print(f"  step {nm}: median P {np.median(s_[fast]):.6f}   frac>=0.5 {100*np.mean(s_[fast]>=0.5):.1f}%"
          f"   frac>=thr {100*np.mean(s_[fast]>=M['thr']):.1f}%   zero {100*np.mean(s_[fast]==0):.1f}%")
szA=sum(f.stat().st_size for f in A_DIR.glob("*.npz")); szB=sum(f.stat().st_size for f in B_DIR.glob("*.npz"))
print("\n"+"="*104); print("6. COST")
print(f"  16 maps: 0.01 = {szA/2**30:.2f} GiB   0.005 = {szB/2**30:.2f} GiB  ({szB/szA:.2f}x)")
print(f"  667 maps projected: {szA/16*667/2**30:.0f} GiB vs {szB/16*667/2**30:.0f} GiB")
print(f"  build time/map: ~9 min (0.01) vs 13-42 min (0.005); peak RSS 13-17 GB at 0.005")
json.dump(dict(metrics_001=MA,metrics_0005=MB,coverage=[100*float(cA.mean()),100*float(cB.mean())],
    flips=int(fl.sum()),median_absd=float(np.median(d)),max_absd=float(d.max()),
    fast_neo=dict(n=int(fast.sum()),median_001=float(np.median(a[fast])),
                  median_0005=float(np.median(b[fast]))),
    storage_gib=[szA/2**30,szB/2**30],projected_667_gib=[szA/16*667/2**30,szB/16*667/2**30]),
    open(OUT/"E0_CLOSURE.json","w"),indent=2)
F.to_csv(OUT/"E0_CLOSURE_RESIDUAL_FLIPS.csv",index=False)
T.to_csv(OUT/"E0_CLOSURE_METRICS.csv",index=False)
print("\nwrote E0_CLOSURE.json / E0_CLOSURE_METRICS.csv / E0_CLOSURE_RESIDUAL_FLIPS.csv")
