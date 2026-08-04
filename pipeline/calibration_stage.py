#!/usr/bin/env python3
"""Calibration stage: 4-candidate OOF comparison on CAL. Manifest sha f8bdf0f399587593...
TEST sealed. No MODEL_SEAL written here.
"""
import glob, json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from scipy.optimize import minimize
from scipy.special import expit, logit
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, precision_recall_curve, brier_score_loss, log_loss
W=Path("/mmfs1/gscratch/dirac/ds2004/sorcha"); OUT=W/"outputs/calibration"; OUT.mkdir(parents=True,exist_ok=True)
D=W/"prob_maps_grid_neomod3_GEN_final"; POPS=["NEO","MBA","TNO","Trojans"]
BINS=["14_16","16_18","18_20","mag20","mag21","mag22","mag23","mag24+"]
EPS=1e-12; NB=15
cal=pd.read_parquet(W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet")
z=np.load(W/"outputs/e1_results/E1_PER_ROW_SCORES.npz")
P_raw=z["C"]; m=z["mask"]; y=z["y"]
print(f"variant C scores loaded: {len(P_raw):,} rows; evaluation mask {int(m.sum()):,}")

# ---- pi_ref: NEO prior implied by the density normalisation ----
def integ(f):
    zz=np.load(f,allow_pickle=True); x=zz["x_grid"]; dA=float((x[1]-x[0])**2)
    out={}
    for p in POPS:
        out[p]=float(sum(np.nan_to_num(np.asarray(zz[f"density_raw__{p}__{b}"],float)).sum()
                         for b in BINS if f"density_raw__{p}__{b}" in zz.files)*dA)
    return out
files=sorted(glob.glob(str(D/"*.npz")))
t0=time.time()
ints=Parallel(n_jobs=32,verbose=1)(delayed(integ)(f) for f in files)
tot={p:sum(d[p] for d in ints) for p in POPS}
pi_ref=tot["NEO"]/sum(tot.values())
per_map=np.array([d["NEO"]/max(sum(d.values()),EPS) for d in ints])
print(f"\n=== pi_ref (density-implied NEO prior) ===")
print(f"  integrated over {len(files)} maps x 8 bins x +-5 grid  ({time.time()-t0:.0f}s)")
print(f"  pi_ref = {pi_ref:.6e}   ({100*pi_ref:.4f}%)")
print(f"  per-map spread: median {np.median(per_map):.3e}  p5 {np.percentile(per_map,5):.3e}  "
      f"p95 {np.percentile(per_map,95):.3e}   (30-deg patches overlap at 10-deg spacing)")
yv=y[m]; pi_cal=float(yv.mean())
print(f"  pi_CAL = {pi_cal:.6e}  ({100*pi_cal:.4f}%)")
print(f"  pi_ref / pi_CAL = {pi_ref/pi_cal:.4f}  -> "
      f"{'IDENTICAL' if abs(pi_ref-pi_cal)<1e-9 else 'DIFFERENT -- pi_CAL must NOT be substituted for pi_ref'}")

# ---- fold assignment: grouped by sky cell, no parent crossing ----
sub=cal[m].reset_index(drop=True); s_raw=logit(np.clip(P_raw[m],EPS,1-EPS))
assert sub.ObjID.is_unique, "ObjID repeats -- grouped-by-cell folds could split a parent object"
print(f"\n  fold assertion: ObjID unique across {len(sub):,} rows -> no parent crosses folds  PASS")
cells=sorted(sub.prob_map_file.unique())
cn=sub.groupby("prob_map_file").apply(lambda t:(t.population=="NEO").sum(),include_groups=False)
order=cn.sort_values(ascending=False).index.tolist()
fold_of={}; load=np.zeros(5)
for c in order:                       # greedy balance of NEO counts across folds
    k=int(np.argmin(load)); fold_of[c]=k; load[k]+=cn[c]
fold=sub.prob_map_file.map(fold_of).to_numpy()
print(f"  5 folds grouped by sky cell; NEO per fold: {[int(v) for v in load]}")

def fit_platt(s,yy):
    def nll(th):
        p=expit(np.exp(th[0])*s+th[1]); p=np.clip(p,EPS,1-EPS)
        return -(yy*np.log(p)+(1-yy)*np.log(1-p)).mean()
    r=minimize(nll,[0.0,0.0],method="Nelder-Mead",options=dict(maxiter=2000,xatol=1e-8,fatol=1e-10))
    return dict(a=float(np.exp(r.x[0])),b=float(r.x[1]))
def fit_temp(s,yy):
    def nll(th):
        p=expit(s/np.exp(th[0])); p=np.clip(p,EPS,1-EPS)
        return -(yy*np.log(p)+(1-yy)*np.log(1-p)).mean()
    r=minimize(nll,[0.0],method="Nelder-Mead",options=dict(maxiter=2000,xatol=1e-8,fatol=1e-10))
    return dict(T=float(np.exp(r.x[0])))

oof={"R":P_raw[m].copy(),"P":np.zeros(len(sub)),"T":np.zeros(len(sub)),"I":np.zeros(len(sub))}
params=[]
for k in range(5):
    tr=fold!=k; te=fold==k
    pp=fit_platt(s_raw[tr],yv[tr]); tt=fit_temp(s_raw[tr],yv[tr])
    iso=IsotonicRegression(out_of_bounds="clip",y_min=0,y_max=1).fit(P_raw[m][tr],yv[tr])
    oof["P"][te]=expit(pp["a"]*s_raw[te]+pp["b"])
    oof["T"][te]=expit(s_raw[te]/tt["T"])
    oof["I"][te]=iso.predict(P_raw[m][te])
    params.append(dict(fold=k,n_train=int(tr.sum()),n_test=int(te.sum()),
                       neo_test=int(yv[te].sum()),**pp,**tt,iso_knots=len(iso.f_.x)))
PA=pd.DataFrame(params)
print("\n=== per-fold fitted parameters ===")
print(PA.to_string(index=False,float_format=lambda v:f"{v:,.6g}"))

def ece(yy,pp,nb=NB):
    e=np.linspace(0,1,nb+1); idx=np.clip(np.digitize(pp,e)-1,0,nb-1); tot_=0.0; rows=[]
    for b in range(nb):
        s_=idx==b
        if not s_.any(): rows.append((b,0,np.nan,np.nan)); continue
        o=yy[s_].mean(); c=pp[s_].mean(); tot_+=s_.mean()*abs(o-c); rows.append((b,int(s_.sum()),c,o))
    return tot_,rows
rng=np.random.default_rng(0); i1=np.where(yv==1)[0]; i0=np.where(yv==0)[0]
rows=[]
for k in ["R","P","T","I"]:
    p_=np.clip(oof[k],EPS,1-EPS); e_,_=ece(yv,p_)
    bs=[]
    for _ in range(500):
        b_=np.concatenate([rng.choice(i1,len(i1),True),rng.choice(i0,len(i0),True)])
        bs.append(ece(yv[b_],p_[b_])[0])
    pr,rc,_=precision_recall_curve(yv,p_)
    f1=np.divide(2*pr*rc,pr+rc,out=np.zeros_like(pr),where=(pr+rc)>0)
    rows.append(dict(cand=k,Brier=brier_score_loss(yv,p_),logloss=log_loss(yv,p_),
        ECE=e_,ECE_lo=np.percentile(bs,2.5),ECE_hi=np.percentile(bs,97.5),
        ROC=roc_auc_score(yv,p_),pAUC=roc_auc_score(yv,p_,max_fpr=0.01),F1=f1[:-1].max(),
        n_unique=len(np.unique(p_))))
R=pd.DataFrame(rows); pd.set_option("display.width",250)
print("\n=== OOF CALIBRATION COMPARISON (primary: Brier, log loss) ===")
print(R[["cand","Brier","logloss","ECE","ECE_lo","ECE_hi"]].to_string(index=False,float_format=lambda v:f"{v:,.6g}"))
print("\n=== INVARIANCE DIAGNOSTIC ONLY (monotonic transforms -- NOT improvement) ===")
print(R[["cand","ROC","pAUC","F1","n_unique"]].to_string(index=False,float_format=lambda v:f"{v:,.6g}"))
print("\n=== reliability, 15 equal-width bins (raw R vs best-by-Brier) ===")
best=R.loc[R.Brier.idxmin(),"cand"]
for k in ["R",best] if best!="R" else ["R"]:
    _,rr=ece(yv,np.clip(oof[k],EPS,1-EPS))
    t=pd.DataFrame(rr,columns=["bin","n","pred","obs"])
    t=t[t.n>0]; t["lo"]=t.apply(lambda r:max(0,r.obs-1.96*np.sqrt(max(r.obs*(1-r.obs),1e-12)/r.n)),axis=1)
    t["hi"]=t.apply(lambda r:min(1,r.obs+1.96*np.sqrt(max(r.obs*(1-r.obs),1e-12)/r.n)),axis=1)
    print(f"\n  {k}:"); print(t.to_string(index=False,float_format=lambda v:f"{v:,.4f}"))
print(f"\n  best by Brier: {best}   best by log loss: {R.loc[R.logloss.idxmin(),'cand']}")
np.savez_compressed(OUT/"CAL_OOF_SCORES.npz",**oof,y=yv,fold=fold)
R.to_csv(OUT/"CALIBRATION_COMPARISON.csv",index=False); PA.to_csv(OUT/"CALIBRATION_FOLD_PARAMS.csv",index=False)
json.dump(dict(pi_ref=pi_ref,pi_cal=pi_cal,pi_ratio=pi_ref/pi_cal,
    pi_ref_identical_to_pi_cal=bool(abs(pi_ref-pi_cal)<1e-9),
    pi_ref_per_map=dict(median=float(np.median(per_map)),p5=float(np.percentile(per_map,5)),
                        p95=float(np.percentile(per_map,95))),
    results=R.to_dict("records"),fold_params=PA.to_dict("records"),
    ece_rule="15 equal-width bins on [0,1]; 500-replicate truth-stratified bootstrap CI",
    label_shift_assumption="P_target = sigmoid(logit(P_cal) - logit(pi_CAL) + logit(pi_target)); "
                           "valid only under label shift (p(X|y) unchanged)"),
    open(OUT/"CALIBRATION_COMPARISON.json","w"),indent=2)
print(f"\nwrote {OUT}/CALIBRATION_COMPARISON.csv/.json, FOLD_PARAMS.csv, CAL_OOF_SCORES.npz")
