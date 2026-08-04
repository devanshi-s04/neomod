#!/usr/bin/env python3
"""THE SINGLE FINAL TEST EVALUATION. Everything is read from MODEL_SEAL.json; nothing is tuned.
Scores TEST once with frozen variant C + Platt and applies the three CAL thresholds UNCHANGED.
"""
import hashlib, json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from scipy.special import expit, logit
from sklearn.metrics import roc_auc_score, precision_recall_curve, brier_score_loss, log_loss
W=Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0,str(W/"neomod"/"src")); os.chdir(W/"neomod")
import velocity_density_pipeline_neomod_clone_only as v
S=W/"outputs/splits"; OUT=W/"outputs/test_final"; OUT.mkdir(parents=True,exist_ok=True)
SEAL=json.loads((S/"MODEL_SEAL.json").read_text())
A_,B_=SEAL["calibration"]["a"],SEAL["calibration"]["b"]
THR=SEAL["operating_thresholds_from_CAL"]
D=W/SEAL["map_build"]["maps_dir"]; POPS=["NEO","MBA","TNO","Trojans"]; EPS=1e-12
BINS=[("14_16",14.,16.),("16_18",16.,18.),("18_20",18.,20.),("mag20",20.,21.),
      ("mag21",21.,22.),("mag22",22.,23.),("mag23",23.,24.),("mag24+",24.,25.)]
CTR=np.array([(lo+hi)/2 for _,lo,hi in BINS]); LABS=[b for b,_,_ in BINS]
print(f"MODEL_SEAL a={A_} b={B_}"); print(f"frozen thresholds: {THR['best_F1_DIAGNOSTIC_ONLY']:.6f} "
      f"(best-F1 diag), {THR['contamination_le_5pct']:.6f} (<=5%), {THR['contamination_le_10pct']:.6f} (<=10%)")
te=pd.read_parquet(W/"outputs/test_tracklets_neomod3/tracklets_benchmark_neomod3.parquet")
y=(te.population=="NEO").astype(int).to_numpy()
vmax=np.maximum(te.vlam.abs(),te.vbeta.abs()).to_numpy()
print(f"TEST rows {len(te):,}  NEO {int(y.sum()):,}  cells {te.prob_map_file.nunique()}")

def bil(a,x,yy,xs,ys):
    ix=np.clip(np.searchsorted(x,xs)-1,0,len(x)-2); iy=np.clip(np.searchsorted(yy,ys)-1,0,len(yy)-2)
    tx=np.clip((xs-x[ix])/(x[ix+1]-x[ix]),0,1); ty=np.clip((ys-yy[iy])/(yy[iy+1]-yy[iy]),0,1)
    return (a[iy,ix]*(1-tx)*(1-ty)+a[iy,ix+1]*tx*(1-ty)+a[iy+1,ix]*(1-tx)*ty+a[iy+1,ix+1]*tx*ty)
def run(cen):
    g=te[te.prob_map_file==cen]; xs=g.vlam.to_numpy(float); ys=g.vbeta.to_numpy(float)
    mags=g.mean_mag.to_numpy(float); z=np.load(D/cen,allow_pickle=True)
    xg,yg=z["x_grid"],z["y_grid"]
    j=np.clip(np.searchsorted(CTR,mags)-1,0,len(CTR)-2); t=np.clip((mags-CTR[j])/(CTR[j+1]-CTR[j]),0,1)
    dens={p:np.zeros(len(g)) for p in POPS}
    for jj in np.unique(j):
        s_=j==jj
        for p in POPS:
            d0=bil(np.nan_to_num(np.asarray(z[f"density_raw__{p}__{LABS[jj]}"],float)),xg,yg,xs[s_],ys[s_])
            d1=bil(np.nan_to_num(np.asarray(z[f"density_raw__{p}__{LABS[jj+1]}"],float)),xg,yg,xs[s_],ys[s_])
            dens[p][s_]=d0*(1-t[s_])+d1*t[s_]
    tot=sum(dens[p] for p in POPS)
    praw=np.where(tot>0,dens["NEO"]/np.where(tot>0,tot,1),0.0)
    return g.index.to_numpy(),praw,tot
t0=time.time()
parts=Parallel(n_jobs=12,verbose=5)(delayed(run)(c) for c in sorted(te.prob_map_file.unique()))
P_raw=np.full(len(te),np.nan); TOT=np.full(len(te),np.nan)
for i,p_,t_ in parts: P_raw[i]=p_; TOT[i]=t_
inb=(te.vlam.abs()<=5)&(te.vbeta.abs()<=5)&(te.mean_mag>=14)&(te.mean_mag<25)
m=inb.to_numpy()&np.isfinite(P_raw)&np.isfinite(TOT)&(TOT>0)
P_cal=np.full(len(te),np.nan); P_cal[m]=expit(A_*logit(np.clip(P_raw[m],EPS,1-EPS))+B_)
print(f"\nscored in {(time.time()-t0)/60:.1f} min")
print("\n"+"="*100); print("1. TECHNICAL COVERAGE")
print(f"  overall {100*m.mean():.4f}%  ({int(m.sum()):,}/{len(te):,})")
cov=pd.DataFrame({"truth":te.population,"covered":m,"vmax":vmax,"mag":te.mean_mag})
print("\n  by truth:"); print(cov.groupby("truth")["covered"].agg(n="size",covered="sum",
      pct=lambda t:100*t.mean()).to_string(float_format=lambda z:f"{z:,.4f}"))
print("\n  by velocity band:")
print(cov.groupby(pd.cut(cov.vmax,[0,0.25,0.5,1,2,5,1e3]),observed=True)["covered"].agg(
      n="size",covered="sum",pct=lambda t:100*t.mean()).to_string(float_format=lambda z:f"{z:,.4f}"))
print("\n  by magnitude:")
print(cov.groupby(pd.cut(cov.mag,[14,20,22,23,24,25]),observed=True)["covered"].agg(
      n="size",covered="sum",pct=lambda t:100*t.mean()).to_string(float_format=lambda z:f"{z:,.4f}"))
yv=y[m]; praw=P_raw[m]; pcal=P_cal[m]; vv=vmax[m]
pi_test=float(yv.mean())
print(f"\n  pi_TEST = {100*pi_test:.4f}%   (pi_CAL was 0.7318%)")
print("\n"+"="*100); print("2. VDP RANKING (raw == calibrated; Platt is strictly monotonic)")
rk=dict(ROC=roc_auc_score(yv,praw),pAUC=roc_auc_score(yv,praw,max_fpr=0.01))
pr,rc,th=precision_recall_curve(yv,praw)
f1=np.divide(2*pr*rc,pr+rc,out=np.zeros_like(pr),where=(pr+rc)>0); rk["F1_max"]=f1[:-1].max()
print(f"  ROC AUC {rk['ROC']:.6f}   std partial AUC (FPR<=0.01) {rk['pAUC']:.6f}   max F1 {rk['F1_max']:.6f}")
print(f"  calibrated ROC {roc_auc_score(yv,pcal):.6f}  (|d| {abs(roc_auc_score(yv,pcal)-rk['ROC']):.2e})")
print("\n"+"="*100); print("3. CALIBRATED VDP: Brier / log loss / ECE / reliability")
def ece(yy,pp,nb=15):
    e=np.linspace(0,1,nb+1); idx=np.clip(np.digitize(pp,e)-1,0,nb-1); tot_=0; rows=[]
    for b in range(nb):
        s_=idx==b
        if not s_.any(): continue
        o=yy[s_].mean(); c=pp[s_].mean(); n=int(s_.sum()); tot_+=s_.mean()*abs(o-c)
        zc=1.959963985; p_=o; dd=1+zc*zc/n; ct=(p_+zc*zc/(2*n))/dd
        h=zc/dd*np.sqrt(p_*(1-p_)/n+zc*zc/(4*n*n))
        rows.append(dict(bin=b,n=n,pred=c,obs=o,wilson_lo=max(0,ct-h),wilson_hi=min(1,ct+h)))
    return tot_,pd.DataFrame(rows)
for nm,pp in [("raw",praw),("calibrated",pcal)]:
    e_,_=ece(yv,np.clip(pp,EPS,1-EPS))
    print(f"  {nm:11s} Brier {brier_score_loss(yv,np.clip(pp,0,1)):.6e}   "
          f"log loss {log_loss(yv,np.clip(pp,EPS,1-EPS)):.6e}   ECE {e_:.6e}")
_,rel=ece(yv,np.clip(pcal,EPS,1-EPS))
print("\n  reliability (calibrated, 15 equal-width bins, Wilson 95%):")
print(rel.to_string(index=False,float_format=lambda z:f"{z:,.4f}"))
print("\n"+"="*100); print("4. FROZEN CAL THRESHOLDS APPLIED TO TEST UNCHANGED")
rows=[]
for nm,t_ in [("best-F1 (diagnostic)",THR["best_F1_DIAGNOSTIC_ONLY"]),
              ("CAL contamination<=5%",THR["contamination_le_5pct"]),
              ("CAL contamination<=10%",THR["contamination_le_10pct"])]:
    sel=pcal>=t_; tp=int((sel&(yv==1)).sum()); fp=int((sel&(yv==0)).sum())
    rows.append(dict(rule=nm,threshold=t_,n_selected=int(sel.sum()),TP=tp,FP=fp,
        completeness=100*tp/max(int(yv.sum()),1),contamination=100*fp/max(tp+fp,1),
        TPR=tp/max(int(yv.sum()),1),FPR=fp/max(int((yv==0).sum()),1)))
OPS=pd.DataFrame(rows); print(OPS.to_string(index=False,float_format=lambda z:f"{z:,.6g}"))
np.savez_compressed(OUT/"TEST_SCORES.npz",P_raw=P_raw,P_cal=P_cal,mask=m,y=y,vmax=vmax)
out=te[["ObjID","population","prob_map_file","mean_mag","vlam","vbeta","ra0","dec0","mjd0_utc",
        "ra1","dec1","mjd1_utc","mag0","mag1"]].copy()
out["P_raw"]=P_raw; out["P_cal"]=P_cal; out["technical_valid"]=m; out["is_neo"]=y
out.to_parquet(OUT/"TEST_SCORED.parquet",index=False)
json.dump(dict(model_seal_sha256=hashlib.sha256((S/"MODEL_SEAL.json").read_bytes()).hexdigest(),
    n_rows=int(len(te)),n_neo=int(y.sum()),n_valid=int(m.sum()),pi_test=pi_test,
    coverage_overall=100*float(m.mean()),ranking=rk,operating=OPS.to_dict("records"),
    brier_cal=float(brier_score_loss(yv,np.clip(pcal,0,1))),
    logloss_cal=float(log_loss(yv,np.clip(pcal,EPS,1-EPS))),ece_cal=float(ece(yv,np.clip(pcal,EPS,1-EPS))[0])),
    open(OUT/"TEST_RESULTS.json","w"),indent=2)
print(f"\nwrote {OUT}/TEST_SCORED.parquet, TEST_SCORES.npz, TEST_RESULTS.json")
