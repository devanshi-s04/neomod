#!/usr/bin/env python3
"""E1 interpolation ablation A-D on the sealed 667-map grid. CAL only; TEST sealed.
Manifest: docs/E1_INTERPOLATION_ABLATION_MANIFEST.md sha 8dbdf474d2d3dbd7...
"""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score, precision_recall_curve, brier_score_loss, log_loss
W=Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0,str(W/"neomod"/"src")); os.chdir(W/"neomod")
import velocity_density_pipeline_neomod_clone_only as v
D=W/"prob_maps_grid_neomod3_GEN_final"; OUT=W/"outputs/e1_results"; OUT.mkdir(parents=True,exist_ok=True)
POPS=["NEO","MBA","TNO","Trojans"]
BINS=[("14_16",14.,16.),("16_18",16.,18.),("18_20",18.,20.),("mag20",20.,21.),
      ("mag21",21.,22.),("mag22",22.,23.),("mag23",23.,24.),("mag24+",24.,25.)]
CTR=np.array([(lo+hi)/2 for _,lo,hi in BINS]); LABS=[b for b,_,_ in BINS]
LAT=sorted({float(x) for x in [0,1,2,3,4,5,8,12,18,25,35,50]}|{-float(x) for x in [0,1,2,3,4,5,8,12,18,25,35,50]})
DLON=[float(d) for d in range(-140,150,10)]
cal=pd.read_parquet(W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet")
assert {"dlon_from_antisun_deg","beta_deg"} <= set(cal.columns)
y=(cal.population=="NEO").astype(int).to_numpy()
vmax=np.maximum(cal.vlam.abs(),cal.vbeta.abs()).to_numpy()
print(f"CAL rows {len(cal):,}  NEO {int(y.sum()):,}  cells {cal.prob_map_file.nunique()}")
def mapname(d,l): return f"prob_maps_grid_dlon{int(round(d)):+04d}_lat{int(round(l)):+03d}.npz"

def bil(a,x,yy,xs,ys):
    ix=np.clip(np.searchsorted(x,xs)-1,0,len(x)-2); iy=np.clip(np.searchsorted(yy,ys)-1,0,len(yy)-2)
    tx=np.clip((xs-x[ix])/(x[ix+1]-x[ix]),0,1); ty=np.clip((ys-yy[iy])/(yy[iy+1]-yy[iy]),0,1)
    return (a[iy,ix]*(1-tx)*(1-ty)+a[iy,ix+1]*tx*(1-ty)+a[iy+1,ix]*(1-tx)*ty+a[iy+1,ix+1]*tx*ty)

_cache={}
def getmap(n):
    if n not in _cache:
        if len(_cache)>6: _cache.pop(next(iter(_cache)))
        _cache[n]=np.load(D/n,allow_pickle=True)
    return _cache[n]

def dens_at(n,lab,xs,ys):
    z=getmap(n); xg,yg=z["x_grid"],z["y_grid"]
    return {p:bil(np.nan_to_num(np.asarray(z[f"density_raw__{p}__{lab}"],float)),xg,yg,xs,ys) for p in POPS}

def dens_maginterp(n,mags,xs,ys):
    """density-first LINEAR interpolation between adjacent magnitude bins"""
    j=np.clip(np.searchsorted(CTR,mags)-1,0,len(CTR)-2)
    t=np.clip((mags-CTR[j])/(CTR[j+1]-CTR[j]),0,1)
    out={p:np.zeros(len(mags)) for p in POPS}
    for jj in np.unique(j):
        s=j==jj
        d0=dens_at(n,LABS[jj],xs[s],ys[s]); d1=dens_at(n,LABS[jj+1],xs[s],ys[s])
        for p in POPS: out[p][s]=d0[p]*(1-t[s])+d1[p]*t[s]
    return out

def run(cen):
    g=cal[cal.prob_map_file==cen]; xs=g.vlam.to_numpy(float); ys=g.vbeta.to_numpy(float)
    mags=g.mean_mag.to_numpy(float); idx=g.index.to_numpy()
    # D must interpolate from each object's TRUE sky position, not the assigned cell centre.
    # Using the cell coords put every object exactly on a grid node, so the bilinear weights
    # collapsed to 1.0 on its own cell and D degenerated to C (identical to every digit).
    obj_dl=g.dlon_from_antisun_deg.to_numpy(float); obj_la=g.beta_deg.to_numpy(float)
    res={}
    pm=v.ProbMapSet.from_npz(str(D/cen),support_mask_min=None,mask_radius_deg_per_day=np.inf)
    r=pm.score_visible(xs,ys,mags); res["A"]=np.asarray(r["NEO"],float)
    totA=np.sum([np.asarray(r[p],float) for p in POPS],axis=0); del pm
    dB=dens_at(cen,None,None,None) if False else None
    # B: nearest bin, density-first velocity interpolation
    lab=np.array(LABS)[np.clip(np.searchsorted([b[2] for b in BINS],mags),0,len(BINS)-1)]
    dB={p:np.zeros(len(g)) for p in POPS}
    for L in np.unique(lab):
        s=lab==L; d=dens_at(cen,L,xs[s],ys[s])
        for p in POPS: dB[p][s]=d[p]
    tB=sum(dB[p] for p in POPS); res["B"]=np.where(tB>0,dB["NEO"]/np.where(tB>0,tB,1),0.0)
    # C: + magnitude interpolation
    dC=dens_maginterp(cen,mags,xs,ys); tC=sum(dC[p] for p in POPS)
    res["C"]=np.where(tC>0,dC["NEO"]/np.where(tC>0,tC,1),0.0)
    # D: + sky-cell interpolation (bilinear over the 4 surrounding cells)
    DA=np.array(DLON); LA=np.array(LAT)
    id_=np.clip(np.searchsorted(DA,obj_dl),1,len(DA)-1); il_=np.clip(np.searchsorted(LA,obj_la),1,len(LA)-1)
    d0,d1=DA[id_-1],DA[id_]; l0,l1=LA[il_-1],LA[il_]
    wd=np.where(d1>d0,(obj_dl-d0)/np.where(d1>d0,d1-d0,1),0.0)
    wl=np.where(l1>l0,(obj_la-l0)/np.where(l1>l0,l1-l0,1),0.0)
    wd=np.clip(wd,0,1); wl=np.clip(wl,0,1)
    # hull membership: an object outside the sky-centre hull is CLIPPED to the boundary map,
    # not genuinely interpolated. Track it so D is not credited with interpolation it did not do.
    out_lat=(obj_la<LA.min())|(obj_la>LA.max())
    out_dlon=(obj_dl<DA.min())|(obj_dl>DA.max())
    outside_hull=out_lat|out_dlon
    n_corners=np.zeros(len(g),int)
    acc={p:np.zeros(len(g)) for p in POPS}; wsum=np.zeros(len(g))
    for dsel,wdd in [(d0,1-wd),(d1,wd)]:
        for lsel,wll in [(l0,1-wl),(l1,wl)]:
            wgt=wdd*wll
            for nm in np.unique([mapname(a_,b_) for a_,b_ in zip(dsel,lsel)]):
                sel=(np.array([mapname(a_,b_) for a_,b_ in zip(dsel,lsel)])==nm)&(wgt>0)
                if not sel.any() or not (D/nm).exists(): continue
                dd_=dens_maginterp(nm,mags[sel],xs[sel],ys[sel])
                for p in POPS: acc[p][sel]+=wgt[sel]*dd_[p]
                wsum[sel]+=wgt[sel]; n_corners[sel]+=1
    ok=wsum>0
    for p in POPS: acc[p]=np.where(ok,acc[p]/np.where(ok,wsum,1),0.0)
    tD=sum(acc[p] for p in POPS)
    res["D"]=np.where(ok&(tD>0),acc["NEO"]/np.where(tD>0,tD,1),res["C"])
    totD=np.where(ok&(tD>0),tD,tC)
    return idx,res,dict(A=totA,B=tB,C=tC,D=totD),dict(outside_hull=outside_hull,
        out_lat=out_lat,out_dlon=out_dlon,n_corners=n_corners,wsum=wsum)

t0=time.time()
cells=sorted(cal.prob_map_file.unique())
parts=Parallel(n_jobs=12,verbose=5)(delayed(run)(c) for c in cells)
S={k:np.full(len(cal),np.nan) for k in "ABCD"}; T={k:np.full(len(cal),np.nan) for k in "ABCD"}
H={k:np.zeros(len(cal),bool) for k in ["outside_hull","out_lat","out_dlon"]}
NC=np.zeros(len(cal),int)
for idx,res,tot,h in parts:
    for k in "ABCD": S[k][idx]=res[k]; T[k][idx]=tot[k]
    for k in H: H[k][idx]=h[k]
    NC[idx]=h["n_corners"]
rt=time.time()-t0
inb=(cal.vlam.abs()<=5)&(cal.vbeta.abs()<=5)&(cal.mean_mag>=14)&(cal.mean_mag<25)
valid={k:inb.to_numpy()&np.isfinite(S[k])&np.isfinite(T[k])&(T[k]>0) for k in "ABCD"}
m=np.logical_and.reduce([valid[k] for k in "ABCD"])
print(f"\ntechnical-valid per variant: {{ {', '.join(f'{k}:{int(valid[k].sum()):,}' for k in 'ABCD')} }}")
print(f"intersection {int(m.sum()):,} of {len(cal):,}")
yv=y[m]; vv=vmax[m]
print(f"  coverage overall {100*m.mean():.4f}%   NEO {100*m[y==1].mean():.4f}% "
      f"({int(m[y==1].sum())}/{int((y==1).sum())})   non-NEO {100*m[y==0].mean():.4f}%")
rows=[]
for k in "ABCD":
    s=S[k][m]; assert np.isfinite(s).all()
    p_,r_,t_=precision_recall_curve(yv,s)
    f1=np.divide(2*p_*r_,p_+r_,out=np.zeros_like(p_),where=(p_+r_)>0); i=int(np.argmax(f1[:-1]))
    fast=(vv>2)&(yv==1)
    rows.append(dict(variant=k,ROC_AUC=roc_auc_score(yv,s),
        stdpAUC_fpr01=roc_auc_score(yv,s,max_fpr=0.01),F1=f1[i],
        Brier=brier_score_loss(yv,np.clip(s,0,1)),
        logloss=log_loss(yv,np.clip(s,1e-15,1-1e-15)),
        thr=float(t_[i]) if i<len(t_) else 1.0,completeness=r_[i]*100,contamination=(1-p_[i])*100,
        fastNEO_median=float(np.median(s[fast])),fastNEO_ge05=100*float(np.mean(s[fast]>=0.5)),
        fastNEO_zero=100*float(np.mean(s[fast]==0))))
R=pd.DataFrame(rows); pd.set_option("display.width",250)
print(f"\n=== A-D on identical rows (raw posterior, NO calibration) ===")
print(R.to_string(index=False,float_format=lambda z:f"{z:,.6g}"))
print(f"\n  runtime {rt/60:.1f} min for 4 variants x {len(cells)} cells")
print("\n=== matched-completeness contamination (frozen comparison rule) ===")
mc=[]
for k in "ABCD":
    s=S[k][m]; p_,r_,_=precision_recall_curve(yv,s); e={"variant":k}
    for tgt in [0.60,0.70,0.80,0.85,0.90]:
        j=np.argmin(np.abs(r_-tgt)); e[f"@{int(tgt*100)}%"]=(1-p_[j])*100
    mc.append(e)
print(pd.DataFrame(mc).to_string(index=False,float_format=lambda z:f"{z:,.3f}"))
print("\n=== reliability (variant A vs D) ===")
for k in ["A","D"]:
    s=S[k][m]; b=pd.cut(s,[0,.1,.3,.5,.7,.9,1.001],include_lowest=True)
    t=pd.DataFrame({"p":s,"y":yv}).groupby(b,observed=True).agg(n=("y","size"),pred=("p","mean"),obs=("y","mean"))
    print(f"\n  {k}:"); print(t.to_string(float_format=lambda z:f"{z:,.4f}"))
print("\n=== SKY INTERPOLATION HULL (variant D) ===")
hull=pd.DataFrame({"truth":cal.population.to_numpy()[m],"outside":H["outside_hull"][m],
                   "out_lat":H["out_lat"][m],"out_dlon":H["out_dlon"][m],"n_corners":NC[m]})
print("  rows OUTSIDE the hull (clipped to a boundary map, NOT interpolated):")
print(hull.groupby("truth").agg(n=("outside","size"),outside=("outside","sum"),
      pct=("outside",lambda t:100*t.mean()),out_lat=("out_lat","sum"),
      out_dlon=("out_dlon","sum")).to_string(float_format=lambda z:f"{z:,.3f}"))
print("\n  corner accounting (how many of the 4 sky corners contributed):")
cc=hull.n_corners.value_counts().sort_index()
tot_=len(hull)
for k_,vv_ in cc.items():
    lab={4:"four valid corners (full bilinear)",3:"3 corners (renormalised)",
         2:"2 corners (renormalised)",1:"1 corner (complete boundary fallback)",
         0:"0 corners (fell back to C)"}.get(k_,f"{k_} corners")
    print(f"    {lab:44s} {vv_:>9,}  ({100*vv_/tot_:5.2f}%)")
print("\n=== PAIRED DIFFERENCES (truth-stratified bootstrap, B=300) ===")
rng=np.random.default_rng(0); ii1=np.where(yv==1)[0]; ii0=np.where(yv==0)[0]
def paired(k1,k2):
    s1=S[k1][m]; s2=S[k2][m]; d={"ROC":[],"pAUC":[],"F1":[],"Brier":[]}
    for _ in range(300):
        b_=np.concatenate([rng.choice(ii1,len(ii1),True),rng.choice(ii0,len(ii0),True)])
        yy=yv[b_]
        if yy.sum()<5: continue
        d["ROC"].append(roc_auc_score(yy,s2[b_])-roc_auc_score(yy,s1[b_]))
        d["pAUC"].append(roc_auc_score(yy,s2[b_],max_fpr=0.01)-roc_auc_score(yy,s1[b_],max_fpr=0.01))
        def f1o(t):
            p_,r_,_=precision_recall_curve(yy,t)
            f=np.divide(2*p_*r_,p_+r_,out=np.zeros_like(p_),where=(p_+r_)>0); return f[:-1].max()
        d["F1"].append(f1o(s2[b_])-f1o(s1[b_]))
        d["Brier"].append(brier_score_loss(yy,np.clip(s2[b_],0,1))-brier_score_loss(yy,np.clip(s1[b_],0,1)))
    return {k_:(float(np.mean(vv_)),float(np.percentile(vv_,2.5)),float(np.percentile(vv_,97.5))) for k_,vv_ in d.items()}
for pair in [("A","C"),("C","D"),("A","D")]:
    pr=paired(*pair)
    print(f"  {pair[1]} minus {pair[0]}:")
    for k_,(mu,lo,hi) in pr.items():
        print(f"    {k_:6s} {mu:+.3e}  95% CI [{lo:+.3e}, {hi:+.3e}]  "
              f"{'excludes 0' if lo>0 or hi<0 else 'includes 0'}")
np.savez_compressed(OUT/"E1_PER_ROW_SCORES.npz",**{k:S[k] for k in "ABCD"},mask=m,y=y,vmax=vmax,
                    outside_hull=H["outside_hull"],n_corners=NC)
R.to_csv(OUT/"E1_ABLATION.csv",index=False)
json.dump(dict(rows=R.to_dict("records"),matched_completeness=mc,runtime_min=rt/60,
    n_eval=int(m.sum()),n_neo=int(yv.sum()),
    coverage=dict(overall=100*float(m.mean()),NEO=100*float(m[y==1].mean()))),
    open(OUT/"E1_ABLATION.json","w"),indent=2)
print(f"\nwrote {OUT}/E1_ABLATION.csv/.json")
