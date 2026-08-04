#!/usr/bin/env python3
"""CORRECTED TNO exact-query test: pure raster-resolution, both confounds removed.

Confound 1 (FIXED): the rasterised TNO term is bilinear over the support-MASKED array, but the
exact-query term was inserted unmasked -- that changed support semantics as well as resolution.
The exact term now carries the identical support-mask decision, applied with the same
interpolation the production scorer uses.
Confound 2 (FIXED): the frozen E0 common row set (15,422 rows / 181 NEOs) is used. All 11 original
abstentions -- every one of them a NEO -- stay abstained.
"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score, precision_recall_curve
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
os.environ.setdefault("NEOMOD3_CACHE_DIR", str(W/"outputs/neomod3_projection_cache/by_pixel"))
import velocity_density_pipeline_neomod_clone_only as v
from astropy.time import Time
from astropy.utils import iers; iers.conf.auto_max_age = None
from astropy.coordinates import get_sun, GeocentricTrueEcliptic
EPOCH="2027-08-25T00:00:00"; MAPS=W/"prob_maps_e0_thr10"; OUT=W/"outputs/e0_results"
BINS=[("14_16",14.,16.),("16_18",16.,18.),("18_20",18.,20.),("mag20",20.,21.),
      ("mag21",21.,22.),("mag22",22.,23.),("mag23",23.,24.),("mag24+",24.,25.)]
CENTERS=[(0,0),(20,-12),(10,2),(-10,-2),(0,8),(50,-1),(50,1),(-50,1),
         (90,0),(140,0),(-140,0),(0,25),(0,-25),(0,50),(0,-50),(-50,-50)]
POPS=["NEO","MBA","TNO","Trojans"]; THRS=[2,3,5,10]
v.NONNEO_SPLIT_FRACTIONS=json.loads((W/"outputs/splits/split_provenance.json").read_text())

def bil(a,x,y,xs,ys):
    ix=np.clip(np.searchsorted(x,xs)-1,0,len(x)-2); iy=np.clip(np.searchsorted(y,ys)-1,0,len(y)-2)
    tx=np.clip((xs-x[ix])/(x[ix+1]-x[ix]),0,1); ty=np.clip((ys-y[iy])/(y[iy+1]-y[iy]),0,1)
    return (a[iy,ix]*(1-tx)*(1-ty)+a[iy,ix+1]*tx*(1-ty)+a[iy+1,ix]*(1-tx)*ty+a[iy+1,ix+1]*tx*ty)

cal=pd.read_parquet(W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet")
names={f"prob_maps_grid_dlon{a:+04d}_lat{b:+03d}.npz" for a,b in CENTERS}
cal=cal[cal.prob_map_file.isin(names)].reset_index(drop=True)

# ---- reproduce the FROZEN common row set exactly as E0_SIMPLE computed it ----
cov={t:np.zeros(len(cal),bool) for t in THRS}
for cen,g in cal.groupby("prob_map_file"):
    idx=g.index.to_numpy()
    for t in THRS:
        pm=v.ProbMapSet.from_npz(str(W/f"prob_maps_e0_thr{t}"/cen),support_mask_min=1,
                                 mask_radius_deg_per_day=np.inf)
        r=pm.score_visible(g.vlam.to_numpy(float),g.vbeta.to_numpy(float),g.mean_mag.to_numpy(float))
        cov[t][idx]=np.sum([np.nan_to_num(np.asarray(r[p],float)) for p in POPS],axis=0)>0
        del pm
common=np.logical_and.reduce([cov[t] for t in THRS])
abst=~common
print(f"frozen common rows: {int(common.sum()):,}  NEO {int((cal.population=='NEO')[common].sum())}")
print(f"preserved abstentions: {int(abst.sum())}  of which NEO {int((cal.population=='NEO')[abst].sum())}")

t0=Time(EPOCH,scale="utc")
sun=get_sun(t0).transform_to(GeocentricTrueEcliptic(obstime=t0))
man=pd.read_parquet(W/"outputs/splits/nonneo_split_manifest.parquet")
gen=set(man.ObjID[man.split=="GEN"])
cache=pd.read_parquet(W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet")
cache=cache[cache.ObjID.isin(gen)|(cache.population=="NEO")].reset_index(drop=True)
_,scorer=v.load_s3m_population("neo",verbose=False)

def one(cen):
    g=cal[(cal.prob_map_file==cen)]
    z=np.load(MAPS/cen,allow_pickle=True); xg,yg=z["x_grid"],z["y_grid"]
    dl=int(cen.split("dlon")[1][:4]); la=int(cen.split("lat")[1][:3])
    pm=v.ProbMapSet.from_npz(str(MAPS/cen),support_mask_min=1,mask_radius_deg_per_day=np.inf)
    prod=np.asarray(pm.score_visible(g.vlam.to_numpy(float),g.vbeta.to_numpy(float),
                                     g.mean_mag.to_numpy(float))["NEO"],float); del pm
    vis=v.build_visible_subset_dataframe(cache[cache.population=="TNO"].reset_index(drop=True),
        obstime_str=EPOCH,scorer=scorer,max_sep_deg=30.0,chunk=200_000,show_progress=False,
        center_mode="custom_ecliptic",
        center_lon_deg=(sun.lon.deg+180.0+dl)%360.0,center_lat_deg=float(la))
    rec=[]
    for lab,lo,hi in BINS:
        m=(g.mean_mag>=lo)&(g.mean_mag<hi)
        if not m.any(): continue
        gg=g[m]; vs=vis[(vis.mag_app>=lo)&(vis.mag_app<hi)]
        xs,ys=gg.vlam.to_numpy(float),gg.vbeta.to_numpy(float)
        dens={}; 
        for p in POPS:
            a=np.nan_to_num(np.asarray(z[f"density_raw__{p}__{lab}"],float))
            s=np.nan_to_num(np.asarray(z[f"support_count__{p}__{lab}"],float))
            if p!="NEO": a=np.where(s>=1,a,0.0)
            dens[p]=bil(a,xg,yg,xs,ys)
        sup_t=np.nan_to_num(np.asarray(z[f"support_count__TNO__{lab}"],float))
        mask_ind=bil((sup_t>=1).astype(float),xg,yg,xs,ys)   # SAME mask decision + interpolation
        sup_at=bil(sup_t,xg,yg,xs,ys)
        tot=sum(dens[p] for p in POPS)
        p_man=np.where(tot>0,dens["NEO"]/np.where(tot>0,tot,1),0.0)
        if len(vs)>=11:
            tree=cKDTree(np.column_stack([vs.vlam.to_numpy(float),vs.vbeta.to_numpy(float)]))
            fs=v._nonneo_split_fraction("TNO",lo,hi)
            raw=np.array([v.estimate_density_full_posterior_2d(tree,a1,b1,k=v.DEFAULT_K_MAP,
                          n_d0_grid=v.DEFAULT_N_D0_GRID_MAP) for a1,b1 in zip(xs,ys)])/fs
            nnd=tree.query(np.column_stack([xs,ys]),k=1)[0]
        else:
            raw=dens["TNO"].copy(); nnd=np.full(len(xs),np.nan)
        ex_m=raw*mask_ind                                     # masked (primary)
        t_m=tot-dens["TNO"]+ex_m;  p_m=np.where(t_m>0,dens["NEO"]/np.where(t_m>0,t_m,1),0.0)
        t_u=tot-dens["TNO"]+raw;   p_u=np.where(t_u>0,dens["NEO"]/np.where(t_u>0,t_u,1),0.0)
        rec.append(pd.DataFrame(dict(idx=gg.index.to_numpy(),magbin=lab,P_prod=prod[m.to_numpy()],
            P_manual=p_man,P_exact_masked=p_m,P_exact_unmasked=p_u,
            rho_TNO_grid=dens["TNO"],rho_TNO_exact=ex_m,rho_TNO_exact_raw=raw,
            rho_NEO=dens["NEO"],rho_MBA=dens["MBA"],rho_Troj=dens["Trojans"],
            tno_support=sup_at,nn_tno_dist=nnd,vlam=xs,vbeta=ys)))
    return pd.concat(rec,ignore_index=True) if rec else None

parts=Parallel(n_jobs=8,verbose=5)(delayed(one)(c) for c in sorted(names))
D=pd.concat([p for p in parts if p is not None],ignore_index=True).set_index("idx").sort_index()
D["truth"]=cal.population.reindex(D.index).values
D["was_abstention"]=abst[D.index.to_numpy()]
D["map"]=cal.prob_map_file.reindex(D.index).values

err=np.abs(D.P_prod-D.P_manual)
print(f"\n=== baseline reproduction (manual bilinear+mask vs production scorer) ===")
print(f"  max |err| = {err.max():.3e}   median = {np.median(err):.3e}   rows = {len(D):,}")

C=D[~D.was_abstention]
y=(C.truth=="NEO").astype(int).to_numpy()
a=C.P_prod.to_numpy(); bm=C.P_exact_masked.to_numpy(); bu=C.P_exact_unmasked.to_numpy()
def mets(s):
    p_,r_,t_=precision_recall_curve(y,s)
    f1=np.divide(2*p_*r_,p_+r_,out=np.zeros_like(p_),where=(p_+r_)>0); i=int(np.argmax(f1[:-1]))
    return roc_auc_score(y,s),roc_auc_score(y,s,max_fpr=0.01),f1[i],(float(t_[i]) if i<len(t_) else 1.0)
A=mets(a); B=mets(bm); U=mets(bu); thr=A[3]
d=np.abs(bm-a)
print(f"\n=== PRIMARY: masked exact-query, frozen {len(C):,} rows / {int(y.sum())} NEOs ===")
print(f"  changed {int((bm!=a).sum()):,} ({100*(bm!=a).mean():.2f}%)")
print(f"  |dP| median {np.median(d):.3e}  p95 {np.percentile(d,95):.3e}  "
      f"p99 {np.percentile(d,99):.3e}  max {d.max():.3e}")
print(f"  ROC AUC {A[0]:.6f} -> {B[0]:.6f}  (d={B[0]-A[0]:+.2e})")
print(f"  pAUC    {A[1]:.6f} -> {B[1]:.6f}  (d={B[1]-A[1]:+.2e})")
print(f"  best F1 {A[2]:.6f} -> {B[2]:.6f}  (d={B[2]-A[2]:+.2e})")
print(f"  flips at production best-F1 threshold {thr:.4g}: {int(((a>=thr)!=(bm>=thr)).sum())}")
S=pd.DataFrame({"truth":C.truth.values,"magbin":C.magbin.values,"absd":d,
                "flip":(a>=thr)!=(bm>=thr)})
print("\n  by truth population:")
print(S.groupby("truth").agg(n=("absd","size"),median=("absd","median"),max=("absd","max"),
      flips=("flip","sum")).to_string(float_format=lambda z:f"{z:,.3e}"))
print("\n  by magnitude bin:")
print(S.groupby("magbin").agg(n=("absd","size"),median=("absd","median"),max=("absd","max"),
      flips=("flip","sum")).to_string(float_format=lambda z:f"{z:,.3e}"))
du=np.abs(bu-a)
print(f"\n=== DIAGNOSTIC ONLY (unmasked exact query -- CHANGES SUPPORT SEMANTICS, not comparable) ===")
print(f"  |dP| median {np.median(du):.3e}  max {du.max():.3e}   "
      f"flips {int(((a>=thr)!=(bu>=thr)).sum())}   F1 {U[2]:.6f}")

# ---- audit of the previous run's flips (unmasked, all rows incl. the 11 abstentions) ----
prev_flip=(D.P_prod>=thr)!=(D.P_exact_unmasked>=thr)
AU=D[prev_flip].copy()
AU["P_before"]=AU.P_prod; AU["P_after_unmasked"]=AU.P_exact_unmasked
cols=["truth","map","magbin","vlam","vbeta","tno_support","nn_tno_dist","rho_TNO_grid",
      "rho_TNO_exact_raw","rho_NEO","rho_MBA","rho_Troj","P_before","P_after_unmasked",
      "P_exact_masked","was_abstention"]
AU[cols].to_csv(OUT/"TNO_FLIP_AUDIT.csv",index=True)
print(f"\n=== previous-run flip audit: {len(AU)} rows -> TNO_FLIP_AUDIT.csv ===")
print(AU[["truth","magbin","tno_support","nn_tno_dist","rho_TNO_grid","rho_TNO_exact_raw",
          "P_before","P_after_unmasked","was_abstention"]].head(12).to_string(
          float_format=lambda z:f"{z:,.4g}"))
print(f"\n  of these {len(AU)}: NEO {int((AU.truth=='NEO').sum())}, "
      f"original abstentions {int(AU.was_abstention.sum())}, "
      f"TNO support<1 at query {int((AU.tno_support<1).sum())}")
json.dump(dict(baseline_max_err=float(err.max()),rows=int(len(C)),neos=int(y.sum()),
    changed=int((bm!=a).sum()),median=float(np.median(d)),p95=float(np.percentile(d,95)),
    p99=float(np.percentile(d,99)),max=float(d.max()),
    roc=[A[0],B[0]],pauc=[A[1],B[1]],f1=[A[2],B[2]],
    flips=int(((a>=thr)!=(bm>=thr)).sum()),thr=float(thr),
    unmasked_diagnostic=dict(median=float(np.median(du)),max=float(du.max()),
        flips=int(((a>=thr)!=(bu>=thr)).sum())),
    prev_flip_audit_rows=int(len(AU))),
    open(OUT/"TNO_EXACT_QUERY_CORRECTED.json","w"),indent=2)
print(f"\nwrote TNO_EXACT_QUERY_CORRECTED.json")
