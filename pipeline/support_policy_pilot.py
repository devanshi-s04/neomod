#!/usr/bin/env python3
"""Support-policy pilot (no rebuild) -- standardized global support score Z.

    z_c(X) = d_k,c(X) / q_c(map, mag_bin)      d_k from DIRECT GEN cKDTree queries
    Z(X)   = min_c z_c(X)                      q_c from GEN leave-one-out kth-NN distances
    Z <= s*  -> full UNMASKED posterior over every population
    Z >  s*  -> abstain (NaN) for the whole observation

No density is ever zeroed per population. Coverage targets are FIXED in advance; no classifier
metric is optimised. Frozen 15,422 CAL rows; abstention rows applied identically to the unmasked
0.01 and 0.005 scores.
"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score, precision_recall_curve, brier_score_loss
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
os.environ.setdefault("NEOMOD3_CACHE_DIR", str(W/"outputs/neomod3_projection_cache/by_pixel"))
import velocity_density_pipeline_neomod_clone_only as v
from astropy.time import Time
from astropy.utils import iers; iers.conf.auto_max_age = None
from astropy.coordinates import get_sun, GeocentricTrueEcliptic
OUT = W/"outputs/e0_results"
A_DIR, B_DIR = W/"prob_maps_e0_thr10", W/"prob_maps_e0_thr10_step005"
CENTERS = [(0,0),(20,-12),(10,2),(-10,-2),(0,8),(50,-1),(50,1),(-50,1),
           (90,0),(140,0),(-140,0),(0,25),(0,-25),(0,50),(0,-50),(-50,-50)]
POPS = ["NEO","MBA","TNO","Trojans"]; THRS=[2,3,5,10]; K = v.DEFAULT_K_MAP
BINS = [("14_16",14.,16.),("16_18",16.,18.),("18_20",18.,20.),("mag20",20.,21.),
        ("mag21",21.,22.),("mag22",22.,23.),("mag23",23.,24.),("mag24+",24.,25.)]
COVS = [99.9, 99.5, 99.0, 98.0]     # FIXED IN ADVANCE
v.NONNEO_SPLIT_FRACTIONS = json.loads((W/"outputs/splits/split_provenance.json").read_text())

cal = pd.read_parquet(W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet")
names = sorted({f"prob_maps_grid_dlon{a:+04d}_lat{b:+03d}.npz" for a,b in CENTERS})
cal = cal[cal.prob_map_file.isin(names)].reset_index(drop=True)
y = (cal.population=="NEO").astype(int).to_numpy()
vmax = np.maximum(cal.vlam.abs(), cal.vbeta.abs()).to_numpy()

def score_dir(dirp, smin):
    s=np.full(len(cal),np.nan); c=np.zeros(len(cal),bool)
    for cen,g in cal.groupby("prob_map_file"):
        pm=v.ProbMapSet.from_npz(str(dirp/cen),support_mask_min=smin,mask_radius_deg_per_day=np.inf)
        r=pm.score_visible(g.vlam.to_numpy(float),g.vbeta.to_numpy(float),g.mean_mag.to_numpy(float))
        i=g.index.to_numpy(); s[i]=np.asarray(r["NEO"],float)
        c[i]=np.sum([np.nan_to_num(np.asarray(r[p],float)) for p in POPS],axis=0)>0
        del pm
    return s,c
cov={}
for t in THRS: _,cov[t]=score_dir(W/f"prob_maps_e0_thr{t}",1)
common=np.logical_and.reduce([cov[t] for t in THRS])
print(f"frozen rows {int(common.sum()):,}  NEO {int(y[common].sum())}")
sA,_=score_dir(A_DIR,0); sB,_=score_dir(B_DIR,0)      # UNMASKED
print("unmasked scores computed for 0.01 and 0.005")

t0=Time("2027-08-25T00:00:00",scale="utc")
sun=get_sun(t0).transform_to(GeocentricTrueEcliptic(obstime=t0))
man=pd.read_parquet(W/"outputs/splits/nonneo_split_manifest.parquet")
gen=set(man.ObjID[man.split=="GEN"])
cache=pd.read_parquet(W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet")
cache=cache[cache.ObjID.isin(gen)|(cache.population=="NEO")].reset_index(drop=True)
_,scorer=v.load_s3m_population("neo",verbose=False)

def per_center(cen):
    g=cal[cal.prob_map_file==cen]
    dl=int(cen.split("dlon")[1][:4]); la=float(cen.split("lat")[1][:3])
    clon=(sun.lon.deg+180.0+dl)%360.0
    vis={}
    for p in POPS:
        if p=="NEO":
            df,_=v._load_neomod3_cache(center_lon_deg=clon,center_lat_deg=la,
                                       obstime_str="2027-08-25T00:00:00",max_sep_deg=30.0)
        else:
            df=cache[cache.population==p].reset_index(drop=True)
        vis[p]=v.build_visible_subset_dataframe(df,obstime_str="2027-08-25T00:00:00",
            scorer=scorer,max_sep_deg=30.0,chunk=200_000,show_progress=False,
            center_mode="custom_ecliptic",center_lon_deg=clon,center_lat_deg=la)
    rec=[]
    for lab,lo,hi in BINS:
        m=(g.mean_mag>=lo)&(g.mean_mag<hi)
        if not m.any(): continue
        gg=g[m]; Q=np.column_stack([gg.vlam.to_numpy(float),gg.vbeta.to_numpy(float)])
        zc={}
        for p in POPS:
            vs=vis[p]; vs=vs[(vs.mag_app>=lo)&(vs.mag_app<hi)]
            if len(vs)<K+1:
                zc[p]=np.full(len(gg),np.inf); continue
            P=np.column_stack([vs.vlam.to_numpy(float),vs.vbeta.to_numpy(float)])
            tr=cKDTree(P)
            dk=tr.query(Q,k=K)[0][:,-1]                     # direct kth-NN distance
            loo=tr.query(P,k=K+1)[0][:,-1]                  # leave-one-out reference spacing
            q=float(np.median(loo)) or 1e-12
            zc[p]=dk/q
        Z=np.min(np.column_stack([zc[p] for p in POPS]),axis=1)
        who=np.array(POPS)[np.argmin(np.column_stack([zc[p] for p in POPS]),axis=1)]
        rec.append(pd.DataFrame(dict(idx=gg.index.to_numpy(),Z=Z,who=who,magbin=lab,
                                     **{f"z_{p}":zc[p] for p in POPS})))
    return pd.concat(rec,ignore_index=True) if rec else None

parts=Parallel(n_jobs=8,verbose=5)(delayed(per_center)(c) for c in names)
Zdf=pd.concat([p for p in parts if p is not None],ignore_index=True).set_index("idx").sort_index()
Z=np.full(len(cal),np.inf); Z[Zdf.index.to_numpy()]=Zdf.Z.to_numpy()

m=common; yv=y[m]; Zc=Z[m]; a=np.nan_to_num(sA[m]); b=np.nan_to_num(sB[m])
print(f"\n=== Z distribution on the frozen rows ===")
S=pd.DataFrame({"truth":cal.population.to_numpy()[m],"Z":Zc,
                "vband":pd.cut(vmax[m],[0,0.25,0.5,1.0,2.0,5.1])})
print(S.groupby("truth").Z.describe(percentiles=[.5,.9,.99])[["count","50%","90%","99%","max"]]
      .to_string(float_format=lambda z:f"{z:,.4g}"))
print("\n  by velocity band:")
print(S.groupby("vband",observed=True).Z.describe(percentiles=[.5,.9,.99])[["count","50%","90%","99%","max"]]
      .to_string(float_format=lambda z:f"{z:,.4g}"))

def mets(s,keep):
    yy,ss=yv[keep],s[keep]
    if yy.sum()<5: return dict(ROC=np.nan,pAUC=np.nan,F1=np.nan,Brier=np.nan)
    p_,r_,_=precision_recall_curve(yy,ss)
    f1=np.divide(2*p_*r_,p_+r_,out=np.zeros_like(p_),where=(p_+r_)>0)
    return dict(ROC=roc_auc_score(yy,ss),pAUC=roc_auc_score(yy,ss,max_fpr=0.01),
                F1=f1[:-1].max(),Brier=brier_score_loss(yy,np.clip(ss,0,1)))
rng=np.random.default_rng(0); fast=(vmax[m]>2)&(yv==1)
rows=[]
for cv in COVS:
    s_star=np.percentile(Zc,cv); keep=Zc<=s_star
    kr=np.zeros_like(keep); kr[rng.choice(len(keep),int(keep.sum()),replace=False)]=True
    MA,MB,MR=mets(a,keep),mets(b,keep),mets(a,kr)
    thrA=0.5
    rows.append(dict(target_cov=cv,s_star=s_star,actual_cov=100*keep.mean(),
        n_kept=int(keep.sum()),n_abst=int((~keep).sum()),
        abst_NEO=int(((~keep)&(yv==1)).sum()),abst_nonNEO=int(((~keep)&(yv==0)).sum()),
        fastNEO_retained=100*float(keep[fast].mean()) if fast.sum() else np.nan,
        ROC=MA["ROC"],pAUC=MA["pAUC"],F1=MA["F1"],Brier=MA["Brier"],
        ROC_rand=MR["ROC"],pAUC_rand=MR["pAUC"],F1_rand=MR["F1"],
        flips_001_vs_0005=int(((a[keep]>=0.5)!=(b[keep]>=0.5)).sum()),
        max_absd_001_0005=float(np.abs(b[keep]-a[keep]).max())))
R=pd.DataFrame(rows)
pd.set_option("display.width",250)
print("\n=== fixed coverage points (NO metric optimised) ===")
print(R[["target_cov","s_star","actual_cov","n_kept","n_abst","abst_NEO","abst_nonNEO",
         "fastNEO_retained","ROC","pAUC","F1","Brier"]].to_string(index=False,float_format=lambda z:f"{z:,.5g}"))
print("\n  random-abstention control (same coverage, informative only):")
print(R[["target_cov","ROC","ROC_rand","pAUC","pAUC_rand","F1","F1_rand"]].to_string(index=False,float_format=lambda z:f"{z:,.5g}"))
print("\n  score stability on RETAINED rows (0.01 vs 0.005, unmasked):")
print(R[["target_cov","flips_001_vs_0005","max_absd_001_0005"]].to_string(index=False))
print("\n=== error rate vs Z (calibration direction) ===")
Zf=np.where(np.isfinite(Zc),Zc,np.nanmax(Zc[np.isfinite(Zc)])*10)
E=pd.DataFrame({"Zbin":pd.qcut(Zf,10,duplicates="drop"),"Z":Zf,"y":yv,"p":a})
G=E.groupby("Zbin",observed=True)
print(pd.DataFrame({"n":G.size(),"meanZ":G.Z.mean(),"NEOfrac":G.y.mean(),
                    "meanP":G.p.mean(),
                    "absErr":G.apply(lambda t:np.abs(t.p-t.y).mean(),include_groups=False)}
                   ).to_string(float_format=lambda z:f"{z:,.4g}"))
print("\n=== abstention composition at 99% coverage ===")
s99=np.percentile(Zc,99.0); k99=Zc<=s99
A2=pd.DataFrame({"truth":cal.population.to_numpy()[m],"mag":cal.mean_mag.to_numpy()[m],
                 "field":cal.prob_map_file.to_numpy()[m],"vband":S.vband.values,"abst":~k99})
for lab,key in [("truth","truth"),("magnitude",pd.cut(A2.mag,[14,20,22,23,24,25])),
                ("velocity band","vband")]:
    print(f"\n  by {lab}:")
    print(A2.groupby(key,observed=True).abst.agg(n="size",abstained="sum",
          pct=lambda t:100*t.mean()).to_string(float_format=lambda z:f"{z:,.3f}"))
print(f"\n  fields with the most abstentions:")
print(A2[A2.abst].field.value_counts().head(5).to_string())
R.to_csv(OUT/"SUPPORT_POLICY_PILOT.csv",index=False)
json.dump(dict(coverage_points=R.to_dict("records"),
               Z_by_truth=S.groupby("truth").Z.median().to_dict()),
          open(OUT/"SUPPORT_POLICY_PILOT.json","w"),indent=2,default=str)
print(f"\nwrote SUPPORT_POLICY_PILOT.csv/json")
