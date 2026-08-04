#!/usr/bin/env python3
"""Bounded TNO numerical-representation check. NOT model selection.

(1) extend the representative convergence to grid step 0.00125
(2) evaluate TNO kNN density at each CAL object's EXACT (vlam, vbeta) -- same source points, k,
    magnitude bin, normalisation and sky map as the rasterised map
(3) recompute P(NEO) replacing ONLY the rasterised TNO denominator term with the exact query
(4) compare current-grid vs exact-query scores
"""
import copy, json, os, sys
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

EPOCH = "2027-08-25T00:00:00"; MAPS = W/"prob_maps_e0_thr10"
BINS = [("14_16",14.,16.),("16_18",16.,18.),("18_20",18.,20.),("mag20",20.,21.),
        ("mag21",21.,22.),("mag22",22.,23.),("mag23",23.,24.),("mag24+",24.,25.)]
CENTERS = [(0,0),(20,-12),(10,2),(-10,-2),(0,8),(50,-1),(50,1),(-50,1),
           (90,0),(140,0),(-140,0),(0,25),(0,-25),(0,50),(0,-50),(-50,-50)]
POPS = ["NEO","MBA","TNO","Trojans"]
prov = json.loads((W/"outputs/splits/split_provenance.json").read_text())
v.NONNEO_SPLIT_FRACTIONS = prov
OUT = W/"outputs/e0_results"

def bilinear(arr, x, y, xs, ys):
    ix = np.clip(np.searchsorted(x, xs)-1, 0, len(x)-2)
    iy = np.clip(np.searchsorted(y, ys)-1, 0, len(y)-2)
    tx = np.clip((xs-x[ix])/(x[ix+1]-x[ix]), 0, 1); ty = np.clip((ys-y[iy])/(y[iy+1]-y[iy]), 0, 1)
    return (arr[iy,ix]*(1-tx)*(1-ty) + arr[iy,ix+1]*tx*(1-ty) +
            arr[iy+1,ix]*(1-tx)*ty + arr[iy+1,ix+1]*tx*ty)

# ---------- (1) convergence to 0.00125 ----------
d0 = pd.read_csv(OUT/"TNO_RESOLUTION_DIAGNOSTIC.csv")
meta0 = json.loads((OUT/"TNO_RESOLUTION_DIAGNOSTIC.json").read_text())
wst = meta0["worst"]; LIM = float(meta0["window"])
mlo, mhi = dict((b,(lo,hi)) for b,lo,hi in BINS)[wst["bin"]]
t = Time(EPOCH, scale="utc")
sun = get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
clon = (sun.lon.deg + 180.0 + float(wst["dlon"])) % 360.0
man = pd.read_parquet(W/"outputs/splits/nonneo_split_manifest.parquet")
gen_ids = set(man.ObjID[man.split=="GEN"])
cache = pd.read_parquet(W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet")
cache = cache[cache.ObjID.isin(gen_ids) | (cache.population=="NEO")].reset_index(drop=True)
_, scorer = v.load_s3m_population("neo", verbose=False)
cs = {}
for pop, cfg in copy.deepcopy(v.DEFAULT_POPULATION_SETTINGS).items():
    sub = cache[cache.population==pop].reset_index(drop=True)
    cs[pop] = dict(df=sub, scorer=scorer, clone_factor=cfg["clone_factor"], use_conditional_cloner=True,
                   scatter_size=cfg["scatter_size"], scatter_alpha=cfg["scatter_alpha"],
                   _mag_app=sub["mag_app"].to_numpy(float))
x,y,X0,Y0,gp = v.make_default_grid(grid_lim=(-LIM,LIM), grid_step=0.00125)
r = v.build_cloned_maps_for_center_magbin(
    center_lon_deg=clon, center_lat_deg=float(wst["lat"]), center_label="tno0.00125",
    clone_sources=cs, obstime_str=EPOCH, max_sep_deg=30.0, rng=np.random.default_rng(1),
    overlay_max=1000, x_grid=x, y_grid=y, grid_points=gp, X0=X0, k_map=v.DEFAULT_K_MAP,
    n_d0_grid_map=v.DEFAULT_N_D0_GRID_MAP, mag_min=mlo, mag_max=mhi,
    smooth_density_maps=True, smooth_population_names=("NEO",), smooth_support_threshold=10.0,
    smooth_sigma_pixels=v.DEFAULT_SMOOTH_SIGMA_PIXELS,
    smooth_truncate_sigma=v.DEFAULT_SMOOTH_TRUNCATE_SIGMA,
    smooth_support_scale_by_clone_factor=False, smooth_presmoothing_passes=None, n_jobs=16)
dA = float((x[1]-x[0])**2)
u = np.nan_to_num(np.asarray(r["density_maps_unsmoothed"]["TNO"], float))
sup = np.nan_to_num(np.asarray(r["support_count_maps"]["TNO"], float)).sum()
Nfull = sup/float(wst["f"]); R = u.sum()*dA/Nfull
b_exp = 10/9
print("=== (1) convergence, extended ===")
for _, row in d0.iterrows():
    print(f"  step {row.step:.5f}  npix {int(row.npix):>8,}  R = {row.R:8.3f}")
print(f"  step 0.00125  npix {u.size:>8,}  R = {R:8.3f}")
print(f"  expected b = k/(k-1) = {b_exp:.4f}   |R-b| = {abs(R-b_exp):.3f}")
seq = list(d0.R) + [R]
print(f"  sequence: {' -> '.join(f'{q:.3f}' for q in seq)}  "
      f"(monotone decreasing: {all(seq[i]>seq[i+1] for i in range(len(seq)-1))})")

# ---------- (2)(3) exact-query TNO density on the pooled CAL set ----------
cal = pd.read_parquet(W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet")
names = {f"prob_maps_grid_dlon{a:+04d}_lat{b:+03d}.npz" for a,b in CENTERS}
cal = cal[cal.prob_map_file.isin(names)].reset_index(drop=True)
print(f"\n=== (2)(3) exact-query TNO on {len(cal):,} pooled CAL rows ===")

def one_center(cen):
    g = cal[cal.prob_map_file == cen]
    z = np.load(MAPS/cen, allow_pickle=True)
    xg, yg = z["x_grid"], z["y_grid"]
    dl = int(cen.split("dlon")[1][:4]); la = int(cen.split("lat")[1][:3])
    clon_i = (sun.lon.deg + 180.0 + dl) % 360.0
    sub_t = cache[cache.population == "TNO"].reset_index(drop=True)
    vis = v.build_visible_subset_dataframe(sub_t, obstime_str=EPOCH, scorer=scorer,
            max_sep_deg=30.0, chunk=200_000, show_progress=False,
            center_mode="custom_ecliptic", center_lon_deg=clon_i, center_lat_deg=float(la))
    idx, cur, exq = [], [], []
    for lab, lo, hi in BINS:
        m = (g.mean_mag >= lo) & (g.mean_mag < hi)
        if not m.any(): continue
        gg = g[m]; vs = vis[(vis.mag_app >= lo) & (vis.mag_app < hi)]
        dens = {}
        for p in POPS:
            a = np.nan_to_num(np.asarray(z[f"density_raw__{p}__{lab}"], float))
            s = np.nan_to_num(np.asarray(z[f"support_count__{p}__{lab}"], float))
            if p != "NEO":
                a = np.where(s >= 1, a, 0.0)          # support_mask_min=1, NEO exempt
            dens[p] = bilinear(a, xg, yg, gg.vlam.to_numpy(float), gg.vbeta.to_numpy(float))
        tot = sum(dens[p] for p in POPS)
        p_cur = np.where(tot > 0, dens["NEO"]/np.where(tot > 0, tot, 1), 0.0)
        if len(vs) >= 11:
            tree = cKDTree(np.column_stack([vs.vlam.to_numpy(float), vs.vbeta.to_numpy(float)]))
            f_s = v._nonneo_split_fraction("TNO", lo, hi)
            ex = np.array([v.estimate_density_full_posterior_2d(
                tree, xx, yy, k=v.DEFAULT_K_MAP, n_d0_grid=v.DEFAULT_N_D0_GRID_MAP)
                for xx, yy in zip(gg.vlam.to_numpy(float), gg.vbeta.to_numpy(float))])/f_s
        else:
            ex = dens["TNO"]
        tot2 = tot - dens["TNO"] + ex
        p_ex = np.where(tot2 > 0, dens["NEO"]/np.where(tot2 > 0, tot2, 1), 0.0)
        idx.append(gg.index.to_numpy()); cur.append(p_cur); exq.append(p_ex)
    return (np.concatenate(idx) if idx else np.array([], int),
            np.concatenate(cur) if cur else np.array([]),
            np.concatenate(exq) if exq else np.array([]))

res = Parallel(n_jobs=8, verbose=5)(delayed(one_center)(c) for c in sorted(names))
P_cur = np.full(len(cal), np.nan); P_ex = np.full(len(cal), np.nan)
for i, c, e in res:
    if len(i): P_cur[i] = c; P_ex[i] = e
ok = np.isfinite(P_cur) & np.isfinite(P_ex)
yv = (cal.population == "NEO").astype(int).to_numpy()[ok]
a, bq = P_cur[ok], P_ex[ok]; dd = bq - a
print(f"\n=== (4) current-grid vs exact-query ===")
print(f"  rows compared      : {ok.sum():,}   NEO {int(yv.sum())}")
print(f"  changed            : {int((dd != 0).sum()):,} ({100*(dd != 0).mean():.2f}%)")
q = np.abs(dd)
print(f"  |dP| median {np.median(q):.3e}  p95 {np.percentile(q,95):.3e}  "
      f"p99 {np.percentile(q,99):.3e}  max {q.max():.3e}")
def mets(s):
    p_,r_,t_ = precision_recall_curve(yv, s)
    f1 = np.divide(2*p_*r_, p_+r_, out=np.zeros_like(p_), where=(p_+r_)>0); i = int(np.argmax(f1[:-1]))
    return roc_auc_score(yv,s), roc_auc_score(yv,s,max_fpr=0.01), f1[i], (float(t_[i]) if i<len(t_) else 1.0)
A = mets(a); B = mets(bq)
print(f"  ROC AUC  {A[0]:.6f} -> {B[0]:.6f}   (d={B[0]-A[0]:+.2e})")
print(f"  pAUC@.01 {A[1]:.6f} -> {B[1]:.6f}   (d={B[1]-A[1]:+.2e})")
print(f"  best F1  {A[2]:.6f} -> {B[2]:.6f}   (d={B[2]-A[2]:+.2e})")
thr = A[3]
flips = int(((a >= thr) != (bq >= thr)).sum())
print(f"  classification flips at the current best-F1 threshold ({thr:.4g}): {flips}")
sub = pd.DataFrame({"pop": cal.population.to_numpy()[ok], "mag": cal.mean_mag.to_numpy()[ok],
                    "absd": q, "flip": (a >= thr) != (bq >= thr)})
print("\n  by population:")
print(sub.groupby("pop").agg(n=("absd","size"), median_absd=("absd","median"),
      max_absd=("absd","max"), flips=("flip","sum")).to_string(float_format=lambda z: f"{z:,.3e}"))
print("\n  by magnitude bin:")
mb = pd.cut(sub.mag, [14,20,22,23,24,25])
print(sub.groupby(mb, observed=True).agg(n=("absd","size"), median_absd=("absd","median"),
      max_absd=("absd","max"), flips=("flip","sum")).to_string(float_format=lambda z: f"{z:,.3e}"))
json.dump(dict(convergence=[{"step":float(s),"R":float(rr)} for s,rr in
                            zip(list(d0.step)+[0.00125], seq)],
               b_expected=b_exp, R_final=float(R),
               rows=int(ok.sum()), changed=int((dd!=0).sum()),
               median_absd=float(np.median(q)), p95=float(np.percentile(q,95)),
               p99=float(np.percentile(q,99)), max_absd=float(q.max()),
               roc_auc=[A[0],B[0]], pauc=[A[1],B[1]], f1=[A[2],B[2]], flips=flips),
          open(OUT/"TNO_EXACT_QUERY_VALIDATION.json","w"), indent=2)
print(f"\nwrote {OUT}/TNO_EXACT_QUERY_VALIDATION.json")
