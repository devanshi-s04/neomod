#!/usr/bin/env python3
"""Acceptance checks for the full 667-map NEOMOD3 grid + full-grid scoring vs production VDP/digest2.

A) file count / integrity
B) the two T1-T4-validated centers must still be bit-identical to the validated rebuild
C) NEO velocity support must reach beyond |v|=2 across the whole grid (the point of the +-5 domain)
D) score ALL tracklets at ALL centers and compare NEOMOD3 vs production VDP vs digest2
"""
from __future__ import annotations
import os, sys, glob
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import precision_recall_curve, roc_auc_score
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
FULL = W/"prob_maps_grid_neomod3_full"; REF = W/"prob_maps_grid_neomod3_vlim5_hpx"
SORCHA = W/"outputs/s3m_linking/case1/sorcha_comparison_case1_nbody_Vband.parquet"
OUT = W/"outputs/neomod3_fullgrid"; OUT.mkdir(parents=True, exist_ok=True)
NIGHT = 61642
import velocity_density_pipeline_neomod_clone_only as vnm

print("="*72); print("A) file count")
files = sorted(glob.glob(str(FULL/"*.npz")))
print(f"   maps present: {len(files)} / 667   total {sum(os.path.getsize(f) for f in files)/2**30:.1f} GiB")
bad = [f for f in files if os.path.getsize(f) < 10*2**20]
print(f"   suspiciously small files: {len(bad)}")

# A2 -- INTEGRITY SCAN. A size check is NOT enough: a task preempted mid-write leaves a truncated
# .npz that is still ~100 MB, still counted by `ls | wc -l`, still reported COMPLETED by Slurm, and
# -- worst -- is SKIPPED by a plain resubmit because skip-existing keys on path existence. The only
# thing that catches it is actually opening the file and forcing a payload array read. (This is how
# prob_maps_grid_dlon+050_lat-25.npz was found; see docs §9.2.)
def _integrity(f):
    try:
        z = np.load(f, allow_pickle=True)
        nk = len(z.files)
        _ = np.asarray(z["density_raw__NEO__mag23"]).sum()
        return (f, None if nk == 151 else f"expected 151 keys, got {nk}")
    except Exception as e:
        return (f, f"{type(e).__name__}: {e}")
scan = Parallel(n_jobs=32)(delayed(_integrity)(f) for f in files)
corrupt = [(f, e) for f, e in scan if e]
print(f"   integrity scan (open + read a payload array): {len(files)-len(corrupt)} OK, "
      f"{len(corrupt)} CORRUPT")
for f, e in corrupt:
    print(f"      !! {os.path.basename(f):45s} {os.path.getsize(f)/2**20:8.1f} MiB  {e[:60]}")
if corrupt:
    print("\n   ABORTING before scoring. Rebuild these WITH --overwrite (a plain resubmit skips\n"
          "   them). Map filename -> array index: neomod/pipeline/slurm/grid_map_manifest.csv")
    sys.exit(2)

print("="*72); print("B) validated centers still bit-identical to the T3b rebuild")
for cen in ["prob_maps_grid_dlon+000_lat+00.npz", "prob_maps_grid_dlon+020_lat-12.npz"]:
    a = np.load(REF/cen, allow_pickle=True); b = np.load(FULL/cen, allow_pickle=True)
    worst, wk = 0.0, ""
    for k in sorted(set(a.files) & set(b.files)):
        x, y = a[k], b[k]
        if x.dtype.kind not in "fi" or x.shape != y.shape: continue
        x = np.nan_to_num(np.asarray(x, float), posinf=0.); y = np.nan_to_num(np.asarray(y, float), posinf=0.)
        d = np.abs(x-y).max()/max(np.abs(x).max(), 1e-300)
        if d > worst: worst, wk = d, k
    print(f"   {cen:42s} worst rel diff {worst:.3e}  ({wk or 'all equal'})")

print("="*72); print("C) NEO velocity support extent across the grid (production ceiling was |v|=2.00)")
def support_extent(f):
    z = np.load(f, allow_pickle=True); x = z["x_grid"]; y = z["y_grid"]
    out = {}
    for b in ["mag22", "mag23", "mag24+"]:
        k = f"support_count__NEO__{b}"
        if k not in z.files: continue
        s = np.asarray(z[k], float)
        if not np.isfinite(s).any() or s.max() <= 0: out[b] = 0.0; continue
        iy, ix = np.nonzero(s > 0)
        out[b] = float(max(np.abs(x[ix]).max(), np.abs(y[iy]).max()))
    return os.path.basename(f), out
samp = files[::max(1, len(files)//60)]
ext = Parallel(n_jobs=24)(delayed(support_extent)(f) for f in samp)
for b in ["mag22", "mag23", "mag24+"]:
    v = np.array([e[1].get(b, 0.0) for e in ext])
    print(f"   {b:>7}: NEO support |v|max over {len(v)} sampled maps -> "
          f"min {v.min():.2f}  median {np.median(v):.2f}  max {v.max():.2f}   "
          f"maps reaching >2.0: {(v>2.0).sum()}/{len(v)}")

print("="*72); print("D) full-grid scoring")
trk = pd.read_parquet(SORCHA, columns=["ObjID","night","population","vlam","vbeta","mean_mag_V",
                                       "prob_map_file","P_NEO_vdp_Vband","P_NEO_d2"])
have = {os.path.basename(f) for f in files}
trk = trk[trk.prob_map_file.isin(have)].reset_index(drop=True)
print(f"   scoring {len(trk):,} tracklets over {trk.prob_map_file.nunique()} centers")

def score_center(cen, sub):
    pm = vnm.ProbMapSet.from_npz(str(FULL/cen), support_mask_min=1, mask_radius_deg_per_day=np.inf)
    s = pm.score_visible(sub.vlam.to_numpy(float), sub.vbeta.to_numpy(float),
                         sub.mean_mag_V.to_numpy(float))["NEO"]
    del pm
    return sub.index.to_numpy(), np.asarray(s, float)

groups = [(c, g) for c, g in trk.groupby("prob_map_file")]
res = Parallel(n_jobs=20, verbose=5)(delayed(score_center)(c, g) for c, g in groups)
trk["P_NEO_neomod3"] = np.nan
for idx, s in res: trk.loc[idx, "P_NEO_neomod3"] = s
trk["is_neo"] = (trk.population == "NEO").astype(int)
trk["maxabs_v"] = np.maximum(trk.vlam.abs(), trk.vbeta.abs())
trk.to_parquet(OUT/"fullgrid_scored_tracklets.parquet", index=False)

def metrics(y, s):
    y = np.asarray(y, int); s = np.asarray(s, float); ok = np.isfinite(s)
    if ok.sum() < 10 or y[ok].sum() < 3: return None
    y, s = y[ok], s[ok]
    p, r, th = precision_recall_curve(y, s)
    f1 = np.divide(2*p*r, p+r, out=np.zeros_like(p), where=(p+r) > 0); bi = int(np.argmax(f1[:-1]))
    return dict(n=len(y), n_neo=int(y.sum()), AUC=roc_auc_score(y, s), bestF1=f1[bi],
                completeness=r[bi]*100, contamination=(1-p[bi])*100)

rows = []
for setname, d in [("full_2yr", trk), (f"night_{NIGHT}", trk[trk.night == NIGHT])]:
    for label, col in [("NEOMOD3 full grid (±5)", "P_NEO_neomod3"),
                       ("production VDP (±2)", "P_NEO_vdp_Vband"), ("digest2", "P_NEO_d2")]:
        m = metrics(d.is_neo.to_numpy(), d[col].to_numpy())
        if m: rows.append(dict(tracklets=setname, classifier=label, **m))
summ = pd.DataFrame(rows); summ.to_csv(OUT/"fullgrid_roc_summary.csv", index=False)
pd.set_option("display.width", 220)
print("\n" + summ.to_string(index=False, float_format=lambda v: f"{v:,.4g}"))

fast = trk[trk.maxabs_v > 2.0]
print(f"\n|v|>2 tracklets across the whole grid: {len(fast):,} "
      f"(NEO {int(fast.is_neo.sum()):,}, non-NEO {int((1-fast.is_neo).sum()):,})")
if len(fast):
    print(f"   production VDP == 0 : {100*(fast.P_NEO_vdp_Vband.fillna(0)==0).mean():.1f}%"
          f"     NEOMOD3 == 0 : {100*(fast.P_NEO_neomod3.fillna(0)==0).mean():.1f}%")
    for lab, g in fast.groupby(fast.is_neo.map({1:"NEO",0:"non-NEO"})):
        print(f"   {lab:>8}: production median {g.P_NEO_vdp_Vband.median():.4f} -> "
              f"NEOMOD3 median {g.P_NEO_neomod3.median():.4f}  (n={len(g):,})")
print(f"\nwrote -> {OUT}")
