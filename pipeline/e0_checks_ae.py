#!/usr/bin/env python3
"""Preregistered Checks A-E on identical common CAL rows (E0_PILOT_PREREGISTRATION + Amendment 1).

Rules are frozen: macro pAUC at FPR<=0.01 over eligible centers, identical common rows for every
candidate, >2pp coverage loss disqualifies, ties -> LARGER threshold. No new metric, no post-hoc
exception. TEST is not touched.
"""
import hashlib, json, sys, os
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
import velocity_density_pipeline_neomod_clone_only as v

THRS = [2, 3, 5, 10]
CENTERS = [(0,0),(20,-12),(10,2),(-10,-2),(0,8),(50,-1),(50,1),(-50,1),
           (90,0),(140,0),(-140,0),(0,25),(0,-25),(0,50),(0,-50),(-50,-50)]
CAL = W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet"
OUT = W/"outputs/e0_results"; OUT.mkdir(parents=True, exist_ok=True)
POPS = ["NEO", "MBA", "TNO", "Trojans"]
BINS = ["14_16","16_18","18_20","mag20","mag21","mag22","mag23","mag24+"]
pd.set_option("display.width", 260)

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()

d = pd.read_parquet(CAL)
names = {f"prob_maps_grid_dlon{dl:+04d}_lat{la:+03d}.npz" for dl, la in CENTERS}
d = d[d.prob_map_file.isin(names)].reset_index(drop=True)
d["is_neo"] = (d.population == "NEO").astype(int)
d["vmax"] = np.maximum(d.vlam.abs(), d.vbeta.abs())
print(f"CAL rows at the 16 pilot centers: {len(d):,}  (NEO {int(d.is_neo.sum()):,})")

# ---- score every candidate on the SAME rows ----
scores = {t: np.full(len(d), np.nan) for t in THRS}
covered = {t: np.zeros(len(d), bool) for t in THRS}
for cen, g in d.groupby("prob_map_file"):
    idx = g.index.to_numpy()
    for t in THRS:
        pm = v.ProbMapSet.from_npz(str(W/f"prob_maps_e0_thr{t}"/cen), support_mask_min=1,
                                   mask_radius_deg_per_day=np.inf)
        res = pm.score_visible(g.vlam.to_numpy(float), g.vbeta.to_numpy(float),
                               g.mean_mag.to_numpy(float))
        tot = np.sum([np.nan_to_num(np.asarray(res[p], float)) for p in POPS], axis=0)
        scores[t][idx] = np.asarray(res["NEO"], float)
        covered[t][idx] = tot > 0          # some population has density -> a real prediction
        del pm
common = np.logical_and.reduce([covered[t] for t in THRS])
print(f"common scorable rows (all 4 candidates): {int(common.sum()):,} "
      f"({100*common.mean():.2f}% of eligible)")

y = d.is_neo.to_numpy()
rows = []
for t in THRS:
    cov = covered[t]
    paucs, elig = [], 0
    for cen, g in d.groupby("prob_map_file"):
        m = common & d.prob_map_file.eq(cen).to_numpy()
        yy, ss = y[m], scores[t][m]
        if yy.sum() >= 30 and (1-yy).sum() >= 30 and np.isfinite(ss).all():
            elig += 1
            paucs.append(roc_auc_score(yy, ss, max_fpr=0.01))
    macro = float(np.mean(paucs)) if paucs else np.nan
    covm = d[cov]
    rec = dict(threshold=t, eligible_centers=elig,
               common_NEO=int(y[common].sum()), common_nonNEO=int((1-y[common]).sum()),
               macro_pAUC_fpr01=macro,
               pAUC_min=float(np.min(paucs)) if paucs else np.nan,
               pAUC_median=float(np.median(paucs)) if paucs else np.nan,
               pAUC_max=float(np.max(paucs)) if paucs else np.nan,
               coverage_pct=100*float(cov.mean()), abstention_pct=100*float(1-cov.mean()),
               cov_NEO=100*float(cov[y == 1].mean()), cov_nonNEO=100*float(cov[y == 0].mean()))
    for b in ["(0,0.25]", "(0.25,0.5]", "(0.5,1]", "(1,2]", ">2"]:
        pass
    vb = pd.cut(d.vmax, [0,0.25,0.5,1.0,2.0,5.1])
    for lab, gg in pd.Series(cov).groupby(vb, observed=True):
        rec[f"cov_v_{str(lab)}"] = 100*float(gg.mean())
    mb = pd.cut(d.mean_mag, [14,20,22,24,25])
    for lab, gg in pd.Series(cov).groupby(mb, observed=True):
        rec[f"cov_mag_{str(lab)}"] = 100*float(gg.mean())
    rows.append(rec)
res = pd.DataFrame(rows)
best_cov = res.coverage_pct.max()
res["coverage_loss_pp"] = best_cov - res.coverage_pct
res["disqualified_gt2pp"] = res.coverage_loss_pp > 2.0

print("\n" + "="*118); print("CHECK D/E -- candidate comparison on identical common CAL rows")
print(res[["threshold","eligible_centers","common_NEO","common_nonNEO","macro_pAUC_fpr01",
           "pAUC_min","pAUC_median","pAUC_max","coverage_pct","abstention_pct",
           "coverage_loss_pp","disqualified_gt2pp"]].to_string(index=False, float_format=lambda x: f"{x:,.5g}"))
print("\ncoverage by population / velocity / magnitude:")
cols = ["threshold","cov_NEO","cov_nonNEO"] + [c for c in res.columns if c.startswith(("cov_v_","cov_mag_"))]
print(res[cols].to_string(index=False, float_format=lambda x: f"{x:,.3f}"))

# ---- Check A: density normalisation on density_unsmoothed ----
print("\n" + "="*118); print("CHECK A -- density normalisation (unsmoothed; expect R ~ b ~ 1.11)")
prov = json.loads((W/"outputs/splits/split_provenance.json").read_text())
A = []
for t in [THRS[0]]:                       # unsmoothed array is identical across thresholds (proven)
    for dl, la in CENTERS[:6]:
        z = np.load(W/f"prob_maps_e0_thr{t}"/f"prob_maps_grid_dlon{dl:+04d}_lat{la:+03d}.npz",
                    allow_pickle=True)
        x, yy2 = z["x_grid"], z["y_grid"]; dA = float((x[1]-x[0])*(yy2[1]-yy2[0]))
        for pop in ["MBA", "TNO", "Trojans"]:
            for b in BINS:
                k = f"density_unsmoothed__{pop}__{b}"
                sk = f"support_count__{pop}__{b}"
                if k not in z.files: continue
                sup = np.nan_to_num(np.asarray(z[sk], float)).sum()
                if sup < 1000: continue
                integ = np.nan_to_num(np.asarray(z[k], float)).sum()*dA
                f = prov["f_by_population_magbin"][pop].get(b) or prov["f_gen_by_population"][pop]
                # N_full = GEN support / f  (GEN-built support represents fraction f of the truth)
                A.append(dict(center=f"{dl:+04d}/{la:+03d}", pop=pop, bin=b,
                              gen_support=sup, f=f, N_full=sup/f, R=integ/(sup/f)))
Ad = pd.DataFrame(A)
if len(Ad):
    s = Ad.groupby("pop").R.agg(["count","min","median","max"]); s["spread"] = s["max"]-s["min"]
    print(s.to_string(float_format=lambda x: f"{x:,.4f}"))
    okA = bool(((Ad.R >= 1.05) & (Ad.R <= 1.20)).all() and (s["spread"] <= 0.10).all())
else:
    okA = False
print(f"  gated bins: {len(Ad)}   in [1.05,1.20]: {int(((Ad.R>=1.05)&(Ad.R<=1.20)).sum()) if len(Ad) else 0}")
print(f"  CHECK A: {'PASS' if okA else 'FAIL'}")

# ---- Check B: support semantics ----
print("\n" + "="*118); print("CHECK B -- support is statistical (integral, unscaled)")
nonint = 0; ratios = []
for dl, la in CENTERS[:6]:
    z = np.load(W/f"prob_maps_e0_thr2"/f"prob_maps_grid_dlon{dl:+04d}_lat{la:+03d}.npz", allow_pickle=True)
    z10 = np.load(W/f"prob_maps_e0_thr10"/f"prob_maps_grid_dlon{dl:+04d}_lat{la:+03d}.npz", allow_pickle=True)
    for pop in POPS:
        for b in BINS:
            k = f"support_count__{pop}__{b}"
            if k not in z.files: continue
            s2 = np.nan_to_num(np.asarray(z[k], float))
            if not np.all(s2 == np.floor(s2)): nonint += 1
            ratios.append(float(np.array_equal(s2, np.nan_to_num(np.asarray(z10[k], float)))))
okB = (nonint == 0) and all(r == 1.0 for r in ratios)
print(f"  all support values integral: {nonint == 0}   ({nonint} violations)")
print(f"  support identical across thresholds 2 vs 10: {all(r == 1.0 for r in ratios)}")
print(f"  CHECK B: {'PASS' if okB else 'FAIL'}")

# ---- Check C: GEN/CAL disjointness ----
man = pd.read_parquet(W/"outputs/splits/nonneo_split_manifest.parquet")
ov = len(set(man.ObjID[man.split == "GEN"]) & set(man.ObjID[man.split == "CAL"]))
okC = ov == 0
print("\n" + "="*118); print(f"CHECK C -- GEN n CAL overlap = {ov}   {'PASS' if okC else 'FAIL'}")

# ---- selection: preregistered rule ----
elig = res[~res.disqualified_gt2pp]
best = elig.macro_pAUC_fpr01.max()
tied = elig[np.isclose(elig.macro_pAUC_fpr01, best, rtol=0, atol=1e-6)]
selected = int(tied.threshold.max())          # frozen tie-break: LARGER threshold
res["selected"] = res.threshold == selected
print("\n" + "="*118)
print(f"SELECTION (preregistered): max macro pAUC@FPR<=0.01 among non-disqualified; ties -> larger threshold")
print(f"  best macro pAUC = {best:.6f}   tied candidates = {list(tied.threshold)}   SELECTED = {selected}")
overall = okA and okB and okC and len(elig) > 0
print(f"\nCHECKS: A={'PASS' if okA else 'FAIL'} B={'PASS' if okB else 'FAIL'} "
      f"C={'PASS' if okC else 'FAIL'} D=reported E=reported -> E0 {'PASS' if overall else 'FAIL'}")

res.to_csv(OUT/"E0_RESULTS.csv", index=False)
meta = {
    "checks": {"A": bool(okA), "B": bool(okB), "C": bool(okC)},
    "selected_threshold": selected, "best_macro_pAUC_fpr01": float(best),
    "tie_break": "larger threshold (frozen)",
    "common_rows": int(common.sum()), "common_NEO": int(y[common].sum()),
    "cal_parquet_sha256": sha(CAL),
    "cal_seal_sha256": sha(W/"outputs/splits/CAL_DATASET_SEAL.json"),
    "gen_manifest_sha256": sha(W/"outputs/splits/GEN_MANIFEST.json"),
    "gen_receipt_sha256": sha(W/"outputs/splits/GEN_VERIFICATION_RECEIPT.json"),
    "frozen_env_seal_sha256": sha(W/"outputs/splits/FROZEN_ENV_SEAL.json"),
    "prereg_sha256": sha(W/"neomod/docs/E0_PILOT_PREREGISTRATION.md"),
    "git_commit": os.popen(f"git -C {W}/neomod rev-parse HEAD").read().strip(),
    "map_dirs": {t: f"prob_maps_e0_thr{t}" for t in THRS},
    "results_csv_sha256": sha(OUT/"E0_RESULTS.csv"),
}
(OUT/"E0_RESULTS.json").write_text(json.dumps(meta, indent=2, sort_keys=True))
print(f"\nwrote {OUT}/E0_RESULTS.csv and .json")
print(f"  E0_RESULTS.json sha256 = {sha(OUT/'E0_RESULTS.json')}")
