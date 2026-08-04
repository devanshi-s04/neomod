#!/usr/bin/env python3
"""E0_SIMPLE_COMPARISON -- pooled map-version comparison on sealed CAL v2.

Replaces the preregistered per-center macro selection (which was undefined: 0/16 centers reached
the >=30 NEO eligibility rule). Pools all 16 pilot centers, uses IDENTICAL rows for all four
candidates, and restricts to rows scorable by all four. Primary comparator: pooled normalized
pAUC at FPR<=0.01. Checks B/C remain QA and already pass; Check A is a density diagnostic and is
NOT used to select (it is identical across candidates by construction -- density_unsmoothed is
bit-identical, proven).
"""
import hashlib, json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); os.chdir(W/"neomod")
import velocity_density_pipeline_neomod_clone_only as v

THRS = [2, 3, 5, 10]
CENTERS = [(0,0),(20,-12),(10,2),(-10,-2),(0,8),(50,-1),(50,1),(-50,1),
           (90,0),(140,0),(-140,0),(0,25),(0,-25),(0,50),(0,-50),(-50,-50)]
CAL = W/"outputs/cal_tracklets_neomod3_v2/tracklets_benchmark_neomod3.parquet"
OUT = W/"outputs/e0_results"; OUT.mkdir(parents=True, exist_ok=True)
POPS = ["NEO", "MBA", "TNO", "Trojans"]
COL = {2: "tab:blue", 3: "tab:orange", 5: "tab:green", 10: "tab:red"}

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()

d = pd.read_parquet(CAL)
names = {f"prob_maps_grid_dlon{a:+04d}_lat{b:+03d}.npz" for a, b in CENTERS}
d = d[d.prob_map_file.isin(names)].reset_index(drop=True)
y_all = (d.population == "NEO").astype(int).to_numpy()
print(f"CAL rows at 16 pilot centers: {len(d):,}  (NEO {int(y_all.sum())})")

scores = {t: np.full(len(d), np.nan) for t in THRS}
cov = {t: np.zeros(len(d), bool) for t in THRS}
for cen, g in d.groupby("prob_map_file"):
    idx = g.index.to_numpy()
    for t in THRS:
        pm = v.ProbMapSet.from_npz(str(W/f"prob_maps_e0_thr{t}"/cen), support_mask_min=1,
                                   mask_radius_deg_per_day=np.inf)
        r = pm.score_visible(g.vlam.to_numpy(float), g.vbeta.to_numpy(float),
                             g.mean_mag.to_numpy(float))
        tot = np.sum([np.nan_to_num(np.asarray(r[p], float)) for p in POPS], axis=0)
        scores[t][idx] = np.asarray(r["NEO"], float); cov[t][idx] = tot > 0
        del pm
common = np.logical_and.reduce([cov[t] for t in THRS])
y = y_all[common]
print(f"common scorable rows: {int(common.sum()):,}  (NEO {int(y.sum())}, non-NEO {int((1-y).sum())})")

def comp_at(prec, rec, c):
    ok = (1 - prec) <= c
    return float(rec[ok].max()*100) if ok.any() else 0.0

rows, curves = [], {}
for t in THRS:
    s = np.nan_to_num(scores[t][common])
    auc = roc_auc_score(y, s)
    pauc = roc_auc_score(y, s, max_fpr=0.01)
    p, r, th = precision_recall_curve(y, s)
    f1 = np.divide(2*p*r, p+r, out=np.zeros_like(p), where=(p+r) > 0)
    i = int(np.argmax(f1[:-1]))
    fpr, tpr, _ = roc_curve(y, s)
    curves[t] = dict(fpr=fpr, tpr=tpr, prec=p, rec=r)
    rows.append(dict(threshold=t, coverage_pct=100*float(cov[t].mean()),
                     n_common=int(common.sum()), n_NEO=int(y.sum()),
                     ROC_AUC=auc, pAUC_fpr01_norm=pauc, best_F1=f1[i],
                     completeness_at_5pct_contam=comp_at(p, r, 0.05),
                     completeness_at_10pct_contam=comp_at(p, r, 0.10),
                     best_F1_threshold=float(th[i]) if i < len(th) else 1.0))
res = pd.DataFrame(rows)
pd.set_option("display.width", 240)
print("\n" + "="*110); print("POOLED METRICS -- identical rows, all four map versions")
print(res.to_string(index=False, float_format=lambda x: f"{x:,.6g}"))
spread = res.pAUC_fpr01_norm.max() - res.pAUC_fpr01_norm.min()
print(f"\n  primary comparator: pooled normalized pAUC @ FPR<=0.01")
print(f"  range {res.pAUC_fpr01_norm.min():.6f} .. {res.pAUC_fpr01_norm.max():.6f}   spread = {spread:.6f}")
print(f"  ROC AUC spread = {res.ROC_AUC.max()-res.ROC_AUC.min():.6f}")

print("\n" + "="*110); print("SCORE DIFFERENCES between thresholds (same objects, common rows)")
diffs = []
for a in range(len(THRS)):
    for b in range(a+1, len(THRS)):
        ta, tb = THRS[a], THRS[b]
        da = np.nan_to_num(scores[ta][common]); db = np.nan_to_num(scores[tb][common])
        dd = db - da
        diffs.append(dict(pair=f"{ta}->{tb}", n_changed=int((dd != 0).sum()),
                          pct_changed=100*float((dd != 0).mean()),
                          max_abs=float(np.abs(dd).max()), mean_abs=float(np.abs(dd).mean()),
                          p99_abs=float(np.percentile(np.abs(dd), 99)),
                          neo_mean_abs=float(np.abs(dd[y == 1]).mean())))
dfd = pd.DataFrame(diffs)
print(dfd.to_string(index=False, float_format=lambda x: f"{x:,.6g}"))

fig, ax = plt.subplots(1, 1, figsize=(7.2, 5.6))
for t in THRS:
    c = curves[t]
    ax.plot(c["rec"]*100, (1-c["prec"])*100, lw=1.8, color=COL[t],
            label=f"threshold {t}  (pAUC {res.loc[res.threshold==t,'pAUC_fpr01_norm'].iloc[0]:.4f})")
ax.set_xlabel("NEO completeness (%)"); ax.set_ylabel("contamination (%)")
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.grid(alpha=.3)
ax.set_title(f"E0 CAL: completeness vs contamination\n{int(common.sum()):,} common rows, {int(y.sum())} NEOs")
ax.legend(fontsize=9); fig.tight_layout()
fig.savefig(OUT/"E0_completeness_vs_contamination.png", dpi=150); plt.close(fig)

fig, ax = plt.subplots(1, 1, figsize=(7.2, 5.6))
for t in THRS:
    c = curves[t]
    ax.plot(c["fpr"], c["tpr"], lw=1.8, color=COL[t],
            label=f"threshold {t}  (AUC {res.loc[res.threshold==t,'ROC_AUC'].iloc[0]:.4f})")
ax.plot([0, 1], [0, 1], "k--", lw=.8, alpha=.5)
ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate")
ax.grid(alpha=.3); ax.set_title(f"E0 CAL: ROC\n{int(common.sum()):,} common rows, {int(y.sum())} NEOs")
ax.legend(fontsize=9); fig.tight_layout()
fig.savefig(OUT/"E0_roc.png", dpi=150); plt.close(fig)

res.to_csv(OUT/"E0_SIMPLE_RESULTS.csv", index=False)
dfd.to_csv(OUT/"E0_SIMPLE_SCORE_DIFFS.csv", index=False)
meta = {
    "test": "E0_SIMPLE_COMPARISON", "primary": "pooled normalized pAUC at FPR<=0.01",
    "pooled": True, "per_center_eligibility_applied": False,
    "n_common_rows": int(common.sum()), "n_NEO": int(y.sum()), "n_nonNEO": int((1-y).sum()),
    "pAUC_spread": float(spread), "ROC_AUC_spread": float(res.ROC_AUC.max()-res.ROC_AUC.min()),
    "indistinguishable": bool(spread < 1e-3),
    "metrics": res.to_dict("records"), "score_diffs": dfd.to_dict("records"),
    "cal_parquet_sha256": sha(CAL),
    "cal_seal_sha256": sha(W/"outputs/splits/CAL_DATASET_SEAL.json"),
    "gen_manifest_sha256": sha(W/"outputs/splits/GEN_MANIFEST.json"),
    "gen_receipt_sha256": sha(W/"outputs/splits/GEN_VERIFICATION_RECEIPT.json"),
    "frozen_env_seal_sha256": sha(W/"outputs/splits/FROZEN_ENV_SEAL.json"),
    "git_commit": os.popen(f"git -C {W}/neomod rev-parse HEAD").read().strip(),
    "results_csv_sha256": sha(OUT/"E0_SIMPLE_RESULTS.csv"),
    "note": "Preregistered per-center macro selection was undefined (0/16 centers reached the >=30 "
            "NEO rule); its outputs are preserved as diagnostics. Check A is a density diagnostic "
            "and cannot select (density_unsmoothed is bit-identical across candidates).",
}
(OUT/"E0_SIMPLE_RESULTS.json").write_text(json.dumps(meta, indent=2, default=str))
print(f"\nwrote {OUT}/E0_SIMPLE_RESULTS.{{csv,json}}, E0_SIMPLE_SCORE_DIFFS.csv, 2 figures")
print(f"  E0_SIMPLE_RESULTS.json sha256 = {sha(OUT/'E0_SIMPLE_RESULTS.json')}")
