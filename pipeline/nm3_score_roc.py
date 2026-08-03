#!/usr/bin/env python3
"""Score Sorcha tracklets with the NEOMOD3 clone-only maps and compare to production VDP + digest2.

For each test center and each tracklet set (full 2-yr and night-61642):
  * re-score VDP against prob_maps_grid_neomod3_vlim5 (production settings:
    support_mask_min=1, nearest-dist mask OFF, Johnson-V magnitudes),
  * compare against the STORED production scores (P_NEO_vdp_Vband from the +-2 S3M/GMM maps)
    and digest2 (P_NEO_d2),
  * report AUC / best-F1 / completeness / contamination, and specifically what happens to the
    tracklets with |v| > 2 that the production maps could not score at all.

Read-only. Output: outputs/neomod3_score_roc/
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, roc_auc_score, roc_curve

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
os.chdir(W / "neomod")
NEWMAPS = Path(os.environ.get("NM3_MAPS_DIR", str(W / "prob_maps_grid_neomod3_vlim5")))
SORCHA = W / "outputs/s3m_linking/case1/sorcha_comparison_case1_nbody_Vband.parquet"
OUT = Path(os.environ.get("NM3_OUT_DIR", str(W / "outputs" / "neomod3_score_roc")))
CENTERS = ["prob_maps_grid_dlon+000_lat+00.npz", "prob_maps_grid_dlon+020_lat-12.npz"]
NIGHT = 61642
import velocity_density_pipeline_neomod_clone_only as vnm


def metrics(y, s):
    s = np.asarray(s, float); y = np.asarray(y, int)
    ok = np.isfinite(s)
    if ok.sum() < 10 or y[ok].sum() < 3 or (1 - y[ok]).sum() < 3:
        return None
    s, y = s[ok], y[ok]
    p, r, th = precision_recall_curve(y, s)
    f1 = np.divide(2 * p * r, p + r, out=np.zeros_like(p), where=(p + r) > 0)
    bi = int(np.argmax(f1[:-1]))
    return dict(n=len(y), n_neo=int(y.sum()), AUC=roc_auc_score(y, s), bestF1=f1[bi],
                completeness=r[bi] * 100, contamination=(1 - p[bi]) * 100,
                thresh=float(th[bi]) if bi < len(th) else 1.0)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    trk = pd.read_parquet(SORCHA, columns=[
        "ObjID", "night", "population", "vlam", "vbeta", "mean_mag_V", "prob_map_file",
        "P_NEO_vdp_Vband", "P_NEO_d2"])
    rows, scored_all = [], []
    for cen in CENTERS:
        sub = trk[trk.prob_map_file == cen].reset_index(drop=True)
        if not len(sub):
            continue
        pm = vnm.ProbMapSet.from_npz(str(NEWMAPS / cen), support_mask_min=1,
                                     mask_radius_deg_per_day=np.inf)
        sub = sub.copy()
        sub["P_NEO_neomod3"] = pm.score_visible(
            sub.vlam.to_numpy(float), sub.vbeta.to_numpy(float), sub.mean_mag_V.to_numpy(float))["NEO"]
        sub["maxabs_v"] = np.maximum(sub.vlam.abs(), sub.vbeta.abs())
        sub["is_neo"] = (sub.population == "NEO").astype(int)
        sub["center"] = cen
        scored_all.append(sub)
        del pm

        for setname, d in [("full_2yr", sub), (f"night_{NIGHT}", sub[sub.night == NIGHT])]:
            for label, col in [("NEOMOD3 (new, ±5)", "P_NEO_neomod3"),
                               ("production VDP (±2)", "P_NEO_vdp_Vband"),
                               ("digest2", "P_NEO_d2")]:
                m = metrics(d.is_neo.to_numpy(), d[col].to_numpy())
                if m:
                    rows.append(dict(center=cen, tracklets=setname, classifier=label, **m))
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "neomod3_roc_summary.csv", index=False)
    allsc = pd.concat(scored_all, ignore_index=True)
    allsc.to_parquet(OUT / "neomod3_scored_tracklets.parquet", index=False)

    pd.set_option("display.width", 240)
    print("\n================ ROC / F1 ================")
    print(res[["center", "tracklets", "classifier", "n", "n_neo", "AUC", "bestF1",
               "completeness", "contamination"]].to_string(index=False, float_format=lambda v: f"{v:,.4g}"))

    # ---- the |v|>2 population: what production structurally could not score ----
    print("\n================ the |v| > 2 tracklets (production grid edge) ================")
    fast = allsc[allsc.maxabs_v > 2.0]
    print(f"tracklets with |v|>2 at these centers: {len(fast):,}  "
          f"(NEO {int(fast.is_neo.sum()):,}, non-NEO {int((1-fast.is_neo).sum()):,})")
    if len(fast):
        pz = (fast.P_NEO_vdp_Vband.fillna(0) == 0).mean()
        nz = (fast.P_NEO_neomod3.fillna(0) == 0).mean()
        print(f"  production VDP == 0 : {100*pz:.1f}%     NEOMOD3 == 0 : {100*nz:.1f}%")
        for lab, g in fast.groupby(fast.is_neo.map({1: "NEO", 0: "non-NEO"})):
            print(f"  {lab:>8}: production median {g.P_NEO_vdp_Vband.median():.4f} -> "
                  f"NEOMOD3 median {g.P_NEO_neomod3.median():.4f}   (n={len(g):,})")

    # ---- figure ----
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    for i, cen in enumerate(CENTERS):
        for j, setname in enumerate(["full_2yr", f"night_{NIGHT}"]):
            ax = axes[i, j]
            d = allsc[allsc.center == cen]
            if setname.startswith("night"):
                d = d[d.night == NIGHT]
            for label, col, c in [("NEOMOD3 (new, ±5)", "P_NEO_neomod3", "tab:green"),
                                  ("production VDP (±2)", "P_NEO_vdp_Vband", "tab:blue"),
                                  ("digest2", "P_NEO_d2", "tab:orange")]:
                m = metrics(d.is_neo.to_numpy(), d[col].to_numpy())
                if not m:
                    continue
                s = np.nan_to_num(d[col].to_numpy(float))
                p, r, _ = precision_recall_curve(d.is_neo.to_numpy(), s)
                ax.plot(r * 100, (1 - p) * 100, lw=2, color=c,
                        label=f"{label}  F1={m['bestF1']:.3f} AUC={m['AUC']:.3f}")
            ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.grid(alpha=0.3)
            ax.set_xlabel("NEO completeness (%)"); ax.set_ylabel("contamination (%)")
            ax.set_title(f"{cen.replace('prob_maps_grid_','').replace('.npz','')}  |  {setname}\n"
                         f"N={len(d):,}  NEO={int(d.is_neo.sum()):,}", fontsize=10)
            ax.legend(fontsize=8)
    fig.suptitle("NEOMOD3 clone-only (±5) vs production VDP (±2) vs digest2", fontsize=13)
    fig.savefig(OUT / "neomod3_roc.png", dpi=150, bbox_inches="tight")
    print(f"\nwrote -> {OUT}")


if __name__ == "__main__":
    main()
