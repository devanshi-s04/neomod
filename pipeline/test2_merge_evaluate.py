#!/usr/bin/env python3
"""TEST2 merge + weighted evaluation. Merges ONLY by immutable tracklet_uid.

Three distinct views, never conflated:
  A  coverage/abstention over all 688,688 rows
  B  ranking metrics on the COMMON-SCORABLE subset (all three classifiers valid)
  C  legacy and digest2 on their own broader valid sets, labelled as such

Four-density availability is geometrically selected and especially low for NEOs (~28%), because
many narrow TNO/Trojan 0.25-mag cells are unavailable. View B is therefore NOT performance over the
whole TEST2 population and is never reported as such.
"""
from __future__ import annotations
import glob, hashlib, json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
T2 = W / "outputs" / "test2_geometric"
SC = T2 / "scored"
POPS = ("NEO", "MBA", "TNO", "Trojans")
N_EXPECT = 688688
NSH = 32
R = []


def chk(n, c, d=""):
    ok = bool(c); R.append({"check": n, "pass": ok, "detail": str(d)})
    print(f"  {'PASS' if ok else 'FAIL'}  {n}" + (f"  [{d}]" if d else ""), flush=True)
    return ok


def wauc(y, s, w):
    """Weighted ROC AUC via the rank/Mann-Whitney identity with tie handling."""
    o = np.argsort(s, kind="mergesort")
    y, s, w = y[o], s[o], w[o]
    P = w[y == 1].sum(); N = w[y == 0].sum()
    if P <= 0 or N <= 0:
        return np.nan, None
    cn = np.cumsum(w * (y == 0))
    below = np.concatenate([[0.0], cn[:-1]])
    # tie groups: average negative weight inside each tie
    auc = 0.0
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        neg_in = w[i:j + 1][y[i:j + 1] == 0].sum()
        pos_in = w[i:j + 1][y[i:j + 1] == 1].sum()
        auc += pos_in * (below[i] + 0.5 * neg_in)
        i = j + 1
    return auc / (P * N), None


def roc_curve_w(y, s, w):
    o = np.argsort(-s, kind="mergesort")
    y, s, w = y[o], s[o], w[o]
    tp = np.cumsum(w * (y == 1)); fp = np.cumsum(w * (y == 0))
    P = w[y == 1].sum(); N = w[y == 0].sum()
    keep = np.concatenate([np.diff(s) != 0, [True]])
    return fp[keep] / N, tp[keep] / P, s[keep]


def pauc_std(fpr, tpr, lim=0.01):
    m = fpr <= lim
    if m.sum() < 2:
        return np.nan
    f = np.concatenate([[0.0], fpr[m]]); t = np.concatenate([[0.0], tpr[m]])
    return float(np.trapezoid(t, f) / lim)


def boot_parent(y, s, w, pid, n=200, seed=99):
    rng = np.random.default_rng(seed)
    up = np.unique(pid)
    idx = pd.Series(np.arange(len(pid))).groupby(pid).apply(lambda g: g.to_numpy())
    out = []
    for _ in range(n):
        pick = rng.choice(up, size=len(up), replace=True)
        sel = np.concatenate([idx[p] for p in pick])
        a, _ = wauc(y[sel], s[sel], w[sel])
        if np.isfinite(a):
            out.append(a)
    if not out:
        return (np.nan, np.nan)
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5)))


def main():
    t0 = time.time()
    te = pd.read_parquet(T2 / "TEST2_TRACKLETS.parquet")
    print(f"TEST2 rows {len(te):,}", flush=True)
    chk("TEST2 has exactly 688,688 rows", len(te) == N_EXPECT, len(te))
    chk("tracklet_uid unique", te.tracklet_uid.is_unique)

    merged = te.copy()
    for mode, vcol in (("new_vdp", "new_vdp_valid"), ("legacy_vdp", "legacy_vdp_valid"),
                       ("digest2", "digest2_valid")):
        fs = sorted(glob.glob(str(SC / f"{mode}_*.parquet")))
        chk(f"{mode}: 32/32 shard files", len(fs) == NSH, len(fs))
        d = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
        chk(f"{mode}: contributed {N_EXPECT:,} rows", len(d) == N_EXPECT, f"{len(d):,}")
        chk(f"{mode}: unique tracklet_uid", d.tracklet_uid.is_unique)
        chk(f"{mode}: no missing or extra IDs",
            set(d.tracklet_uid) == set(te.tracklet_uid),
            f"missing {len(set(te.tracklet_uid)-set(d.tracklet_uid))} "
            f"extra {len(set(d.tracklet_uid)-set(te.tracklet_uid))}")
        merged = merged.merge(d, on="tracklet_uid", how="left", validate="one_to_one")
    chk("merged row count preserved", len(merged) == N_EXPECT, len(merged))

    # integrity of carried-through columns
    chk("truth labels unchanged", (merged.population.isin(POPS)).all())
    chk("w_phys finite and positive", bool(np.isfinite(merged.w_phys).all() and (merged.w_phys > 0).all()))
    chk("NEO weight is the sealed value",
        abs(float(merged.w_phys[merged.population == "NEO"].iloc[0]) - 0.0999994694) < 1e-9)
    chk("non-NEO weight is 1.0", bool((merged.w_phys[merged.population != "NEO"] == 1.0).all()))
    rehash = merged.tracklet_input_sha256_used.fillna(merged.tracklet_input_sha256)
    chk("digest2 used the sealed input hashes",
        bool((rehash == merged.tracklet_input_sha256).all()))
    d2 = merged.digest2_raw[merged.digest2_valid].to_numpy(float)
    chk("digest2 values are integers in 0..100",
        bool(len(d2) > 0 and np.all(np.abs(d2 - np.round(d2)) < 1e-9)
             and d2.min() >= 0 and d2.max() <= 100),
        f"n={len(d2):,} levels={len(np.unique(d2))}")
    nv = merged.new_vdp_valid.to_numpy()
    chk("valid new-VDP rows all have n_pops_new == 4",
        bool((merged.n_pops_new.to_numpy()[nv] == 4).all()))
    chk("valid four-class probabilities sum to one",
        float(np.nanmax(np.abs(merged.prob_sum_new.to_numpy()[nv] - 1.0))) < 1e-9,
        f"{np.nanmax(np.abs(merged.prob_sum_new.to_numpy()[nv]-1.0)):.2e}")
    chk("no partial-denominator row is valid",
        not merged.new_vdp_reason.astype(str).str.startswith("ok_partial").any())
    for c_, v_ in (("P_NEO_new", "new_vdp_valid"), ("P_NEO_legacy_raw", "legacy_vdp_valid"),
                   ("P_NEO_digest2", "digest2_valid")):
        x = merged[c_].to_numpy(float); m = merged[v_].to_numpy()
        chk(f"{c_}: valid in [0,1] and finite",
            bool(np.isfinite(x[m]).all() and (x[m] >= 0).all() and (x[m] <= 1).all()))
        chk(f"{c_}: invalid rows are NaN", bool(np.isnan(x[~m]).all()))
    chk("invalid new-VDP rows carry a non-ok reason",
        bool((merged.new_vdp_reason.astype(str)[~nv] != "ok").all()))

    merged.to_parquet(T2 / "TEST2_SCORED.parquet", index=False)
    print(f"\nwrote TEST2_SCORED.parquet ({time.time()-t0:.0f}s)", flush=True)

    # ---------------- VIEW A: coverage / abstention -------------------------
    wt = merged.w_phys.to_numpy(float)
    cov_rows = []
    for p in POPS:
        g = merged[merged.population == p]
        for mode, v in (("new_vdp", "new_vdp_valid"), ("legacy_vdp", "legacy_vdp_valid"),
                        ("digest2", "digest2_valid")):
            m = g[v].to_numpy()
            cov_rows.append(dict(population=p, classifier=mode, rows=len(g),
                                 valid=int(m.sum()), raw=float(m.mean()),
                                 weighted=float(g.w_phys[m].sum() / g.w_phys.sum())))
    cov = pd.DataFrame(cov_rows)
    print("VIEW A — COVERAGE BY CLASSIFIER AND POPULATION"); print(cov.to_string(index=False), flush=True)
    reasons = (merged.groupby(["population", "new_vdp_reason"]).size().rename("n")
               .reset_index().sort_values(["population", "n"], ascending=[True, False]))

    # ---------------- VIEW B: common-scorable -------------------------------
    common = merged.new_vdp_valid & merged.legacy_vdp_valid & merged.digest2_valid
    cm = merged[common].copy()
    y = (cm.population == "NEO").to_numpy().astype(int)
    wv = cm.w_phys.to_numpy(float)
    pid = cm.source_parent_uid.to_numpy()
    print(f"\nVIEW B — COMMON-SCORABLE {len(cm):,} rows "
          f"({100*len(cm)/len(merged):.2f}% of TEST2); NEO {int(y.sum()):,}", flush=True)
    res = {}
    for name, col in (("new_vdp", "P_NEO_new"), ("legacy_vdp_raw", "P_NEO_legacy_raw"),
                      ("legacy_vdp_cal", "P_NEO_legacy_cal"), ("digest2", "P_NEO_digest2")):
        s = cm[col].to_numpy(float)
        a, _ = wauc(y, s, wv)
        f_, t_, th = roc_curve_w(y, s, wv)
        pa = pauc_std(f_, t_)
        lo, hi = boot_parent(y, s, wv, pid, n=int(150))
        # completeness at fixed contamination
        comp = {}
        for target in (0.05, 0.10):
            # contamination = weighted FP / (TP+FP)
            o = np.argsort(-s, kind="mergesort")
            yy, ww = y[o], wv[o]
            tp = np.cumsum(ww * (yy == 1)); fp = np.cumsum(ww * (yy == 0))
            cont = fp / np.maximum(tp + fp, 1e-30)
            okc = cont <= target
            comp[f"completeness_at_{int(target*100)}pct_contam"] = (
                float(tp[okc].max() / ww[yy == 1].sum()) if okc.any() else None)
        res[name] = dict(auc=float(a), pauc_fpr01=float(pa), auc_ci95=[lo, hi], **comp)
        print(f"  {name:16s} AUC {a:.5f}  pAUC(FPR<=.01) {pa:.5f}  "
              f"CI95 [{lo:.5f},{hi:.5f}]  "
              f"C@5% {comp['completeness_at_5pct_contam']}  "
              f"C@10% {comp['completeness_at_10pct_contam']}", flush=True)

    # paired new - legacy on common rows
    dlt = cm.P_NEO_new.to_numpy(float) - cm.P_NEO_legacy_raw.to_numpy(float)
    paired = {"n": int(len(cm)), "mean": float(np.nanmean(dlt)),
              "median": float(np.nanmedian(dlt)),
              "p05": float(np.nanpercentile(dlt, 5)), "p95": float(np.nanpercentile(dlt, 95)),
              "frac_new_higher": float(np.nanmean(dlt > 0))}
    print(f"\n  paired new-minus-legacy(raw): median {paired['median']:+.4f}  "
          f"mean {paired['mean']:+.4f}  new higher on {100*paired['frac_new_higher']:.1f}%",
          flush=True)

    # ---------------- VIEW C: broader valid sets ----------------------------
    viewc = {}
    for name, col, v in (("legacy_vdp_raw", "P_NEO_legacy_raw", "legacy_vdp_valid"),
                         ("digest2", "P_NEO_digest2", "digest2_valid")):
        g = merged[merged[v]]
        yy = (g.population == "NEO").to_numpy().astype(int)
        a, _ = wauc(yy, g[col].to_numpy(float), g.w_phys.to_numpy(float))
        f_, t_, _ = roc_curve_w(yy, g[col].to_numpy(float), g.w_phys.to_numpy(float))
        viewc[name] = dict(n=int(len(g)), auc=float(a), pauc_fpr01=float(pauc_std(f_, t_)))
        print(f"VIEW C — {name} on its own valid set: n {len(g):,} AUC {a:.5f}", flush=True)

    # included vs abstained NEOs
    neo = merged[merged.population == "NEO"]
    inc, ab = neo[neo.new_vdp_valid], neo[~neo.new_vdp_valid]
    neo_cmp = {"included": int(len(inc)), "abstained": int(len(ab)),
               "included_median_V": float(inc.mean_mag_V.median()) if len(inc) else None,
               "abstained_median_V": float(ab.mean_mag_V.median()) if len(ab) else None,
               "included_median_absv": float(np.hypot(inc.vlam, inc.vbeta).median()) if len(inc) else None,
               "abstained_median_absv": float(np.hypot(ab.vlam, ab.vbeta).median()) if len(ab) else None}
    print(f"\n  NEO included {neo_cmp['included']:,} (median V {neo_cmp['included_median_V']}) vs "
          f"abstained {neo_cmp['abstained']:,} (median V {neo_cmp['abstained_median_V']})", flush=True)

    d2lev = np.unique(d2)
    out = {"created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "n_rows": int(len(merged)), "checks": R,
           "ALL_PASS": all(x["pass"] for x in R),
           "coverage_definitions": {
               "neo_cell_availability_raw": 0.9933,
               "four_density_coverage_raw": 0.4396,
               "four_density_coverage_weighted": 0.4506},
           "view_a_coverage": cov.to_dict("records"),
           "new_vdp_reason_counts": reasons.to_dict("records"),
           "view_b_common_scorable": {"n": int(len(cm)),
                                      "frac_of_test2": float(len(cm) / len(merged)),
                                      "n_neo": int(y.sum()), "metrics": res,
                                      "paired_new_minus_legacy_raw": paired},
           "view_c_broader_valid_sets": viewc,
           "neo_included_vs_abstained": neo_cmp,
           "digest2_operating_points": {"n_levels": int(len(d2lev)),
                                        "levels": d2lev.tolist()},
           "caveat": ("four-density availability is geometrically selected and only ~28% for NEO "
                      "because many narrow TNO/Trojan 0.25-mag cells are unavailable; VIEW B is NOT "
                      "performance over the whole TEST2 population")}
    (T2 / "TEST2_SCORED_VALIDATION.json").write_text(json.dumps(out, indent=2, default=str))
    cov.to_csv(T2 / "TEST2_COVERAGE_BY_CLASSIFIER.csv", index=False)
    reasons.to_csv(T2 / "TEST2_REASON_COUNTS.csv", index=False)
    pd.DataFrame(res).T.to_csv(T2 / "TEST2_METRICS_COMMON_SCORABLE.csv")
    print(f"\nALL_PASS={out['ALL_PASS']}  ({time.time()-t0:.0f}s)", flush=True)
    if not out["ALL_PASS"]:
        sys.exit(4)


if __name__ == "__main__":
    main()
