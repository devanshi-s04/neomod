#!/usr/bin/env python3
"""Pilot scorer: p_NEO_vs_OTHER = rho_NEO / (rho_NEO + rho_OTHER) on the selected-centre pilot maps.

SEPARATE from the frozen Rubin v1 scorer. It does not import, wrap, modify or re-seal it. The
frozen scorer and its sealed maps are untouched.

Geometry is deliberately identical to pipeline/score_test2.score_new_vdp: the same bilinear-in-
velocity-only interpolation helper is imported from that module rather than reimplemented, the same
single-containing-0.25-mag-slice rule applies, and the same +-5 deg/day bound is enforced. The only
difference is the denominator: two densities instead of four.

Missing densities are never replaced by zero or an epsilon. A row whose cell lacks either rho_NEO or
rho_OTHER is invalid with an explicit reason, exactly as the frozen scorer treats a missing
population.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"pipeline"))
from score_test2 import bil                      # identical interpolation, imported not copied

PILOT_MAPS = W/"outputs"/"neomod3_pneo_other_pilot"/"pilot_maps"


def score_pneo_other(te, pilot_maps=PILOT_MAPS):
    n = len(te)
    rho_neo   = np.full(n, np.nan)
    rho_other = np.full(n, np.nan)
    P         = np.full(n, np.nan)
    tot       = np.full(n, np.nan)
    reason    = np.array(["ok"]*n, dtype=object)

    for cen, g in te.groupby("center_label", sort=True):
        f = pilot_maps/f"pilot_{cen}.npz"
        idx = g.index.to_numpy()
        if not f.exists():
            reason[idx] = "missing_pilot_map"; continue
        z = np.load(f, allow_pickle=True)
        keys = set(z.keys())
        xg = np.asarray(z["x_grid"], float); yg = np.asarray(z["y_grid"], float)
        for lab, gg in g.groupby("magnitude_bin", sort=True):
            ii = gg.index.to_numpy()
            if lab is None or (isinstance(lab, float) and np.isnan(lab)):
                reason[ii] = "v_out_of_range"; continue
            xs = gg.vlam.to_numpy(float); ys = gg.vbeta.to_numpy(float)
            inb = (np.abs(xs) <= 5.0) & (np.abs(ys) <= 5.0)

            kn, ko = f"rho_NEO__{lab}", f"rho_OTHER__{lab}"
            missing = [nm for nm, k in (("NEO", kn), ("OTHER", ko)) if k not in keys]
            if missing:
                # absent means the cell was INVALID under the support rule; never zero-filled
                reason[ii] = "missing_population_density:" + ",".join(missing)
                continue
            dn = bil(np.asarray(z[kn], np.float64), xg, yg, xs, ys)
            do = bil(np.asarray(z[ko], np.float64), xg, yg, xs, ys)
            rho_neo[ii] = dn; rho_other[ii] = do
            tt = dn + do
            good = inb & np.isfinite(tt) & (tt > 0) & np.isfinite(dn) & np.isfinite(do)
            pn = np.full(len(ii), np.nan)
            pn[good] = dn[good] / tt[good]
            P[ii] = pn; tot[ii] = tt
            reason[ii] = np.where(~inb, "outside_velocity_grid",
                          np.where(~np.isfinite(tt), "nonfinite_density",
                          np.where(tt <= 0, "zero_total_density", "ok")))

    out = pd.DataFrame({"p_NEO_vs_OTHER": P, "rho_NEO_pilot": rho_neo,
                        "rho_OTHER_pilot": rho_other, "total_density_pilot": tot,
                        "pilot_reason": reason}, index=te.index)
    out["pilot_valid"] = np.isfinite(P)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracklets", required=True, help="parquet with the TEST2 rows to score")
    ap.add_argument("--output", required=True)
    ap.add_argument("--pilot-maps", default=str(PILOT_MAPS))
    a = ap.parse_args()
    t0 = time.time()
    te = pd.read_parquet(a.tracklets).reset_index(drop=True)
    print(f"rows {len(te):,}  centres {te.center_label.nunique()}", flush=True)
    res = score_pneo_other(te, Path(a.pilot_maps))
    out = pd.concat([te, res], axis=1)
    out.to_parquet(a.output, index=False)
    print(f"valid {int(res.pilot_valid.sum()):,}/{len(res):,} "
          f"({100*res.pilot_valid.mean():.2f}%)  {time.time()-t0:.0f}s")
    print(res.pilot_reason.value_counts().to_string())


if __name__ == "__main__":
    main()
