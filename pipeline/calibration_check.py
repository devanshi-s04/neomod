#!/usr/bin/env python3
"""Calibration check: does a probability map's IMPLIED population mix match what Sorcha detected?

THE TEST
--------
A density map factorises as  rho_pop(v) = N_pop * f_pop(v)  with  int f dv = 1.
So integrating the map over velocity recovers the number of objects the map thinks are in this
(patch, magnitude bin):

    int rho_pop dv  =  N_pop      ->  predicted NEO fraction = N_NEO / sum_pop N_pop

That is an absolute, falsifiable prediction. Sorcha's tracklets assigned to the same map + magnitude
bin give the OBSERVED fraction (we know each object's true population). If the cross-population
normalisation is right, the two should agree.

CAVEAT (stated, not hidden): the map describes the VISIBLE population (sky cut + mag bin) while Sorcha
tracklets are the DETECTED population (after footprint, cadence, SNR, trailing losses, linking).
Trailing loss preferentially removes fast movers, i.e. NEOs, so we expect
    observed NEO fraction  <=  predicted NEO fraction
by a modest factor. This makes the test STRONG for gross (many-x) normalisation errors and WEAK for
fine calibration. The sharper signal is the TREND ACROSS MAGNITUDE BINS -- detection efficiency varies
smoothly with magnitude, so a correct per-bin normalisation should track the observed fractions up to
a smooth factor. Per-bin normalisation is exactly what drives AUC (a single global factor leaves the
ranking, and therefore AUC, unchanged).

Usage:
    python calibration_check.py --maps-dir prob_maps_grid_s3m_nbody \
        --centers prob_maps_grid_dlon+000_lat+00.npz,prob_maps_grid_dlon+020_lat-12.npz
"""
from __future__ import annotations
import argparse, os, sys
from pathlib import Path
import numpy as np, pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
SORCHA = W / "outputs/s3m_linking/case1/sorcha_comparison_case1_nbody_Vband.parquet"
OUT = W / "outputs" / "calibration_check"
# tracklet population label -> map population name
POP_MAP = {"NEO": "NEO", "MBA": "MBA", "TNO": "TNO", "Trojan": "Trojans"}


def check_map(npz_path: Path, trk: pd.DataFrame, night: int | None = None):
    z = np.load(npz_path, allow_pickle=True)
    pops = [str(p) for p in z["population_names"]]
    labels = [str(l) for l in z["mag_bin_labels"]]
    xg, yg = z["x_grid"], z["y_grid"]
    cell = float((xg[1] - xg[0]) * (yg[1] - yg[0]))

    sub = trk if night is None else trk[trk.night == night]
    rows = []
    for lab in labels:
        # ---- PREDICTED: integrate each population's density map ----
        integ = {}
        for p in pops:
            key = f"density_raw__{p}__{lab}"
            integ[p] = float(z[key].sum()) * cell if key in z else 0.0
        tot_pred = sum(integ.values())
        pred_neo = integ.get("NEO", 0.0) / tot_pred if tot_pred > 0 else np.nan

        # ---- OBSERVED: Sorcha tracklets in this map + V-band mag bin ----
        obs = sub[sub.mag_bin_label_Vband == lab]
        n_obs = len(obs)
        obs_counts = {p: int((obs.population.map(POP_MAP) == p).sum()) for p in pops}
        n_class = sum(obs_counts.values())          # excludes 'other'
        obs_neo = obs_counts.get("NEO", 0) / n_class if n_class > 0 else np.nan

        rows.append(dict(
            mag_bin=lab,
            pred_N_NEO=integ.get("NEO", 0.0), pred_N_MBA=integ.get("MBA", 0.0),
            pred_NEO_frac=pred_neo,
            obs_n_tracklets=n_obs, obs_n_classified=n_class,
            obs_n_NEO=obs_counts.get("NEO", 0), obs_n_MBA=obs_counts.get("MBA", 0),
            obs_NEO_frac=obs_neo,
            ratio_pred_over_obs=(pred_neo / obs_neo) if (obs_neo and obs_neo > 0 and np.isfinite(pred_neo)) else np.nan,
        ))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maps-dir", default="prob_maps_grid_s3m_nbody")
    ap.add_argument("--centers", required=True, help="comma-separated npz filenames")
    ap.add_argument("--night", type=int, default=None, help="restrict tracklets to one night")
    ap.add_argument("--tag", default="production")
    a = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    mapdir = W / a.maps_dir
    trk = pd.read_parquet(SORCHA, columns=["population", "prob_map_file", "night", "mag_bin_label_Vband"])
    pd.set_option("display.width", 240)

    all_rows = []
    for cen in a.centers.split(","):
        cen = cen.strip()
        p = mapdir / cen
        if not p.exists():
            print(f"!! missing {p}"); continue
        sub = trk[trk.prob_map_file == cen]
        df = check_map(p, sub, night=a.night)
        df.insert(0, "center", cen)
        all_rows.append(df)
        print(f"\n{'='*100}\n{cen}   (maps: {a.maps_dir}"
              + (f", night {a.night}" if a.night else ", full 2-yr") + ")")
        show = df[["mag_bin", "pred_N_NEO", "pred_N_MBA", "pred_NEO_frac",
                   "obs_n_NEO", "obs_n_MBA", "obs_NEO_frac", "ratio_pred_over_obs"]]
        print(show.to_string(index=False, float_format=lambda v: f"{v:,.4g}"))
        ok = df.dropna(subset=["ratio_pred_over_obs"])
        if len(ok):
            print(f"  --> predicted/observed NEO fraction: median {ok.ratio_pred_over_obs.median():.2f}x  "
                  f"range {ok.ratio_pred_over_obs.min():.2f}-{ok.ratio_pred_over_obs.max():.2f}")
            print(f"      (expect slightly >1 from trailing losses; a smooth, near-constant ratio across "
                  f"bins = correct per-bin normalisation)")
    if all_rows:
        res = pd.concat(all_rows, ignore_index=True)
        suffix = "" if a.night is None else f"_n{a.night}"
        f = OUT / f"calibration_{a.tag}{suffix}.csv"
        res.to_csv(f, index=False)
        print(f"\nwrote {f}")


if __name__ == "__main__":
    main()
