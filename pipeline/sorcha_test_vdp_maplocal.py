#!/usr/bin/env python3
"""VDP scoring sharded by MAP CELL, so every .npz map file is loaded exactly once.

The replaced layout sharded by ROW over a uid-sorted table. Each row-shard therefore spanned
nearly all 667 map cells and reloaded ~667 map files (~100 MB each) to score ~22 rows per cell
-- roughly 53,000 map loads and ~5 TB of redundant I/O across 80 shards, which is why VDP took
16.8 min per shard instead of seconds.

Here the work is partitioned by prob_map_file: a shard owns a disjoint set of map cells, loads
each of its maps once, and scores every row belonging to them. Total map loads across all shards
equals the number of distinct cells (<= 667), not shards x cells.

Model parameters are read from MODEL_SEAL and applied unchanged. Output merges on tracklet_uid.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.special import expit, logit

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
SEAL = json.loads((W/"outputs"/"splits"/"MODEL_SEAL.json").read_text())
A_, B_ = SEAL["calibration"]["a"], SEAL["calibration"]["b"]
D = W/SEAL["map_build"]["maps_dir"]
POPS = ["NEO", "MBA", "TNO", "Trojans"]; EPS = 1e-12
BINS = [("14_16",14.,16.),("16_18",16.,18.),("18_20",18.,20.),("mag20",20.,21.),
        ("mag21",21.,22.),("mag22",22.,23.),("mag23",23.,24.),("mag24+",24.,25.)]
CTR = np.array([(lo+hi)/2 for _,lo,hi in BINS]); LABS = [b for b,_,_ in BINS]


def bil(arr, x, yy, xs, ys):
    ix = np.clip(np.searchsorted(x, xs)-1, 0, len(x)-2); iy = np.clip(np.searchsorted(yy, ys)-1, 0, len(yy)-2)
    tx = np.clip((xs-x[ix])/(x[ix+1]-x[ix]), 0, 1); ty = np.clip((ys-yy[iy])/(yy[iy+1]-yy[iy]), 0, 1)
    return (arr[iy,ix]*(1-tx)*(1-ty)+arr[iy,ix+1]*tx*(1-ty)
            + arr[iy+1,ix]*(1-tx)*ty+arr[iy+1,ix+1]*tx*ty)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--first-tracklet", required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--out-dir", required=True)
    a_ = ap.parse_args()
    OUT = Path(a_.out_dir); OUT.mkdir(parents=True, exist_ok=True)
    dst = OUT/f"vdp_shard_{a_.shard:04d}.parquet"
    if dst.exists():
        print(f"shard {a_.shard} already complete"); return

    cols = ["tracklet_uid", "prob_map_file", "vlam", "vbeta", "mean_mag_V",
            "in_domain_mag245", "in_domain_map", "sorcha_object_linked", "population"]
    te = pd.read_parquet(a_.first_tracklet, columns=cols)
    cells = sorted(te.prob_map_file.unique())
    mine = cells[a_.shard::a_.nshards]                 # DISJOINT map cells per shard
    te = te[te.prob_map_file.isin(mine)].reset_index(drop=True)
    print(f"shard {a_.shard}/{a_.nshards}: {len(mine)} map cells, {len(te):,} rows "
          f"(each map loaded exactly once)", flush=True)

    P_raw = np.full(len(te), np.nan); TOT = np.full(len(te), np.nan)
    t0 = time.time()
    for ci, cen in enumerate(mine):
        g = te.index[te.prob_map_file.to_numpy() == cen].to_numpy()
        if not len(g):
            continue
        f = D/cen
        if not f.exists():
            continue                                    # stays NaN -> explicit invalid
        z = np.load(f, allow_pickle=True); xg, yg = z["x_grid"], z["y_grid"]
        xs = te.vlam.to_numpy(float)[g]; ys = te.vbeta.to_numpy(float)[g]
        mags = te.mean_mag_V.to_numpy(float)[g]
        j = np.clip(np.searchsorted(CTR, mags)-1, 0, len(CTR)-2)
        t = np.clip((mags-CTR[j])/(CTR[j+1]-CTR[j]), 0, 1)
        dens = {p: np.zeros(len(g)) for p in POPS}
        for jj in np.unique(j):
            s_ = j == jj
            for p in POPS:
                d0 = bil(np.nan_to_num(np.asarray(z[f"density_raw__{p}__{LABS[jj]}"], float)), xg, yg, xs[s_], ys[s_])
                d1 = bil(np.nan_to_num(np.asarray(z[f"density_raw__{p}__{LABS[jj+1]}"], float)), xg, yg, xs[s_], ys[s_])
                dens[p][s_] = d0*(1-t[s_])+d1*t[s_]
        tot = sum(dens[p] for p in POPS)
        P_raw[g] = np.where(tot > 0, dens["NEO"]/np.where(tot > 0, tot, 1), np.nan)
        TOT[g] = tot
        del z
        if (ci+1) % 25 == 0:
            print(f"  {ci+1}/{len(mine)} cells  {time.time()-t0:.0f}s", flush=True)

    inb = ((te.vlam.abs() <= 5) & (te.vbeta.abs() <= 5)
           & (te.mean_mag_V >= 14) & (te.mean_mag_V < 25)).to_numpy()
    m = inb & np.isfinite(P_raw) & np.isfinite(TOT) & (TOT > 0)
    P_cal = np.full(len(te), np.nan)
    P_cal[m] = expit(A_*logit(np.clip(P_raw[m], EPS, 1-EPS))+B_)
    out = te[["tracklet_uid"]].copy()
    out["P_raw"], out["P_cal"], out["total_density"], out["vdp_valid"] = P_raw, P_cal, TOT, m
    assert out.tracklet_uid.is_unique
    out.to_parquet(dst, index=False)
    h = hashlib.sha256(pd.util.hash_pandas_object(out, index=False).values.tobytes()).hexdigest()
    print(f"DONE {len(out):,} rows in {(time.time()-t0)/60:.2f} min  valid {int(m.sum()):,}")
    print(f"CONTENT_HASH {h}")


if __name__ == "__main__":
    main()
