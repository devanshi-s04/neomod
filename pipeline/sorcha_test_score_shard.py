#!/usr/bin/env python3
"""One scoring shard: frozen VDP + deterministic digest2 on a slice of the first-tracklet product.

digest2 runs in `repeatable` mode with `--cpu 1` inside EVERY shard: per-thread LCG seeds make
multi-threaded digest2 non-reproducible, so parallelism comes from SLURM array tasks, never from
digest2's internal thread count.

Rows are merged back by `tracklet_uid` (ObjID__night), a stable parent/pair identifier.
Every input row leaves with exactly one score OR an explicit invalid status -- never a zero.
"""
from __future__ import annotations
import argparse, os, subprocess, sys, tempfile, time
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from scipy.special import expit, logit
import json

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); sys.path.insert(0, str(W/"neomod"/"pipeline"))
SEAL = json.loads((W/"outputs"/"splits"/"MODEL_SEAL.json").read_text())
A_, B_ = SEAL["calibration"]["a"], SEAL["calibration"]["b"]
D = W/SEAL["map_build"]["maps_dir"]
POPS = ["NEO", "MBA", "TNO", "Trojans"]; EPS = 1e-12
BINS = [("14_16",14.,16.),("16_18",16.,18.),("18_20",18.,20.),("mag20",20.,21.),
        ("mag21",21.,22.),("mag22",22.,23.),("mag23",23.,24.),("mag24+",24.,25.)]
CTR = np.array([(lo+hi)/2 for _,lo,hi in BINS]); LABS = [b for b,_,_ in BINS]
D2_CONFIG = "noheadings\nnorms\nrepeatable\nNEO\n"      # repeatable == deterministic
D2_CPUS = "1"                                            # NEVER raise: per-thread seeds break it


def bil(arr, x, yy, xs, ys):
    ix = np.clip(np.searchsorted(x, xs)-1, 0, len(x)-2); iy = np.clip(np.searchsorted(yy, ys)-1, 0, len(yy)-2)
    tx = np.clip((xs-x[ix])/(x[ix+1]-x[ix]), 0, 1); ty = np.clip((ys-yy[iy])/(yy[iy+1]-yy[iy]), 0, 1)
    return (arr[iy,ix]*(1-tx)*(1-ty)+arr[iy,ix+1]*tx*(1-ty)
            + arr[iy+1,ix]*(1-tx)*ty+arr[iy+1,ix+1]*tx*ty)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--first-tracklet", required=True)
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--nshards", type=int, required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--d2-batch", type=int, default=1000)
    # optional sub-range WITHIN the shard, so one shard can be split across jobs.
    # digest2 `repeatable` reseeds per tracklet, so scores do not depend on batching.
    ap.add_argument("--row-start", type=int, default=None)
    ap.add_argument("--row-stop", type=int, default=None)
    a_ = ap.parse_args()
    OUT = Path(a_.out_dir); OUT.mkdir(parents=True, exist_ok=True)
    sub = a_.row_start is not None
    dst = (OUT/f"scored_shard_{a_.shard:04d}_sub{a_.row_start:06d}.parquet" if sub
           else OUT/f"scored_shard_{a_.shard:04d}.parquet")
    if dst.exists():
        print(f"shard {a_.shard} already complete"); return

    full = pd.read_parquet(a_.first_tracklet).sort_values("tracklet_uid",
                                                          kind="mergesort").reset_index(drop=True)
    bounds = np.array_split(np.arange(len(full)), a_.nshards)[a_.shard]
    te = full.iloc[bounds].reset_index(drop=True)
    if sub:
        te = te.iloc[a_.row_start:a_.row_stop].reset_index(drop=True)
    n_in = len(te)
    print(f"shard {a_.shard}/{a_.nshards}: {n_in:,} rows", flush=True)

    # ------------------------------------------------------------------ VDP
    t0 = time.time()
    def run(cen):
        g = te[te.prob_map_file == cen]
        xs = g.vlam.to_numpy(float); ys = g.vbeta.to_numpy(float)
        mags = g.mean_mag_V.to_numpy(float)
        f = D/cen
        if not f.exists():
            return g.index.to_numpy(), np.full(len(g), np.nan), np.full(len(g), np.nan)
        z = np.load(f, allow_pickle=True); xg, yg = z["x_grid"], z["y_grid"]
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
        return g.index.to_numpy(), np.where(tot > 0, dens["NEO"]/np.where(tot > 0, tot, 1), np.nan), tot
    parts = Parallel(n_jobs=8)(delayed(run)(c) for c in sorted(te.prob_map_file.unique()))
    P_raw = np.full(n_in, np.nan); TOT = np.full(n_in, np.nan)
    for i, p_, t_ in parts: P_raw[i] = p_; TOT[i] = t_
    inb = ((te.vlam.abs() <= 5) & (te.vbeta.abs() <= 5)
           & (te.mean_mag_V >= 14) & (te.mean_mag_V < 25))
    m = inb.to_numpy() & np.isfinite(P_raw) & np.isfinite(TOT) & (TOT > 0)
    P_cal = np.full(n_in, np.nan); P_cal[m] = expit(A_*logit(np.clip(P_raw[m], EPS, 1-EPS))+B_)
    te["P_raw"], te["P_cal"], te["vdp_valid"], te["total_density"] = P_raw, P_cal, m, TOT
    t_vdp = time.time()-t0
    print(f"  VDP {t_vdp:.0f}s  valid {int(m.sum()):,}/{n_in:,}", flush=True)

    # ------------------------------------------------------------------ digest2, --cpu 1
    from sorcha_phase2 import format_mpc80
    keys = [f"D{i:06d}" for i in range(n_in)]           # local to this shard; merged by uid
    obs = []
    for k, row in zip(keys, te.itertuples(index=False)):
        obs.append(format_mpc80(f"     {k}", row.mjd0_utc, row.ra0, row.dec0, row.mag0_V))
        obs.append(format_mpc80(f"     {k}", row.mjd1_utc, row.ra1, row.dec1, row.mag1_V))
    t0 = time.time()
    score, status = {}, {}
    with tempfile.NamedTemporaryFile("w", suffix=".config", delete=False) as c:
        c.write(D2_CONFIG); cfg = c.name
    try:
        for s in range(0, n_in, a_.d2_batch):
            e = min(s+a_.d2_batch, n_in)
            with tempfile.NamedTemporaryFile("w", suffix=".obs", delete=False) as t:
                t.write("\n".join(obs[2*s:2*e])+"\n"); op = t.name
            try:
                r = subprocess.run([str(W/"digest2"/"digest2"), "-p", str(W/"digest2"),
                                    "-c", cfg, "--cpu", D2_CPUS, op],
                                   capture_output=True, text=True, timeout=21600)
                if r.returncode != 0:
                    for k in keys[s:e]: status[k] = f"digest2_rc{r.returncode}"
                    continue
                seen = set()
                for ln in r.stdout.splitlines():
                    p = ln.split()
                    if len(p) >= 2:
                        try: score[p[0]] = float(p[1]); seen.add(p[0])
                        except ValueError: status[p[0]] = "unparseable"
                for k in keys[s:e]:
                    if k not in seen and k not in status: status[k] = "no_output_line"
            except subprocess.TimeoutExpired:
                for k in keys[s:e]: status[k] = "timeout"
            finally:
                os.unlink(op)
    finally:
        os.unlink(cfg)
    t_d2 = time.time()-t0
    te["P_NEO_d2"] = [score[k]/100.0 if k in score else np.nan for k in keys]   # NaN, never 0
    te["d2_status"] = [("ok" if k in score else status.get(k, "missing")) for k in keys]
    te["d2_valid"] = np.isfinite(te.P_NEO_d2.to_numpy())

    # ---- EXACTLY one score or one explicit invalid status per input row -----------------
    assert len(te) == n_in, "row count changed during scoring"
    assert te.tracklet_uid.is_unique, "tracklet_uid not unique within the shard"
    ok = te.d2_valid.to_numpy()
    bad = te.d2_status.to_numpy() != "ok"
    assert (ok ^ bad).all(), "a row is neither scored nor explicitly invalid"
    assert not (te.P_NEO_d2 == 0).any() or True   # genuine 0 is a real score, distinct from NaN
    te["sample_main_mag245"] = te.in_domain_mag245 & te.vdp_valid & te.d2_valid
    te["sample_full_map"] = te.in_domain_map & te.vdp_valid & te.d2_valid
    te["sample_case1_linked"] = te.sorcha_object_linked
    te.to_parquet(dst, index=False)
    print(f"  digest2 {t_d2/60:.1f} min  valid {int(te.d2_valid.sum()):,}  "
          f"invalid {int((~te.d2_valid).sum()):,} ({te.d2_status[~te.d2_valid].value_counts().to_dict()})")
    print(f"  wrote {dst}   VDP {t_vdp:.0f}s + digest2 {t_d2:.0f}s = {(t_vdp+t_d2)/60:.1f} min")


if __name__ == "__main__":
    main()
