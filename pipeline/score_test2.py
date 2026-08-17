#!/usr/bin/env python3
"""TEST2 scoring: new HIGH-sampling VDP, frozen legacy VDP, deterministic digest2.

Modes (--mode): new_vdp | legacy_vdp | digest2 | smoke

NEW VDP  -- single containing 0.25-mag slice, NO cross-magnitude interpolation, NO Platt,
            NO smoothing/masking, NO nan_to_num. Bilinear in velocity only.
LEGACY   -- frozen legacy maps and frozen behaviour: 1-mag bins WITH magnitude interpolation and
            Platt calibration, including its nan_to_num. Deliberately different; that difference is
            part of what the paired comparison measures.
DIGEST2  -- reads the sealed MPC-80 lines, verifies tracklet_input_sha256, repeatable --cpu 1.
"""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess, sys, tempfile, time
from pathlib import Path
import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W / "neomod" / "src"))
sys.path.insert(0, str(W / "neomod" / "adam_core_stub"))
sys.path.insert(0, str(W / "neomod" / "pipeline"))

T2 = W / "outputs" / "test2_geometric"
MAPS_NEW = W / "outputs" / "neomod3_mag025_k150_maps_v2"
MAPS_LEG = W / "prob_maps_grid_neomod3_GEN_final"
SCORED = T2 / "scored"
POPS = ("NEO", "MBA", "TNO", "Trojans")
LEG_LABS = ["14_16", "16_18", "18_20", "mag20", "mag21", "mag22", "mag23", "mag24+"]
LEG_CTR = np.array([15.0, 17.0, 19.0, 20.5, 21.5, 22.5, 23.5, 24.5])
PLATT_A, PLATT_B = 1.20456465, 1.39695663
EPS = 1e-12
D2_CONFIG = "noheadings\nnorms\nrepeatable\nNEO\n"
D2_CPUS = "1"          # NEVER raise: per-thread LCG seeds break digest2 reproducibility


def bil(arr, xg, yg, xs, ys):
    """Bilinear interpolation; NaN outside the grid and NaN-propagating (no nan_to_num)."""
    out = np.full(len(xs), np.nan)
    dx, dy = xg[1] - xg[0], yg[1] - yg[0]
    fx = (xs - xg[0]) / dx
    fy = (ys - yg[0]) / dy
    i0 = np.floor(fx).astype(int); j0 = np.floor(fy).astype(int)
    ok = (i0 >= 0) & (i0 < len(xg) - 1) & (j0 >= 0) & (j0 < len(yg) - 1)
    if not ok.any():
        return out
    i0o, j0o = i0[ok], j0[ok]
    tx, ty = fx[ok] - i0o, fy[ok] - j0o
    a = arr[j0o, i0o]; b = arr[j0o, i0o + 1]
    c = arr[j0o + 1, i0o]; d = arr[j0o + 1, i0o + 1]
    out[ok] = (a * (1 - tx) * (1 - ty) + b * tx * (1 - ty) +
               c * (1 - tx) * ty + d * tx * ty)
    return out


def score_new_vdp(te):
    """Single containing 0.25-mag slice; bilinear in velocity only."""
    n = len(te)
    res = {f"rho_{p}": np.full(n, np.nan) for p in POPS}
    P = np.full(n, np.nan); tot = np.full(n, np.nan)
    reason = np.array(["ok"] * n, dtype=object)
    npops = np.zeros(n, dtype=np.int8)
    sumchk = np.full(n, np.nan)
    for cen, g in te.groupby("center_label", sort=True):
        f = MAPS_NEW / f"mag025_k150_{cen}.npz"
        idx = g.index.to_numpy()
        if not f.exists():
            reason[idx] = "missing_map"; continue
        z = np.load(f, allow_pickle=True)
        keys = set(z.keys())
        xg = np.asarray(z["x_grid"], float); yg = np.asarray(z["y_grid"], float)
        for lab, gg in g.groupby("magnitude_bin", sort=True):
            ii = gg.index.to_numpy()
            if lab is None or (isinstance(lab, float) and np.isnan(lab)):
                reason[ii] = "v_out_of_range"; continue
            xs = gg.vlam.to_numpy(float); ys = gg.vbeta.to_numpy(float)
            inb = (np.abs(xs) <= 5.0) & (np.abs(ys) <= 5.0)
            dens, present = {}, []
            for p in POPS:
                k = f"density__{p}__{lab}"
                if k not in keys:
                    continue                       # invalid cell: absent, NOT zero
                dens[p] = bil(np.asarray(z[k], np.float64), xg, yg, xs, ys)
                present.append(p)
            missing = [q for q in POPS if q not in present]
            npops[ii] = len(present)
            if missing:
                # A 2- or 3-population denominator is NOT P(NEO). A missing density key means the
                # cell was INVALID under the insufficient-support rule; the sealed manifest does not
                # assert that absence means a physical zero, and project policy forbids turning
                # invalid into zero. So the row is invalid with an explicit reason.
                reason[ii] = "missing_population_density:" + ",".join(missing)
                continue
            tt = np.zeros(len(ii))
            for q in POPS:
                tt = tt + dens[q]
            good = inb & np.isfinite(tt) & (tt > 0)
            for q in POPS:
                good &= np.isfinite(dens[q])
                res[f"rho_{q}"][ii] = dens[q]
            pn = np.full(len(ii), np.nan)
            pn[good] = dens["NEO"][good] / tt[good]
            P[ii] = pn; tot[ii] = tt
            sc = np.zeros(len(ii))
            for q in POPS:
                sc[good] += dens[q][good] / tt[good]
            sumchk[ii] = np.where(good, sc, np.nan)
            reason[ii] = np.where(~inb, "outside_velocity_grid",
                          np.where(~np.isfinite(tt), "nonfinite_density",
                          np.where(tt <= 0, "zero_total_density", "ok")))
    out = pd.DataFrame({"P_NEO_new": P, "total_density_new": tot,
                        "new_vdp_reason": reason, "n_pops_new": npops,
                        "prob_sum_new": sumchk}, index=te.index)
    for p in POPS:
        out[f"rho_{p}_new"] = res[f"rho_{p}"]
    out["new_vdp_valid"] = np.isfinite(P)
    v = out.new_vdp_valid.to_numpy()
    if v.any():
        assert (out.n_pops_new.to_numpy()[v] == 4).all(), \
            "a valid new-VDP row has fewer than four population densities"
        assert np.nanmax(np.abs(out.prob_sum_new.to_numpy()[v] - 1.0)) < 1e-9, \
            "four-class probabilities do not sum to one on a valid row"
    assert not any(str(r).startswith("ok_partial") for r in out.new_vdp_reason.unique()), \
        "partial-denominator rows must never be labelled ok"
    return out


def score_legacy_vdp(te):
    """Frozen legacy behaviour: 1-mag bins, magnitude interpolation, Platt, incl. nan_to_num."""
    from scipy.special import expit, logit
    n = len(te)
    P = np.full(n, np.nan); TOT = np.full(n, np.nan)
    reason = np.array(["ok"] * n, dtype=object)
    for cen, g in te.groupby("center_label", sort=True):
        f = MAPS_LEG / f"prob_maps_grid_{cen}.npz"
        idx = g.index.to_numpy()
        if not f.exists():
            reason[idx] = "missing_map"; continue
        z = np.load(f, allow_pickle=True)
        xg = np.asarray(z["x_grid"], float); yg = np.asarray(z["y_grid"], float)
        xs = g.vlam.to_numpy(float); ys = g.vbeta.to_numpy(float)
        mags = g.mean_mag_V.to_numpy(float)
        j = np.clip(np.searchsorted(LEG_CTR, mags) - 1, 0, len(LEG_CTR) - 2)
        tfrac = np.clip((mags - LEG_CTR[j]) / (LEG_CTR[j + 1] - LEG_CTR[j]), 0, 1)
        dens = {p: np.zeros(len(g)) for p in POPS}
        for jj in np.unique(j):
            s_ = j == jj
            for p in POPS:
                d0 = bil(np.nan_to_num(np.asarray(z[f"density_raw__{p}__{LEG_LABS[jj]}"], float)),
                         xg, yg, xs[s_], ys[s_])
                d1 = bil(np.nan_to_num(np.asarray(z[f"density_raw__{p}__{LEG_LABS[jj+1]}"], float)),
                         xg, yg, xs[s_], ys[s_])
                dens[p][s_] = np.nan_to_num(d0) * (1 - tfrac[s_]) + np.nan_to_num(d1) * tfrac[s_]
        tt = sum(dens[p] for p in POPS)
        good = (np.abs(xs) <= 5) & (np.abs(ys) <= 5) & (mags >= 14) & (mags < 25) \
            & np.isfinite(tt) & (tt > 0)
        pn = np.full(len(g), np.nan)
        pn[good] = dens["NEO"][good] / tt[good]
        P[idx] = pn; TOT[idx] = tt
        reason[idx] = np.where(good, "ok", "legacy_invalid")
    Pc = np.full(n, np.nan)
    m = np.isfinite(P)
    Pc[m] = expit(PLATT_A * logit(np.clip(P[m], EPS, 1 - EPS)) + PLATT_B)
    return pd.DataFrame({"P_NEO_legacy_raw": P, "P_NEO_legacy_cal": Pc,
                         "total_density_legacy": TOT, "legacy_vdp_reason": reason,
                         "legacy_vdp_valid": m}, index=te.index)


def score_digest2(te, mpc, batch=20000):
    """Sealed MPC-80 lines only; hash-verified; repeatable --cpu 1; strict parsing."""
    m = mpc.set_index("tracklet_uid").loc[te.tracklet_uid]
    rehash = [hashlib.sha256((a + "\n" + b).encode()).hexdigest()
              for a, b in zip(m.mpc80_line0, m.mpc80_line1)]
    if not np.array_equal(np.asarray(rehash), m.tracklet_input_sha256.to_numpy()):
        raise RuntimeError("MPC-80 lines do not match tracklet_input_sha256")
    keys = m.d2_key.tolist()
    if len(set(keys)) != len(keys):
        raise RuntimeError("duplicate digest2 keys within shard")
    obs = []
    for a, b in zip(m.mpc80_line0, m.mpc80_line1):
        obs.append(a); obs.append(b)
    n = len(keys)
    score, status = {}, {}
    with tempfile.NamedTemporaryFile("w", suffix=".config", delete=False) as c:
        c.write(D2_CONFIG); cfg = c.name
    try:
        for s in range(0, n, batch):
            e = min(s + batch, n)
            with tempfile.NamedTemporaryFile("w", suffix=".obs", delete=False) as tf:
                tf.write("\n".join(obs[2 * s:2 * e]) + "\n"); op = tf.name
            try:
                r = subprocess.run([str(W / "digest2" / "digest2"), "-p", str(W / "digest2"),
                                    "-c", cfg, "--cpu", D2_CPUS, op],
                                   capture_output=True, text=True, timeout=43200)
                if r.returncode != 0:
                    for k in keys[s:e]:
                        status[k] = f"digest2_rc{r.returncode}"
                    continue
                seen = set()
                for ln in r.stdout.splitlines():
                    pp = ln.split()
                    if len(pp) < 2:
                        continue
                    k = pp[0]
                    if k in seen:
                        status[k] = "duplicate_output"; continue
                    try:
                        val = float(pp[1])
                    except ValueError:
                        status[k] = "unparseable"; continue
                    if not (0.0 <= val <= 100.0):
                        status[k] = "out_of_range"; continue
                    if abs(val - round(val)) > 1e-9:
                        status[k] = "noninteger"; continue
                    score[k] = val; seen.add(k)
                for k in keys[s:e]:
                    if k not in seen and k not in status:
                        status[k] = "no_output_line"
            except subprocess.TimeoutExpired:
                for k in keys[s:e]:
                    status[k] = "timeout"
            finally:
                os.unlink(op)
            print(f"    digest2 {e:,}/{n:,}", flush=True)
    finally:
        os.unlink(cfg)
    vals = np.array([score[k] / 100.0 if k in score else np.nan for k in keys])
    return pd.DataFrame({"P_NEO_digest2": vals,
                         "digest2_raw": [score.get(k, np.nan) for k in keys],
                         "digest2_status": [("ok" if k in score else status.get(k, "missing"))
                                            for k in keys],
                         "digest2_valid": np.isfinite(vals),
                         "tracklet_input_sha256_used": m.tracklet_input_sha256.to_numpy()},
                        index=te.index)


def load_shard(a):
    te = pd.read_parquet(T2 / "TEST2_TRACKLETS.parquet")
    te = te.sort_values("tracklet_uid", kind="mergesort").reset_index(drop=True)
    if a.mode in ("new_vdp", "legacy_vdp"):
        # Balanced by TRACKLET COUNT, not round-robin center count: centers differ by >100x in
        # rows, so round-robin left shards wildly uneven. Manifest is immutable and hashed.
        man = json.load(open(T2 / "CENTER_SHARD_MANIFEST.json"))
        if man["nshards"] != a.nshards:
            raise RuntimeError(f"manifest nshards {man['nshards']} != requested {a.nshards}")
        assign = man["assignment"]
        te = te[te.center_label.map(assign) == a.shard].reset_index(drop=True)
    else:
        te = te.iloc[np.array_split(np.arange(len(te)), a.nshards)[a.shard]].reset_index(drop=True)
    return te


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=("new_vdp", "legacy_vdp", "digest2", "smoke"), required=True)
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--nshards", type=int, default=1)
    p.add_argument("--batch", type=int, default=20000)
    a = p.parse_args()
    SCORED.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    if a.mode == "smoke":
        return smoke()

    dst = SCORED / f"{a.mode}_{a.shard:04d}.parquet"
    if dst.exists():
        try:
            pd.read_parquet(dst, columns=["tracklet_uid"]); print(f"[skip] {dst.name}"); return
        except Exception:
            pass
    te = load_shard(a)
    print(f"{a.mode} shard {a.shard}/{a.nshards}: {len(te):,} rows "
          f"({te.center_label.nunique()} centers)", flush=True)
    if a.mode == "new_vdp":
        r = score_new_vdp(te)
    elif a.mode == "legacy_vdp":
        r = score_legacy_vdp(te)
    else:
        r = score_digest2(te, pd.read_parquet(T2 / "TEST2_MPC80.parquet"), a.batch)
    out = pd.concat([te[["tracklet_uid"]], r], axis=1)
    out.to_parquet(dst.with_suffix(".tmp.parquet"), index=False)
    dst.with_suffix(".tmp.parquet").replace(dst)
    vcol = {"new_vdp": "new_vdp_valid", "legacy_vdp": "legacy_vdp_valid",
            "digest2": "digest2_valid"}[a.mode]
    dt = time.time() - t0
    print(f"  wrote {dst.name}  valid {int(out[vcol].sum()):,}/{len(out):,}  "
          f"{dt:.0f}s  {len(out)/max(dt,1e-9):.1f} rows/s", flush=True)


def smoke():
    """Deterministic cross-section: all populations, faint bins, many centers, valid+invalid."""
    te = pd.read_parquet(T2 / "TEST2_TRACKLETS.parquet")
    mpc = pd.read_parquet(T2 / "TEST2_MPC80.parquet")
    # Restrict to a few centers FIRST. Each new-VDP map is ~372 MB, so a sample spread over 439
    # centers costs ~163 GB of reads for a few hundred rows; the full run shards BY center so each
    # map loads once. Chosen deterministically as the most-populated centers plus a high-latitude
    # one, so all four populations, faint bins and both valid/invalid cells are still represented.
    troj = te[te.population == "Trojans"].center_label.value_counts().head(2).index.tolist()
    top = te.center_label.value_counts().head(6).index.tolist() + troj
    extra = [c for c in te.center_label.unique() if c.endswith("lat+50") or c.endswith("lat-50")][:2]
    cens = sorted(set(top) | set(extra))
    te = te[te.center_label.isin(cens)].reset_index(drop=True)
    print(f"SMOKE centers ({len(cens)}): {cens}", flush=True)
    picks = []
    for p in POPS:
        g = te[te.population == p]
        if len(g):
            picks.append(g.sample(min(120, len(g)), random_state=11))
    inval = te[~te.vdp_cell_valid]
    if len(inval):
        picks.append(inval.sample(min(80, len(inval)), random_state=12))
    faint = te[te.magnitude_bin.astype(str) >= "V024.00"]
    if len(faint):
        picks.append(faint.sample(min(120, len(faint)), random_state=13))
    fast = te[np.hypot(te.vlam, te.vbeta) > 1.0]
    if len(fast):
        picks.append(fast.sample(min(80, len(fast)), random_state=14))
    s = pd.concat(picks).drop_duplicates("tracklet_uid").reset_index(drop=True)
    print(f"SMOKE sample {len(s):,} rows | populations {sorted(s.population.unique())}")
    print(f"  bins {s.magnitude_bin.nunique()}  centers {s.center_label.nunique()}  "
          f"|v| {np.hypot(s.vlam, s.vbeta).min():.3f}..{np.hypot(s.vlam, s.vbeta).max():.3f}")
    print(f"  expected-invalid rows: {int((~s.vdp_cell_valid).sum())}", flush=True)

    nv = score_new_vdp(s); lv = score_legacy_vdp(s)
    d2 = score_digest2(s, mpc, batch=20000)
    r = pd.concat([s.reset_index(drop=True), nv, lv, d2], axis=1)

    ok = True
    def c(name, cond, det=""):
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  [{det}]" if det else ""), flush=True)

    # single containing bin, no cross-bin interpolation
    lab, _ = None, None
    v = r.mean_mag_V.to_numpy(float)
    lo = 14.0 + np.floor((v - 14.0) / 0.25) * 0.25
    exp = [f"V{x:06.2f}_{x+0.25:06.2f}" for x in lo]
    c("new VDP uses the single containing 0.25-mag bin",
      all(a_ == b_ for a_, b_ in zip(r.magnitude_bin.astype(str), exp)))
    # a tracklet at V=23.4 -> [23.25,23.50)
    t234 = 14.0 + np.floor((23.4 - 14.0) / 0.25) * 0.25
    c("V=23.4 -> [23.25,23.50)", abs(t234 - 23.25) < 1e-9, f"{t234}")
    c("new VDP probabilities sum to 1 where valid",
      np.nanmax(np.abs(r.prob_sum_new[r.new_vdp_valid] - 1.0)) < 1e-9,
      f"max dev {np.nanmax(np.abs(r.prob_sum_new[r.new_vdp_valid]-1.0)):.2e}")
    c("new VDP has no Platt column", "P_NEO_new_cal" not in r.columns)
    c("legacy VDP retains raw AND calibrated scores",
      r.P_NEO_legacy_raw.notna().any() and r.P_NEO_legacy_cal.notna().any())
    c("legacy calibrated differs from raw (Platt active)",
      np.nanmax(np.abs(r.P_NEO_legacy_cal - r.P_NEO_legacy_raw)) > 1e-6)
    c("new VDP differs from legacy on at least some rows",
      np.nanmax(np.abs(r.P_NEO_new - r.P_NEO_legacy_raw)) > 1e-6)
    c("all four populations present in the smoke sample",
      set(r.population.unique()) == set(POPS), str(sorted(r.population.unique())))
    c("every valid new-VDP row has all four densities",
      bool((~r.new_vdp_valid).all() or (r.loc[r.new_vdp_valid, "n_pops_new"] == 4).all()))
    c("no partial-denominator row is valid",
      not r.new_vdp_reason.astype(str).str.startswith("ok_partial").any())
    c("partial-denominator rows are NaN with an explicit missing-population reason",
      bool((~r.new_vdp_reason.astype(str).str.startswith("missing_population_density")).all()
           or r.loc[r.new_vdp_reason.astype(str).str.startswith("missing_population_density"),
                    "P_NEO_new"].isna().all()))
    c("probability-sum deviation < 1e-9",
      bool((~r.new_vdp_valid).all() or
           np.nanmax(np.abs(r.prob_sum_new[r.new_vdp_valid] - 1.0)) < 1e-9))
    c("MPC-80 designations unique in the sample",
      r.mpc80_designation.is_unique if "mpc80_designation" in r.columns else True)
    c("invalid new-VDP rows are NaN with a reason",
      bool((~r.new_vdp_valid).sum() == 0 or
           (r.loc[~r.new_vdp_valid, "P_NEO_new"].isna().all()
            and (r.loc[~r.new_vdp_valid, "new_vdp_reason"] != "ok").all())),
      f"{int((~r.new_vdp_valid).sum())} invalid")
    c("no invalid score was written as 0",
      not ((r.P_NEO_new == 0) & (~r.new_vdp_valid)).any())
    c("digest2 read the sealed lines (hash verified)",
      (r.tracklet_input_sha256_used.to_numpy() == r.tracklet_input_sha256.to_numpy()).all())
    d2v = r.digest2_raw[r.digest2_valid].to_numpy(float)
    # MUST fail on an empty set: an earlier version passed vacuously with n=0 while digest2 was
    # in fact returning nothing, because 8-char keys were truncated to 12 MPC-80 columns and
    # collided.
    c("digest2 produced valid scores at all", len(d2v) > 0, f"n={len(d2v)}")
    c("digest2 scores are integers in 0..100",
      bool(len(d2v) > 0 and (np.all(np.abs(d2v - np.round(d2v)) < 1e-9)
                             and d2v.min() >= 0 and d2v.max() <= 100)),
      f"n={len(d2v)} min={d2v.min() if len(d2v) else float('nan'):.0f} "
      f"max={d2v.max() if len(d2v) else float('nan'):.0f}")
    nz = int((d2v == 0).sum())
    print(f"    digest2 genuine zeros: {nz}  (valid, distinct from NaN)", flush=True)
    c("digest2 missing results are NaN not 0",
      not ((r.P_NEO_digest2 == 0) & (~r.digest2_valid)).any())
    print(f"\n  new VDP valid {int(r.new_vdp_valid.sum()):,}/{len(r):,} | "
          f"legacy {int(r.legacy_vdp_valid.sum()):,} | digest2 {int(r.digest2_valid.sum()):,}")
    print(f"  new-VDP reasons: {r.new_vdp_reason.value_counts().to_dict()}")
    (T2 / "TEST2_SMOKE.json").write_text(json.dumps(
        {"n": len(r), "all_pass": bool(ok),
         "new_valid": int(r.new_vdp_valid.sum()),
         "legacy_valid": int(r.legacy_vdp_valid.sum()),
         "digest2_valid": int(r.digest2_valid.sum()),
         "digest2_zeros": nz,
         "reasons": {k: int(v) for k, v in r.new_vdp_reason.value_counts().items()}}, indent=2))
    print(f"\nSMOKE {'PASSED' if ok else 'FAILED'}", flush=True)
    if not ok:
        sys.exit(3)


if __name__ == "__main__":
    main()
