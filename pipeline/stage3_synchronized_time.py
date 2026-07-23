#!/usr/bin/env python3
"""STAGE 3 — synchronized-time benchmark/Sorcha test (explanatory).

Explains the benchmark-vs-Sorcha digest2/VDP score differences by removing the observation-time
mismatch: propagate each matched NEO's n-body orbit to Sorcha's EXACT t0/t1 and re-score. This does
NOT touch the Stage-2 conclusion (digest2 > VDP on identical inputs); it decomposes the benchmark↔Sorcha
gap into time / measurement / band / residual-position pieces.

Four controlled cases per object (794-NEO pilot; ObjID==s3m_objid, Sorcha night 61642):
  A  benchmark original   — benchmark tracklet at its fixed epoch (MJD 61642.0, +30 min), native V
  B  benchmark synced     — SAME n-body orbit propagated (ASSIST) to Sorcha's t0/t1, V mags
  C  Sorcha raw           — Sorcha's detected tracklet, raw LSST mags
  D  Sorcha V-corrected   — Sorcha's detected tracklet, V-corrected mags
All scored with deterministic digest2 (repeatable, --cpu 1) and VDP (n-body maps).

Decomposition (both classifiers): B-A time/epoch; C-B & D-B Sorcha measurement/position; D-C band;
D-B remaining benchmark-vs-Sorcha (band-matched, both V).

Commands
  propagate   Case-B n-body propagation to Sorcha t0/t1 (heavy) -> caseB.parquet
  analyze     assemble A/B/C/D, score both classifiers, decompose, figures, tables

Outputs: outputs/mag245_nbody_synchronized_time_test/  (production + Stage 1/2 outputs untouched)
Unmappable / unpropagatable objects are reported explicitly, never dropped or zeroed.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "src"))
import audit_digest2 as ad          # noqa: E402
import rescore_vdp_Vband as rv       # noqa: E402

# audit_digest2 sets iers.conf.auto_download=False; the obs are in 2027 (future), so astropy would
# refuse to interpolate the predictive IERS table on age grounds. The cached table DOES cover 2027
# (Case-B smoke test confirmed), so lift the age gate to allow predictive extrapolation offline.
from astropy.utils import iers  # noqa: E402
iers.conf.auto_max_age = None

WORK = ad.WORKDIR
OUT_DIR = WORK / "outputs" / "mag245_nbody_synchronized_time_test"
MAPDIR = WORK / "prob_maps_grid_s3m_nbody"
S1_DIR = WORK / "outputs" / "mag245_nbody_deterministic_rescore"
SORCHA_P = WORK / "outputs/s3m_linking/case1/sorcha_comparison_case1_nbody_Vband.parquet"
BENCH_P = WORK / "outputs/phase2_benchmark_s3m_nbody_mag245/benchmark_comparison_s3m_nbody_mag245.parquet"
S0_S3M = WORK / "neomod" / "S3Mdata" / "S0.s3m"
CASEB = OUT_DIR / "caseB_propagated.parquet"
NIGHT = 61642


def _hash(l0, l1):
    return hashlib.sha256((l0 + "\n" + l1 + "\n").encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------- base 794 table
def load_base() -> pd.DataFrame:
    """The 794 matched NEOs with Sorcha obs (C/D) + benchmark obs (A). Verifies unique cache mapping."""
    s = pd.read_parquet(SORCHA_P, columns=[
        "ObjID", "night", "population", "ra0", "dec0", "mjd0_utc", "mag0", "filter0",
        "ra1", "dec1", "mjd1_utc", "mag1", "filter1", "mean_mag", "mean_mag_V",
        "vlam", "vbeta", "prob_map_file"])
    # per-detection Sorcha V-mags (Stage-1 sidecar, row-aligned to the FULL sorcha parquet) —
    # attach BEFORE any filtering so they travel with the rows.
    vm = pd.read_parquet(S1_DIR / "sorcha_vmags.parquet")
    assert len(vm) == len(s), f"vmags {len(vm)} != sorcha {len(s)}"
    s["mag0_V"] = vm.mag0_V.to_numpy(); s["mag1_V"] = vm.mag1_V.to_numpy()
    s = s[s.night == NIGHT].copy(); s["ObjID"] = s.ObjID.astype(str)
    b = pd.read_parquet(BENCH_P, columns=[
        "s3m_objid", "population", "ra0", "dec0", "mjd0_utc", "mag0", "ra1", "dec1", "mjd1_utc",
        "mag1", "mean_mag", "vlam", "vbeta", "prob_map_file", "H"])
    b["s3m_objid"] = b.s3m_objid.astype(str)
    b = b[b.population == "NEO"].drop_duplicates("s3m_objid").set_index("s3m_objid")
    m = s[s.ObjID.isin(b.index)].reset_index(drop=True)      # benchmark-NEO matched
    bm = b.reindex(m.ObjID.to_numpy())
    for c in ["ra0", "dec0", "mjd0_utc", "mag0", "ra1", "dec1", "mjd1_utc", "mag1", "mean_mag",
              "vlam", "vbeta", "prob_map_file", "H"]:
        m[f"b_{c}"] = bm[c].to_numpy()
    return m


# ---------------------------------------------------------------------------- propagate Case B
def propagate(args) -> None:
    from astropy.time import Time
    import neoscore as nsc
    import s3m_loader
    import velocity_density_pipeline_gmm as vdp

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    base = load_base()
    ids = base.ObjID.tolist()
    print(f"Case B: propagating {len(ids)} NEOs to Sorcha t0/t1 (ASSIST n-body)", flush=True)

    src = s3m_loader.define_s3m(pattern=str(S0_S3M), verbose=False)
    src["OID"] = src["OID"].astype(str).str.strip(); src["a"] = src["q"] / (1.0 - src["e"])
    src = src.set_index("OID")
    missing = [o for o in ids if o not in src.index]
    if missing:
        print(f"  !! {len(missing)} objects missing S3M elements: {missing[:10]}", flush=True)
    scorer = nsc.NEOMODScorer(None, None, None, None, None)

    def predict(obj_row, mjd):
        t = Time(float(mjd), format="mjd", scale="utc")
        v = vdp.build_visible_subset_dataframe(obj_row, obstime_str=t, scorer=scorer, max_sep_deg=180.0,
            chunk=1, show_progress=False, center_mode="custom_ecliptic",
            center_lon_deg=0.0, center_lat_deg=0.0)
        return v.iloc[0] if len(v) else None

    rows = []
    t0 = time.time()
    for k, (_, r) in enumerate(base.iterrows()):
        if r.ObjID in missing:
            rows.append(dict(ObjID=r.ObjID, prop_status="no_s3m_elements")); continue
        obj = src.loc[[r.ObjID]].reset_index().rename(columns={"index": "OID"})
        try:
            p0 = predict(obj, r.mjd0_utc); p1 = predict(obj, r.mjd1_utc)
            if p0 is None or p1 is None:
                rows.append(dict(ObjID=r.ObjID, prop_status="empty_visible")); continue
            rows.append(dict(
                ObjID=r.ObjID, prop_status="ok",
                B_ra0=p0.ra_deg, B_dec0=p0.dec_deg, B_ra1=p1.ra_deg, B_dec1=p1.dec_deg,
                B_mag0=p0.mag_app, B_mag1=p1.mag_app,
                B_vlam=p0.vlam, B_vbeta=p0.vbeta))
        except Exception as e:  # explicit, not dropped
            rows.append(dict(ObjID=r.ObjID, prop_status=f"error:{type(e).__name__}:{str(e)[:60]}"))
        if (k + 1) % 100 == 0:
            print(f"  {k+1}/{len(base)}  ({time.time()-t0:.0f}s)", flush=True)
    out = pd.DataFrame(rows)
    out.to_parquet(CASEB, index=False)
    vc = out.prop_status.value_counts().to_dict()
    print(f"wrote {CASEB}  status={vc}  ({time.time()-t0:.0f}s)", flush=True)


# ---------------------------------------------------------------------------- digest2 (all cases)
# digest2's statistical ranging can INFINITE-LOOP on a degenerate short-arc tracklet. We therefore
# run with a per-call timeout and BISECT on timeout to isolate the exact hanging tracklet(s), which
# are marked `digest2_hang` (NaN — reported explicitly, never zeroed). `repeatable` scores are
# cpu-count-independent (verified cpu1==cpu4), so we use several cpus: identical scores, and the
# timeout cleanly distinguishes a hang from mere slowness.
D2_CPU = os.environ.get("STAGE3_D2_CPU", "8")
D2_TIMEOUT = int(os.environ.get("STAGE3_D2_TIMEOUT", "90"))   # >> a BATCH's completion time
D2_BATCH = int(os.environ.get("STAGE3_D2_BATCH", "200"))      # ~11s at cpu8, never false-times-out


def _run_d2(obs_lines):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".obs", delete=False) as o:
        o.write("\n".join(obs_lines) + "\n"); op = o.name
    with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as c:
        c.write(ad.D2_CONFIG); cp = c.name          # repeatable
    try:
        pr = subprocess.run([str(ad.DIGEST2_EXEC), "-p", str(ad.DIGEST2_DIR), "-c", cp,
                             "--cpu", D2_CPU, op], capture_output=True, text=True, timeout=D2_TIMEOUT)
        parsed = {}
        for ln in pr.stdout.splitlines():
            q = ln.strip().split()
            if len(q) < 2:
                continue
            try:
                parsed.setdefault(q[0], []).append(int(q[1]))
            except ValueError:
                parsed.setdefault(q[0], [])
        return parsed, False
    except subprocess.TimeoutExpired:
        return {}, True
    finally:
        os.unlink(op); os.unlink(cp)


def _score_batch(triples):
    """Recursive bisection to isolate hanging tracklets within a batch that timed out."""
    obs = []
    for _, l0, l1 in triples:
        obs += [l0, l1]
    parsed, timed_out = _run_d2(obs)
    if not timed_out:
        return {d: _strict(parsed.get(d)) for d, _, _ in triples}
    if len(triples) == 1:
        print(f"    digest2 HANG on {triples[0][0]}", flush=True)
        return {triples[0][0]: (np.nan, "digest2_hang")}
    mid = len(triples) // 2
    r = _score_batch(triples[:mid]); r.update(_score_batch(triples[mid:]))
    return r


def _score_digest2(triples):
    """Score in fixed D2_BATCH-sized batches (each completes well under D2_TIMEOUT, so only a real
    infinite-loop tracklet times out); bisect only the batch that hangs. {desig: (norm, status)}."""
    out = {}
    for j in range(0, len(triples), D2_BATCH):
        out.update(_score_batch(triples[j:j + D2_BATCH]))
        if (j // D2_BATCH) % 4 == 0:
            print(f"    digest2 {min(j + D2_BATCH, len(triples))}/{len(triples)}", flush=True)
    return out


def _strict(vals):
    if vals is None:
        return np.nan, "missing"
    if len(vals) == 0:
        return np.nan, "malformed"
    if len(vals) > 1:
        return np.nan, "duplicate"
    raw = vals[0]
    if not (0 <= raw <= 100):
        return np.nan, "out_of_range"
    return raw / 100.0, "ok"


# ---------------------------------------------------------------------------- analyze
def analyze(args) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.time import Time
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    base = load_base()
    B = pd.read_parquet(CASEB)
    base = base.merge(B, on="ObjID", how="left")
    okB = (base.prop_status == "ok")
    print(f"=== STAGE 3 pilot: {len(base)} NEOs; Case-B propagated ok: {int(okB.sum())} "
          f"({base.prop_status.value_counts().to_dict()}) ===", flush=True)

    # ---- assemble the four cases as one long table (one row per ObjID x case) ----
    def case_rows(letter, t0, t1, ra0, dec0, ra1, dec1, mag0, mag1, f0, f1, vlam, vbeta, pmf):
        d = pd.DataFrame(dict(ObjID=base.ObjID, population="NEO", case=letter,
            t0=t0, t1=t1, RA0=ra0, Dec0=dec0, RA1=ra1, Dec1=dec1,
            mag0=mag0, mag1=mag1, filter0=f0, filter1=f1, vlam=vlam, vbeta=vbeta, prob_map_file=pmf))
        return d
    cases = pd.concat([
        case_rows("A", base.b_mjd0_utc, base.b_mjd1_utc, base.b_ra0, base.b_dec0, base.b_ra1, base.b_dec1,
                  base.b_mag0, base.b_mag1, "V", "V", base.b_vlam, base.b_vbeta, base.b_prob_map_file),
        case_rows("B", base.mjd0_utc, base.mjd1_utc, base.B_ra0, base.B_dec0, base.B_ra1, base.B_dec1,
                  base.B_mag0, base.B_mag1, "V", "V", base.B_vlam, base.B_vbeta, base.b_prob_map_file),
        case_rows("C", base.mjd0_utc, base.mjd1_utc, base.ra0, base.dec0, base.ra1, base.dec1,
                  base.mag0, base.mag1, base.filter0, base.filter1, base.vlam, base.vbeta, base.prob_map_file),
        case_rows("D", base.mjd0_utc, base.mjd1_utc, base.ra0, base.dec0, base.ra1, base.dec1,
                  base.mag0_V, base.mag1_V, "V", "V", base.vlam, base.vbeta, base.prob_map_file),
    ], ignore_index=True)
    cases["dt_days"] = cases.t1 - cases.t0

    # ---- MPC lines + hashes (vectorised Time) ----
    valid = cases[["t0", "t1", "RA0", "Dec0", "RA1", "Dec1", "mag0", "mag1"]].notna().all(axis=1)
    dt0 = np.full(len(cases), None, object); dt1 = np.full(len(cases), None, object)
    vi = valid.to_numpy()
    dt0[vi] = Time(cases.t0.to_numpy(float)[vi], format="mjd", scale="utc").utc.datetime
    dt1[vi] = Time(cases.t1.to_numpy(float)[vi], format="mjd", scale="utc").utc.datetime
    l0s, l1s, hs, desigs = [], [], [], []
    ci = {c: i for i, c in enumerate(["A", "B", "C", "D"])}
    for i, r in enumerate(cases.itertuples(index=False)):
        if not vi[i]:
            l0s.append(None); l1s.append(None); hs.append(None); desigs.append(None); continue
        key = f"{r.case}{i % 100000:05d}"; d12 = f"     {key}"
        # inline format using precomputed dt (audit_digest2.format_mpc80 recomputes Time; replicate)
        l0 = _fmt(d12, dt0[i], r.RA0, r.Dec0, r.mag0)
        l1 = _fmt(d12, dt1[i], r.RA1, r.Dec1, r.mag1)
        l0s.append(l0); l1s.append(l1); hs.append(_hash(l0, l1)); desigs.append(key)
    cases["desig"] = desigs; cases["mpc_line0"] = l0s; cases["mpc_line1"] = l1s; cases["sha256_input"] = hs

    # ---- deterministic digest2 for every valid case-tracklet (bisecting on hangs) ----
    valid_idx = [i for i in range(len(cases)) if vi[i]]
    (OUT_DIR / "raw").mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "raw" / "stage3_all.obs").write_text(
        "\n".join(sum(([cases.mpc_line0.iloc[i], cases.mpc_line1.iloc[i]] for i in valid_idx), [])) + "\n")
    triples = [(cases.desig.iloc[i], cases.mpc_line0.iloc[i], cases.mpc_line1.iloc[i]) for i in valid_idx]
    print(f"  scoring {len(triples)} tracklets with deterministic digest2 (cpu={D2_CPU}, "
          f"timeout={D2_TIMEOUT}s, bisecting on hangs) ...", flush=True)
    scored = _score_digest2(triples)
    d2 = np.full(len(cases), np.nan); d2st = np.array(["no_tracklet"] * len(cases), object)
    for i in valid_idx:
        sc, st = scored.get(cases.desig.iloc[i], (np.nan, "missing"))
        d2[i] = sc; d2st[i] = st
    cases["digest2_det"] = d2; cases["digest2_status"] = d2st
    print(f"  digest2 status={pd.Series(d2st).value_counts().to_dict()}", flush=True)

    # ---- VDP (n-body maps) for every case (uses vlam/vbeta + mean mag + prob_map_file) ----
    cases["mean_mag_case"] = 0.5 * (cases.mag0 + cases.mag1)
    vdp_scores = np.full(len(cases), np.nan)
    scv = cases[cases[["vlam", "vbeta", "prob_map_file", "mean_mag_case"]].notna().all(axis=1)].copy()
    P = rv.rescore(scv.rename(columns={}), scv.mean_mag_case.to_numpy(float), mapdir=str(MAPDIR))
    vdp_scores[scv.index.to_numpy()] = P
    cases["vdp"] = vdp_scores

    # ---- sky separation of benchmark cases vs Sorcha detection (Case C positions) ----
    Cpos = cases[cases.case == "C"].set_index("ObjID")
    def sep_to_C(sub):
        c = Cpos.reindex(sub.ObjID.to_numpy())
        s0 = SkyCoord(sub.RA0.to_numpy() * u.deg, sub.Dec0.to_numpy() * u.deg).separation(
            SkyCoord(c.RA0.to_numpy() * u.deg, c.Dec0.to_numpy() * u.deg)).arcsec
        return s0
    cases["sky_sep_vs_sorcha_arcsec"] = np.nan
    for L in ("A", "B"):
        mask = cases.case == L
        cases.loc[mask, "sky_sep_vs_sorcha_arcsec"] = sep_to_C(cases[mask])

    cases.to_parquet(OUT_DIR / "stage3_cases_794NEO.parquet", index=False)

    # ---- decomposition ----
    piv_d2 = cases.pivot_table(index="ObjID", columns="case", values="digest2_det")
    piv_vd = cases.pivot_table(index="ObjID", columns="case", values="vdp")

    def stat_block(piv, a, b, clf):
        d = (piv[a] - piv[b]).dropna()
        ad_ = d.abs()
        return dict(classifier=clf, diff=f"{a}-{b}", n=len(d),
                    median=round(d.median(), 4), mean=round(d.mean(), 4),
                    median_abs=round(ad_.median(), 4), p95_abs=round(np.percentile(ad_, 95), 4),
                    p99_abs=round(np.percentile(ad_, 99), 4), max_abs=round(ad_.max(), 4),
                    frac_gt_0_05=round((ad_ > 0.05).mean(), 4), frac_gt_0_10=round((ad_ > 0.10).mean(), 4))
    comps = [("B", "A"), ("C", "B"), ("D", "B"), ("D", "C")]
    rows = []
    for a, b in comps:
        rows.append(stat_block(piv_d2, a, b, "digest2"))
        rows.append(stat_block(piv_vd, a, b, "VDP"))
    summ = pd.DataFrame(rows)
    summ.to_csv(OUT_DIR / "STAGE3_decomposition.csv", index=False)

    # ---- per-object table ----
    per = piv_d2.rename(columns={c: f"d2_{c}" for c in "ABCD"}).join(
          piv_vd.rename(columns={c: f"vdp_{c}" for c in "ABCD"}))
    per["d2_B_minus_A"] = per.d2_B - per.d2_A
    per["d2_D_minus_B"] = per.d2_D - per.d2_B
    per["d2_D_minus_C"] = per.d2_D - per.d2_C
    per = per.join(cases[cases.case == "B"].set_index("ObjID")[["sky_sep_vs_sorcha_arcsec", "dt_days"]])
    per.to_parquet(OUT_DIR / "stage3_per_object_794NEO.parquet")

    # ---- figures ----
    _figs(cases, piv_d2, piv_vd, per)

    pd.set_option("display.width", 240)
    print("\n=== STAGE 3 decomposition (abs score differences) ===", flush=True)
    print(summ.to_string(index=False), flush=True)
    print(f"\n  median sky_sep A vs Sorcha: {cases[cases.case=='A'].sky_sep_vs_sorcha_arcsec.median():.0f}\"  "
          f"B vs Sorcha: {cases[cases.case=='B'].sky_sep_vs_sorcha_arcsec.median():.0f}\"", flush=True)
    print(f"  wrote stage3_cases_794NEO.parquet, stage3_per_object_794NEO.parquet, "
          f"STAGE3_decomposition.csv, figures", flush=True)


def _fmt(d12, dt, ra_deg, dec_deg, mag):
    day_frac = (dt.day + dt.hour / 24.0 + dt.minute / 1440.0
                + (dt.second + dt.microsecond / 1e6) / 86400.0)
    ra_h = (float(ra_deg) % 360.0) / 15.0
    rah = int(ra_h); ram_f = (ra_h - rah) * 60.0; ram = int(ram_f)
    ras = min((ram_f - ram) * 60.0, 59.99)
    dec_a = abs(float(np.clip(dec_deg, -89.99, 89.99)))
    sign = "+" if dec_deg >= 0 else "-"
    decd = int(dec_a); decm_f = (dec_a - decd) * 60.0; decm = int(decm_f)
    decs = min((decm_f - decm) * 60.0, 59.9)
    mag_v = max(0.0, min(99.9, float(mag) if np.isfinite(mag) else 21.0))
    line = (f"{d12[:12]:12s}  C{dt.year:04d} {dt.month:02d} {day_frac:08.5f}"
            f" {rah:02d} {ram:02d} {ras:05.2f} {sign}{decd:02d} {decm:02d} {decs:04.1f}"
            f"          {mag_v:4.1f} V      {ad.OBSCODE:3s}")
    assert len(line) == 80, f"len {len(line)}"
    return line


def _figs(cases, piv_d2, piv_vd, per):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # histograms of the decomposition (digest2)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for a, (aa, bb, ttl) in zip(ax, [("B", "A", "time/epoch  B-A"), ("D", "B", "remaining (V)  D-B"),
                                     ("D", "C", "band  D-C")]):
        d = (piv_d2[aa] - piv_d2[bb]).dropna()
        a.hist(d, bins=np.linspace(-1, 1, 41), color="tab:orange", alpha=0.8)
        a.axvline(0, color="k", lw=1, ls="--")
        a.set_title(f"digest2 {ttl}\nmed {d.median():+.3f}, |Δ|>0.1 {100*(d.abs()>0.1).mean():.0f}%", fontsize=10)
        a.set_xlabel("Δ score"); a.grid(alpha=0.25)
    fig.suptitle("Stage 3 digest2 score-difference decomposition (794 NEOs)", fontsize=12)
    fig.tight_layout(); fig.savefig(OUT_DIR / "stage3_digest2_decomposition_hist.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    # scatter benchmark(A,B) vs Sorcha(D) digest2
    fig, ax = plt.subplots(1, 2, figsize=(11, 5.2))
    for a, L, ttl in zip(ax, ["A", "B"], ["A original", "B synced"]):
        d = piv_d2[[L, "D"]].dropna()
        a.scatter(d[L], d["D"], s=14, alpha=0.5, color="tab:purple")
        a.plot([0, 1], [0, 1], "k--", lw=1); a.set_aspect("equal")
        a.set_xlim(-.03, 1.03); a.set_ylim(-.03, 1.03); a.grid(alpha=0.25)
        a.set_xlabel(f"benchmark {ttl} digest2"); a.set_ylabel("Sorcha D digest2")
        a.set_title(f"{ttl} vs Sorcha  (med|Δ|={ (d['D']-d[L]).abs().median():.3f})", fontsize=10)
    fig.suptitle("Stage 3 benchmark vs Sorcha digest2 — before/after time-sync", fontsize=12)
    fig.tight_layout(); fig.savefig(OUT_DIR / "stage3_bench_vs_sorcha_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("propagate").set_defaults(func=propagate)
    sub.add_parser("analyze").set_defaults(func=analyze)
    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
