#!/usr/bin/env python3
"""STAGE 1 — deterministic digest2 production rescore of the mag245 n-body sets.

Reruns digest2 for the benchmark and Sorcha mag245 n-body sets with the bugs the audit found fixed:
  * `repeatable` mode (production ran `random` -> ±1-3 pt jitter; repeatable is deterministic and
    cpu-count-independent, verified) so scores are exactly reproducible;
  * strict parsing — a missing designation is NaN (never 0); a raw score outside 0-100 is clamped and
    counted; duplicates/malformed counted;
  * full provenance — every MPC input line + SHA256 hash + raw stdout/stderr saved.

Three deterministic outputs (production parquets are NOT modified; new column P_NEO_d2_det is added
to a COPY):
  bench          benchmark-native synthetic V   (mag0, mag1)
  sorcha_rawfix  Sorcha raw LSST mags labelled V (mag0, mag1)               -- the reproduce-but-fixed set
  sorcha_vcorr   Sorcha V-corrected mags         (mag0_V, mag1_V)           -- consistent V-band vs benchmark

Commands
  prep-sorcha-vmags   compute per-detection mag0_V/mag1_V for the full Sorcha set (once) -> sidecar parquet
  run --source S --chunk-index i --chunk-size N   score one chunk (Slurm array task)
  combine --source S                              merge chunks -> replacement parquet + diff table
  analyze --source S                              ROC/F1/AUC + figures + summary vs stored P_NEO_d2

Determinism: digest2 `repeatable` + any --cpu (score is cpu-count-independent). Uses DIGEST2_AUDIT_CPUS.

Nothing here overwrites production columns; P_NEO_d2 stays, P_NEO_d2_det is added.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_digest2 as ad  # noqa: E402  format_mpc80, DIGEST2_*, D2_CONFIG, constants

WORKDIR = ad.WORKDIR
OUT_DIR = WORKDIR / "outputs" / "mag245_nbody_deterministic_rescore"
BENCH_PARQUET = WORKDIR / "outputs/phase2_benchmark_s3m_nbody_mag245" / \
    "benchmark_comparison_s3m_nbody_mag245.parquet"
SORCHA_PARQUET = WORKDIR / "outputs/s3m_linking/case1/sorcha_comparison_case1_nbody_Vband.parquet"
SORCHA_VMAGS = OUT_DIR / "sorcha_vmags.parquet"
PHYS_CSV = WORKDIR / "inputs" / "s3m_sorcha_phys.csv"
COLOR_COLS = {"u": "u-r", "g": "g-r", "i": "i-r", "z": "z-r", "y": "y-r"}  # r -> 0
CPUS = os.environ.get("DIGEST2_AUDIT_CPUS", "8")
CLAMP_MAX = 100

# source key -> (parquet, mag0 col, mag1 col, uses vmags sidecar)
SOURCES = {
    "bench":         (BENCH_PARQUET,  "mag0",   "mag1",   False),
    "sorcha_rawfix": (SORCHA_PARQUET, "mag0",   "mag1",   False),
    "sorcha_vcorr":  (SORCHA_PARQUET, "mag0_V", "mag1_V", True),
}


def _hash(l0, l1):
    return hashlib.sha256((l0 + "\n" + l1 + "\n").encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------- prep V-mags
def prep_sorcha_vmags(args) -> None:
    """Per-detection Johnson-V mags for the full Sorcha set, from stored mean_mag_V + phys colours.

    mean_mag == 0.5*(mag0+mag1) exactly (verified), and
        mean_mag_V = mean_mag - 0.5*(c0+c1) + (H_V-H_r),   c=colour(filter), c_r=0
    => mag_i_V = mag_i - mean_mag + mean_mag_V + 0.5*(c_other - c_i).  No census scan needed.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    s = pd.read_parquet(SORCHA_PARQUET,
                        columns=["ObjID", "mag0", "mag1", "mean_mag", "mean_mag_V", "filter0", "filter1"])
    s["ObjID"] = s.ObjID.astype(str)
    ids = set(s.ObjID)
    print(f"reading phys colours for {len(ids):,} unique objids ...", flush=True)
    usecols = ["ObjID"] + list(COLOR_COLS.values())
    parts = []
    for ch in pd.read_csv(PHYS_CSV, usecols=usecols, dtype={"ObjID": str}, chunksize=2_000_000):
        sub = ch[ch.ObjID.isin(ids)]
        if len(sub):
            parts.append(sub)
    phys = pd.concat(parts, ignore_index=True).drop_duplicates("ObjID").set_index("ObjID")

    def colour(filters):
        out = np.zeros(len(s))
        f = filters.to_numpy(); oid = s.ObjID.to_numpy()
        for band, col in COLOR_COLS.items():
            m = f == band
            if m.any():
                out[m] = phys.reindex(oid[m])[col].to_numpy(float)
        return out

    c0 = colour(s.filter0); c1 = colour(s.filter1)
    base = s.mean_mag_V - s.mean_mag
    s["mag0_V"] = s.mag0 + base + 0.5 * (c1 - c0)
    s["mag1_V"] = s.mag1 + base + 0.5 * (c0 - c1)
    # sanity: (mag0_V+mag1_V)/2 == mean_mag_V
    err = (0.5 * (s.mag0_V + s.mag1_V) - s.mean_mag_V).abs().max()
    nmiss = int(s.mag0_V.isna().sum())
    s[["mag0_V", "mag1_V"]].to_parquet(SORCHA_VMAGS, index=False)
    print(f"wrote {SORCHA_VMAGS}  rows={len(s):,}  max|mean check|={err:.2e}  NaN={nmiss}", flush=True)


# ---------------------------------------------------------------------------- run one chunk
def _digest2(obs_lines, tag):
    raw_dir = OUT_DIR / "raw"; raw_dir.mkdir(parents=True, exist_ok=True)
    obs_path = raw_dir / f"{tag}.obs"; obs_path.write_text("\n".join(obs_lines) + "\n")
    with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as cfg:
        cfg.write(ad.D2_CONFIG); cfg_path = cfg.name          # includes `repeatable`
    t0 = time.time()
    try:
        p = subprocess.run([str(ad.DIGEST2_EXEC), "-p", str(ad.DIGEST2_DIR), "-c", cfg_path,
                            "--cpu", CPUS, str(obs_path)], capture_output=True, text=True, timeout=7200)
        out, err, rc = p.stdout, p.stderr, p.returncode
    finally:
        os.unlink(cfg_path)
    (raw_dir / f"{tag}.out").write_text(out); (raw_dir / f"{tag}.err").write_text(err)
    parsed = {}
    for ln in out.splitlines():
        q = ln.strip().split()
        if len(q) < 2:
            continue
        try:
            parsed.setdefault(q[0], []).append(int(q[1]))
        except ValueError:
            parsed.setdefault(q[0], [])
    return parsed, rc, time.time() - t0


def run(args) -> None:
    parquet, m0, m1, use_vmags = SOURCES[args.source]
    cols = ["ObjID", "population", "ra0", "dec0", "mjd0_utc", "ra1", "dec1", "mjd1_utc"]
    base_m0, base_m1 = ("mag0", "mag1") if use_vmags else (m0, m1)
    cols += [base_m0, base_m1]
    lo = args.chunk_index * args.chunk_size
    n_rows = ad_num_rows(parquet)
    hi = min(lo + args.chunk_size, n_rows)
    if lo >= n_rows:
        print(f"chunk {args.chunk_index} beyond {n_rows} rows; nothing to do", flush=True); return
    df = pd.read_parquet(parquet, columns=cols).iloc[lo:hi].reset_index(drop=True)
    if use_vmags:
        vm = pd.read_parquet(SORCHA_VMAGS).iloc[lo:hi].reset_index(drop=True)
        df["mag0_V"] = vm["mag0_V"].to_numpy(); df["mag1_V"] = vm["mag1_V"].to_numpy()
    df["ObjID"] = df.ObjID.astype(str)

    obs_lines, desigs, hashes = [], [], []
    for i, r in enumerate(df.itertuples(index=False)):
        g = lo + i
        key = f"D{g:06d}"; d12 = f"     {key}"
        l0 = ad.format_mpc80(d12, getattr(r, "mjd0_utc"), r.ra0, r.dec0, getattr(r, m0))
        l1 = ad.format_mpc80(d12, getattr(r, "mjd1_utc"), r.ra1, r.dec1, getattr(r, m1))
        obs_lines += [l0, l1]; desigs.append(key); hashes.append(_hash(l0, l1))

    tag = f"{args.source}_chunk{args.chunk_index:04d}"
    parsed, rc, dt = _digest2(obs_lines, tag)
    print(f"{tag}: rows {lo}-{hi} rc={rc} runtime={dt:.1f}s out_desigs={len(parsed)}", flush=True)

    rows = []
    for i, (key, h) in enumerate(zip(desigs, hashes)):
        vals = parsed.get(key)
        if vals is None:
            status, raw, norm, oor = "missing", np.nan, np.nan, False
        elif len(vals) == 0:
            status, raw, norm, oor = "malformed", np.nan, np.nan, False
        elif len(vals) > 1:
            status, raw, norm, oor = "duplicate", np.nan, np.nan, False
        else:
            raw = vals[0]; oor = not (0 <= raw <= CLAMP_MAX)
            status, norm = "ok", min(max(raw, 0), CLAMP_MAX) / 100.0
        rows.append({"global_row": lo + i, "ObjID": df.ObjID.iloc[i], "population": df.population.iloc[i],
                     "digest2_id_det": key, "input_hash": h, "P_NEO_d2_det": norm,
                     "raw_int": raw, "parse_status": status, "out_of_range": oor})
    res_dir = OUT_DIR / "results" / args.source; res_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_parquet(res_dir / f"chunk{args.chunk_index:04d}.parquet", index=False)
    print(f"  wrote {len(rows)} rows", flush=True)


def ad_num_rows(parquet):
    import pyarrow.parquet as pq
    return pq.read_metadata(parquet).num_rows


# ---------------------------------------------------------------------------- combine
def combine(args) -> None:
    parquet = SOURCES[args.source][0]
    res_dir = OUT_DIR / "results" / args.source
    shards = sorted(res_dir.glob("chunk*.parquet"))
    if not shards:
        raise FileNotFoundError(f"no chunks in {res_dir}")
    det = pd.concat([pd.read_parquet(p) for p in shards], ignore_index=True).sort_values("global_row")
    n_rows = ad_num_rows(parquet)
    assert len(det) == n_rows, f"expected {n_rows} rows, got {len(det)} (missing chunks?)"
    assert det.global_row.is_monotonic_increasing and det.global_row.iloc[-1] == n_rows - 1

    stored = pd.read_parquet(parquet)
    out = stored.copy()
    out["P_NEO_d2_det"] = det.P_NEO_d2_det.to_numpy()
    out["d2_det_raw_int"] = det.raw_int.to_numpy()
    out["d2_det_parse_status"] = det.parse_status.to_numpy()
    out["d2_det_input_hash"] = det.input_hash.to_numpy()
    repl = OUT_DIR / f"{args.source}_deterministic.parquet"
    out.to_parquet(repl, index=False)

    ps = det.parse_status.value_counts().to_dict()
    oor = int(det.out_of_range.sum())
    d2 = det.P_NEO_d2_det
    stored_d2 = stored["P_NEO_d2"].to_numpy()
    change = np.abs(d2.to_numpy() - stored_d2)
    valid = ~np.isnan(change)
    print(f"\n=== combine {args.source}  ({len(det):,} rows) ===", flush=True)
    print(f"  parse: {ps}  out_of_range={oor}", flush=True)
    print(f"  det score range: [{np.nanmin(d2):.3f}, {np.nanmax(d2):.3f}]  NaN={int(d2.isna().sum())}", flush=True)
    print(f"  vs stored P_NEO_d2: changed(|Δ|>0)={int((change[valid] > 0).sum()):,}  "
          f"median|Δ|={np.nanmedian(change):.4f}  max|Δ|={np.nanmax(change):.4f}", flush=True)
    print(f"  wrote {repl}", flush=True)
    # score-diff table
    diff = pd.DataFrame({"ObjID": det.ObjID, "population": det.population,
                         "stored": stored_d2, "det": d2.to_numpy(), "abs_change": change})
    diff.to_parquet(OUT_DIR / f"{args.source}_scorediff.parquet", index=False)


# ---------------------------------------------------------------------------- analyze
def analyze(args) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import precision_recall_curve, roc_auc_score, roc_curve

    repl = OUT_DIR / f"{args.source}_deterministic.parquet"
    df = pd.read_parquet(repl)
    det_all = df["P_NEO_d2_det"].to_numpy()
    n_unscorable = int(np.isnan(det_all).sum())
    # Fair comparison on the SCORABLE set: rows digest2 could actually score deterministically.
    # (production silently zeroed the unscorable ones; we exclude them from BOTH curves instead.)
    scorable = ~np.isnan(det_all)
    df = df[scorable].reset_index(drop=True)
    is_neo = (df.population == "NEO").to_numpy().astype(int)
    stored = df["P_NEO_d2"].to_numpy()
    det = df["P_NEO_d2_det"].to_numpy()

    def metrics(y, s):
        p, r, th = precision_recall_curve(y, s)
        f1 = np.divide(2 * p * r, p + r, out=np.zeros_like(p), where=(p + r) > 0)
        bi = int(np.argmax(f1[:-1]))
        thr = th[bi] if bi < len(th) else 1.0
        return dict(AUC=roc_auc_score(y, s), bestF1=f1[bi], completeness=r[bi] * 100,
                    contamination=(1 - p[bi]) * 100, thresh=float(thr))

    ms, md = metrics(is_neo, stored), metrics(is_neo, det)
    print(f"  scorable rows: {len(df):,}  unscorable (det NaN): {n_unscorable}", flush=True)
    # threshold-classification changes at the deterministic best threshold
    thr = md["thresh"]
    cls_s = (np.nan_to_num(stored) >= thr).astype(int)
    cls_d = (np.nan_to_num(det) >= thr).astype(int)
    n_flip = int((cls_s != cls_d).sum())

    fig, ax = plt.subplots(1, 2, figsize=(12, 5.4))
    for a, (name, s, m) in zip(ax, [("stored (random)", stored, ms), ("deterministic", det, md)]):
        y = is_neo; sc = np.nan_to_num(s)
        p, r, _ = precision_recall_curve(y, sc)
        a.plot(r * 100, (1 - p) * 100, lw=2, color="tab:red" if "sorcha" in args.source else "tab:blue")
        a.set_xlim(0, 100); a.set_ylim(0, 100); a.grid(alpha=0.3)
        a.set_xlabel("NEO completeness (%)"); a.set_ylabel("contamination (%)")
        a.set_title(f"{name}\nF1={m['bestF1']:.3f} AUC={m['AUC']:.3f} "
                    f"compl={m['completeness']:.0f}% contam={m['contamination']:.0f}%", fontsize=10)
    fig.suptitle(f"Stage 1 digest2 ROC — {args.source} (n-body mag245)", fontsize=12)
    fig.tight_layout(); fig.savefig(OUT_DIR / f"{args.source}_roc.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    summ = pd.DataFrame([
        {"source": args.source, "which": "stored", **ms},
        {"source": args.source, "which": "deterministic", **md},
    ])
    summ["n_threshold_flips_vs_stored"] = [np.nan, n_flip]
    summ.to_csv(OUT_DIR / f"{args.source}_metrics.csv", index=False)
    print(f"\n=== analyze {args.source} ===", flush=True)
    print(summ.to_string(index=False), flush=True)
    print(f"  threshold flips (@det thr {thr:.3f}) vs stored: {n_flip:,}", flush=True)
    print(f"  wrote {args.source}_roc.png, {args.source}_metrics.csv", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prep-sorcha-vmags").set_defaults(func=prep_sorcha_vmags)
    for name, fn in [("run", run), ("combine", combine), ("analyze", analyze)]:
        sp = sub.add_parser(name)
        sp.add_argument("--source", choices=list(SOURCES), required=True)
        if name == "run":
            sp.add_argument("--chunk-index", type=int, required=True)
            sp.add_argument("--chunk-size", type=int, default=20000)
        sp.set_defaults(func=fn)
    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
