#!/usr/bin/env python3
"""STAGE 2 — identical-input VDP-vs-digest2 comparison on a canonical tracklet manifest.

For each of the three Stage-1 variants (bench native V, sorcha_rawfix raw mags, sorcha_vcorr
V-corrected), build a canonical tracklet manifest and score it with BOTH classifiers from the SAME
tracklet:

  digest2 : the Stage-1 deterministic score (P_NEO_d2_det), computed from the exact MPC lines whose
            SHA256 we recompute here and VERIFY against the Stage-1 input hash (proves same lines).
  VDP     : re-scored against the n-body maps (prob_maps_grid_s3m_nbody) from the SAME tracklet's
            stored kinematics (vlam, vbeta) and the SAME per-variant magnitude digest2 used
            (bench mean_mag / sorcha mean_mag / sorcha mean_mag_V). Only the magnitude differs
            between variants; vlam/vbeta are fixed — exactly parallel to digest2.

Both classifiers are therefore evaluated on the identical canonical tracklet (same parquet row,
same photometry), and the hash is the machine-checkable proof for the digest2 side.

Invalid digest2 cases (missing/duplicate/malformed from Stage 1) are kept EXPLICIT — never zeroed.
Metrics are computed on the digest2-scorable set and the excluded counts are reported.

Commands
  score-vdp --source S   re-score VDP for one source against the n-body maps -> sidecar parquet
  build     --source S   build manifest, verify hashes, compute ROC/F1/metrics, figures

Outputs: outputs/mag245_nbody_exact_classifier_test/
Nothing here modifies production parquets.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
import audit_digest2 as ad          # noqa: E402  format_mpc80, OBSCODE
import rescore_vdp_Vband as rv       # noqa: E402  rescore() against a given map dir

WORKDIR = ad.WORKDIR
S1_DIR = WORKDIR / "outputs" / "mag245_nbody_deterministic_rescore"
OUT_DIR = WORKDIR / "outputs" / "mag245_nbody_exact_classifier_test"
MAPDIR_NBODY = WORKDIR / "prob_maps_grid_s3m_nbody"
SORCHA_VMAGS = S1_DIR / "sorcha_vmags.parquet"

# source -> (det parquet, VDP magnitude column, per-detection mag0/mag1 cols, filters?)
SOURCES = {
    "bench":         dict(vdp_mag="mean_mag",   m0="mag0",   m1="mag1",   filters=False),
    "sorcha_rawfix": dict(vdp_mag="mean_mag",   m0="mag0",   m1="mag1",   filters=True),
    "sorcha_vcorr":  dict(vdp_mag="mean_mag_V", m0="mag0_V", m1="mag1_V", filters=True),
}


def _det_parquet(source):
    return S1_DIR / f"{source}_deterministic.parquet"


def _hash(l0, l1):
    return hashlib.sha256((l0 + "\n" + l1 + "\n").encode()).hexdigest()[:16]


def _fmt_line(d12: str, dt, ra_deg: float, dec_deg: float, mag: float) -> str:
    """Byte-identical to audit_digest2.format_mpc80 but takes a precomputed UTC datetime `dt`
    (so the astropy Time conversion can be vectorised once over the whole column)."""
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
    if len(line) != 80:
        raise AssertionError(f"MPC line length {len(line)}: {line!r}")
    return line


def _utc_datetimes(mjd):
    from astropy.time import Time
    return Time(np.asarray(mjd, float), format="mjd", scale="utc").utc.datetime


# ---------------------------------------------------------------------------- score-vdp
def score_vdp(args) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cfg = SOURCES[args.source]
    df = pd.read_parquet(_det_parquet(args.source),
                         columns=["vlam", "vbeta", "prob_map_file", cfg["vdp_mag"]])
    mag = df[cfg["vdp_mag"]].to_numpy(float)
    print(f"[{args.source}] VDP re-score {len(df):,} rows vs n-body maps, mag={cfg['vdp_mag']}", flush=True)
    P = rv.rescore(df, mag, mapdir=str(MAPDIR_NBODY))
    out = OUT_DIR / f"vdp_{args.source}.parquet"
    pd.DataFrame({"P_NEO_vdp_canon": P}).to_parquet(out, index=False)
    print(f"[{args.source}] wrote {out}  non-null {100*np.isfinite(P).mean():.1f}%", flush=True)


# ---------------------------------------------------------------------------- build manifest + metrics
def _metrics(y, s):
    from sklearn.metrics import precision_recall_curve, roc_auc_score
    p, r, th = precision_recall_curve(y, s)
    f1 = np.divide(2 * p * r, p + r, out=np.zeros_like(p), where=(p + r) > 0)
    bi = int(np.argmax(f1[:-1]))
    return dict(AUC=roc_auc_score(y, s), bestF1=f1[bi], thresh=float(th[bi] if bi < len(th) else 1.0),
                completeness=r[bi] * 100, contamination=(1 - p[bi]) * 100)


def build(args) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cfg = SOURCES[args.source]
    cols = ["ObjID", "population", "ra0", "dec0", "mjd0_utc", "ra1", "dec1", "mjd1_utc",
            "mag0", "mag1", "P_NEO_d2_det", "d2_det_parse_status", "d2_det_input_hash"]
    if cfg["filters"]:
        cols += ["filter0", "filter1"]
    det = pd.read_parquet(_det_parquet(args.source), columns=cols).reset_index(drop=True)
    det["ObjID"] = det.ObjID.astype(str)
    if cfg["m0"] not in det.columns:  # vcorr per-detection V-mags come from the sidecar
        vm = pd.read_parquet(SORCHA_VMAGS).reset_index(drop=True)
        det["mag0_V"] = vm["mag0_V"].to_numpy(); det["mag1_V"] = vm["mag1_V"].to_numpy()
    vdp = pd.read_parquet(OUT_DIR / f"vdp_{args.source}.parquet").reset_index(drop=True)
    det["P_NEO_vdp_canon"] = vdp["P_NEO_vdp_canon"].to_numpy()

    m0col, m1col = cfg["m0"], cfg["m1"]
    n = len(det)
    # vectorise the (slow) astropy Time conversion once per column, then format per row
    dt0 = _utc_datetimes(det.mjd0_utc.to_numpy()); dt1 = _utc_datetimes(det.mjd1_utc.to_numpy())
    ra0 = det.ra0.to_numpy(); dec0 = det.dec0.to_numpy(); m0 = det[m0col].to_numpy(float)
    ra1 = det.ra1.to_numpy(); dec1 = det.dec1.to_numpy(); m1 = det[m1col].to_numpy(float)
    stored_hash = det.d2_det_input_hash.to_numpy()
    # spot-check the fast path is byte-identical to audit_digest2.format_mpc80 on a few rows
    for i in (0, n // 2, n - 1):
        d12 = f"     D{i:06d}"
        assert _fmt_line(d12, dt0[i], ra0[i], dec0[i], m0[i]) == \
            ad.format_mpc80(d12, det.mjd0_utc.iloc[i], ra0[i], dec0[i], m0[i]), "fast MPC path drift"

    line0 = [None] * n; line1 = [None] * n; sha = [None] * n; match = np.zeros(n, bool)
    for i in range(n):
        d12 = f"     D{i:06d}"
        l0 = _fmt_line(d12, dt0[i], ra0[i], dec0[i], m0[i])
        l1 = _fmt_line(d12, dt1[i], ra1[i], dec1[i], m1[i])
        h = _hash(l0, l1)
        line0[i] = l0; line1[i] = l1; sha[i] = h; match[i] = (h == stored_hash[i])
    hash_ok = int(match.sum())

    man = pd.DataFrame({
        "ObjID": det.ObjID.to_numpy(), "population": det.population.to_numpy(), "source": args.source,
        "t0": det.mjd0_utc.to_numpy(), "t1": det.mjd1_utc.to_numpy(),
        "RA0": ra0, "Dec0": dec0, "RA1": ra1, "Dec1": dec1, "mag0": m0, "mag1": m1,
        "filter0": (det.filter0.to_numpy() if cfg["filters"] else "V"),
        "filter1": (det.filter1.to_numpy() if cfg["filters"] else "V"),
        "obscode": ad.OBSCODE, "mpc_line0": line0, "mpc_line1": line1,
        "sha256_input": sha, "hash_matches_stage1": match,
        "digest2_det": det.P_NEO_d2_det.to_numpy(), "digest2_parse_status": det.d2_det_parse_status.to_numpy(),
        "vdp_canon": det.P_NEO_vdp_canon.to_numpy(),
    })
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    man.to_parquet(OUT_DIR / f"manifest_{args.source}.parquet", index=False)
    print(f"\n=== STAGE 2 build [{args.source}] — {n:,} canonical tracklets ===", flush=True)
    print(f"  hash verification (recomputed == Stage-1 digest2 input hash): {hash_ok:,}/{n:,} "
          f"({100*hash_ok/n:.4f}%)", flush=True)
    if hash_ok != n:
        print("  !! WARNING: some hashes do not match — digest2 and manifest lines differ", flush=True)

    # ---- classifier metrics ----
    is_neo = (man.population == "NEO").to_numpy().astype(int)
    d2_ok = (man.digest2_parse_status == "ok").to_numpy()
    vdp_ok = np.isfinite(man.vdp_canon.to_numpy())
    d2 = man.digest2_det.to_numpy()
    vd = man.vdp_canon.to_numpy()
    # common scorable set: digest2 ok AND VDP finite (same rows for a fair head-to-head)
    common = d2_ok & vdp_ok
    n_excl_d2 = int((~d2_ok).sum()); n_excl_vdp = int((~vdp_ok).sum())
    by_status = man.digest2_parse_status.value_counts().to_dict()
    print(f"  digest2 parse: {by_status}  | VDP non-finite: {n_excl_vdp}", flush=True)
    print(f"  metrics computed on the COMMON scorable set ({int(common.sum()):,} rows); "
          f"excluded digest2-invalid={n_excl_d2}, VDP-nonfinite={n_excl_vdp} (NOT zeroed)", flush=True)

    md2 = _metrics(is_neo[common], d2[common])
    mvd = _metrics(is_neo[common], vd[common])
    summ = pd.DataFrame([
        {"source": args.source, "classifier": "digest2_det", "n": int(common.sum()), **md2},
        {"source": args.source, "classifier": "VDP_canon", "n": int(common.sum()), **mvd},
    ])
    summ["n_excl_digest2_invalid"] = n_excl_d2
    summ["n_excl_vdp_nonfinite"] = n_excl_vdp
    summ.to_csv(OUT_DIR / f"metrics_{args.source}.csv", index=False)
    pd.set_option("display.width", 220)
    print(summ.to_string(index=False), flush=True)

    # ---- ROC figure (contamination vs completeness) ----
    from sklearn.metrics import precision_recall_curve
    fig, ax = plt.subplots(figsize=(6.4, 5.8))
    for name, s, c in [("digest2 (det)", d2, "tab:orange"), ("VDP (canonical)", vd, "tab:blue")]:
        p, rr, _ = precision_recall_curve(is_neo[common], s[common])
        m = _metrics(is_neo[common], s[common])
        ax.plot(rr * 100, (1 - p) * 100, lw=2, color=c,
                label=f"{name}  F1={m['bestF1']:.3f} AUC={m['AUC']:.3f}")
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.grid(alpha=0.3)
    ax.set_xlabel("NEO completeness (%)"); ax.set_ylabel("contamination (%)")
    ax.set_title(f"Stage 2 identical-input VDP vs digest2 — {args.source}\n"
                 f"(common scorable N={int(common.sum()):,}; hash-verified)", fontsize=10)
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(OUT_DIR / f"roc_{args.source}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote manifest_{args.source}.parquet, metrics_{args.source}.csv, roc_{args.source}.png", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, fn in [("score-vdp", score_vdp), ("build", build)]:
        sp = sub.add_parser(name)
        sp.add_argument("--source", choices=list(SOURCES), required=True)
        sp.set_defaults(func=fn)
    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
