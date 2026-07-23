#!/usr/bin/env python3
"""Controlled digest2 EXACT-INPUT audit — separate Monte-Carlo noise from input differences.

For the 794 identity-matched night-61642 NEOs (ObjID == s3m_objid), and the three input variants
already built by the digest2 audit:
    B = benchmark-native V input
    S = Sorcha raw LSST magnitudes labelled V
    V = Sorcha V-corrected magnitudes
this runs digest2 DETERMINISTICALLY (repeatable mode, --cpu 1) on the exact same two MPC lines
TWICE per object and measures `delta_same_input = |run2 - run1|`. In repeatable mode that delta
must be ~0 — which proves any B/S/V score difference is caused by the INPUT tracklet, not the RNG.

It also reports the cross-input comparisons S-B, V-B, V-S using the deterministic (run1) score.

Strict parsing: a missing designation is NaN (never 0); a value > 100 is clamped to 100 and
counted; duplicates and malformed lines are counted. Nothing is silently defaulted.

Reuses the exact MPC-line generator (format_mpc80) and constants from audit_digest2.py, so the
lines are byte-identical to the main audit. Reads the pairs table it already built.

Outputs (new dir, nothing production touched):
    outputs/digest2_exact_input_audit_night61642/
        per_object_scores.parquet      one row per (ObjID, variant): both runs + delta + obs fields + hash
        same_input_delta_summary.csv   per-variant reproducibility stats
        cross_input_summary.csv        S-B, V-B, V-S distributions (deterministic scores)
        parse_issues.csv               missing / duplicate / malformed / out-of-range counts
        raw/<variant>_run{1,2}.{obs,out,err}

Usage:
    PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
    $PY neomod/pipeline/audit_digest2_exact_input.py
"""
from __future__ import annotations

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
import audit_digest2 as ad  # noqa: E402  (format_mpc80, constants, pairs path)

OUT_DIR = ad.WORKDIR / "outputs" / "digest2_exact_input_audit_night61642"
RAW_DIR = OUT_DIR / "raw"
VARIANTS = ("B", "S", "V")
CLAMP_MAX = 100  # digest2 NEO score is a 0-100 percentage; anything above is rejected/clamped
# `repeatable` is thread-count-independent and reproducible (verified: cpu4 == cpu1, run1 == run2),
# so we parallelise for speed without changing any score. Override with DIGEST2_AUDIT_CPUS.
CPUS = os.environ.get("DIGEST2_AUDIT_CPUS", "8")


def line_hash(l0: str, l1: str) -> str:
    return hashlib.sha256((l0 + "\n" + l1 + "\n").encode()).hexdigest()[:16]


def run_digest2(obs_lines: list[str], tag: str) -> tuple[dict[str, list[int]], int, float]:
    """Deterministic digest2 (repeatable + --cpu 1). Saves obs/out/err. Returns
    (desig -> list of raw int scores, returncode, runtime_s)."""
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    obs_path = RAW_DIR / f"{tag}.obs"
    obs_path.write_text("\n".join(obs_lines) + "\n")
    with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as cfg:
        cfg.write(ad.D2_CONFIG)      # includes `repeatable`
        cfg_path = cfg.name
    t0 = time.time()
    try:
        proc = subprocess.run([str(ad.DIGEST2_EXEC), "-p", str(ad.DIGEST2_DIR), "-c", cfg_path,
                               "--cpu", CPUS, str(obs_path)],
                              capture_output=True, text=True, timeout=3600)
        out, err, rc = proc.stdout, proc.stderr, proc.returncode
    finally:
        os.unlink(cfg_path)
    dt = time.time() - t0
    (RAW_DIR / f"{tag}.out").write_text(out)
    (RAW_DIR / f"{tag}.err").write_text(err)
    parsed: dict[str, list[int]] = {}
    for ln in out.splitlines():
        p = ln.strip().split()
        if len(p) < 2:
            continue
        try:
            parsed.setdefault(p[0], []).append(int(p[1]))
        except ValueError:
            parsed.setdefault(p[0], [])   # seen but unparseable -> empty list -> 'malformed'
    return parsed, rc, dt


def strict_score(vals):
    """(norm, status, out_of_range_flag) from the raw-int occurrence list for one designation."""
    if vals is None:
        return np.nan, "missing", False
    if len(vals) == 0:
        return np.nan, "malformed", False
    if len(vals) > 1:
        return np.nan, "duplicate", False
    raw = vals[0]
    oor = not (0 <= raw <= CLAMP_MAX)
    raw_c = min(max(raw, 0), CLAMP_MAX)     # clamp out-of-range
    return raw_c / 100.0, "ok", oor


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_parquet(ad.PAIRS_FILE)
    neo = pairs[(pairs.is_neo == 1) & (pairs.pair_index < ad.SPECIAL_PI0)].copy()
    n_obj = neo.pair_index.nunique()
    print(f"exact-input audit: {n_obj} identity-matched night-61642 NEOs x {len(VARIANTS)} variants",
          flush=True)

    rows = []
    parse_counts = {v: dict(missing=0, duplicate=0, malformed=0, out_of_range=0) for v in VARIANTS}
    for v in VARIANTS:
        vs = neo[neo.variant == v].sort_values("pair_index")
        obs_lines, blocks = [], []
        for _, r in vs.iterrows():
            l0 = ad.format_mpc80(r.desig12, r.mjd0_utc, r.ra0, r.dec0, r.mag0)
            l1 = ad.format_mpc80(r.desig12, r.mjd1_utc, r.ra1, r.dec1, r.mag1)
            obs_lines += [l0, l1]
            blocks.append((r, l0, l1))
        # two independent invocations of the identical input
        p1, rc1, dt1 = run_digest2(obs_lines, f"{v}_run1")
        p2, rc2, dt2 = run_digest2(obs_lines, f"{v}_run2")
        print(f"  variant={v}: n={len(vs)} rc=({rc1},{rc2}) runtime=({dt1:.1f},{dt2:.1f})s "
              f"desigs_out=({len(p1)},{len(p2)})", flush=True)
        for r, l0, l1 in blocks:
            n1, s1, o1 = strict_score(p1.get(r.desig))
            n2, s2, o2 = strict_score(p2.get(r.desig))
            for st in (s1, s2):
                if st != "ok":
                    parse_counts[v][st] += 1
            parse_counts[v]["out_of_range"] += int(o1) + int(o2)
            rows.append({
                "ObjID": r.ObjID, "variant": v, "input_hash": line_hash(l0, l1),
                "desig": r.desig,
                "mjd0": r.mjd0_utc, "mjd1": r.mjd1_utc,
                "ra0": r.ra0, "dec0": r.dec0, "ra1": r.ra1, "dec1": r.dec1,
                "mag0": r.mag0, "mag1": r.mag1, "filter0": r.filter0, "filter1": r.filter1,
                "mag_source": r.mag_source,
                "digest2_score_run1": n1, "digest2_score_run2": n2,
                "parse_status_run1": s1, "parse_status_run2": s2,
                "delta_same_input": (abs(n2 - n1) if (s1 == "ok" and s2 == "ok") else np.nan),
                "mpc_line0": l0, "mpc_line1": l1,
            })

    tbl = pd.DataFrame(rows)
    tbl.to_parquet(OUT_DIR / "per_object_scores.parquet", index=False)

    # ---- per-variant same-input reproducibility ----
    def q(x, p): return float(np.nanpercentile(x, p)) if len(x) else float("nan")
    summ = []
    for v in VARIANTS:
        d = tbl[tbl.variant == v].delta_same_input.to_numpy(float)
        d = d[~np.isnan(d)]
        pc = parse_counts[v]
        summ.append({
            "variant": v, "n": int((tbl.variant == v).sum()), "n_delta": len(d),
            "median_abs_delta": float(np.median(d)) if len(d) else float("nan"),
            "p95_abs_delta": q(d, 95), "p99_abs_delta": q(d, 99),
            "max_abs_delta": float(d.max()) if len(d) else float("nan"),
            "frac_gt_0.05": float((d > 0.05).mean()) if len(d) else float("nan"),
            "frac_gt_0.10": float((d > 0.10).mean()) if len(d) else float("nan"),
            "missing": pc["missing"], "duplicate": pc["duplicate"],
            "malformed": pc["malformed"], "out_of_range": pc["out_of_range"],
        })
    same = pd.DataFrame(summ)
    same.to_csv(OUT_DIR / "same_input_delta_summary.csv", index=False)
    pd.DataFrame([{"variant": v, **parse_counts[v]} for v in VARIANTS]).to_csv(
        OUT_DIR / "parse_issues.csv", index=False)

    # ---- cross-input comparisons (deterministic run1 score per input) ----
    wide = tbl.pivot_table(index="ObjID", columns="variant", values="digest2_score_run1")
    cross = []
    for name, a, b in [("S_minus_B", "S", "B"), ("V_minus_B", "V", "B"), ("V_minus_S", "V", "S")]:
        if not {a, b}.issubset(wide.columns):
            continue
        d = (wide[a] - wide[b]).dropna()
        ad_ = d.abs()
        cross.append({
            "comparison": name, "n": len(d),
            "median": float(d.median()), "mean": float(d.mean()),
            "median_abs": float(ad_.median()), "p95_abs": q(ad_.to_numpy(), 95),
            "p99_abs": q(ad_.to_numpy(), 99), "max_abs": float(ad_.max()),
            "frac_abs_gt_0.05": float((ad_ > 0.05).mean()),
            "frac_abs_gt_0.10": float((ad_ > 0.10).mean()),
        })
    crossdf = pd.DataFrame(cross)
    crossdf.to_csv(OUT_DIR / "cross_input_summary.csv", index=False)

    pd.set_option("display.width", 200)
    print(f"\n=== SAME-INPUT reproducibility (repeatable, --cpu {CPUS} [score is cpu-count-independent]; "
          f"rerun of IDENTICAL MPC lines) ===", flush=True)
    print(same.to_string(index=False), flush=True)
    print("\n=== CROSS-INPUT differences (deterministic score; real input changes) ===", flush=True)
    print(crossdf.to_string(index=False), flush=True)
    print(f"\nwrote per_object_scores.parquet, same_input_delta_summary.csv, cross_input_summary.csv, "
          f"parse_issues.csv -> {OUT_DIR}", flush=True)

    same_ok = float(np.nanmax(same.max_abs_delta.to_numpy())) if len(same) else float("nan")
    print(f"\nHEADLINE: max |delta_same_input| across all variants/objects = {same_ok:.4f}", flush=True)
    print("  -> if ~0, identical MPC lines reproduce EXACTLY; all B/S/V differences are input-driven.",
          flush=True)


if __name__ == "__main__":
    main()
