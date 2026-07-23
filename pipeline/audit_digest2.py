#!/usr/bin/env python3
"""Auditable digest2 re-scoring for the night-61642 benchmark/Sorcha comparison.

Purpose (see the plan in the 2026-07-22 chat / mag245_nbody_benchmark_vs_sorcha_roc.ipynb):
separate **digest2 implementation errors** (parsing, silent-zero defaults, the P_NEO_d2 = 1.17
out-of-range value, designation-echo mismatches) from **real score changes** caused by different
observations — most importantly the magnitude-band mismatch (production labels every MPC line "V"
but feeds Sorcha's raw r/i/g/z/y/u magnitudes; the benchmark feeds synthetic Johnson V).

This script NEVER touches the production parquets. Everything lands in
``outputs/digest2_audit_night61642/``.

Three input variants are run per matched object, each with its own unique designation:
  B  benchmark-native  — the exact two benchmark observations (mags already synthetic Johnson V)
  S  sorcha-native     — the exact Sorcha observations, raw mags labelled "V" (reproduces the
                         stored P_NEO_d2 — this is the reproduction gate)
  V  sorcha-Vcorrected — same Sorcha positions/times, but each r/i/... magnitude converted to
                         Johnson V before writing the MPC line (isolates the band effect: S->V
                         changes ONLY the magnitude band)
Benchmark and Sorcha pairs are NOT forced to share times/positions — that difference is part of
what we measure. S->V holds Sorcha position/time fixed and changes only the band.

Subcommands
    build-pairs   join the two parquets on identity + night 61642, compute per-detection V mags,
                  emit one row per (matched object x variant) with the exact obs to send to digest2.
    run           score one Slurm chunk of objects (all requested variants) with full auditing:
                  exact MPC lines, raw digest2 stdout+stderr, runtime, return code, strict parse.
    combine       merge chunk results; reproduction check; parse log; outliers; summary.

Digest2 is invoked ONCE per (chunk x variant): scores are per-tracklet independent, so batching
cannot change a score, and it amortises the ~0.66 s model load over the whole chunk.

Usage
    PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
    $PY neomod/pipeline/audit_digest2.py build-pairs
    $PY neomod/pipeline/audit_digest2.py run --select neo --chunk-index 0 --chunk-size 200
    $PY neomod/pipeline/audit_digest2.py run --pilot                 # the fixed pilot cohort
    $PY neomod/pipeline/audit_digest2.py combine --select neo
"""
from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
Path("/tmp/astropy").mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.utils import iers

iers.conf.auto_download = False

# ---------------------------------------------------------------------------- paths / constants
WORKDIR = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
NEOMD = WORKDIR / "neomod"
SORCHA_PARQUET = WORKDIR / "outputs/s3m_linking/case1/sorcha_comparison_case1_nbody_Vband.parquet"
BENCH_PARQUET = (WORKDIR / "outputs/phase2_benchmark_s3m_nbody_mag245"
                 / "benchmark_comparison_s3m_nbody_mag245.parquet")
DIGEST2_DIR = WORKDIR / "digest2"
DIGEST2_EXEC = DIGEST2_DIR / "digest2"
OBSCODE = "X05"
NIGHT = 61642

AUDIT_DIR = WORKDIR / "outputs/digest2_audit_night61642"
PAIRS_FILE = AUDIT_DIR / "digest2_audit_pairs.parquet"
RESULTS_DIR = AUDIT_DIR / "results"
RAW_DIR = AUDIT_DIR / "raw"

# V-magnitude conversion (same convention as pipeline/rescore_vdp_Vband.py):
#   m_V = m_f - (f - r color) + (H_V - H_r);   color_r = 0
PHYS_CSV = WORKDIR / "inputs" / "s3m_sorcha_phys.csv"
S3M_GLOB = str(NEOMD / "S3Mdata" / "S*.s3m")
COLOR_COLS = {"u": "u-r", "g": "g-r", "i": "i-r", "z": "z-r", "y": "y-r"}  # r -> 0

# digest2 config: no heading row, no RMS column, single NEO score column.
# `repeatable` fixes the RNG seed (digest2 default is `random` == srand(time(0)), which is what
# PRODUCTION used -> stored P_NEO_d2 carry ±1-3 pt Monte-Carlo jitter). We add `repeatable` and pin
# `--cpu 1` (below) so the audit's own scores are 100% deterministic/reproducible; the residual
# gap to the stored random-mode values is then either MC noise (small) or a real bug (large).
D2_CONFIG = "noheadings\nnorms\nrepeatable\nNEO\n"
D2_CPUS = "1"                      # pin threads: per-thread LCG seeds make multi-thread non-repeatable
MC_TOL = 0.03                      # |Δnorm| within this vs stored == digest2 MC noise, not a bug

VARIANTS = ("B", "S", "V")  # benchmark-native, sorcha-native, sorcha-Vcorrected


# ---------------------------------------------------------------------------- MPC formatting
def format_mpc80(desig12: str, mjd_utc: float, ra_deg: float, dec_deg: float, mag: float) -> str:
    """One 80-character MPC observation line.

    Copied VERBATIM from pipeline/sorcha_phase2.py::format_mpc80 so the audit's MPC lines are
    byte-identical to production (guarantees the S variant reproduces the stored P_NEO_d2). The
    hardcoded ' V ' band label is the very thing under audit and is preserved deliberately.
    """
    dt = Time(float(mjd_utc), format="mjd", scale="utc").utc.datetime
    day_frac = (
        dt.day
        + dt.hour / 24.0
        + dt.minute / 1440.0
        + (dt.second + dt.microsecond / 1e6) / 86400.0
    )
    ra_h = (float(ra_deg) % 360.0) / 15.0
    rah = int(ra_h)
    ram_f = (ra_h - rah) * 60.0
    ram = int(ram_f)
    ras = min((ram_f - ram) * 60.0, 59.99)
    dec_a = abs(float(np.clip(dec_deg, -89.99, 89.99)))
    sign = "+" if dec_deg >= 0 else "-"
    decd = int(dec_a)
    decm_f = (dec_a - decd) * 60.0
    decm = int(decm_f)
    decs = min((decm_f - decm) * 60.0, 59.9)
    mag_v = max(0.0, min(99.9, float(mag) if np.isfinite(mag) else 21.0))
    line = (
        f"{desig12[:12]:12s}"
        f"  C"
        f"{dt.year:04d} {dt.month:02d} {day_frac:08.5f}"
        f" {rah:02d} {ram:02d} {ras:05.2f} "
        f"{sign}{decd:02d} {decm:02d} {decs:04.1f}"
        f"          "
        f"{mag_v:4.1f} V      "
        f"{OBSCODE:3s}"
    )
    if len(line) != 80:
        raise AssertionError(f"MPC line has length {len(line)}: {line!r}")
    return line


def desig_for(variant: str, pair_index: int) -> tuple[str, str]:
    """(7-char echoed key, 12-char MPC field) for a variant + pair index."""
    key = f"{variant}{pair_index:06d}"        # e.g. 'S000123' — what digest2 echoes back
    return key, f"     {key}"                  # 5 blanks + 7 chars = MPC cols 1-12


# ---------------------------------------------------------------------------- V-mag conversion
# H_V (Johnson-V absolute mag): for matched objects we take it from the benchmark parquet's `H`
# column, which IS the S3M Johnson-V H (fixing_integrator.md §9.9: "H_v = the S3M H, which is
# Johnson V") — no census scan needed. Only the unmatched 1.17 pilot (a NEO) reads S0.s3m.

def _read_phys_colors(objids: set[str]) -> pd.DataFrame:
    """ObjID -> (u-r,g-r,i-r,z-r,y-r) LSST colours, from the 1.9 GB phys CSV.
    Chunked + filtered to `objids` so it never materialises the whole file."""
    usecols = ["ObjID"] + list(COLOR_COLS.values())
    parts = []
    for chunk in pd.read_csv(PHYS_CSV, usecols=usecols, dtype={"ObjID": str}, chunksize=2_000_000):
        sub = chunk[chunk.ObjID.isin(objids)]
        if len(sub):
            parts.append(sub)
    if not parts:
        return pd.DataFrame(columns=usecols).set_index("ObjID")
    return pd.concat(parts, ignore_index=True).drop_duplicates("ObjID").set_index("ObjID")


def _hv_from_s0(objids: set[str]) -> pd.Series:
    """ObjID -> H_V from the S3M NEO census only (S0.s3m, col0=OID, col8=H_V)."""
    parts = []
    for chunk in pd.read_csv(NEOMD / "S3Mdata" / "S0.s3m", sep=r"\s+", header=None,
                             usecols=[0, 8], names=["OID", "H_V"], dtype={0: str},
                             chunksize=500_000):
        sub = chunk[chunk.OID.isin(objids)]
        if len(sub):
            parts.append(sub)
    if not parts:
        return pd.Series(dtype=float)
    return pd.concat(parts, ignore_index=True).drop_duplicates("OID").set_index("OID")["H_V"].astype(float)


def compute_Vmags(df: pd.DataFrame, phys_colors: pd.DataFrame, hv: pd.Series) -> pd.DataFrame:
    """Add H_V, mag0_V, mag1_V: mag_i_V = mag_i - colour(filter_i) + (H_V - H_r), colour_r = 0."""
    oid = df.ObjID.astype(str).to_numpy()

    def colour_of(filters: pd.Series) -> np.ndarray:
        out = np.zeros(len(df), dtype=float)         # r -> 0
        f = filters.to_numpy()
        for band, col in COLOR_COLS.items():
            m = f == band
            if m.any():
                out[m] = phys_colors.reindex(oid[m])[col].to_numpy(float)
        return out

    H_V = hv.reindex(oid).to_numpy(float)
    cVr = H_V - df.H_r.to_numpy(float)               # (V - r) per object
    df = df.copy()
    df["H_V"] = H_V
    df["mag0_V"] = df.mag0.to_numpy(float) - colour_of(df.filter0) + cVr
    df["mag1_V"] = df.mag1.to_numpy(float) - colour_of(df.filter1) + cVr
    n_missing = int(np.isnan(df.mag0_V).sum())
    if n_missing:
        print(f"  WARNING: {n_missing:,} rows have no H_V/colour -> mag*_V NaN", flush=True)
    return df


# ---------------------------------------------------------------------------- build-pairs
def build_pairs(args: argparse.Namespace) -> None:
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)

    scols = ["ObjID", "night", "population", "ra0", "dec0", "mjd0_utc", "mag0", "filter0", "snr0",
             "ra1", "dec1", "mjd1_utc", "mag1", "filter1", "snr1", "mean_mag", "mean_mag_V", "H_r",
             "night_span_min", "mean_ra", "mean_dec", "vlam", "vbeta", "prob_map_file",
             "P_NEO_d2", "P_NEO_vdp_Vband", "mag_bin_label_Vband"]
    s = pd.read_parquet(SORCHA_PARQUET, columns=scols)
    s = s[s.night == NIGHT].copy()
    s["ObjID"] = s.ObjID.astype(str)

    bcols = ["s3m_objid", "population", "ra0", "dec0", "mjd0_utc", "mag0",
             "ra1", "dec1", "mjd1_utc", "mag1", "mean_mag", "H", "e", "vlam", "vbeta",
             "mean_ra", "mean_dec", "prob_map_file", "P_NEO_d2", "P_NEO_vdp"]
    b = pd.read_parquet(BENCH_PARQUET, columns=bcols)
    b["s3m_objid"] = b.s3m_objid.astype(str)
    b = b.drop_duplicates("s3m_objid").set_index("s3m_objid")

    s = s[s.ObjID.isin(b.index)].reset_index(drop=True)   # only matched objects need V-mags
    special = _load_special_hits()                        # 1.17 rows (full set, not night 61642)

    # ONE phys scan over matched ∪ special. H_V from benchmark `H` (= S3M Johnson-V H) for matched;
    # from the NEO census (S0.s3m) for the unmatched special object.
    phys_colors = _read_phys_colors(set(s.ObjID) | set(special.ObjID))
    hv = pd.concat([b["H"].astype(float), _hv_from_s0(set(special.ObjID))])
    hv = hv[~hv.index.duplicated()]
    s = compute_Vmags(s, phys_colors, hv)
    special = compute_Vmags(special, phys_colors, hv) if len(special) else special

    matched = s.reset_index(drop=True)
    matched["pair_index"] = np.arange(len(matched))
    bm = b.reindex(matched.ObjID.to_numpy())
    # NEO cohort = benchmark-population NEO (matches the notebook's 794); keep both labels since
    # 49 objects are labelled differently by the two pipelines (benchmark MBA vs Sorcha NEO/other).
    matched["population_bench"] = bm.population.to_numpy()
    matched["population_sorcha"] = matched.population.to_numpy()
    matched["is_neo"] = (matched.population_bench == "NEO").astype(int)
    print(f"matched night-{NIGHT} objects: {len(matched):,}  "
          f"(benchmark-NEO {int(matched.is_neo.sum()):,}, other {int((1 - matched.is_neo).sum()):,}; "
          f"Sorcha-NEO {int((matched.population_sorcha == 'NEO').sum()):,})", flush=True)

    rows = []
    for pi, (_, r) in enumerate(matched.iterrows()):
        br = bm.iloc[pi]
        common = dict(
            pair_index=pi, ObjID=r.ObjID, population=r.population_bench, is_neo=int(r.is_neo),
            population_bench=r.population_bench, population_sorcha=r.population_sorcha,
            night=NIGHT,
            # metadata for step-7 correlation (variant-independent)
            span_days_sorcha=r.night_span_min / 1440.0,
            span_days_bench=float(br.mjd1_utc - br.mjd0_utc),
            filter0=r.filter0, filter1=r.filter1, snr0=r.snr0, snr1=r.snr1,
            sky_sep_deg=float(np.hypot(
                (r.mean_ra - br.mean_ra) * np.cos(np.radians(r.mean_dec)),
                r.mean_dec - br.mean_dec)),
            dv_deg_day=float(np.hypot(r.vlam - br.vlam, r.vbeta - br.vbeta)),
            dmag_V=float(r.mean_mag_V - br.mean_mag),
            H_V=r.H_V, e=float(br.e),
            prob_map_sorcha=r.prob_map_file, prob_map_bench=br.prob_map_file,
            stored_vdp_sorcha=r.P_NEO_vdp_Vband, stored_vdp_bench=br.P_NEO_vdp,
        )
        specs = {
            # variant : (ra0,dec0,mjd0,mag0, ra1,dec1,mjd1,mag1, stored_expected, mag_source)
            "B": (br.ra0, br.dec0, br.mjd0_utc, br.mag0, br.ra1, br.dec1, br.mjd1_utc, br.mag1,
                  br.P_NEO_d2, "benchmark_synthetic_V"),
            "S": (r.ra0, r.dec0, r.mjd0_utc, r.mag0, r.ra1, r.dec1, r.mjd1_utc, r.mag1,
                  r.P_NEO_d2, "sorcha_raw_labelledV"),
            "V": (r.ra0, r.dec0, r.mjd0_utc, r.mag0_V, r.ra1, r.dec1, r.mjd1_utc, r.mag1_V,
                  np.nan, "sorcha_converted_V"),
        }
        for v, (ra0, dec0, mjd0, mag0, ra1, dec1, mjd1, mag1, exp, msrc) in specs.items():
            key, d12 = desig_for(v, pi)
            rows.append({**common, "variant": v, "desig": key, "desig12": d12,
                         "ra0": ra0, "dec0": dec0, "mjd0_utc": mjd0, "mag0": mag0,
                         "ra1": ra1, "dec1": dec1, "mjd1_utc": mjd1, "mag1": mag1,
                         "mag_source": msrc, "stored_score_expected": exp,
                         "pilot": False, "pilot_reason": ""})

    pairs = pd.DataFrame(rows)
    _flag_pilot(pairs, matched, bm)
    _append_special_pilots(pairs, special)  # the 1.17 object (not in the matched night-61642 set)

    pairs.to_parquet(PAIRS_FILE, index=False)
    n_obj = pairs.pair_index.nunique()
    print(f"wrote {len(pairs):,} variant-rows for {n_obj:,} objects -> {PAIRS_FILE}", flush=True)
    print(f"  pilot objects flagged: {pairs[pairs.pilot].pair_index.nunique()}", flush=True)
    for sel, n in [("neo", pairs[pairs.is_neo == 1].pair_index.nunique()),
                   ("all", pairs.pair_index.nunique())]:
        print(f"  select={sel}: {n} objects -> "
              f"{int(np.ceil(n / args.chunk_size))} chunks at chunk-size {args.chunk_size}", flush=True)


def _flag_pilot(pairs: pd.DataFrame, matched: pd.DataFrame, bm: pd.DataFrame) -> None:
    """Flag a small deterministic pilot cohort (stored-agree, biggest disagreements, 1->0 outlier)."""
    m = matched.copy()
    m["d2_b"] = bm.P_NEO_d2.to_numpy()
    m["d2_s"] = m.P_NEO_d2.to_numpy()
    m["dP_d2"] = m.d2_s - m.d2_b
    picks: dict[int, str] = {}
    # biggest |dP_d2| disagreements (includes the digest2 1->0 outlier S00002XMa)
    for pi in m.reindex(m.dP_d2.abs().sort_values(ascending=False).index).pair_index.head(4):
        picks[int(pi)] = "big_d2_disagreement"
    # stored-agree at 1 and near 0
    for pi in m[m.dP_d2.abs() < 0.01].pair_index.head(3):
        picks.setdefault(int(pi), "stored_agree")
    # exact-zero Sorcha d2 (silent-zero suspects)
    for pi in m[m.d2_s == 0].pair_index.head(3):
        picks.setdefault(int(pi), "sorcha_d2_exact_zero")
    for pi, reason in picks.items():
        sel = pairs.pair_index == pi
        pairs.loc[sel, "pilot"] = True
        pairs.loc[sel, "pilot_reason"] = reason


def _load_special_hits() -> pd.DataFrame:
    """Full-set Sorcha rows with P_NEO_d2 > 1.0 (the 1.17 out-of-range object). ObjID as str.
    V-mags are added later by the single shared phys scan in build_pairs."""
    s = pd.read_parquet(SORCHA_PARQUET,
                        columns=["ObjID", "night", "population", "ra0", "dec0", "mjd0_utc", "mag0",
                                 "filter0", "snr0", "ra1", "dec1", "mjd1_utc", "mag1", "filter1",
                                 "snr1", "mean_mag", "mean_mag_V", "H_r", "night_span_min",
                                 "prob_map_file", "P_NEO_d2", "P_NEO_vdp_Vband"])
    hit = s[s.P_NEO_d2 > 1.0].copy()
    hit["ObjID"] = hit.ObjID.astype(str)
    return hit


def _append_special_pilots(pairs: pd.DataFrame, special: pd.DataFrame) -> None:
    """Add the P_NEO_d2 = 1.17 object as an S/V pilot pair (V-mags already computed in build_pairs).

    No benchmark match required — the point is to capture digest2's raw output and confirm it
    emits 117 for this tracklet (the source of the 1.17 = 117/100 out-of-range value).
    """
    if special is None or special.empty:
        return
    pi0 = 900000
    extra = []
    for k, (_, r) in enumerate(special.iterrows()):
        pi = pi0 + k
        common = dict(pair_index=pi, ObjID=r.ObjID, population=r.population, is_neo=0,
                      population_bench=None, population_sorcha=r.population,
                      night=int(r.night), span_days_sorcha=r.night_span_min / 1440.0, span_days_bench=np.nan,
                      filter0=r.filter0, filter1=r.filter1, snr0=r.snr0, snr1=r.snr1,
                      sky_sep_deg=np.nan, dv_deg_day=np.nan, dmag_V=np.nan, H_V=r.H_V, e=np.nan,
                      prob_map_sorcha=r.prob_map_file, prob_map_bench=None,
                      stored_vdp_sorcha=r.P_NEO_vdp_Vband, stored_vdp_bench=np.nan)
        for v, mag0, mag1, exp, msrc in [
                ("S", r.mag0, r.mag1, r.P_NEO_d2, "sorcha_raw_labelledV"),
                ("V", r.mag0_V, r.mag1_V, np.nan, "sorcha_converted_V")]:
            key, d12 = desig_for(v, pi)
            extra.append({**common, "variant": v, "desig": key, "desig12": d12,
                          "ra0": r.ra0, "dec0": r.dec0, "mjd0_utc": r.mjd0_utc, "mag0": mag0,
                          "ra1": r.ra1, "dec1": r.dec1, "mjd1_utc": r.mjd1_utc, "mag1": mag1,
                          "mag_source": msrc, "stored_score_expected": exp,
                          "pilot": True, "pilot_reason": "d2_gt_1_(1.17)"})
    if extra:
        pairs2 = pd.concat([pairs, pd.DataFrame(extra)], ignore_index=True)
        pairs.drop(pairs.index, inplace=True)
        for c in pairs2.columns:
            pairs[c] = pairs2[c].values


# ---------------------------------------------------------------------------- run (one chunk)
def _run_digest2(obs_lines: list[str], tag: str) -> tuple[list[str], str, int, float]:
    """Invoke digest2 once on obs_lines. Returns (stdout_lines, stderr, returncode, runtime_s).
    Saves the exact .obs sent and raw .out/.err under RAW_DIR/. Does NOT raise on nonzero rc."""
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    obs_path = RAW_DIR / f"{tag}.obs"
    obs_path.write_text("\n".join(obs_lines) + "\n")
    with tempfile.NamedTemporaryFile(mode="w", suffix=".config", delete=False) as cfg:
        cfg.write(D2_CONFIG)
        cfg_path = cfg.name
    t0 = time.time()
    try:
        proc = subprocess.run([str(DIGEST2_EXEC), "-p", str(DIGEST2_DIR), "-c", cfg_path,
                               "--cpu", D2_CPUS, str(obs_path)],
                              capture_output=True, text=True, timeout=3600)
        rc, out, err = proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired as exc:
        rc, out, err = 124, exc.stdout or "", f"TIMEOUT\n{exc.stderr or ''}"
    finally:
        os.unlink(cfg_path)
    dt = time.time() - t0
    (RAW_DIR / f"{tag}.out").write_text(out)
    (RAW_DIR / f"{tag}.err").write_text(err)
    return out.splitlines(), err, rc, dt


def _strict_parse(out_lines: list[str]) -> dict[str, list[int]]:
    """desig -> list of raw int scores. Keeps ALL occurrences so duplicates are detectable.
    Malformed second tokens are dropped from the value list (surfaced as parse_status)."""
    d: dict[str, list[int]] = {}
    for line in out_lines:
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        try:
            val = int(parts[1])
        except ValueError:
            d.setdefault(parts[0], [])   # seen but unparseable -> empty list => 'malformed'
            continue
        d.setdefault(parts[0], []).append(val)
    return d


SPECIAL_PI0 = 900000   # special (non-cohort) pilots like the 1.17 live at pair_index >= this


def _select_chunk(pairs: pd.DataFrame, select: str, chunk_index: int, chunk_size: int) -> pd.DataFrame:
    # Cohort = all matched night-61642 objects (pair_index < SPECIAL_PI0). The `pilot` flag is only
    # a marker for the pilot PREVIEW — flagged cohort objects must still appear in the full run;
    # only the special out-of-cohort pilots (the 1.17, night 61081) are excluded here.
    reg = pairs[pairs.pair_index < SPECIAL_PI0]
    if select == "neo":
        obj = reg[reg.is_neo == 1]
    else:
        obj = reg
    order = np.sort(obj.pair_index.unique())
    lo, hi = chunk_index * chunk_size, (chunk_index + 1) * chunk_size
    keep = set(order[lo:hi].tolist())
    return pairs[pairs.pair_index.isin(keep)].copy()


def run(args: argparse.Namespace) -> None:
    pairs = pd.read_parquet(PAIRS_FILE)
    if args.pilot:
        sub = pairs[pairs.pilot].copy()
        tagbase, select = "pilot", "pilot"
    else:
        sub = _select_chunk(pairs, args.select, args.chunk_index, args.chunk_size)
        tagbase, select = f"{args.select}_chunk{args.chunk_index:04d}", args.select
    if sub.empty:
        print(f"no pairs for {tagbase}; nothing to do", flush=True)
        return

    variants = tuple(args.variants.split(",")) if args.variants else VARIANTS
    results = []
    for v in variants:
        vs = sub[sub.variant == v].sort_values("pair_index")
        if vs.empty:
            continue
        obs_lines, order = [], []
        for _, r in vs.iterrows():
            obs_lines.append(format_mpc80(r.desig12, r.mjd0_utc, r.ra0, r.dec0, r.mag0))
            obs_lines.append(format_mpc80(r.desig12, r.mjd1_utc, r.ra1, r.dec1, r.mag1))
            order.append(r.desig)
        out_lines, err, rc, dt = _run_digest2(obs_lines, f"{tagbase}_{v}")
        parsed = _strict_parse(out_lines)
        print(f"  {tagbase} variant={v}: n={len(vs)} rc={rc} runtime={dt:.2f}s "
              f"out_desigs={len(parsed)}", flush=True)
        for _, r in vs.iterrows():
            vals = parsed.get(r.desig, None)
            if vals is None:
                status, raw = "missing", np.nan
            elif len(vals) == 0:
                status, raw = "malformed", np.nan
            elif len(vals) > 1:
                status, raw = "duplicate", np.nan
            else:
                status, raw = "ok", vals[0]
            norm = raw / 100.0 if status == "ok" else np.nan
            out_of_range = bool(status == "ok" and not (0 <= raw <= 100))
            results.append({
                "pair_index": int(r.pair_index), "ObjID": r.ObjID, "variant": v,
                "population": r.population, "is_neo": int(r.is_neo), "night": int(r.night),
                "desig": r.desig, "mag_source": r.mag_source,
                "raw_int": raw, "norm": norm, "parse_status": status,
                "out_of_range": out_of_range, "n_output_occurrences": 0 if vals is None else len(vals),
                "returncode": rc, "runtime_s": round(dt, 3),
                "mpc_line0": obs_lines[2 * order.index(r.desig)],
                "mpc_line1": obs_lines[2 * order.index(r.desig) + 1],
                "stored_score_expected": r.stored_score_expected,
                "filter0": r.filter0, "filter1": r.filter1, "snr0": r.snr0, "snr1": r.snr1,
                "span_days_sorcha": r.span_days_sorcha, "span_days_bench": r.span_days_bench,
                "sky_sep_deg": r.sky_sep_deg, "dv_deg_day": r.dv_deg_day, "dmag_V": r.dmag_V,
                "H_V": r.H_V, "e": r.e, "prob_map_sorcha": r.prob_map_sorcha,
                "prob_map_bench": r.prob_map_bench, "stored_vdp_sorcha": r.stored_vdp_sorcha,
                "pilot": bool(r.pilot), "pilot_reason": r.pilot_reason,
            })

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = RESULTS_DIR / f"{tagbase}.parquet"
    pd.DataFrame(results).to_parquet(out, index=False)
    print(f"wrote {len(results):,} result rows -> {out}", flush=True)

    if args.verify_determinism:
        _verify_determinism(sub, variants)


def _verify_determinism(sub: pd.DataFrame, variants) -> None:
    """Re-run each variant a second time; assert identical raw scores (digest2 is deterministic)."""
    ok = True
    for v in variants:
        vs = sub[sub.variant == v].sort_values("pair_index")
        if vs.empty:
            continue
        runs = []
        for rep in range(2):
            obs = []
            for _, r in vs.iterrows():
                obs.append(format_mpc80(r.desig12, r.mjd0_utc, r.ra0, r.dec0, r.mag0))
                obs.append(format_mpc80(r.desig12, r.mjd1_utc, r.ra1, r.dec1, r.mag1))
            out_lines, _, _, _ = _run_digest2(obs, f"determinism_{v}_rep{rep}")
            runs.append(_strict_parse(out_lines))
        same = all(runs[0].get(d) == runs[1].get(d) for d in vs.desig)
        print(f"  determinism variant={v}: {'IDENTICAL' if same else 'MISMATCH !!'}", flush=True)
        ok = ok and same
    print(f"determinism gate: {'PASS' if ok else 'FAIL'}", flush=True)


# ---------------------------------------------------------------------------- combine
def combine(args: argparse.Namespace) -> None:
    if args.pilot:
        paths = [RESULTS_DIR / "pilot.parquet"]
    else:
        paths = sorted(RESULTS_DIR.glob(f"{args.select}_chunk*.parquet"))
    paths = [p for p in paths if p.exists()]
    if not paths:
        raise FileNotFoundError(f"no result shards for select={args.select} in {RESULTS_DIR}")
    res = pd.concat([pd.read_parquet(p) for p in paths], ignore_index=True)
    tag = "pilot" if args.pilot else args.select

    # ---- parse log (per variant) ----
    log_rows = []
    for v, g in res.groupby("variant"):
        log_rows.append({
            "select": tag, "variant": v, "n": len(g),
            "ok": int((g.parse_status == "ok").sum()),
            "missing": int((g.parse_status == "missing").sum()),
            "duplicate": int((g.parse_status == "duplicate").sum()),
            "malformed": int((g.parse_status == "malformed").sum()),
            "out_of_range": int(g.out_of_range.sum()),
            "returncodes": sorted(g.returncode.unique().tolist()),
        })
    parse_log = pd.DataFrame(log_rows)
    parse_log.to_csv(AUDIT_DIR / f"digest2_audit_parse_log_{tag}.csv", index=False)

    # ---- reproduction gate: repeatable-mode rerun vs stored (random-mode production) ----
    # Production ran digest2 in `random` mode, so exact reproduction is IMPOSSIBLE (MC jitter).
    # Classify: exact | within MC noise (|Δ|<=MC_TOL) | real discrepancy (a parser/silent-zero bug).
    print(f"\n=== reproduction: repeatable rerun vs stored random-mode score, select={tag} "
          f"(MC_TOL={MC_TOL}) ===", flush=True)
    for v in ("B", "S"):
        g = res[(res.variant == v) & (res.parse_status == "ok") & res.stored_score_expected.notna()]
        if g.empty:
            continue
        diff = (g.norm - g.stored_score_expected).abs()
        exact = int((diff < 1e-9).sum())
        within = int((diff <= MC_TOL).sum())
        real = g[diff > MC_TOL]
        print(f"  {v}: n={len(g)}  exact={exact}  within_MC_noise(<= {MC_TOL})={within}  "
              f"real_discrepancy(> {MC_TOL})={len(real)}  max|Δ|={diff.max():.3f}", flush=True)
        if len(real):
            show = real.assign(dabs=diff[diff > MC_TOL]).sort_values("dabs", ascending=False)
            for _, rr in show.head(8).iterrows():
                print(f"      {rr.ObjID}: stored={rr.stored_score_expected:.2f} rerun={rr.norm:.2f} "
                      f"Δ={rr.norm - rr.stored_score_expected:+.2f}", flush=True)

    # ---- impossible stored values (norm > 1, e.g. the 1.17): can't come from one digest2 draw ----
    imposs = res[(res.variant == "S") & (res.stored_score_expected > 1.0)].drop_duplicates("pair_index")
    if len(imposs):
        print(f"\n=== impossible stored values (stored norm > 1) ===", flush=True)
        for _, rr in imposs.iterrows():
            print(f"  {rr.ObjID}: stored={rr.stored_score_expected:.2f}  rerun={rr.norm:.2f} "
                  f"(digest2 max is 1.00 -> stored is a production parse corruption)", flush=True)

    # ---- silent-zero unmasking: stored==0 but rerun>0 ----
    s_ok = res[(res.variant == "S") & (res.parse_status == "ok")]
    silent = s_ok[(s_ok.stored_score_expected == 0) & (s_ok.norm > 0)]
    genuine0 = s_ok[(s_ok.stored_score_expected == 0) & (s_ok.norm == 0)]
    print(f"\n=== silent-zero analysis (S variant) ===", flush=True)
    print(f"  stored==0 total: {int((s_ok.stored_score_expected == 0).sum())}  "
          f"-> genuine (rerun 0): {len(genuine0)}   silently-defaulted (rerun>0): {len(silent)}", flush=True)

    # ---- outliers ----
    outliers = res[(res.out_of_range) | (res.parse_status != "ok")].copy()
    if len(silent):
        silent2 = silent.copy(); silent2["flag"] = "silent_zero_unmasked"; outliers = pd.concat([outliers, silent2])
    outliers.to_csv(AUDIT_DIR / f"digest2_audit_outliers_{tag}.csv", index=False)

    # ---- wide per-object table + band effect (S vs V) ----
    wide = res[res.parse_status == "ok"].pivot_table(index="pair_index", columns="variant",
                                                     values="norm").rename(
        columns={"B": "P_d2_B", "S": "P_d2_S", "V": "P_d2_V"})
    meta = res.drop_duplicates("pair_index").set_index("pair_index")[
        ["ObjID", "population", "is_neo", "stored_vdp_sorcha", "dmag_V", "sky_sep_deg", "dv_deg_day"]]
    wide = meta.join(wide)
    wide.to_parquet(AUDIT_DIR / f"digest2_audit_results_{tag}.parquet")

    if {"P_d2_S", "P_d2_V"}.issubset(wide.columns):
        band = (wide.P_d2_V - wide.P_d2_S).dropna()
        print(f"\n=== band effect S->V (Sorcha raw-labelledV -> converted-V), select={tag} ===", flush=True)
        print(f"  N={len(band)}  median dP={band.median():+.4f}  mean={band.mean():+.4f}  "
              f"|dP|>0.1: {int((band.abs() > 0.1).sum())} ({100 * (band.abs() > 0.1).mean():.1f}%)", flush=True)
    if {"P_d2_B", "P_d2_S"}.issubset(wide.columns):
        bs = (wide.P_d2_S - wide.P_d2_B).dropna()
        print(f"=== benchmark-native vs Sorcha-native (B->S), select={tag} ===", flush=True)
        print(f"  N={len(bs)}  median dP={bs.median():+.4f}  |dP|>0.1: "
              f"{int((bs.abs() > 0.1).sum())} ({100 * (bs.abs() > 0.1).mean():.1f}%)", flush=True)

    print(f"\nparse log:\n{parse_log.to_string(index=False)}", flush=True)
    print(f"\nwrote: digest2_audit_results_{tag}.parquet, "
          f"digest2_audit_parse_log_{tag}.csv, digest2_audit_outliers_{tag}.csv", flush=True)


# ---------------------------------------------------------------------------- report (plots)
def report(args: argparse.Namespace) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tag = args.select
    paths = sorted(RESULTS_DIR.glob(f"{tag}_chunk*.parquet"))
    res = pd.concat([pd.read_parquet(p) for p in paths], ignore_index=True)
    ok = res[res.parse_status == "ok"]
    norm = ok.pivot_table(index="pair_index", columns="variant", values="norm")
    stored = ok.pivot_table(index="pair_index", columns="variant", values="stored_score_expected")
    meta = ok.drop_duplicates("pair_index").set_index("pair_index")[["ObjID", "is_neo", "dmag_V"]]
    w = meta.join(norm.rename(columns={"B": "S_B", "S": "S_S", "V": "S_V"}))
    w["stored_S"] = stored["S"]; w["stored_B"] = stored["B"]

    # 1) reproduction scatter: stored (random-mode production) vs audit rerun (repeatable)
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.6))
    for a, (col, stcol, name) in zip(ax, [("S_B", "stored_B", "benchmark-native (B)"),
                                          ("S_S", "stored_S", "Sorcha-native (S)")]):
        d = w.dropna(subset=[col, stcol])
        diff = (d[col] - d[stcol]).abs()
        big = diff > 0.05
        a.scatter(d[stcol][~big], d[col][~big], s=14, alpha=0.5, color="tab:blue", label="within noise")
        a.scatter(d[stcol][big], d[col][big], s=34, color="tab:red", edgecolors="k",
                  linewidths=0.4, label=f"|Δ|>0.05 (n={int(big.sum())})", zorder=5)
        a.plot([0, 1], [0, 1], "k--", lw=1)
        a.set_xlabel(f"stored P_NEO_d2 (production, random mode)")
        a.set_ylabel("audit rerun (repeatable)")
        a.set_title(name, fontsize=10); a.set_aspect("equal")
        a.set_xlim(-0.03, 1.03); a.set_ylim(-0.03, 1.03); a.grid(alpha=0.25); a.legend(fontsize=8)
    fig.suptitle(f"digest2 audit — stored vs faithful rerun ({tag}, N={len(w)})", fontsize=12)
    fig.tight_layout(); fig.savefig(AUDIT_DIR / f"digest2_audit_reproduction_{tag}.png", dpi=150,
                                    bbox_inches="tight"); plt.close(fig)

    # 2) native vs V-corrected (the band effect), coloured by |Δmag_V|
    fig, a = plt.subplots(figsize=(6.2, 5.8))
    d = w.dropna(subset=["S_S", "S_V"])
    sc = a.scatter(d.S_S, d.S_V, s=18, c=d.dmag_V.abs(), cmap="viridis", alpha=0.8, vmin=0, vmax=1.2)
    a.plot([0, 1], [0, 1], "k--", lw=1)
    a.set_xlabel("Sorcha-native P_d2 (raw mag, labelled V)")
    a.set_ylabel("Sorcha V-corrected P_d2")
    a.set_title(f"Band effect S→V ({tag}, N={len(d)})\nmedian Δ={ (d.S_V-d.S_S).median():+.3f}, "
                f"{100*((d.S_V-d.S_S).abs()>0.1).mean():.1f}% shift >0.1", fontsize=10)
    a.set_aspect("equal"); a.set_xlim(-0.03, 1.03); a.set_ylim(-0.03, 1.03); a.grid(alpha=0.25)
    fig.colorbar(sc, label="|Δmag_V| (band error, mag)")
    fig.tight_layout(); fig.savefig(AUDIT_DIR / f"digest2_audit_native_vs_Vcorrected_{tag}.png",
                                    dpi=150, bbox_inches="tight"); plt.close(fig)

    # 3) summary table
    def frac(x, t): return f"{int((x > t).sum())} ({100 * (x > t).mean():.1f}%)"
    dS = (w.S_S - w.stored_S).abs().dropna()
    dB = (w.S_B - w.stored_B).abs().dropna()
    band = (w.S_V - w.S_S).abs().dropna()
    bs = (w.S_S - w.S_B).abs().dropna()
    summ = pd.DataFrame([
        {"metric": "N objects", "value": len(w)},
        {"metric": "S reproduces stored within 0.03", "value": frac(-dS, -0.03).replace("(", "→ mismatch (")},
        {"metric": "S |Δ vs stored|>0.05 (candidate bugs)", "value": frac(dS, 0.05)},
        {"metric": "B |Δ vs stored|>0.05", "value": frac(dB, 0.05)},
        {"metric": "band S→V |Δ|>0.1", "value": frac(band, 0.1)},
        {"metric": "benchmark-vs-Sorcha native |Δ|>0.1", "value": frac(bs, 0.1)},
    ])
    summ.to_csv(AUDIT_DIR / f"digest2_audit_summary_{tag}.csv", index=False)
    print(f"wrote plots + digest2_audit_summary_{tag}.csv to {AUDIT_DIR}", flush=True)
    print(summ.to_string(index=False), flush=True)


# ---------------------------------------------------------------------------- cli
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    bp = sub.add_parser("build-pairs")
    bp.add_argument("--chunk-size", type=int, default=200)
    bp.set_defaults(func=build_pairs)

    rp = sub.add_parser("run")
    rp.add_argument("--select", choices=["neo", "all"], default="neo")
    rp.add_argument("--chunk-index", type=int, default=0)
    rp.add_argument("--chunk-size", type=int, default=200)
    rp.add_argument("--variants", default="", help="comma list, default B,S,V")
    rp.add_argument("--pilot", action="store_true", help="run the flagged pilot cohort only")
    rp.add_argument("--verify-determinism", action="store_true")
    rp.set_defaults(func=run)

    cp = sub.add_parser("combine")
    cp.add_argument("--select", choices=["neo", "all"], default="neo")
    cp.add_argument("--pilot", action="store_true")
    cp.set_defaults(func=combine)

    pp = sub.add_parser("report")
    pp.add_argument("--select", choices=["neo", "all"], default="neo")
    pp.set_defaults(func=report)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
