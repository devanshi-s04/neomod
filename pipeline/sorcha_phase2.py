#!/usr/bin/env python3
"""Phase 2 scoring for Sorcha nightly tracklets.

This script consumes the Phase 1 ``outputs/tracklets/tracklets_*.parquet``
files, scores each tracklet with VDP probability maps, optionally runs digest2
on the same two observed positions, and combines scored shards.
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
Path("/tmp/astropy").mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from astropy.time import Time
from astropy.utils import iers


iers.conf.auto_download = False

_NEOMOD = Path(__file__).resolve().parent.parent   # neomod/pipeline -> neomod/
ROOT = _NEOMOD.parent                               # neomod/ -> sorcha/ (data root)
NEOMOD_SRC = _NEOMOD / "src"
ADAM_STUB = _NEOMOD / "adam_core_stub"
if ADAM_STUB.exists():
    sys.path.insert(0, str(ADAM_STUB))
sys.path.insert(0, str(NEOMOD_SRC))

import velocity_density_pipeline_gmm as vdp  # noqa: E402
# ^ Use the GMM pipeline module: it generated the v5 grid maps AND carries the
#   support-count mask (support_mask_min). The original velocity_density_pipeline
#   has no such param, so score-vdp --support-mask-min crashed against it. Its
#   ProbMapSet.score_observation/from_npz are drop-in for scoring (unmasked AUC
#   matches the original module; verified on the 707k subsample).


OBSCODE = "X05"
DEFAULT_TRACKLET_DIR = Path("outputs/tracklets")
DEFAULT_PROB_MAPS_DIR = Path("prob_maps")
DEFAULT_WORK_DIR = Path("outputs/phase2")
VDP_COLUMNS = [
    "mean_ra",
    "mean_dec",
    "mean_dra",
    "mean_ddec",
    "mean_mag",
    "prob_map_file",
]
DIGEST2_COLUMNS = [
    "ra0",
    "dec0",
    "mjd0_utc",
    "mag0",
    "ra1",
    "dec1",
    "mjd1_utc",
    "mag1",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Score Sorcha Phase 1 tracklets.")
    parser.add_argument(
        "command",
        choices=["audit", "score-vdp", "sample", "run-digest2", "combine"],
        help="Phase 2 step to run.",
    )
    parser.add_argument("--tracklet-dir", type=Path, default=DEFAULT_TRACKLET_DIR)
    parser.add_argument("--prob-maps-dir", type=Path, default=DEFAULT_PROB_MAPS_DIR)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--limit-files", type=int, help="Debug limit for input files.")
    parser.add_argument("--shard-index", type=int, help="Run only one computed shard index.")
    parser.add_argument("--shard-start", type=int, help="First computed shard index to run.")
    parser.add_argument("--shard-stop", type=int, help="Stop before this computed shard index.")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--no-nearest-dist-mask",
        action="store_true",
        help="Disable the nearest-clone-distance mask when loading probability maps. "
             "GMM maps do not need this mask (density naturally → 0 far from data). "
             "Remove this flag to re-enable the mask.",
    )
    parser.add_argument(
        "--support-mask-min",
        type=float,
        default=None,
        help="If set, zero each non-smoothed population's density where its in-cell "
             "clone support_count < this value, before computing P(NEO). Cuts the kNN "
             "estimator's bleed of the dense MBA core into zero-support NEO-wing cells "
             "(diagnosed 2026-06-22). Try 1.",
    )
    parser.add_argument(
        "--subsample-file",
        type=Path,
        default=DEFAULT_WORK_DIR / "sorcha_subsample.parquet",
        help="Subsample parquet written by 'sample' and read by 'run-digest2'.",
    )
    parser.add_argument(
        "--n-sample-nonneo",
        type=int,
        default=500_000,
        help="Number of non-NEO tracklets to include in the digest2 subsample.",
    )
    parser.add_argument(
        "--sample-seed",
        type=int,
        default=42,
        help="Random seed for subsample reproducibility.",
    )
    parser.add_argument("--digest2-dir", type=Path, help="Directory containing digest2.")
    parser.add_argument("--digest2-exec", type=Path, help="Path to digest2 executable.")
    parser.add_argument("--digest2-chunk-tracklets", type=int, default=5_000)
    parser.add_argument("--digest2-timeout-sec", type=int, default=3_600)
    parser.add_argument("--d2-row-start", type=int, default=None,
                        help="First row index within the source parquet (for Slurm array tasks).")
    parser.add_argument("--d2-row-stop", type=int, default=None,
                        help="Exclusive last row index within the source parquet.")
    parser.add_argument(
        "--outfile",
        type=Path,
        default=DEFAULT_WORK_DIR / "sorcha_comparison.parquet",
    )
    return parser.parse_args()


def tracklet_paths(tracklet_dir: Path, limit_files: int | None = None) -> list[Path]:
    paths = sorted(tracklet_dir.glob("tracklets_*.parquet"))
    if limit_files is not None:
        paths = paths[:limit_files]
    if not paths:
        raise FileNotFoundError(f"No tracklets_*.parquet files found in {tracklet_dir}")
    return paths


def vdp_shard_dir(work_dir: Path) -> Path:
    return work_dir / "vdp_shards"


def digest2_shard_dir(work_dir: Path) -> Path:
    return work_dir / "digest2_shards"


def batched(seq: list[Path], n: int):
    if n <= 0:
        raise ValueError("--batch-size must be positive")
    for i in range(0, len(seq), n):
        yield i // n, seq[i : i + n]


def select_shards(shards: list[tuple[int, list[Path]]], args: argparse.Namespace):
    if args.shard_index is not None:
        return [(idx, batch) for idx, batch in shards if idx == args.shard_index]
    start = 0 if args.shard_start is None else args.shard_start
    stop = math.inf if args.shard_stop is None else args.shard_stop
    return [(idx, batch) for idx, batch in shards if idx >= start and idx < stop]


def read_parquet_many(paths: list[Path], columns: list[str] | None = None) -> pd.DataFrame:
    frames = [pd.read_parquet(path, columns=columns) for path in paths]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def audit(args: argparse.Namespace) -> None:
    paths = tracklet_paths(args.tracklet_dir, args.limit_files)
    total_rows = 0
    bad = []
    pop_counts: Counter[str] = Counter()
    map_counts: Counter[str] = Counter()

    for i, path in enumerate(paths, start=1):
        try:
            meta = pq.ParquetFile(path).metadata
            total_rows += meta.num_rows
            df = pd.read_parquet(path, columns=["population", "prob_map"])
            pop_counts.update(df["population"].dropna().astype(str).value_counts().to_dict())
            map_counts.update(df["prob_map"].dropna().astype(str).value_counts().to_dict())
        except Exception as exc:
            bad.append((str(path), repr(exc)))
        if i % 1000 == 0:
            print(f"checked={i:,} rows={total_rows:,}", flush=True)

    print(f"files={len(paths):,}")
    print(f"rows={total_rows:,}")
    print(f"bad_files={len(bad):,}")
    if bad:
        for item in bad[:20]:
            print(f"BAD {item[0]} {item[1]}")
    print("population_counts:")
    for key, val in sorted(pop_counts.items()):
        print(f"  {key}: {val:,}")
    print("prob_map_counts:")
    for key, val in sorted(map_counts.items()):
        print(f"  {key}: {val:,}")


def load_prob_map(prob_maps_dir: Path, name: str, cache: dict[str, vdp.ProbMapSet],
                  mask_radius_deg_per_day: float = vdp.DEFAULT_MASK_RADIUS_DEG_PER_DAY,
                  support_mask_min=None):
    if name not in cache:
        path = prob_maps_dir / name
        if not path.exists():
            cache[name] = None  # map absent from this maps dir; tracklets will be NaN
        else:
            print(f"  loading {path}", flush=True)
            cache[name] = vdp.ProbMapSet.from_npz(path,
                mask_radius_deg_per_day=mask_radius_deg_per_day,
                support_mask_min=support_mask_min)
    return cache[name]


def score_vdp_frame(df: pd.DataFrame, prob_maps_dir: Path, cache: dict[str, vdp.ProbMapSet],
                    mask_radius_deg_per_day: float = vdp.DEFAULT_MASK_RADIUS_DEG_PER_DAY,
                    support_mask_min=None) -> pd.DataFrame:
    df = df.copy()
    df["P_NEO_vdp"] = np.nan
    df["vlam"] = np.nan
    df["vbeta"] = np.nan
    df["mag_bin_label"] = None

    for map_file, idx in df.groupby("prob_map_file", dropna=True).groups.items():
        # Skip the legacy fixed-position maps (off-axis or single-epoch).
        # Maps we score against: the monthly "antisun" maps (v3.3/GMM run) OR the
        # antisun-relative "grid" maps (v5.0 full-sky grid). Legacy fixed-position
        # maps (e.g. lon229 NEOCP) have neither token and are skipped.
        mf = str(map_file)
        if "antisun" not in mf and "grid" not in mf:
            continue
        pms = load_prob_map(prob_maps_dir, str(map_file), cache,
                            mask_radius_deg_per_day=mask_radius_deg_per_day,
                            support_mask_min=support_mask_min)
        if pms is None:
            cache.pop(mf, None)
            continue
        sub = df.loc[idx]
        scored = pms.score_observation(
            ra_deg=sub["mean_ra"].to_numpy(float),
            dec_deg=sub["mean_dec"].to_numpy(float),
            dra_deg_day=sub["mean_dra"].to_numpy(float),
            ddec_deg_day=sub["mean_ddec"].to_numpy(float),
            mag_app=sub["mean_mag"].to_numpy(float),
            scorer=None,
            return_intermediates=True,
        )
        probs = scored["probs"]
        df.loc[idx, "P_NEO_vdp"] = probs.get("NEO", np.zeros(len(sub)))
        df.loc[idx, "vlam"] = scored["vlam"]
        df.loc[idx, "vbeta"] = scored["vbeta"]
        df.loc[idx, "mag_bin_label"] = scored["bin_labels"]
        # Free the map immediately: groupby visits each prob_map_file exactly
        # once per frame, so we never need to keep it. Critical for the 667-map
        # grid — a shard touches hundreds of maps; caching them all OOMs (the
        # 24-monthly-map run never hit this). Bounds RAM to ~1 map + the frame.
        cache.pop(mf, None)
        del pms, scored, sub
    return df


def score_vdp(args: argparse.Namespace) -> None:
    paths = tracklet_paths(args.tracklet_dir, args.limit_files)
    outdir = vdp_shard_dir(args.work_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    cache: dict[str, vdp.ProbMapSet] = {}
    mask_radius = np.inf if getattr(args, "no_nearest_dist_mask", False) else vdp.DEFAULT_MASK_RADIUS_DEG_PER_DAY
    if mask_radius == np.inf:
        print("nearest-dist mask DISABLED (--no-nearest-dist-mask)", flush=True)
    support_mask_min = getattr(args, "support_mask_min", None)
    if support_mask_min is not None:
        print(f"support-count mask ENABLED: min={support_mask_min} clones/cell "
              f"(non-smoothed populations)", flush=True)

    needed = None
    shards = select_shards(list(batched(paths, args.batch_size)), args)
    if not shards:
        raise IndexError("No VDP shards selected")
    for shard_idx, batch in shards:
        outfile = outdir / f"vdp_{shard_idx:05d}.parquet"
        if outfile.exists() and not args.overwrite:
            print(f"skip existing {outfile}", flush=True)
            continue
        print(f"VDP shard {shard_idx:05d}: {len(batch)} files", flush=True)
        df = read_parquet_many(batch, columns=needed)
        missing = [col for col in VDP_COLUMNS if col not in df.columns]
        if missing:
            raise KeyError(f"Missing VDP columns in shard {shard_idx}: {missing}")
        df = score_vdp_frame(df, args.prob_maps_dir, cache,
                             mask_radius_deg_per_day=mask_radius,
                             support_mask_min=support_mask_min)
        tmp = outfile.with_name(f".{outfile.name}.tmp_{os.getpid()}")
        df.to_parquet(tmp, index=False)
        tmp.replace(outfile)
        print(f"  wrote {len(df):,} rows to {outfile}", flush=True)


def format_mpc80(desig12: str, mjd_utc: float, ra_deg: float, dec_deg: float, mag: float) -> str:
    """Return one 80-character MPC observation line."""
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


def run_digest2_lines(
    obs_lines: list[str],
    digest2_exec: Path,
    digest2_dir: Path,
    chunk_tracklets: int,
    timeout_sec: int,
) -> list[str]:
    if len(obs_lines) % 2:
        raise ValueError("Expected exactly two observation lines per tracklet")

    cfg_text = "noheadings\nnorms\nNEO\n"
    n_tracklets = len(obs_lines) // 2
    all_lines: list[str] = []

    with tempfile.NamedTemporaryFile(mode="w", suffix=".config", delete=False) as cfg:
        cfg.write(cfg_text)
        cfg_path = cfg.name

    try:
        for start in range(0, n_tracklets, chunk_tracklets):
            stop = min(start + chunk_tracklets, n_tracklets)
            with tempfile.NamedTemporaryFile(mode="w", suffix=".obs", delete=False) as tmp:
                tmp.write("\n".join(obs_lines[2 * start : 2 * stop]) + "\n")
                obs_path = tmp.name
            try:
                result = subprocess.run(
                    [str(digest2_exec), "-p", str(digest2_dir), "-c", cfg_path, obs_path],
                    capture_output=True,
                    text=True,
                    timeout=timeout_sec,
                )
            finally:
                os.unlink(obs_path)
            if result.returncode != 0:
                raise RuntimeError(
                    f"digest2 failed for tracklets {start}:{stop} "
                    f"with code {result.returncode}\n{result.stderr[:2000]}"
                )
            all_lines.extend(result.stdout.splitlines())
            print(f"  digest2 {stop:,}/{n_tracklets:,}", flush=True)
    finally:
        os.unlink(cfg_path)
    return all_lines


def parse_digest2_scores(out_lines: list[str]) -> dict[str, float]:
    scores: dict[str, float] = {}
    for line in out_lines:
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        try:
            scores[parts[0]] = int(parts[1]) / 100.0
        except ValueError:
            continue
    return scores


def sample(args: argparse.Namespace) -> None:
    """Build a stratified subsample from VDP shards for digest2 scoring.

    Keeps all NEO tracklets and samples --n-sample-nonneo non-NEOs
    proportionally by population across shards. Writes a single shuffled
    parquet that run-digest2 reads instead of the full VDP shard set.
    """
    indir = vdp_shard_dir(args.work_dir)
    paths = sorted(indir.glob("vdp_*.parquet"))
    if args.limit_files is not None:
        paths = paths[: args.limit_files]
    if not paths:
        raise FileNotFoundError(f"No vdp_*.parquet found in {indir}")

    outfile = args.subsample_file
    if outfile.exists() and not args.overwrite:
        raise FileExistsError(f"{outfile} exists; pass --overwrite")

    # Pass 1: count total non-NEOs (read only population column — fast).
    print("Pass 1: counting non-NEO rows …", flush=True)
    total_nonneo = 0
    for i, p in enumerate(paths, 1):
        df_pop = pd.read_parquet(p, columns=["population"])
        total_nonneo += int((df_pop["population"] != "NEO").sum())
        if i % 20 == 0:
            print(f"  counted {i}/{len(paths)} shards, non-NEO so far: {total_nonneo:,}", flush=True)
    print(f"  total non-NEO: {total_nonneo:,}", flush=True)

    rate = min(1.0, args.n_sample_nonneo / max(total_nonneo, 1))
    print(f"  sampling rate: {rate:.5f}  (target non-NEO: {args.n_sample_nonneo:,})", flush=True)

    # Pass 2: stream shards, keep all NEOs, sample non-NEOs at computed rate.
    print("Pass 2: sampling …", flush=True)
    rng = np.random.default_rng(args.sample_seed)
    neo_parts: list[pd.DataFrame] = []
    nonneo_parts: list[pd.DataFrame] = []

    for i, p in enumerate(paths, 1):
        df = pd.read_parquet(p)
        neo_parts.append(df[df["population"] == "NEO"])
        nonneo = df[df["population"] != "NEO"]
        n_take = round(len(nonneo) * rate)
        if 0 < n_take < len(nonneo):
            nonneo = nonneo.sample(n=n_take, random_state=int(rng.integers(2**31)))
        nonneo_parts.append(nonneo)
        if i % 20 == 0:
            print(f"  sampled {i}/{len(paths)} shards", flush=True)

    neo_df = pd.concat(neo_parts, ignore_index=True)
    nonneo_df = pd.concat(nonneo_parts, ignore_index=True)
    print(f"NEO: {len(neo_df):,}  non-NEO sampled: {len(nonneo_df):,}", flush=True)

    result = (
        pd.concat([neo_df, nonneo_df], ignore_index=True)
        .sample(frac=1, random_state=int(rng.integers(2**31)))
        .reset_index(drop=True)
    )

    print("Population breakdown in subsample:")
    for pop, n in result["population"].value_counts().sort_index().items():
        print(f"  {pop}: {n:,}")

    outfile.parent.mkdir(parents=True, exist_ok=True)
    tmp = outfile.with_name(f".{outfile.name}.tmp_{os.getpid()}")
    result.to_parquet(tmp, index=False)
    tmp.replace(outfile)
    print(f"Wrote {len(result):,} rows to {outfile}", flush=True)


def run_digest2(args: argparse.Namespace) -> None:
    # If a subsample file exists (written by 'sample'), use it as a single shard.
    if args.subsample_file.exists():
        print(f"Using subsample file: {args.subsample_file}", flush=True)
        paths = [args.subsample_file]
    else:
        indir = vdp_shard_dir(args.work_dir)
        paths = sorted(indir.glob("vdp_*.parquet"))
        if args.limit_files is not None:
            paths = paths[: args.limit_files]
        if not paths:
            raise FileNotFoundError(f"No vdp_*.parquet files found in {indir}")

    digest2_exec = args.digest2_exec
    digest2_dir = args.digest2_dir
    if digest2_exec is None and digest2_dir is not None:
        digest2_exec = digest2_dir / "digest2"
    if digest2_exec is None:
        raise ValueError("Provide --digest2-dir or --digest2-exec")
    if digest2_dir is None:
        digest2_dir = digest2_exec.parent
    if not digest2_exec.exists():
        raise FileNotFoundError(digest2_exec)

    outdir = digest2_shard_dir(args.work_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    global_offset = 0
    selected = []
    for shard_idx, path in enumerate(paths):
        rows = pq.ParquetFile(path).metadata.num_rows
        if args.shard_index is not None and shard_idx != args.shard_index:
            global_offset += rows
            continue
        if args.shard_start is not None and shard_idx < args.shard_start:
            global_offset += rows
            continue
        if args.shard_stop is not None and shard_idx >= args.shard_stop:
            global_offset += rows
            continue
        selected.append((shard_idx, path, global_offset, rows))
        global_offset += rows
    if not selected:
        raise IndexError("No digest2 shards selected")

    for shard_idx, path, row_offset, rows in selected:
        row_start = args.d2_row_start if args.d2_row_start is not None else 0
        row_stop  = args.d2_row_stop  if args.d2_row_stop  is not None else rows

        if args.d2_row_start is not None or args.d2_row_stop is not None:
            stem = path.stem.replace("vdp_", "", 1)
            outfile = outdir / f"d2_{stem}_r{row_start:07d}.parquet"
        else:
            stem = path.name
            if stem.startswith("vdp_"):
                stem = "d2_" + stem[4:]
            elif not stem.startswith("d2_"):
                stem = "d2_" + stem
            outfile = outdir / stem

        if outfile.exists() and not args.overwrite:
            print(f"skip existing {outfile}", flush=True)
            continue

        print(f"digest2 shard {shard_idx:05d} rows {row_start}-{row_stop}: {path}", flush=True)
        df = pd.read_parquet(path).iloc[row_start:row_stop].reset_index(drop=True)
        missing = [col for col in DIGEST2_COLUMNS if col not in df.columns]
        if missing:
            raise KeyError(f"Missing digest2 columns in {path}: {missing}")

        # d2_key: 7-char "D######" — what digest2 echoes in output (strips 5-char packed field)
        # mpc_desig: 12-char "     D######" — correct MPC 80-col unnumbered-object format
        # Global ID = row_start + local index (unique across parallel chunks)
        global_offset = row_offset + row_start
        d2_keys = [f"D{global_offset + i:06d}" for i in range(len(df))]
        mpc_desigs = [f"     {k}" for k in d2_keys]
        obs_lines = []
        for mpc_desig, row in zip(mpc_desigs, df.itertuples(index=False)):
            obs_lines.append(format_mpc80(mpc_desig, row.mjd0_utc, row.ra0, row.dec0, row.mag0))
            obs_lines.append(format_mpc80(mpc_desig, row.mjd1_utc, row.ra1, row.dec1, row.mag1))

        out_lines = run_digest2_lines(
            obs_lines=obs_lines,
            digest2_exec=digest2_exec,
            digest2_dir=digest2_dir,
            chunk_tracklets=args.digest2_chunk_tracklets,
            timeout_sec=args.digest2_timeout_sec,
        )
        scores = parse_digest2_scores(out_lines)
        df["digest2_id"] = d2_keys
        df["P_NEO_d2"] = [scores.get(k, 0.0) for k in d2_keys]
        print(f"  parsed={len(scores):,} missing={len(df) - len(scores):,}", flush=True)

        tmp = outfile.with_name(f".{outfile.name}.tmp_{os.getpid()}")
        df.to_parquet(tmp, index=False)
        tmp.replace(outfile)
        print(f"  wrote {len(df):,} rows to {outfile}", flush=True)


def combine(args: argparse.Namespace) -> None:
    indir = digest2_shard_dir(args.work_dir)
    paths = sorted(indir.glob("d2_*.parquet"))
    if not paths:
        indir = vdp_shard_dir(args.work_dir)
        paths = sorted(indir.glob("vdp_*.parquet"))
    if args.limit_files is not None:
        paths = paths[: args.limit_files]
    if not paths:
        raise FileNotFoundError(f"No scored shards found under {args.work_dir}")
    if args.outfile.exists() and not args.overwrite:
        raise FileExistsError(f"{args.outfile} exists; pass --overwrite")

    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    frames = []
    total = 0
    for i, path in enumerate(paths, start=1):
        df = pd.read_parquet(path)
        total += len(df)
        frames.append(df)
        print(f"read {i:,}/{len(paths):,} total={total:,}", flush=True)
    out = pd.concat(frames, ignore_index=True)
    tmp = args.outfile.with_name(f".{args.outfile.name}.tmp_{os.getpid()}")
    out.to_parquet(tmp, index=False)
    tmp.replace(args.outfile)
    print(f"wrote {len(out):,} rows to {args.outfile}", flush=True)


def main() -> None:
    args = parse_args()
    if args.command == "audit":
        audit(args)
    elif args.command == "score-vdp":
        score_vdp(args)
    elif args.command == "sample":
        sample(args)
    elif args.command == "run-digest2":
        run_digest2(args)
    elif args.command == "combine":
        combine(args)
    else:
        raise ValueError(args.command)


if __name__ == "__main__":
    main()
