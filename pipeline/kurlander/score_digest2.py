#!/usr/bin/env python3
"""score_digest2.py — run digest2 on the Kurlander referee tracklets (HYAK ONLY: needs the binary).

Builds the digest2 input EXACTLY as the v5 baseline did (`run_digest2_comparison_gmm.py`):
two synthetic MPC 80-column observations per tracklet, the second obtained from the first by
propagating with the RAW alpha_dot: ra1 = ra0 + dra*dt (no cos(dec) divide). This keeps the
comparison apples-to-apples with the stored v5 `P_NEO_d2` and gives digest2 exactly the same
information our ranging engine gets.

Usage (on klone):
  python pipeline/kurlander/score_digest2.py \
      --in  outputs/kurlander/referee_eval.parquet \
      --out outputs/kurlander/referee_eval_d2.parquet \
      --d2-dir /mmfs1/gscratch/dirac/ds2004/digest2
"""
import argparse, os, subprocess, tempfile, time
import numpy as np, pandas as pd

DT_DAYS = 30.0 / 1440.0          # 30-minute synthetic arc, as in the v5 baseline


def format_mpc80(desig12, year, month, day_frac, ra_deg, dec_deg, mag):
    ra_h = ra_deg / 15.0
    rh = int(ra_h); rm = int((ra_h - rh) * 60); rs = ((ra_h - rh) * 60 - rm) * 60
    sign = "-" if dec_deg < 0 else "+"
    ad = abs(dec_deg); dd = int(ad); dm = int((ad - dd) * 60); ds = ((ad - dd) * 60 - dm) * 60
    return (f"     {desig12:<7s}  C{year:4d} {month:02d} {day_frac:08.5f} "
            f"{rh:02d} {rm:02d} {rs:06.3f}{sign}{dd:02d} {dm:02d} {ds:05.2f}"
            f"         {mag:04.1f} V      X05")


def run_digest2(obs_lines, d2_dir, cfg_text="noheadings\nnorms\nNEO\n"):
    with tempfile.TemporaryDirectory() as td:
        cfgp = os.path.join(td, "d2.cfg"); obsp = os.path.join(td, "obs.txt")
        open(cfgp, "w").write(cfg_text)
        open(obsp, "w").write("\n".join(obs_lines) + "\n")
        r = subprocess.run([os.path.join(d2_dir, "digest2"), "-c", cfgp, obsp],
                           capture_output=True, text=True, cwd=d2_dir)
        if r.returncode != 0:
            raise RuntimeError(r.stderr[:2000])
        return r.stdout


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--d2-dir", required=True)
    ap.add_argument("--chunk", type=int, default=20000)
    ap.add_argument("--row-start", type=int, default=None,
                     help="slice input to [row-start:row-stop) before scoring (for parallel Slurm workers)")
    ap.add_argument("--row-stop", type=int, default=None)
    a = ap.parse_args()

    df = pd.read_parquet(a.inp).reset_index(drop=True)
    if a.row_start is not None or a.row_stop is not None:
        rs = a.row_start or 0
        re_ = a.row_stop if a.row_stop is not None else len(df)
        df = df.iloc[rs:re_].reset_index(drop=True)
        print(f"sliced to rows [{rs}:{re_}) -> {len(df):,} tracklets")
    print(f"{len(df):,} tracklets -> digest2")
    scores = np.full(len(df), np.nan)

    t0 = time.time()
    for s in range(0, len(df), a.chunk):
        e = min(s + a.chunk, len(df)); sub = df.iloc[s:e]
        lines, keys = [], []
        for j, (_, row) in enumerate(sub.iterrows()):
            desig = f"D{s+j:06d}"
            ra0, dec0 = float(row.mean_ra), float(row.mean_dec)
            ra1 = (ra0 + float(row.mean_dra) * DT_DAYS) % 360.0     # RAW alpha_dot, as v5 baseline
            dec1 = dec0 + float(row.mean_ddec) * DT_DAYS
            # place the 30-min arc at a fixed epoch (digest2 uses geometry+rates, not absolute date,
            # for classification; the v5 baseline did the same)
            from astropy.time import Time
            t = Time(float(row.mjd_mean), format="mjd").datetime
            d0 = t.day + (t.hour*3600 + t.minute*60 + t.second)/86400.0
            lines.append(format_mpc80(desig, t.year, t.month, d0, ra0, dec0, float(row.mean_mag)))
            lines.append(format_mpc80(desig, t.year, t.month, d0 + DT_DAYS, ra1, dec1, float(row.mean_mag)))
            keys.append(desig)
        out = run_digest2(lines, a.d2_dir)
        m = {}
        for ln in out.splitlines():
            p = ln.split()
            if len(p) >= 2 and p[0].startswith("D"):
                try: m[p[0]] = int(p[1]) / 100.0
                except ValueError: pass
        for j, k in enumerate(keys):
            scores[s + j] = m.get(k, np.nan)
        print(f"  {e:,}/{len(df):,}  ({time.time()-t0:.0f}s)", flush=True)

    df["P_NEO_d2"] = scores
    df.to_parquet(a.out, index=False)
    print(f"wrote {a.out}; digest2 nan frac {np.mean(~np.isfinite(scores)):.4f}")
