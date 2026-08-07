#!/usr/bin/env python3
"""digest2-ONLY shard. Determinism is preserved EXACTLY: repeatable config + --cpu 1.

Replaces the combined VDP+digest2 scorer, whose layout was the throughput problem:
  * VDP was sharded by ROW over a uid-sorted table, so every shard touched ~all 667 map
    cells and reloaded ~667 .npz files to score ~22 rows each (~5 TB redundant I/O).
  * 8 CPU / 64 GB was reserved for a process pinned to ONE cpu.
  * No batch progress was logged, so throughput could not be read off a log.

This script fixes only the architecture, never the digest2 invocation.
VDP is scored separately (map-local partitioning); the two merge on tracklet_uid.
"""
from __future__ import annotations
import argparse, hashlib, os, subprocess, sys, tempfile, time
from pathlib import Path
import numpy as np, pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"pipeline"))
D2_CONFIG = "noheadings\nnorms\nrepeatable\nNEO\n"    # FROZEN -- verified on the live run
D2_CPUS = "1"                                          # FROZEN -- never raise


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--first-tracklet", required=True)
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--nshards", type=int, required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--batch", type=int, default=250)
    a_ = ap.parse_args()
    from sorcha_phase2 import format_mpc80
    OUT = Path(a_.out_dir); OUT.mkdir(parents=True, exist_ok=True)
    dst = OUT/f"d2_shard_{a_.shard:04d}.parquet"
    if dst.exists():
        print(f"shard {a_.shard} already complete", flush=True); return

    cols = ["tracklet_uid", "ObjID", "population", "mjd0_utc", "ra0", "dec0", "mag0_V",
            "mjd1_utc", "ra1", "dec1", "mag1_V"]
    full = pd.read_parquet(a_.first_tracklet, columns=cols).sort_values(
        "tracklet_uid", kind="mergesort").reset_index(drop=True)
    te = full.iloc[np.array_split(np.arange(len(full)), a_.nshards)[a_.shard]].reset_index(drop=True)
    n = len(te)
    assert te.tracklet_uid.is_unique, "tracklet_uid not unique within shard"
    print(f"shard {a_.shard}/{a_.nshards}: {n:,} rows", flush=True)

    keys = [f"D{i:06d}" for i in range(n)]
    obs = []
    for k, r in zip(keys, te.itertuples(index=False)):
        obs.append(format_mpc80(f"     {k}", r.mjd0_utc, r.ra0, r.dec0, r.mag0_V))
        obs.append(format_mpc80(f"     {k}", r.mjd1_utc, r.ra1, r.dec1, r.mag1_V))

    score, status = {}, {}
    with tempfile.NamedTemporaryFile("w", suffix=".config", delete=False) as c:
        c.write(D2_CONFIG); cfg = c.name
    t_start = time.time()
    try:
        for s in range(0, n, a_.batch):
            e = min(s+a_.batch, n)
            tb = time.time()
            with tempfile.NamedTemporaryFile("w", suffix=".obs", delete=False) as t:
                t.write("\n".join(obs[2*s:2*e])+"\n"); op = t.name
            try:
                r = subprocess.run([str(W/"digest2"/"digest2"), "-p", str(W/"digest2"),
                                    "-c", cfg, "--cpu", D2_CPUS, op],
                                   capture_output=True, text=True, timeout=21600)
                if r.returncode != 0:
                    for k in keys[s:e]: status[k] = f"digest2_rc{r.returncode}"
                else:
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
            el = time.time()-t_start
            rate = e/el if el > 0 else 0
            # BATCH PROGRESS -- readable from the log, so throughput never needs forensics again
            print(f"  {e:>7,}/{n:,}  batch {time.time()-tb:6.1f}s  cum {el/60:6.1f}min  "
                  f"{rate:5.2f} rows/s  ETA {(n-e)/max(rate,1e-9)/60:6.1f}min", flush=True)
    finally:
        os.unlink(cfg)

    te["P_NEO_d2"] = [score[k]/100.0 if k in score else np.nan for k in keys]   # NaN, never 0
    te["d2_status"] = [("ok" if k in score else status.get(k, "missing")) for k in keys]
    te["d2_valid"] = np.isfinite(te.P_NEO_d2.to_numpy())
    assert len(te) == n and te.tracklet_uid.is_unique
    assert ((te.d2_status.to_numpy() == "ok") ^ (~te.d2_valid.to_numpy())).all(), \
        "row neither scored nor explicitly invalid"
    out = te[["tracklet_uid", "P_NEO_d2", "d2_status", "d2_valid"]]
    out.to_parquet(dst, index=False)
    h = hashlib.sha256(pd.util.hash_pandas_object(out, index=False).values.tobytes()).hexdigest()
    print(f"DONE {n:,} rows in {(time.time()-t_start)/60:.1f} min  "
          f"valid {int(te.d2_valid.sum()):,}  invalid {int((~te.d2_valid).sum()):,}")
    print(f"CONTENT_HASH {h}")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
