#!/usr/bin/env python3
"""PROVE (not infer) that the NEOMOD3 projection cache was drawn with seed 42.

E0 Amendment 1 A1.7 / verification 4: cache_metadata.json records n_draws but NOT the seed, path or
hash, so it cannot by itself establish GEN provenance. Inference from the submit line
(--export=ALL,N_TOTAL=100000000,NSHARDS=100, no --seed override => argparse default 42) is suggestive
but not proof.

This regenerates shard 0 under the claimed parameters and requires it to reproduce the stored shard
EXACTLY. If it does, provenance is established and a signed GEN_PROVENANCE.json is written. If it
does not, the cache must NOT be relabelled GEN.
"""
import argparse, hashlib, json, sys, time, types
from pathlib import Path
import numpy as np, pandas as pd
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
sys.path.insert(0, str(W/"neomod"/"src")); sys.path.insert(0, str(W/"neomod"/"pipeline"))
import os; os.chdir(W/"neomod")
import neomod3_projection_cache as npc

CLAIM = dict(seed=42, nshards=100, n_orbits_total=100_000_000, shard=0,
             epoch="2027-08-25T00:00:00", chunk=200_000)
STORED = W/"outputs/neomod3_projection_cache/shards/nm3_proj_0000.parquet"
TMP = W/"outputs/neomod3_projection_cache/_provenance_check"
TMP.mkdir(parents=True, exist_ok=True)

npc._shard_path = lambda shard: TMP/f"repro_{shard:04d}.parquet"
args = types.SimpleNamespace(**CLAIM, overwrite=True)
t0 = time.time()
print(f"regenerating shard {CLAIM['shard']} with seed={CLAIM['seed']} "
      f"({CLAIM['n_orbits_total']//CLAIM['nshards']:,} orbits) ...", flush=True)
npc.build(args)
repro = TMP/f"repro_{CLAIM['shard']:04d}.parquet"

a = pd.read_parquet(STORED); b = pd.read_parquet(repro)
print(f"\nstored : {len(a):,} rows  {list(a.columns)}")
print(f"repro  : {len(b):,} rows  {list(b.columns)}")
same_shape = a.shape == b.shape and list(a.columns) == list(b.columns)
exact = same_shape and all(np.array_equal(a[c].to_numpy(), b[c].to_numpy(), equal_nan=True)
                           for c in a.columns)
print(f"\nshape/columns match : {same_shape}")
print(f"BIT-IDENTICAL       : {exact}")
if not exact and same_shape:
    for c in a.columns:
        d = np.abs(np.nan_to_num(a[c].to_numpy(float)) - np.nan_to_num(b[c].to_numpy(float))).max()
        print(f"   {c:12s} max|diff| = {d:.6g}")

if exact:
    meta_p = W/"outputs/neomod3_projection_cache/cache_metadata.json"
    prov = {
        "role": "GEN",
        "proven": True,
        "method": "shard-0 regeneration reproduced the stored shard bit-identically",
        "seed": CLAIM["seed"], "per_shard_seed_rule": "seed + shard_index",
        "nshards": CLAIM["nshards"], "n_draws": CLAIM["n_orbits_total"],
        "epoch": CLAIM["epoch"],
        "cache_parquet": str(W/"outputs/neomod3_projection_cache/neomod3_projection_20270825T000000.parquet"),
        "by_pixel_dir": str(W/"outputs/neomod3_projection_cache/by_pixel"),
        "cache_metadata_sha256": hashlib.sha256(meta_p.read_bytes()).hexdigest(),
        "shard0_sha256": hashlib.sha256(STORED.read_bytes()).hexdigest(),
        "verified_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "distinct_from": {"CAL_seed": 20270825, "TEST_seed": 31415926},
    }
    out = W/"outputs/splits/GEN_PROVENANCE.json"
    out.write_text(json.dumps(prov, indent=2))
    print(f"\nPROVEN -> {out}")
    print(f"  GEN_PROVENANCE sha256[:16] = {hashlib.sha256(out.read_bytes()).hexdigest()[:16]}")
else:
    print("\nNOT PROVEN -- the cache must NOT be relabelled GEN. Regenerate a GEN cache with "
          "recorded provenance instead.")
print(f"\nelapsed {time.time()-t0:.0f}s")
sys.exit(0 if exact else 1)
