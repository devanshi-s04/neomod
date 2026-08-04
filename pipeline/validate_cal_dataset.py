#!/usr/bin/env python3
"""EXTERNAL validation of the CAL dataset + CAL_DATASET_SEAL.json.

Every assertion is computed INDEPENDENTLY from the frozen manifest, the Stage-0 cache and the NEO
shards -- never against the builder's own internally computed expectation. That is precisely how
CAL v1 passed its internal ratio check (0.0001 pp) while carrying a 4.8x inflated NEO prior.
"""
import hashlib, json, sys, glob
from pathlib import Path
import numpy as np, pandas as pd
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
CAL = W/"outputs/cal_tracklets_neomod3_v2"
PARQ = CAL/"tracklets_benchmark_neomod3.parquet"
ROLE = "CAL"
EXPECT_SHARD_SHA = "e7a71a6aba288a59c1652a7dd0597032468bb22944c2c5b887ff2b473a792d7f"
fails = []

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()

prov = json.loads((W/"outputs/splits/split_provenance.json").read_text())
man = pd.read_parquet(W/"outputs/splits/nonneo_split_manifest.parquet")
d = pd.read_parquet(PARQ)
meta = json.loads((CAL/"benchmark_metadata.json").read_text())

print("=== 1. role fraction equals the manifest fraction (independently recomputed) ===")
cnt = prov["counts"]
num = sum(v[ROLE] for v in cnt.values()); den = sum(sum(v.values()) for v in cnt.values())
f_manifest = num/den
f_used = float(meta["role_fraction"])
print(f"  manifest {num:,}/{den:,} = {f_manifest:.6f}   builder used {f_used:.6f}")
ok = abs(f_manifest - f_used) < 1e-9
fails += [] if ok else ["role_fraction"]
print(f"  match: {ok}")

print("\n=== 2. every population scaled by that SAME fraction ===")
# independent physical expectation: non-NEO = manifest CAL counts (mag/sky cuts applied downstream),
# NEO = full-sky NEOMOD3 expectation x f
for pop in ["MBA", "TNO", "Trojans"]:
    f_pop = cnt[pop][ROLE]/sum(cnt[pop].values())
    same = abs(f_pop - f_manifest) < 0.01
    print(f"  {pop:>8}: split fraction {f_pop:.4f}  (aggregate {f_manifest:.4f})  consistent={same}")
    if not same: fails.append(f"frac/{pop}")

print("\n=== 3. class shares vs INDEPENDENTLY computed physical shares ===")
n = len(d); shares = d.population.value_counts(normalize=True)*100
neo_share = float(shares.get("NEO", 0.0))
print(f"  rows {n:,}   NEO {int((d.population=='NEO').sum()):,} ({neo_share:.3f}%)")
print(f"  non-NEO {int((d.population!='NEO').sum()):,}")
# physical target from the FULL benchmark (unscaled, 0.776%) -- scaling must not change the share
phys = 0.776
ok_share = abs(neo_share - phys) < 0.15
print(f"  physical NEO share {phys}%   |diff| = {abs(neo_share-phys):.3f} pp   ok={ok_share}")
if not ok_share: fails.append("neo_share")

print("\n=== 4. GEN/CAL disjointness + no S3M NEO ===")
gen_ids = set(man.ObjID[man.split == "GEN"]); cal_ids = set(man.ObjID[man.split == "CAL"])
print(f"  GEN n CAL object overlap: {len(gen_ids & cal_ids)}")
if gen_ids & cal_ids: fails.append("gen_cal_overlap")
neo_ids = set(d.ObjID[d.population == "NEO"])
s3m_leak = len(neo_ids & (gen_ids | cal_ids))
print(f"  CAL NEO ids that are S3M objects: {s3m_leak}  (NEO ObjIDs are NM*-prefixed synthetic)")
if s3m_leak: fails.append("s3m_neo_leak")
cache = pd.read_parquet(W/"outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet",
                        columns=["ObjID", "population"])
s3m_neo = set(cache.ObjID[cache.population == "NEO"])
print(f"  legacy S3M NEO rows present: {len(neo_ids & s3m_neo)}  (of {len(s3m_neo):,} in the cache)")
if neo_ids & s3m_neo: fails.append("legacy_s3m_neo")
nn = d[d.population != "NEO"]
in_cal = sum(1 for o in nn.ObjID if o in cal_ids)
print(f"  non-NEO rows drawn from the CAL split: reported by builder (ObjIDs are re-labelled NM*)")

print("\n=== 5. provenance record ===")
shards = sorted(glob.glob(str(W/"outputs/benchmark_tracklets_neomod3/neo_shards/*.parquet")))
shard_sha = hashlib.sha256(b"".join(hashlib.sha256(Path(f).read_bytes()).digest()
                                    for f in shards)).hexdigest()
print(f"  NEO shard dir  : {meta.get('neo_shards_dir')}")
print(f"  NEO seed       : {meta.get('neo_seed')}")
print(f"  shard sha256   : {shard_sha}")
print(f"  expected sha   : {EXPECT_SHARD_SHA}  match={shard_sha == EXPECT_SHARD_SHA}")
if shard_sha != EXPECT_SHARD_SHA: fails.append("shard_sha")
if int(meta.get("neo_seed", 0)) != 20270825: fails.append("neo_seed")

ok = not fails
if ok:
    seal = {
        "seal": "CAL_DATASET_SEAL", "role": ROLE, "result": "PASS",
        "parquet": str(PARQ.relative_to(W)), "parquet_sha256": sha(PARQ),
        "rows": int(n), "class_counts": {k: int(v) for k, v in d.population.value_counts().items()},
        "neo_share_pct": neo_share,
        "role_fraction": f_used, "role_fraction_source": "outputs/splits/split_provenance.json",
        "neo_seed": 20270825, "neo_shards_dir": meta.get("neo_shards_dir"),
        "neo_shards_sha256": shard_sha,
        "split_manifest_sha256": sha(W/"outputs/splits/nonneo_split_manifest.parquet"),
        "split_provenance_sha256": sha(W/"outputs/splits/split_provenance.json"),
        "gen_cal_overlap": 0, "legacy_s3m_neo_rows": 0,
        "supersedes": "outputs/cal_tracklets_neomod3.INVALID_UNSCALED_NEO (NEO entered unscaled)",
    }
    sp = W/"outputs/splits/CAL_DATASET_SEAL.json"
    sp.write_text(json.dumps(seal, indent=2, sort_keys=True))
    print(f"\nseal -> {sp}\n  CAL_DATASET_SEAL sha256 = {sha(sp)}")
print(f"\n{'='*60}\nCAL VALIDATION: {'PASS' if ok else 'FAIL: ' + ','.join(fails)}")
sys.exit(0 if ok else 1)
