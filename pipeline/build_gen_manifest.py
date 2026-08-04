#!/usr/bin/env python3
"""READ-ONLY, OFFLINE manifest of the frozen GEN artifact.

This tool may ONLY inspect, hash and report. It must never download, copy, overwrite, import into a
cache, or otherwise mutate anything. An earlier version violated that -- it picked "the newest large
file in the astropy download cache" as the IERS table, which selected de432s.bsp once the kernel was
cached, overwrote the sealed IERS table with a binary SPK, and propagated that into astropy's global
cache. freeze_environment.py is now the ONLY writer of sealed artifacts.

Separation of concerns recorded here:
  frozen_cache_artifact  -- the bytes that ARE the GEN cache (datacube, metadata, monolithic,
                            shards, HEALPix partitions)
  builder_provenance     -- CURRENT source state; the cache was built 2026-07-31 and every current
                            commit postdates it, so these hashes are NOT proof of what built it.
                            The behavioural evidence is job 38078793 (bit-identical elements).
"""
from __future__ import annotations
import hashlib, json, subprocess, sys, time, os
from pathlib import Path
from joblib import Parallel, delayed

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
OUT = W/"outputs/splits"
ARTIFACT_GROUPS = ("datacube", "cache_meta", "monolithic", "shards", "by_pixel")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return (str(p), h.hexdigest(), os.path.getsize(p))


groups = {
    "datacube":   [W/"neomod/NEOMOD3/input_neomod3.dat"],
    "cache_meta": [W/"outputs/neomod3_projection_cache/cache_metadata.json"],
    "monolithic": [W/"outputs/neomod3_projection_cache/neomod3_projection_20270825T000000.parquet"],
    "shards":     sorted(Path(W/"outputs/neomod3_projection_cache/shards").glob("*.parquet")),
    "by_pixel":   sorted(Path(W/"outputs/neomod3_projection_cache/by_pixel").rglob("*.parquet")),
    "source":     [W/"neomod/src/neomod3_sampler.py",
                   W/"neomod/pipeline/neomod3_projection_cache.py",
                   W/"neomod/src/velocity_density_pipeline_gmm.py",
                   W/"neomod/src/velocity_density_pipeline_neomod_clone_only.py",
                   W/"neomod/pipeline/neomod3_cache_healpix.py"],
}

t0 = time.time()
dirty = subprocess.run(["git", "-C", str(W/"neomod"), "status", "--porcelain"],
                       capture_output=True, text=True).stdout.strip()
if dirty:
    print("REFUSING: worktree is dirty; commit first:\n" + dirty)
    sys.exit(2)

manifest = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "tool": "build_gen_manifest.py (READ-ONLY, OFFLINE)", "groups": {}}
for g, files in groups.items():
    files = [f for f in files if Path(f).exists()]
    res = Parallel(n_jobs=32)(delayed(sha)(f) for f in files)
    manifest["groups"][g] = {"n_files": len(res), "total_bytes": sum(r[2] for r in res),
                             "files": {str(Path(p).relative_to(W)): h for p, h, _ in sorted(res)}}
    print(f"  {g:12s} {len(res):>5} files  {sum(r[2] for r in res)/2**30:8.2f} GiB  ({time.time()-t0:.0f}s)", flush=True)

manifest["git_commit"] = subprocess.run(["git", "-C", str(W/"neomod"), "rev-parse", "HEAD"],
                                        capture_output=True, text=True).stdout.strip()
manifest["git_dirty"] = False

# reference the sealed environment; verify it, never create or modify it
seal_p = OUT/"FROZEN_ENV_SEAL.json"
if seal_p.exists():
    seal = json.loads(seal_p.read_text())
    for n in ("iers", "ephemeris"):
        a = seal["artifacts"][n]
        if sha(W/a["path"])[1] != a["sha256"]:
            sys.exit(f"REFUSING: sealed {n} artifact changed since FROZEN_ENV_SEAL was written")
    manifest["frozen_env_seal"] = {"path": str(seal_p.relative_to(W)), "sha256": sha(seal_p)[1],
                                   "iers": seal["artifacts"]["iers"],
                                   "ephemeris": seal["artifacts"]["ephemeris"],
                                   "packages": seal["packages"]}
else:
    sys.exit("REFUSING: FROZEN_ENV_SEAL.json missing -- run freeze_environment.py --install")

import astropy, erfa, numpy, scipy, pandas
manifest["environment"] = {"python": sys.version.split()[0], "astropy": astropy.__version__,
                           "erfa": erfa.__version__, "numpy": numpy.__version__,
                           "scipy": scipy.__version__, "pandas": pandas.__version__}
manifest["artifact_vs_builder"] = {
    "readonly_tool": True,
    "frozen_cache_artifact_groups": list(ARTIFACT_GROUPS),
    "builder_provenance_groups": ["source"],
    "caveat": "cache built 2026-07-31 16:09; all current commits postdate it. The `source` group is "
              "the CURRENT builder state, not what produced the cache. Behavioural equivalence of "
              "the sampler comes from job 38078793 (bit-identical orbital elements).",
}
art = {k: v for g in ARTIFACT_GROUPS for k, v in manifest["groups"][g]["files"].items()}
manifest["frozen_artifact_sha256"] = hashlib.sha256(json.dumps(art, sort_keys=True).encode()).hexdigest()

mp = OUT/"GEN_MANIFEST.json"
mp.write_text(json.dumps(manifest, indent=2, sort_keys=True))
print(f"\nmanifest -> {mp}")
print(f"  files hashed : {sum(v['n_files'] for v in manifest['groups'].values()):,}")
print(f"  MANIFEST SHA256      = {sha(mp)[1]}")
print(f"  frozen_artifact_sha  = {manifest['frozen_artifact_sha256']}")
print(f"  FROZEN_ENV_SEAL sha  = {manifest['frozen_env_seal']['sha256']}")
print(f"  git commit = {manifest['git_commit'][:12]} dirty=False")
