#!/usr/bin/env python3
"""Prove build_gen_manifest.py is read-only: sealed artifact hashes AND mtimes unchanged by it."""
import hashlib, os, subprocess, sys
from pathlib import Path
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
WATCH = [W/"outputs/splits/frozen_env/iers/finals2000A.all.GEN-v1",
         W/"outputs/splits/frozen_env/ephemeris/de432s.bsp",
         W/"outputs/splits/FROZEN_ENV_SEAL.json",
         W/"outputs/neomod3_projection_cache/cache_metadata.json"]
GLOBAL = Path.home()/".astropy/cache/download/url"
def snap():
    s = {}
    for p in WATCH:
        s[str(p)] = (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns) if p.exists() else None
    s["__global_cache__"] = sorted(f"{q.name}:{q.stat().st_size}" for q in GLOBAL.rglob("*") if q.is_file()) if GLOBAL.exists() else "absent"
    return s
before = snap()
r = subprocess.run([sys.executable, str(W/"neomod/pipeline/build_gen_manifest.py")],
                   capture_output=True, text=True)
print(r.stdout[-700:]); print("exit", r.returncode)
after = snap()
fails = [k for k in before if before[k] != after[k]]
for k in fails:
    print(f"  CHANGED: {k}\n    before {before[k]}\n    after  {after[k]}")
print(f"\nwatched artifacts unchanged (hash AND mtime): {not fails}")
print(f"global astropy cache unchanged: {before['__global_cache__'] == after['__global_cache__']}")
print("PASS" if not fails and r.returncode == 0 else "FAIL")
sys.exit(0 if (not fails and r.returncode == 0) else 1)
