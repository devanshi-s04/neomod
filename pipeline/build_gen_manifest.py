#!/usr/bin/env python3
"""Complete SHA256 manifest of the frozen GEN artifact + frozen projection environment.

GEN is accepted as frozen (seed-42 draw proven by bit-identical orbital elements; ~1-ULP projection
drift documented, not regenerated). This pins every byte it depends on so any later change is
detectable, and freezes the projection environment CAL and TEST must reuse.
"""
import hashlib, json, subprocess, sys, time, glob, os
from pathlib import Path
from joblib import Parallel, delayed
W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
OUT = W/"outputs/splits"

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return (str(p), h.hexdigest(), os.path.getsize(p))

# The cache was built 2026-07-31; every current commit postdates it. So the source files hashed
# today are NOT proof of what produced the cache -- they are the CURRENT builder state. The only
# behavioural evidence is the seed-42 reproduction (job 38078793), which showed the sampler still
# emits bit-identical orbital elements. Kept in a separate group so the two are never conflated.
ARTIFACT_GROUPS = ("datacube", "cache_meta", "monolithic", "shards", "by_pixel")
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
manifest = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "groups": {}}
for g, files in groups.items():
    files = [f for f in files if Path(f).exists()]
    res = Parallel(n_jobs=32)(delayed(sha)(f) for f in files)
    manifest["groups"][g] = {
        "n_files": len(res),
        "total_bytes": sum(r[2] for r in res),
        "files": {str(Path(p).relative_to(W)): h for p, h, _ in sorted(res)},
    }
    print(f"  {g:12s} {len(res):>5} files  {sum(r[2] for r in res)/2**30:8.2f} GiB  ({time.time()-t0:.0f}s)", flush=True)

# --- REQUIRE a clean worktree: a dirty tree makes the recorded commit meaningless ---
_dirty = subprocess.run(["git", "-C", str(W/"neomod"), "status", "--porcelain"],
                        capture_output=True, text=True).stdout.strip()
if _dirty:
    print("REFUSING: worktree is dirty; commit before building the frozen manifest:\n" + _dirty)
    sys.exit(2)
manifest["git_commit"] = subprocess.run(
    ["git", "-C", str(W/"neomod"), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
manifest["git_dirty"] = bool(subprocess.run(
    ["git", "-C", str(W/"neomod"), "status", "--porcelain"], capture_output=True, text=True).stdout.strip())

# ---- frozen projection environment (must be reused by CAL and TEST) ----
import astropy, erfa, numpy, scipy, pandas
from astropy.utils import iers
from astropy.utils.data import _get_download_cache_loc
env = {"python": sys.version.split()[0], "astropy": astropy.__version__, "erfa": erfa.__version__,
       "numpy": numpy.__version__, "scipy": scipy.__version__, "pandas": pandas.__version__}
iers_info = {"policy": "REQUIRED for CAL/TEST: load the preserved table with "
                       "iers.IERS_Auto.read(path) (NOT .open(), which silently substitutes the "
                       "table bundled in astropy_iers_data -- 75 ms different UT1-UTC at the map "
                       "epoch, ~1.1 arcsec of Earth rotation)."}
try:
    from astropy.utils.data import download_file
    import shutil
    # Resolve the table BY URL. An earlier version picked 'the newest large file in the astropy
    # download cache', which silently selected de432s.bsp once the kernel was cached and
    # OVERWROTE the preserved IERS table with a binary SPK. Never guess; resolve by URL.
    url = iers.conf.iers_auto_url
    src_path = Path(download_file(url, cache=True, show_progress=False))
    head = src_path.read_bytes()[:8]
    if head.startswith(b"DAF/SPK"):
        raise RuntimeError(f"resolved IERS path {src_path} is a binary SPK kernel, not an IERS table")
    frozen = W/"outputs/splits/frozen_iers"; frozen.mkdir(parents=True, exist_ok=True)
    dst = frozen/"finals2000A.all"
    shutil.copy2(src_path, dst)
    iers_info.update({
        "preserved_table": str(dst.relative_to(W)),
        "preserved_table_sha256": sha(dst)[1],
        "preserved_table_bytes": dst.stat().st_size,
        "source_cache_path": str(src_path),
        "iers_auto_url": url,
        "format": "finals2000A.all (IERS-A ASCII)",
        "loader": "iers.IERS_Auto.read(path)",
    })
except Exception as e:
    iers_info["error"] = f"{type(e).__name__}: {e}"
manifest["environment"] = env
manifest["iers"] = iers_info

mp = OUT/"GEN_MANIFEST.json"
mp.write_text(json.dumps(manifest, indent=2, sort_keys=True))
mh = sha(mp)[1]
print(f"\nmanifest -> {mp}")
print(f"  files hashed: {sum(v['n_files'] for v in manifest['groups'].values()):,}")
print(f"  total bytes : {sum(v['total_bytes'] for v in manifest['groups'].values())/2**30:.2f} GiB")
print(f"  MANIFEST SHA256 = {mh}")
print(f"  git commit = {manifest['git_commit'][:12]}  dirty={manifest['git_dirty']}")
print(f"  env = {env}")
print(f"  iers = {json.dumps(iers_info)[:200]}")

gp = OUT/"GEN_PROVENANCE.json"
prov = json.loads(gp.read_text())
prov["manifest"] = {"path": "outputs/splits/GEN_MANIFEST.json", "sha256": mh,
                    "n_files": sum(v["n_files"] for v in manifest["groups"].values())}
prov["full_bit_reproduction_test"] = {
    "job": 38078793, "result": "FAILED (by design of its own criterion)",
    "seed_draw_provenance": "PASSED -- orbital elements bit-identical across 970,429 rows",
    "projected_column_byte_reproduction": "FAILED -- ~1 ULP float32 drift",
    "cause": "environment drift (IERS_Auto table updated); NOT run-to-run non-determinism "
             "(same-environment reruns agree with worst |diff| exactly 0)",
    "disposition": "accepted; GEN is frozen as-is. Regeneration is NOT required and would itself "
                   "produce a different projection."}
prov["frozen_environment"] = env
prov["iers"] = iers_info
gp.write_text(json.dumps(prov, indent=2, sort_keys=True))
print(f"\nGEN_PROVENANCE updated  sha256[:16] = {sha(gp)[1][:16]}")
