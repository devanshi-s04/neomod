#!/usr/bin/env python3
"""ONE-TIME environment freezer. The ONLY tool permitted to write sealed artifacts.

Explicit URLs, explicit destination filenames, staging-then-atomic-install, format verification,
and a refusal to overwrite anything already sealed. build_gen_manifest.py is read-only and must
never do any of this (an earlier version's "newest large cache file" heuristic overwrote the sealed
IERS table with de432s.bsp and propagated it into astropy's global cache).

    freeze_environment.py --verify-only     # check the sealed environment, write nothing
    freeze_environment.py --install         # first-time install (refuses to clobber)
"""
from __future__ import annotations
import argparse, hashlib, json, os, shutil, sys, tempfile, time
from pathlib import Path
import urllib.request

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
ENV = W/"outputs/splits/frozen_env"
IERS_DIR, EPH_DIR = ENV/"iers", ENV/"ephemeris"
SEAL = W/"outputs/splits/FROZEN_ENV_SEAL.json"

ARTIFACTS = {
    "iers": {
        "dest": IERS_DIR/"finals2000A.all.GEN-v1",
        "urls": ["https://datacenter.iers.org/data/9/finals2000A.all",
                 "https://maia.usno.navy.mil/ser7/finals2000A.all"],
        "sha256": "4b828090fc94114168014b61439fa5e6ec0bdfda518075a32baffea90110954d",
        "bytes": 3758308,
        "format": "IERS-A ASCII (finals2000A.all)",
        "magic_not": b"DAF/SPK",
        "loader": "astropy.utils.iers.IERS_Auto.read(path)",
    },
    "ephemeris": {
        "dest": EPH_DIR/"de432s.bsp",
        "urls": ["https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp"],
        "sha256": "363f32e14f5255359ac32c4d38080cf28ab55564a5e16696a75f63394b666e9b",
        "bytes": 10895360,
        "format": "DAF/SPK binary kernel",
        "magic_is": b"DAF/SPK",
        "loader": "solar_system_ephemeris.set(<path>)",
    },
}


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def check_format(p: Path, spec: dict) -> None:
    head = p.read_bytes()[:8]
    if "magic_is" in spec and not head.startswith(spec["magic_is"]):
        raise RuntimeError(f"{p}: expected {spec['magic_is']!r} magic, got {head!r}")
    if "magic_not" in spec and head.startswith(spec["magic_not"]):
        raise RuntimeError(f"{p}: has forbidden {spec['magic_not']!r} magic -- wrong artifact type")


def verify(name: str, spec: dict) -> dict:
    p = spec["dest"]
    if not p.exists():
        return {"name": name, "ok": False, "reason": "missing", "path": str(p)}
    got_sha, got_sz = sha(p), p.stat().st_size
    check_format(p, spec)
    ok = got_sha == spec["sha256"] and got_sz == spec["bytes"]
    return {"name": name, "ok": ok, "path": str(p.relative_to(W)), "sha256": got_sha,
            "bytes": got_sz, "expected_sha256": spec["sha256"], "expected_bytes": spec["bytes"],
            "format": spec["format"], "loader": spec["loader"],
            "mode": oct(p.stat().st_mode & 0o777)}


def install(name: str, spec: dict) -> dict:
    p = spec["dest"]
    if p.exists():
        raise RuntimeError(f"REFUSING to overwrite sealed artifact {p} (delete deliberately first)")
    p.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=str(p.parent)) as td:   # empty staging dir
        stage = Path(td)/f"{name}.staged"
        last = None
        for url in spec["urls"]:
            try:
                with urllib.request.urlopen(url, timeout=120) as r, open(stage, "wb") as f:
                    shutil.copyfileobj(r, f)
                if sha(stage) == spec["sha256"] and stage.stat().st_size == spec["bytes"]:
                    check_format(stage, spec)
                    os.chmod(stage, 0o444)
                    os.replace(stage, p)                          # atomic install
                    return {**verify(name, spec), "source_url": url}
                last = f"hash/size mismatch from {url}"
            except Exception as e:
                last = f"{type(e).__name__}: {e}"
        raise RuntimeError(f"could not install {name}: {last}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--install", action="store_true")
    a = ap.parse_args()
    if not (a.verify_only or a.install):
        ap.error("choose --verify-only or --install")
    res = {}
    for name, spec in ARTIFACTS.items():
        res[name] = verify(name, spec) if (a.verify_only or spec["dest"].exists()) \
            else install(name, spec)
        st = "OK" if res[name]["ok"] else "FAIL"
        print(f"  {name:10s} {st:4s}  {res[name].get('path')}  sha={str(res[name].get('sha256'))[:24]} "
              f"bytes={res[name].get('bytes')} mode={res[name].get('mode')}")
    ok = all(r["ok"] for r in res.values())
    if a.install or not SEAL.exists():
        import astropy, erfa, numpy
        seal = {
            "seal": "FROZEN_ENV_SEAL", "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "artifacts": res,
            "packages": {"python": sys.version.split()[0], "astropy": astropy.__version__,
                         "erfa": erfa.__version__, "numpy": numpy.__version__},
            "provenance": {
                "iers": "GEN-v1 table. The original cache copy was destroyed by a manifest-tool bug "
                        "on 2026-08-03 and RECOVERED byte-exactly the same day from both the IERS "
                        "datacenter and the USNO mirror (both byte-identical to each other and to "
                        "the original sha 4b828090...). GEN-v1 remains authoritative.",
                "ephemeris": "de432s.bsp preserved from the astropy download cache.",
            },
            "rules": ["build_gen_manifest.py is READ-ONLY and must never write these",
                      "IERS loads via IERS_Auto.read(path) -- .open() substitutes the bundled table "
                      "(75 ms UT1-UTC, ~1.1 arcsec)",
                      "ephemeris is passed BY PATH to solar_system_ephemeris.set(...)"],
        }
        SEAL.write_text(json.dumps(seal, indent=2, sort_keys=True))
        print(f"\nseal -> {SEAL}\n  FROZEN_ENV_SEAL sha256 = {sha(SEAL)}")
    print(f"\nRESULT: {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
