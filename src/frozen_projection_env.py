"""Frozen projection environment for sealed evaluation (CAL / TEST).

Rules enforced here:
  * IERS is loaded DIRECTLY from the sealed file with IERS_Auto.read(path).
    NEVER IERS_Auto.open() -- with auto_download=False it silently substitutes the table bundled in
    astropy_iers_data, which differs by 75 ms of UT1-UTC at the map epoch (~1.1 arcsec of Earth
    rotation) from the GEN-v1 table.
  * The ephemeris is passed BY PATH, so no name->URL->cache lookup occurs. Hard-coded
    `solar_system_ephemeris.set('de432s')` sites read SORCHA_EPHEMERIS, which is set here.
  * A private astropy cache is used for the whole process lifetime; the global cache is never
    touched. `global_cache_fingerprint()` lets a caller prove that.
"""
from __future__ import annotations
import hashlib, json, os
from pathlib import Path

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
SEAL_PATH = W/"outputs/splits/FROZEN_ENV_SEAL.json"
PRIVATE_CACHE = W/"outputs/splits/frozen_env/astropy_cache_private"
GLOBAL_CACHE = Path.home()/".astropy/cache/download/url"
_ACTIVE = None


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def global_cache_fingerprint() -> str:
    """Hash of the global cache's directory listing -- to prove we never modified it."""
    if not GLOBAL_CACHE.exists():
        return "absent"
    ents = sorted(f"{p.name}:{p.stat().st_size}" for p in GLOBAL_CACHE.rglob("*") if p.is_file())
    return hashlib.sha256("\n".join(ents).encode()).hexdigest()


def activate(strict: bool = True) -> dict:
    global _ACTIVE
    if _ACTIVE is not None:
        return _ACTIVE
    seal = json.loads(SEAL_PATH.read_text())
    info = {"seal_sha256": _sha(SEAL_PATH), "global_cache_before": global_cache_fingerprint()}

    PRIVATE_CACHE.mkdir(parents=True, exist_ok=True)
    os.environ["XDG_CACHE_HOME"] = str(PRIVATE_CACHE)      # isolate before astropy imports
    from astropy.config import set_temp_cache
    _ctx = set_temp_cache(str(PRIVATE_CACHE))              # process-lifetime, never exited
    _ctx.__enter__()
    info["private_cache"] = str(PRIVATE_CACHE)

    from astropy.utils import iers
    from astropy.coordinates import solar_system_ephemeris
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    info["auto_download"] = iers.conf.auto_download

    a = seal["artifacts"]["iers"]
    tp = W/a["path"]
    got = _sha(tp)
    if strict and got != a["sha256"]:
        raise RuntimeError(f"sealed IERS hash mismatch: {got[:16]} != {a['sha256'][:16]}")
    iers.earth_orientation_table.set(iers.IERS_Auto.read(str(tp)))
    info["iers_path"] = str(tp); info["iers_sha256"] = got

    e = seal["artifacts"]["ephemeris"]
    kp = W/e["path"]
    kgot = _sha(kp)
    if strict and kgot != e["sha256"]:
        raise RuntimeError(f"sealed kernel hash mismatch: {kgot[:16]} != {e['sha256'][:16]}")
    os.environ["SORCHA_EPHEMERIS"] = str(kp)               # parameterized call sites use this
    solar_system_ephemeris.set(str(kp))                    # BY PATH: no cache lookup
    info["kernel_path"] = str(kp); info["kernel_sha256"] = kgot

    # ---- report what is ACTIVE, not what was requested ----
    from astropy.time import Time
    _ = Time("2027-08-25T00:00:00", scale="utc").ut1
    tab = iers.earth_orientation_table.get()
    info["active_iers_class"] = type(tab).__name__
    info["active_iers_rows"] = int(len(tab))
    info["active_iers_sha256"] = got
    import astropy.coordinates.solar_system as _ss
    with solar_system_ephemeris.set(os.environ["SORCHA_EPHEMERIS"]):
        info["active_kernel_path"] = str(Path(_ss._get_kernel(os.environ["SORCHA_EPHEMERIS"]).daf.file.name))
    info["active_kernel_sha256"] = _sha(Path(info["active_kernel_path"]))
    if strict:
        if info["active_kernel_sha256"] != kgot:
            raise RuntimeError("ACTIVE ephemeris kernel is not the sealed kernel")
    info["global_cache_after"] = global_cache_fingerprint()
    info["global_cache_unchanged"] = info["global_cache_before"] == info["global_cache_after"]
    if strict and not info["global_cache_unchanged"]:
        raise RuntimeError("the global astropy cache was modified -- sealed evaluation must not touch it")
    _ACTIVE = info
    return info
