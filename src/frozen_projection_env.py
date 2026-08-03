"""Frozen projection environment for CAL and TEST (E0 governance item 5).

GEN was built with IERS auto-download ON, which is why its projection cannot be byte-reproduced
today (~1 ULP float32 drift; job 38078793). CAL and TEST must NOT inherit that: import this module
FIRST, before any astropy coordinate work, so the IERS table and ephemeris kernel are pinned.

    import frozen_projection_env as fpe
    fpe.activate()          # raises if the frozen files are missing or hash-mismatched
"""
from __future__ import annotations
import hashlib, json
from pathlib import Path

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
MANIFEST = W/"outputs/splits/GEN_MANIFEST.json"
_ACTIVE = None


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def activate(strict: bool = True) -> dict:
    """Disable IERS downloads, install the preserved table, pin the ephemeris kernel."""
    global _ACTIVE
    if _ACTIVE is not None:
        return _ACTIVE
    man = json.loads(MANIFEST.read_text())
    info = {"manifest_sha256": _sha(MANIFEST)}

    from astropy.utils import iers
    from astropy.utils.data import import_file_to_cache
    from astropy.coordinates import solar_system_ephemeris

    # 1. no downloads, ever -- a silent refresh is exactly what desynchronised GEN
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    info["auto_download"] = iers.conf.auto_download

    # 2. install the preserved table into astropy's cache UNDER THE AUTO URL, then open it.
    # It is finals2000A.all (IERS-A). Opening it with IERS_B.open() raises
    # "Column year failed to convert" -- IERS_B cannot parse IERS-A. Verified to be the exact file
    # the live IERS_Auto was using (same cache path, same sha256).
    t = man.get("iers", {})
    tp = W/t["preserved_table"]
    got = _sha(tp)
    if strict and got != t.get("preserved_table_sha256"):
        raise RuntimeError(f"frozen IERS table hash mismatch: {got[:16]} != "
                           f"{str(t.get('preserved_table_sha256'))[:16]}")
    url = t.get("iers_auto_url") or iers.conf.iers_auto_url
    import_file_to_cache(url, tp, replace=True)
    iers.earth_orientation_table.set(iers.IERS_Auto.open())
    info["iers_table"] = str(tp); info["iers_sha256"] = got; info["iers_url"] = url

    # 3. pin the ephemeris. neoscore.py/NEO_H.py/neoom.py all do
    #    `with solar_system_ephemeris.set("de432s")`, which resolves the NAME -> URL ->
    #    download_file(cache=True). Setting a path globally would NOT survive those contexts, so
    #    the frozen bytes are imported UNDER THE ORIGINAL URL KEY. After this, the name itself
    #    resolves to the frozen kernel and those contexts cannot override it.
    e = man.get("ephemeris", {})
    name = e.get("solar_system_ephemeris", "de432s")
    DE_URL = ("https://naif.jpl.nasa.gov/pub/naif/generic_kernels"
              f"/spk/planets/{name.lower()}.bsp")
    if "preserved_kernel" in e:
        kp = W/e["preserved_kernel"]
        kgot = _sha(kp)
        if strict and kgot != e.get("kernel_sha256"):
            raise RuntimeError(f"frozen ephemeris kernel hash mismatch: {kgot[:16]} != "
                               f"{str(e.get('kernel_sha256'))[:16]}")
        import_file_to_cache(DE_URL, kp, replace=True)
        info["kernel_preserved"] = str(kp); info["kernel_sha256_frozen"] = kgot
        info["kernel_url_key"] = DE_URL
    solar_system_ephemeris.set(name)
    info["solar_system_ephemeris"] = name

    # --- report what is ACTIVE, not what was requested ---
    import astropy.coordinates.solar_system as _ss
    from astropy.time import Time as _T
    _ = _T("2027-08-25T00:00:00", scale="utc").ut1          # force IERS load
    _tab = iers.earth_orientation_table.get()
    _dp = Path(str(dict(getattr(_tab, "meta", {})).get("data_path", "")))/"contents"
    info["active_iers_class"] = type(_tab).__name__
    info["active_iers_path"] = str(_dp)
    info["active_iers_sha256"] = _sha(_dp) if _dp.exists() else None
    with solar_system_ephemeris.set(name):                   # the exact context neoscore.py uses
        _k = _ss._get_kernel(name)
        _kp = Path(_k.daf.file.name)
    info["active_kernel_path"] = str(_kp)
    info["active_kernel_sha256"] = _sha(_kp) if _kp.exists() else None
    if strict:
        if info["active_iers_sha256"] != info.get("iers_sha256"):
            raise RuntimeError(f"ACTIVE IERS table {str(info['active_iers_sha256'])[:16]} != frozen "
                               f"{str(info.get('iers_sha256'))[:16]}")
        if info.get("kernel_sha256_frozen") and \
                info["active_kernel_sha256"] != info["kernel_sha256_frozen"]:
            raise RuntimeError(
                f"ACTIVE ephemeris kernel {str(info['active_kernel_sha256'])[:16]} != frozen "
                f"{info['kernel_sha256_frozen'][:16]} -- a set('{name}') context escaped the freeze")

    _ACTIVE = info
    return info
