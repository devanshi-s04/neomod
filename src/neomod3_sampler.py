"""
neomod3_sampler.py — sample training orbits from the NEOMOD3 debiased NEO model.

NEOMOD3 (Nesvorny et al.) is a 4D probability table array4D[iH, ia, ie, ii]
representing the TRUE debiased NEO orbital distribution. We use it to augment
the visible source NEO training set for the GMM cloner, covering unusual orbital
configurations (high-i, Aten-class, high-e) that the observed population misses.

This module does NOT use NEOMOD3 as a scoring function (that approach failed — see
old_ones/antisunneo.ipynb). We only sample from it to build a richer training set.

File: NEOMOD3/input_neomod3.dat  (62 MB, downloaded from Boulder SWRI on 2026-06-01)
Bins: H(52) 15–28 mag, a(42) 0–4.2 AU, e(25) 0–1, i(22) 0–88°
"""

import re
import numpy as np
import pandas as pd
from pathlib import Path

# --- constants (must match velocity_density_pipeline_gmm.py) ---
AU_KM = 149_597_870.7            # km / AU
MU_SUN = 1.32712440018e11        # km^3 / s^2  (GM_sun)

# --- bin definitions (verified against old_ones/neomod-scratch.ipynb) ---
N_H, H_MIN, H_MAX = 52, 15.0, 28.0   # dH = 0.25 mag
N_A, A_MIN, A_MAX = 42, 0.0,  4.2    # da = 0.10 AU
N_E, E_MIN, E_MAX = 25, 0.0,  1.0    # de = 0.04
N_I, I_MIN, I_MAX = 22, 0.0, 88.0    # di = 4.0°

NEOMOD3_LOCAL = (
    Path(__file__).parent.parent / "NEOMOD3" / "input_neomod3.dat"
)
NEOMOD3_URL = (
    "https://www2.boulder.swri.edu/~davidn/NEOMD_Simulator/input_neomod3.dat"
)

# Module-level cache — loaded once per process
_NEOMOD3_CACHE = None  # (array4D, centers, widths)


def load_neomod3_array(force_reload=False):
    """Load the NEOMOD3 4D probability table.

    Returns (array4D, centers, widths) where:
      array4D  — shape (52, 42, 25, 22), raw model weights
      centers  — dict of bin-center arrays keyed 'H','a','e','i'
      widths   — dict of bin widths keyed 'H','a','e','i'

    Downloads from Boulder SWRI if local copy is absent.
    Uses a module-level cache so the 62 MB file is only parsed once per process.
    """
    global _NEOMOD3_CACHE
    if _NEOMOD3_CACHE is not None and not force_reload:
        return _NEOMOD3_CACHE

    if NEOMOD3_LOCAL.exists():
        text = NEOMOD3_LOCAL.read_text(errors="ignore")
    else:
        print(f"  [neomod3_sampler] Downloading NEOMOD3 from {NEOMOD3_URL}")
        from urllib.request import urlopen
        import ssl
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        text = urlopen(NEOMOD3_URL, context=ctx).read().decode("utf-8", errors="ignore")
        NEOMOD3_LOCAL.parent.mkdir(parents=True, exist_ok=True)
        NEOMOD3_LOCAL.write_text(text)
        print(f"  [neomod3_sampler] Saved to {NEOMOD3_LOCAL}")

    lines = text.splitlines()
    # Header line: "  iH  ia  ie  ii  binned mod.    nu6    3:1  ..."
    header_idx = next(
        i for i, ln in enumerate(lines)
        if re.match(r"^\s*iH\s+ia\s+ie\s+ii\s+binned\s+mod\.", ln)
    )

    # Parse data rows: iH  ia  ie  ii  weight  (columns 4 onward ignored)
    rows = []
    for ln in lines[header_idx + 1:]:
        parts = ln.split()
        if len(parts) >= 5:
            try:
                rows.append([
                    int(parts[0]), int(parts[1]),
                    int(parts[2]), int(parts[3]),
                    float(parts[4]),
                ])
            except ValueError:
                continue

    df = pd.DataFrame(rows, columns=["iH", "ia", "ie", "ii", "weight"])

    dH = (H_MAX - H_MIN) / N_H
    da = (A_MAX - A_MIN) / N_A
    de = (E_MAX - E_MIN) / N_E
    di = (I_MAX - I_MIN) / N_I

    array4D = np.zeros((N_H, N_A, N_E, N_I))
    idx = df[["iH", "ia", "ie", "ii"]].astype(int).values - 1  # 0-based
    array4D[idx[:, 0], idx[:, 1], idx[:, 2], idx[:, 3]] = df["weight"].values

    centers = {
        "H": H_MIN + (np.arange(N_H) + 0.5) * dH,
        "a": A_MIN + (np.arange(N_A) + 0.5) * da,
        "e": E_MIN + (np.arange(N_E) + 0.5) * de,
        "i": I_MIN + (np.arange(N_I) + 0.5) * di,
    }
    widths = {"H": dH, "a": da, "e": de, "i": di}

    _NEOMOD3_CACHE = (array4D, centers, widths)
    return _NEOMOD3_CACHE


def sample_neomod3_orbits(n_samples, obstime_str, rng=None,
                          array4D=None, centers=None, widths=None):
    """Sample n_samples NEO orbital elements from the NEOMOD3 debiased distribution.

    Uses lazy-loaded NEOMOD3 cache if array4D/centers/widths not supplied.
    Phase angles (node, argperi, M) are drawn uniformly — the debiased model
    has no preferred angular orientation.
    t_p is reconstructed from M at obstime_str using two-body mechanics.

    Returns a DataFrame with columns compatible with _gmm_feature_matrix:
        a [AU], e, i [deg], node [deg], argperi [deg], t_p [MJD TDB], H, q [AU]

    Only physically valid NEO orbits are returned (q < 1.3 AU, e < 1, a > 0.05).
    The returned length may be slightly less than n_samples due to validity filtering.
    Caller should oversample if an exact count is needed.
    """
    from astropy.time import Time

    if array4D is None:
        array4D, centers, widths = load_neomod3_array()
    if rng is None:
        rng = np.random.default_rng(42)

    # --- sample bin indices proportional to model weight ---
    flat = array4D.flatten()
    flat = np.maximum(flat, 0.0)
    flat /= flat.sum()
    flat_idx = rng.choice(len(flat), size=n_samples, p=flat)
    iH, ia, ie, ii = np.unravel_index(flat_idx, array4D.shape)

    # --- perturb uniformly within each bin ---
    H = centers["H"][iH] + rng.uniform(-widths["H"] / 2, widths["H"] / 2, n_samples)
    a = centers["a"][ia] + rng.uniform(-widths["a"] / 2, widths["a"] / 2, n_samples)
    e = centers["e"][ie] + rng.uniform(-widths["e"] / 2, widths["e"] / 2, n_samples)
    i = centers["i"][ii] + rng.uniform(-widths["i"] / 2, widths["i"] / 2, n_samples)

    # --- phase angles: uniform ---
    node    = rng.uniform(0.0, 360.0, n_samples)
    argperi = rng.uniform(0.0, 360.0, n_samples)
    M_deg   = rng.uniform(0.0, 360.0, n_samples)

    # --- reconstruct t_p from mean anomaly at obstime ---
    # n [rad/day] = sqrt(GM / a^3) with GM in km^3/s^2, a in km
    t_obs = Time(obstime_str, scale="tdb").mjd
    a_km = a * AU_KM
    n_rad_day = np.sqrt(MU_SUN / (a_km ** 3)) * 86400.0   # sqrt(km^3/s^2 / km^3) * s/day
    M_rad = np.deg2rad(M_deg)
    t_p = t_obs - M_rad / n_rad_day

    # --- derived quantity ---
    q = a * (1.0 - e)

    # --- validity filter: keep only physical NEO orbits ---
    valid = (
        (a > 0.05) & (a < 4.2) &
        (e >= 0.0) & (e < 1.0) &
        (q >= 0.0) & (q < 1.3) &  # NEO definition
        (i >= 0.0) & (i <= 90.0) &
        np.isfinite(H) & np.isfinite(a) & np.isfinite(e) & np.isfinite(i)
    )

    return pd.DataFrame({
        "a":       a[valid],
        "e":       e[valid],
        "i":       i[valid],
        "node":    node[valid],
        "argperi": argperi[valid],
        "t_p":     t_p[valid],
        "H":       H[valid],
        "q":       q[valid],
    })
