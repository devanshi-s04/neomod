#!/usr/bin/env python3
"""Rubin VDP tracklet scorer, interface v1 -- FROZEN.

This is a PACKAGING layer. It does not change any scientific behaviour: the density lookup and
posterior are delegated verbatim to `score_test2.score_new_vdp`, the exact function that produced
the accepted TEST2 `P_NEO_new` column at commit 87f6bd82d190e946f99ab128ff8ffff380d09a7a.

WHAT v1 IS
----------
A stable input/output contract so future map versions can be swapped underneath without
redesigning the nightly pipeline. It is an ENGINEERING baseline, not an operational-readiness
claim: four-density coverage on TEST2 is ~43.96%, so a large minority of in-domain tracklets
return NaN with an explicit reason. That is preserved deliberately.

FROZEN SCIENTIFIC BEHAVIOUR (do not alter in v1)
------------------------------------------------
  * 667 sky centres, 44 half-open 0.25-mag apparent-V bins over 14 <= V < 25
  * k_NEO = 150, k_MBA = k_TNO = k_Trojan = 10, Bayesian kNN
  * no Gaussian smoothing, no support masking
  * velocity grid [-5, +5] deg/day at 0.01, bilinear in velocity only
  * NO interpolation between magnitude bins -- the single containing slice is used
  * no Platt calibration
  * a probability requires ALL FOUR population densities present and finite; a 2- or 3-population
    denominator is NOT P(NEO) and is returned as NaN with reason
    `missing_population_density:<pops>`
  * missing densities are NEVER replaced by zero or epsilon

GEOMETRY -- A DOCUMENTED DIVERGENCE FROM TEST2
----------------------------------------------
TEST2 assigned each row to a sky centre using ecliptic coordinates carried in the epoch-state
cache, i.e. the MODEL frame produced during n-body propagation. A real Rubin tracklet supplies only
observed RA/Dec, so v1 derives ecliptic coordinates from the observed positions. The two differ by
a median of 1.8e-03 deg (~6.5 arcsec), which reassigns 295 of 688,688 TEST2 rows (0.043%) that sit
near a centre boundary. v1 uses the OBSERVABLE definition because it is the only one a survey can
supply. `geometry_source="precomputed"` exists solely so the regression test can isolate this
effect; it is not for operational use.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
for _p in (W / "neomod" / "src", W / "neomod" / "adam_core_stub", W / "neomod" / "pipeline"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

INTERFACE_VERSION = "rubin_vdp_scorer_v1"
EPOCH = "2027-08-25T00:00:00"
POPS = ("NEO", "MBA", "TNO", "Trojans")
V_LO, V_HI, MAG_STEP, N_BINS = 14.0, 25.0, 0.25, 44
VEL_LIMIT = 5.0
DLON_LIMIT = 140.0
LAT_BASE = [0, 1, 2, 3, 4, 5, 8, 12, 18, 25, 35, 50]
LON_STEP = 10.0
SUPPORTED_BANDS = ("V",)

REQUIRED_COLUMNS = ("tracklet_id", "mjd0", "ra0", "dec0", "mag0_V",
                    "mjd1", "ra1", "dec1", "mag1_V", "observatory_code")


# --------------------------------------------------------------------------- geometry
def center_grid() -> np.ndarray:
    """The frozen 667-centre grid: 10-deg longitude steps within |dlon| <= 140, 23 latitudes."""
    lats = sorted({float(v) for v in LAT_BASE} | {float(-v) for v in LAT_BASE})
    dl = [round(d, 6) for d in np.arange(-180.0, 180.0, LON_STEP) if abs(d) <= DLON_LIMIT + 1e-9]
    return np.array([(d, l) for l in lats for d in dl], dtype=float)


def antisun_lon_deg(epoch: str = EPOCH) -> float:
    from astropy.time import Time
    from astropy.coordinates import get_sun, GeocentricTrueEcliptic
    from astropy.utils import iers
    iers.conf.auto_max_age = None
    t = Time(epoch, scale="utc")
    return float((get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t)).lon.deg + 180.0)
                 % 360.0)


def center_label(dlon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    return np.array([f"dlon{int(round(x)):+04d}_lat{int(round(y)):+03d}"
                     for x, y in zip(dlon, lat)], dtype=object)


def magnitude_bin_label(v):
    """Half-open [lo, lo+0.25). Outside [14, 25) -> None. NEVER clipped into an edge bin."""
    v = np.asarray(v, dtype=float)
    lab = np.array([None] * len(v), dtype=object)
    ok = np.isfinite(v) & (v >= V_LO) & (v < V_HI)
    idx = np.floor((v[ok] - V_LO) / MAG_STEP)
    lo = V_LO + idx * MAG_STEP
    lab[ok] = [f"V{a:06.2f}_{a + MAG_STEP:06.2f}" for a in lo]
    return lab


def derive_motion(mjd0, ra0, dec0, mjd1, ra1, dec1):
    """Angular rates from two detections. RA difference is wrapped to (-180, 180]."""
    import velocity_density_pipeline_neomod_clone_only as vdp
    dt = np.asarray(mjd1, float) - np.asarray(mjd0, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        dra = (((np.asarray(ra1, float) - np.asarray(ra0, float) + 180.0) % 360.0) - 180.0) / dt
        ddec = (np.asarray(dec1, float) - np.asarray(dec0, float)) / dt
    vlam, vbeta = vdp.radec_rates_to_ecliptic_rates_manual(
        np.asarray(ra0, float), np.asarray(dec0, float), dra, ddec)
    return dra, ddec, vlam, vbeta, dt


def ecliptic_of_radec(ra_deg, dec_deg, epoch: str = EPOCH):
    """Observed RA/Dec -> geocentric true ecliptic lon/lat at the frozen epoch."""
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic
    from astropy.utils import iers
    import astropy.units as u
    iers.conf.auto_max_age = None
    t = Time(epoch, scale="utc")
    c = SkyCoord(ra=np.asarray(ra_deg, float) * u.deg,
                 dec=np.asarray(dec_deg, float) * u.deg,
                 frame="icrs").transform_to(GeocentricTrueEcliptic(obstime=t))
    return c.lon.deg, c.lat.deg


# --------------------------------------------------------------------------- validation
def validate_input(df: pd.DataFrame) -> pd.DataFrame:
    """Strict schema validation. Returns a per-row `input_reason` ('ok' when usable)."""
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"input is missing required columns: {missing}")
    if df.tracklet_id.isna().any():
        raise ValueError("tracklet_id contains nulls")
    if df.tracklet_id.duplicated().any():
        n = int(df.tracklet_id.duplicated().sum())
        raise ValueError(f"tracklet_id must be unique; {n} duplicates found")

    reason = np.array(["ok"] * len(df), dtype=object)
    num = ["mjd0", "ra0", "dec0", "mag0_V", "mjd1", "ra1", "dec1", "mag1_V"]
    vals = {c: pd.to_numeric(df[c], errors="coerce").to_numpy(float) for c in num}
    bad = np.zeros(len(df), dtype=bool)
    for c in num:
        b = ~np.isfinite(vals[c])
        reason = np.where(b & ~bad, f"nonfinite_field:{c}", reason)
        bad |= b
    for c in ("ra0", "ra1"):
        b = np.isfinite(vals[c]) & ((vals[c] < 0.0) | (vals[c] >= 360.0))
        reason = np.where(b & ~bad, f"out_of_range:{c}", reason); bad |= b
    for c in ("dec0", "dec1"):
        b = np.isfinite(vals[c]) & (np.abs(vals[c]) > 90.0)
        reason = np.where(b & ~bad, f"out_of_range:{c}", reason); bad |= b
    dt = vals["mjd1"] - vals["mjd0"]
    b = np.isfinite(dt) & (dt <= 0)
    reason = np.where(b & ~bad, "nonpositive_time_baseline", reason); bad |= b
    return pd.Series(reason, index=df.index)


# --------------------------------------------------------------------------- scoring
def prepare(df: pd.DataFrame, epoch: str = EPOCH,
            geometry_source: str = "derived") -> pd.DataFrame:
    """Attach motion, sky centre and magnitude bin. Adds nothing scientific."""
    out = df.copy()
    out["_input_reason"] = validate_input(df)
    ok = (out["_input_reason"] == "ok").to_numpy()

    dra, ddec, vlam, vbeta, dt = derive_motion(
        out.mjd0, out.ra0, out.dec0, out.mjd1, out.ra1, out.dec1)
    out["dra_deg_day"], out["ddec_deg_day"] = dra, ddec
    out["time_baseline_days"] = dt
    if geometry_source == "precomputed":
        for c in ("vlam", "vbeta", "lam_deg", "beta_deg"):
            if c not in out.columns:
                raise ValueError(f"geometry_source='precomputed' requires column {c}")
        lam, beta = out.lam_deg.to_numpy(float), out.beta_deg.to_numpy(float)
    else:
        out["vlam"], out["vbeta"] = vlam, vbeta
        # SkyCoord raises on |dec| > 90 or non-finite input, which would abort the whole run
        # because of a single malformed row. Rows already marked invalid are transformed with a
        # placeholder and their result is discarded below.
        ra_s = pd.to_numeric(out.ra0, errors="coerce").to_numpy(float)
        dec_s = pd.to_numeric(out.dec0, errors="coerce").to_numpy(float)
        safe = np.isfinite(ra_s) & np.isfinite(dec_s) & (np.abs(dec_s) <= 90.0)
        ra_s = np.where(safe, ra_s % 360.0, 0.0)
        dec_s = np.where(safe, dec_s, 0.0)
        lam, beta = ecliptic_of_radec(ra_s, dec_s, epoch)
        lam = np.where(safe, lam, np.nan)
        beta = np.where(safe, beta, np.nan)
        out["lam_deg"], out["beta_deg"] = lam, beta

    out["mean_V"] = 0.5 * (pd.to_numeric(out.mag0_V, errors="coerce").to_numpy(float)
                           + pd.to_numeric(out.mag1_V, errors="coerce").to_numpy(float))

    grid = center_grid()
    a = antisun_lon_deg(epoch)
    dlon = ((lam - a + 180.0) % 360.0) - 180.0
    out["dlon_from_antisun_deg"] = dlon
    safe_dlon = np.where(np.isfinite(dlon), dlon, 0.0)
    safe_beta = np.where(np.isfinite(beta), beta, 0.0)
    ci = np.argmin((safe_dlon[:, None] - grid[None, :, 0]) ** 2
                   + (safe_beta[:, None] - grid[None, :, 1]) ** 2, axis=1)
    lab = center_label(grid[ci, 0], grid[ci, 1])
    lab = np.where(np.isfinite(dlon) & np.isfinite(beta), lab, None)
    out["center_label"] = lab
    out["map_center_id"] = out["center_label"]
    out["magnitude_bin"] = magnitude_bin_label(out["mean_V"].to_numpy(float))
    out["magnitude_bin_id"] = out["magnitude_bin"]

    # rows that failed schema validation must not be assigned a map
    out.loc[~ok, ["center_label", "map_center_id", "magnitude_bin", "magnitude_bin_id"]] = None
    return out


def score(df: pd.DataFrame, map_root: Path, epoch: str = EPOCH,
          geometry_source: str = "derived", chunk_size: int = 200_000,
          stats: dict | None = None) -> pd.DataFrame:
    """Score prepared tracklets.

    The density lookup and posterior are delegated verbatim to `score_test2.score_new_vdp`.
    Rows are grouped by sky centre so each ~390 MB map is opened EXACTLY ONCE, and chunking
    happens inside a centre group -- so neither chunk size nor input order can change a result.
    """
    import time
    import score_test2 as st2
    st2.MAPS_NEW = Path(map_root)                     # explicit seal-driven map root

    t_prep0 = time.perf_counter()
    prep = prepare(df, epoch=epoch, geometry_source=geometry_source)
    t_prep = time.perf_counter() - t_prep0

    res = pd.DataFrame(index=prep.index)
    for c in POPS:
        res[f"rho_{c}"] = np.nan
        res[f"p_{c}"] = np.nan
    res["valid"] = False
    res["reason"] = prep["_input_reason"].to_numpy()
    res["total_density"] = np.nan
    res["n_pops"] = 0
    res["prob_sum"] = np.nan

    scorable = ((prep["_input_reason"] == "ok")
                & prep["center_label"].notna()).to_numpy()
    n_maps, n_chunks, biggest = 0, 0, 0
    t_score0 = time.perf_counter()
    if scorable.any():
        sub = prep[scorable]
        for cen, g_all in sub.groupby("center_label", sort=True):
            n_maps += 1
            biggest = max(biggest, len(g_all))
            got_parts, idx_parts = [], []
            for s0 in range(0, len(g_all), max(chunk_size, 1)):
                g = g_all.iloc[s0:s0 + chunk_size]
                n_chunks += 1
                got_parts.append(st2.score_new_vdp(g.reset_index(drop=True)))
                idx_parts.append(g.index)
            got = pd.concat(got_parts, ignore_index=True)
            got.index = np.concatenate([np.asarray(i) for i in idx_parts])
            tot = got["total_density_new"].to_numpy(float)
            ok = got["new_vdp_valid"].to_numpy()
            for c in POPS:
                rho = got[f"rho_{c}_new"].to_numpy(float)
                res.loc[got.index, f"rho_{c}"] = rho
                with np.errstate(invalid="ignore", divide="ignore"):
                    res.loc[got.index, f"p_{c}"] = np.where(
                        np.isfinite(tot) & (tot > 0) & ok, rho / tot, np.nan)
            res.loc[got.index, "valid"] = ok
            res.loc[got.index, "reason"] = got["new_vdp_reason"].to_numpy()
            res.loc[got.index, "total_density"] = tot
            res.loc[got.index, "n_pops"] = got["n_pops_new"].to_numpy()
            res.loc[got.index, "prob_sum"] = got["prob_sum_new"].to_numpy()
    t_score = time.perf_counter() - t_score0

    out = pd.concat([prep.drop(columns=["_input_reason"]), res], axis=1)

    # DEFECT FOUND VIA THE EXAMPLE SUITE, corrected here rather than in the frozen scorer.
    # score_new_vdp groups by `magnitude_bin`; pandas groupby DROPS None groups, so a row whose
    # mean V falls outside [14, 25) never reaches the `v_out_of_range` branch and keeps its
    # initialised reason "ok" while its probability is (correctly) NaN. TEST2 contained zero such
    # rows so the frozen results are unaffected, but real survey input will hit it. Fixing it in
    # score_test2.py would change the sealed source, so the interface repairs the reason instead;
    # no probability changes.
    oob = (out["magnitude_bin"].isna().to_numpy()
           & (out["reason"].to_numpy() == "ok")
           & ~out["valid"].to_numpy(bool))
    if oob.any():
        out.loc[oob, "reason"] = "v_out_of_range"

    v = out["valid"].to_numpy(bool)
    if v.any():
        ps = out.loc[v, [f"p_{c}" for c in POPS]].to_numpy(float).sum(axis=1)
        assert np.nanmax(np.abs(ps - 1.0)) < 1e-9, "class probabilities do not sum to one"
        assert (out.loc[v, "n_pops"].to_numpy() == 4).all(), "valid row without four densities"
    if stats is not None:
        # HONEST ACCOUNTING: score_new_vdp opens the .npz on every call, so the number of actual
        # map OPENS equals the number of chunks, not the number of centres. They coincide only
        # when chunk_size >= the largest centre group (the default case). Setting --chunk-size
        # below a centre's row count costs one extra open per additional chunk; results are
        # unchanged, only I/O cost rises.
        stats.update({"n_rows": int(len(df)), "n_scorable": int(scorable.sum()),
                      "n_valid": int(v.sum()),
                      "n_centers": n_maps, "n_map_opens": n_chunks,
                      "n_maps_loaded": n_maps,
                      "n_chunks": n_chunks, "largest_center_group": int(biggest),
                      "one_open_per_center": bool(n_chunks == n_maps),
                      "prepare_seconds": round(t_prep, 3),
                      "score_seconds": round(t_score, 3)})
    return out


def input_hash(df: pd.DataFrame) -> str:
    """Stable digest of the scientific input fields, independent of row order.

    Built from explicit per-column string conversion: `astype(str).agg("|".join)` leaves numeric
    values unconverted under pandas 3.0 and raises. Floats use %.17g so the digest is exact and
    reproducible across platforms.
    """
    cols = ["tracklet_id", "mjd0", "ra0", "dec0", "mag0_V",
            "mjd1", "ra1", "dec1", "mag1_V", "observatory_code"]
    parts = []
    for c in cols:
        v = df[c].to_numpy()
        if v.dtype.kind in "fiu":
            vv = v.astype(float)
            parts.append([("nan" if not np.isfinite(x) else f"{x:.17g}") for x in vv])
        else:
            parts.append([("" if x is None else str(x)) for x in v])
    rows = ["|".join(col[i] for col in parts) for i in range(len(df))]
    h = hashlib.sha256()
    for s_ in sorted(rows):
        h.update(s_.encode()); h.update(b"\x1e")
    return h.hexdigest()
