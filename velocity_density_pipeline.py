"""
velocity_density_pipeline.py
============================

Refactored solar system population pipeline for building, saving, and
querying ecliptic velocity-space (v_lambda, v_beta) probability maps
from the S3M synthetic populations (MBA, NEO, TNO, Trojan).

This module is the single entry point for:

  1. Loading S3M populations (delegates to user's `s3m_loader` module)
  2. Generating per-population, per-magnitude-bin density maps via
     conditional K|M cloning + Bayesian kNN density estimation
     (Appendix B 2D posterior)
  3. Saving the resulting maps as compressed .npz files (one file per
     obstime + sky center)
  4. Loading saved maps and scoring arbitrary objects against them
     (single observations or bulk S3M arrays for ROC analysis)
  5. Plotting the per-population probability maps with the 0.2 deg/day
     nearest-object mask

Typical usage from a notebook:

    >>> import velocity_density_pipeline as vdp
    >>> # one-time generation for a chosen obstime
    >>> vdp.generate_probability_maps(
    ...     obstime_str="2025-03-21T00:00:00",
    ...     output_path="prob_maps_2025-03-21.npz",
    ... )
    >>> # later, in any other notebook
    >>> maps = vdp.ProbMapSet.from_npz("prob_maps_2025-03-21.npz")
    >>> probs = maps.score_visible(vlam_arr, vbeta_arr, mag_app_arr)
    >>> probs["NEO"]   # array of P(NEO) for each input object

Individual functions (cloning, propagation, density estimation, etc.)
remain importable for users who want to drop down a level.

NOTES ON FIDELITY TO THE ORIGINAL NOTEBOOK
------------------------------------------
- The conditional K|M cloner with sky-cut reapplication is the active
  cloning path for all four populations (per advisor decision).
- A uniform mean-anomaly fallback path is preserved but expects
  `nsc.clone_population_uniform_mean_anomaly` to exist; it is dead code
  in the default configuration.
- `build_visible_subset_dataframe` performs the equatorial->ecliptic
  rotation manually (the astropy GeocentricTrueEcliptic differential
  bug fix carried over from the notebook).
- `elements_to_vlam_vbeta_d_h` retains the legacy
  `nsc.radec_rates_to_ecliptic_rates_at_obstime` call; this function
  is only invoked on the uniform-cloner path which is not used by the
  default pipeline.
"""

# =============================================================================
# Imports
# =============================================================================
from __future__ import annotations

import gc
from math import lgamma
from typing import Optional, Sequence

import numpy as np
import pandas as pd

import astropy.units as u
from astropy.coordinates import (
    GCRS,
    GeocentricTrueEcliptic,
    SkyCoord,
    get_sun,
)
from astropy.time import Time
from astroML.density_estimation import EmpiricalDistribution
from scipy.spatial import cKDTree
from tqdm import tqdm

import neoscore as nsc
import s3m_loader as load_s3m


# =============================================================================
# Module-level constants
# =============================================================================
AU_KM = 149_597_870.7
MU_SUN = 1.32712440018e11   # km^3 / s^2
ECLIPTIC_OBLIQUITY_DEG = 23.439291

# default density-map grid (matches the notebook)
DEFAULT_GRID_LIM = (-0.8, 0.8)
DEFAULT_GRID_STEP = 0.01

# default density-estimator tuning
DEFAULT_K_MAP = 5
DEFAULT_N_D0_GRID_MAP = 400

# default sky-cut geometry
DEFAULT_MAX_SEP_DEG = 30.0
DEFAULT_CENTER_LON_DEG = 180.0
DEFAULT_CENTER_LAT_DEG = 0.0

# default mag-bin layout from the notebook
DEFAULT_MAG_BINS = [
    {"label": "14_16", "mag_min": 14.0, "mag_max": 16.0},
    {"label": "16_18", "mag_min": 16.0, "mag_max": 18.0},
    {"label": "18_20", "mag_min": 18.0, "mag_max": 20.0},
    {"label": "mag20", "mag_min": 20.0, "mag_max": 21.0},
    {"label": "mag21", "mag_min": 21.0, "mag_max": 22.0},
    {"label": "mag22", "mag_min": 22.0, "mag_max": 23.0},
    {"label": "mag23", "mag_min": 23.0, "mag_max": 24.0},
    {"label": "mag24+", "mag_min": 24.0, "mag_max": 26.0},
]

# default per-population settings (clone_factor and overlay style)
DEFAULT_POPULATION_SETTINGS = {
    "MBA":     {"s3m_pop": "mba",    "clone_factor": 100, "scatter_size": 1, "scatter_alpha": 0.05},
    "NEO":     {"s3m_pop": "neo",    "clone_factor": 300, "scatter_size": 4, "scatter_alpha": 0.08},
    "TNO":     {"s3m_pop": "tno",    "clone_factor": 100, "scatter_size": 6, "scatter_alpha": 0.10},
    "Trojans": {"s3m_pop": "trojan", "clone_factor": 100, "scatter_size": 6, "scatter_alpha": 0.10},
}

# default mask radius for probability map plotting (deg/day)
DEFAULT_MASK_RADIUS_DEG_PER_DAY = 0.2


# =============================================================================
# Grid helpers
# =============================================================================
def make_default_grid(grid_lim=DEFAULT_GRID_LIM, grid_step=DEFAULT_GRID_STEP):
    """Build the (vlam, vbeta) evaluation grid used by the density estimator.

    Returns
    -------
    x_grid, y_grid : 1D ndarrays
    X0, Y0 : 2D meshgrids in 'xy' indexing
    grid_points : ndarray, shape (N, 2)
    """
    x_grid = np.arange(grid_lim[0], grid_lim[1] + grid_step, grid_step)
    y_grid = np.arange(grid_lim[0], grid_lim[1] + grid_step, grid_step)
    X0, Y0 = np.meshgrid(x_grid, y_grid, indexing="xy")
    grid_points = np.column_stack((X0.ravel(), Y0.ravel()))
    return x_grid, y_grid, X0, Y0, grid_points


# =============================================================================
# S3M loading
# =============================================================================
def load_s3m_population(s3m_pop_name, verbose=False):
    """Load a single S3M population and build its scorer.

    Parameters
    ----------
    s3m_pop_name : str
        One of "mba", "neo", "tno", "trojan" (whatever
        `s3m_loader.define_s3m` accepts).
    verbose : bool

    Returns
    -------
    df : pandas.DataFrame
        Filtered, augmented S3M DataFrame.
    scorer : nsc.NEOMODScorer
        Population-specific scorer (used for observer-state queries
        and apparent-magnitude calculations).
    """
    df = load_s3m.define_s3m(pop=s3m_pop_name, verbose=verbose)
    df, array4D, H_center, a_center, e_center, i_center = load_s3m.s3m_array(df)
    scorer = nsc.NEOMODScorer(array4D, H_center, a_center, e_center, i_center)
    return df, scorer


def load_all_populations(
    population_settings=None,
    verbose=False,
):
    """Load every S3M population listed in `population_settings`.

    Parameters
    ----------
    population_settings : dict or None
        Mapping pop_label -> {"s3m_pop": ..., "clone_factor": ..., ...}.
        Defaults to DEFAULT_POPULATION_SETTINGS.
    verbose : bool

    Returns
    -------
    clone_sources : dict
        pop_label -> {"df", "scorer", "clone_factor", "use_conditional_cloner",
                      "scatter_size", "scatter_alpha", ...}
        Drop-in replacement for the notebook's `clone_sources` dict.
    """
    if population_settings is None:
        population_settings = DEFAULT_POPULATION_SETTINGS

    clone_sources = {}
    for pop_label, cfg in population_settings.items():
        if verbose:
            print(f"[load_all_populations] loading {pop_label} ({cfg['s3m_pop']})...")
        df, scorer = load_s3m_population(cfg["s3m_pop"], verbose=verbose)
        clone_sources[pop_label] = {
            "df": df,
            "scorer": scorer,
            "clone_factor": cfg["clone_factor"],
            "use_conditional_cloner": cfg.get("use_conditional_cloner", True),
            "scatter_size": cfg.get("scatter_size", 4),
            "scatter_alpha": cfg.get("scatter_alpha", 0.1),
        }
        if verbose:
            print(f"  loaded {len(df)} objects.")

    return clone_sources


# =============================================================================
# Geometry helpers
# =============================================================================
def wrap360_deg(x):
    return np.mod(x, 360.0)


def mean_anomaly_deg_at_obstime(a_AU, tp_mjd, obstime_str):
    """Mean anomaly (deg) at obstime for a Keplerian orbit specified
    by semi-major axis and time of perihelion."""
    t_obs = Time(obstime_str, scale="tdb")
    a_AU = np.asarray(a_AU, dtype=float)
    tp_mjd = np.asarray(tp_mjd, dtype=float)
    a_km = a_AU * AU_KM
    n_rad_s = np.sqrt(MU_SUN / a_km**3)
    n_rad_day = n_rad_s * 86400.0
    M_rad = n_rad_day * (t_obs.mjd - tp_mjd)
    return wrap360_deg(np.rad2deg(M_rad))


def mean_longitude_deg_at_obstime(a_AU, node_deg, argperi_deg, tp_mjd, obstime_str):
    """Mean longitude L = node + argperi + M (deg) at obstime."""
    M_deg = mean_anomaly_deg_at_obstime(a_AU, tp_mjd, obstime_str)
    return wrap360_deg(node_deg + argperi_deg + M_deg)


# =============================================================================
# Apparent magnitude (whole sky, no positional cut)
# =============================================================================
def compute_apparent_magnitude_for_population(
    df,
    obstime_str,
    scorer,
    chunk=100_000,
    show_progress=False,
):
    """Compute apparent V magnitude for every object in df at obstime,
    BEFORE any sky-position cut, using the HG phase function (G=0.15).

    The returned array has the same length as df.
    Non-finite or infinite results indicate propagation failure.
    """
    r_ecl, v_ecl = nsc.elements_to_helio_ecliptic_state(
        a_AU=df["a"].to_numpy(dtype=float),
        e=df["e"].to_numpy(dtype=float),
        inc_deg=df["i"].to_numpy(dtype=float),
        raan_deg=df["node"].to_numpy(dtype=float),
        argp_deg=df["argperi"].to_numpy(dtype=float),
        tp_mjd=df["t_p"].to_numpy(dtype=float),
        obstime_str=obstime_str,
        method="newton",
        n_iter=10,
        chunk=chunk,
        show_progress=show_progress,
    )

    eps = np.deg2rad(ECLIPTIC_OBLIQUITY_DEG)
    c, s = np.cos(eps), np.sin(eps)
    r_eq = np.column_stack([
        r_ecl[:, 0],
        c * r_ecl[:, 1] - s * r_ecl[:, 2],
        s * r_ecl[:, 1] + c * r_ecl[:, 2],
    ])

    _, rE, vE, r_tele = scorer._get_earth_and_observer(obstime_str)

    r_helio_vec = r_eq
    r_obs_vec = r_eq - (rE + r_tele)

    r_sun = np.linalg.norm(r_helio_vec, axis=1)
    Delta = np.linalg.norm(r_obs_vec, axis=1)

    # convert km -> AU if pipeline returned km
    if np.nanmedian(r_sun) > 1000:
        r_sun = r_sun / AU_KM
        Delta = Delta / AU_KM

    u_sun = -r_helio_vec / np.linalg.norm(r_helio_vec, axis=1, keepdims=True)
    u_obs = -r_obs_vec / np.linalg.norm(r_obs_vec, axis=1, keepdims=True)

    cos_alpha = np.sum(u_sun * u_obs, axis=1)
    cos_alpha = np.clip(cos_alpha, -1.0, 1.0)
    alpha = np.arccos(cos_alpha)

    G = 0.15
    A1, A2, B1, B2 = 3.33, 1.87, 0.63, 1.22
    phi1 = np.exp(-A1 * np.tan(0.5 * alpha)**B1)
    phi2 = np.exp(-A2 * np.tan(0.5 * alpha)**B2)
    Phi = (1 - G) * phi1 + G * phi2
    Phi = np.maximum(Phi, 1e-30)

    H = df["H"].to_numpy(dtype=float)
    mag_app = H + 5.0 * np.log10(r_sun * Delta) - 2.5 * np.log10(Phi)
    return mag_app


def select_df_by_mag_bin(df, mag_app, mag_min=None, mag_max=None):
    """Boolean-select rows of df by apparent magnitude.

    Returns
    -------
    df_sel : DataFrame (copy of selected rows, reset index)
    mag_sel : ndarray (apparent mags for selected rows, aligned with df_sel)
    """
    sel = np.isfinite(mag_app)
    if mag_min is not None:
        sel &= (mag_app >= mag_min)
    if mag_max is not None:
        sel &= (mag_app < mag_max)
    df_sel = df.loc[sel].copy()
    mag_sel = mag_app[sel]
    return df_sel, mag_sel


# =============================================================================
# Sky-patch cut + orbital-phase visible subset
# =============================================================================
def build_visible_subset_dataframe(
    df,
    obstime_str,
    scorer,
    max_sep_deg=DEFAULT_MAX_SEP_DEG,
    chunk=100_000,
    show_progress=True,
    center_mode="custom_ecliptic",
    center_lon_deg=DEFAULT_CENTER_LON_DEG,
    center_lat_deg=DEFAULT_CENTER_LAT_DEG,
):
    """Propagate df to obstime, apply the sky cut, and return only the
    visible objects with vlam, vbeta, ra, dec, dra, ddec, etc. attached.

    Uses a manual equatorial->ecliptic rotation to compute vlam, vbeta
    so it does not lose differentials to the astropy
    GeocentricTrueEcliptic frame transform.
    """
    t_obs = Time(obstime_str, scale="tdb")

    orig_idx = np.arange(len(df), dtype=int)

    a_AU = df["a"].to_numpy(dtype=float)
    e = df["e"].to_numpy(dtype=float)
    inc = df["i"].to_numpy(dtype=float)
    raan = df["node"].to_numpy(dtype=float)
    argp = df["argperi"].to_numpy(dtype=float)
    tp_mjd = df["t_p"].to_numpy(dtype=float)

    r_obj_ecl, v_obj_ecl = nsc.elements_to_helio_ecliptic_state(
        a_AU=a_AU, e=e, inc_deg=inc,
        raan_deg=raan, argp_deg=argp, tp_mjd=tp_mjd,
        obstime_str=obstime_str,
        method="newton", n_iter=10,
        chunk=chunk, show_progress=show_progress,
    )

    eps = np.deg2rad(ECLIPTIC_OBLIQUITY_DEG)
    c, s = np.cos(eps), np.sin(eps)

    r_obj_helio = np.column_stack([
        r_obj_ecl[:, 0],
        c * r_obj_ecl[:, 1] - s * r_obj_ecl[:, 2],
        s * r_obj_ecl[:, 1] + c * r_obj_ecl[:, 2],
    ])
    v_obj_helio = np.column_stack([
        v_obj_ecl[:, 0],
        c * v_obj_ecl[:, 1] - s * v_obj_ecl[:, 2],
        s * v_obj_ecl[:, 1] + c * v_obj_ecl[:, 2],
    ])

    _, rE, vE, r_tele = scorer._get_earth_and_observer(obstime_str)
    r_rel = r_obj_helio - (rE + r_tele)
    v_rel = v_obj_helio - vE

    x, y, z = r_rel[:, 0], r_rel[:, 1], r_rel[:, 2]
    vx, vy, vz = v_rel[:, 0], v_rel[:, 1], v_rel[:, 2]

    rho2 = x*x + y*y
    r2 = rho2 + z*z
    rho = np.sqrt(rho2)

    ra = np.arctan2(y, x)
    dec = np.arctan2(z, rho)

    dra = (x*vy - y*vx) / rho2
    ddec = (vz*rho2 - z*(x*vx + y*vy)) / (r2 * rho)

    good = (
        np.isfinite(ra) & np.isfinite(dec) &
        np.isfinite(dra) & np.isfinite(ddec) &
        np.isfinite(rho2) & (rho2 > 0) &
        np.isfinite(r2) & (r2 > 0) &
        np.isfinite(rho) & (rho > 0)
    )

    orig_idx_good = orig_idx[good]

    ra_deg = np.rad2deg(ra[good]) % 360.0
    dec_deg = np.rad2deg(dec[good])
    dra_deg_day = np.rad2deg(dra[good]) * 86400.0
    ddec_deg_day = np.rad2deg(ddec[good]) * 86400.0

    sc = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame=GCRS(obstime=t_obs))

    if center_mode == "antisun":
        sun_gcrs = get_sun(t_obs).transform_to(GCRS(obstime=t_obs))
        center_sc = SkyCoord(
            x=-sun_gcrs.cartesian.x,
            y=-sun_gcrs.cartesian.y,
            z=-sun_gcrs.cartesian.z,
            representation_type="cartesian",
            frame=GCRS(obstime=t_obs),
        )
    elif center_mode == "custom_ecliptic":
        center_ecl = SkyCoord(
            lon=center_lon_deg*u.deg, lat=center_lat_deg*u.deg,
            distance=1.0*u.AU,
            frame=GeocentricTrueEcliptic(obstime=t_obs),
        )
        center_sc = center_ecl.transform_to(GCRS(obstime=t_obs))
    else:
        raise ValueError("center_mode must be 'antisun' or 'custom_ecliptic'")

    sky_sep_deg = sc.separation(center_sc).deg
    m = sky_sep_deg <= max_sep_deg

    sel_idx = orig_idx_good[m]
    out = df.iloc[sel_idx].copy().reset_index(drop=True)

    out["ra_deg"] = ra_deg[m]
    out["dec_deg"] = dec_deg[m]
    out["dra_deg_day"] = dra_deg_day[m]
    out["ddec_deg_day"] = ddec_deg_day[m]
    out["sky_sep_deg"] = sky_sep_deg[m]

    # equatorial -> ecliptic rotation (manual; preserves differentials)
    r_rel_sel = r_rel[good][m]
    v_rel_sel = v_rel[good][m]

    x_ecl = r_rel_sel[:, 0]
    y_ecl = c * r_rel_sel[:, 1] + s * r_rel_sel[:, 2]
    z_ecl = -s * r_rel_sel[:, 1] + c * r_rel_sel[:, 2]

    vx_ecl = v_rel_sel[:, 0]
    vy_ecl = c * v_rel_sel[:, 1] + s * v_rel_sel[:, 2]
    vz_ecl = -s * v_rel_sel[:, 1] + c * v_rel_sel[:, 2]

    rho2_ecl = x_ecl*x_ecl + y_ecl*y_ecl
    r2_ecl = rho2_ecl + z_ecl*z_ecl
    rho_ecl = np.sqrt(rho2_ecl)

    lam = np.arctan2(y_ecl, x_ecl)
    beta = np.arctan2(z_ecl, rho_ecl)

    dlam = (x_ecl*vy_ecl - y_ecl*vx_ecl) / rho2_ecl
    dbeta = (vz_ecl*rho2_ecl - z_ecl*(x_ecl*vx_ecl + y_ecl*vy_ecl)) / (r2_ecl*rho_ecl)

    out["lam_deg"] = np.rad2deg(lam) % 360.0
    out["beta_deg"] = np.rad2deg(beta)
    out["vlam"] = np.rad2deg(dlam) * 86400.0
    out["vbeta"] = np.rad2deg(dbeta) * 86400.0

    out["M_obs_deg"] = mean_anomaly_deg_at_obstime(
        out["a"].to_numpy(), out["t_p"].to_numpy(), obstime_str,
    )
    out["L_obs_deg"] = mean_longitude_deg_at_obstime(
        out["a"].to_numpy(), out["node"].to_numpy(),
        out["argperi"].to_numpy(), out["t_p"].to_numpy(),
        obstime_str,
    )

    return out


# =============================================================================
# Conditional K|M cloner (Trojan / TNO friendly)
# =============================================================================
def clone_population_conditional_K_from_M(
    df,
    clone_factor,
    obstime_str,
    bin_width_deg=10.0,
    unwrap_threshold_deg=180.0,
    rng=None,
    return_diagnostics=False,
):
    """Clone each orbit clone_factor times using:

      1. an empirical distribution of mean anomaly M(obstime), and
      2. a conditional empirical distribution of K = node + argperi
         given M.

    Designed for populations where preserving the M -- (node+argperi)
    relation is important (Trojans, TNOs). Output orbital elements
    keep node fixed and place new argperi = K_clone - node_parent.

    Returns raw cloned orbital elements; reapply the sky cut afterward
    via build_visible_subset_dataframe to get the visible cloned sample.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    t_obs = Time(obstime_str, scale="tdb")

    a_AU = df["a"].to_numpy(dtype=np.float64)
    e = df["e"].to_numpy(dtype=np.float32)
    inc = df["i"].to_numpy(dtype=np.float32)
    raan = df["node"].to_numpy(dtype=np.float32)
    argp = df["argperi"].to_numpy(dtype=np.float32)
    tp = df["t_p"].to_numpy(dtype=np.float64)

    a_rep = np.repeat(a_AU, clone_factor).astype(np.float64)
    e_rep = np.repeat(e, clone_factor).astype(np.float32)
    inc_rep = np.repeat(inc, clone_factor).astype(np.float32)
    raan_rep = np.repeat(raan, clone_factor).astype(np.float32)

    a_km = a_AU * AU_KM
    n_rad_s = np.sqrt(MU_SUN / (a_km ** 3))
    n_rad_day = n_rad_s * 86400.0

    M_obs = np.mod((t_obs.mjd - tp) * n_rad_day, 2.0 * np.pi)

    # circular shift so M is not split at 0 / 2pi
    mean_angle = np.mod(np.angle(np.mean(np.exp(1j * M_obs))), 2.0 * np.pi)
    M_shift = np.mod(M_obs - mean_angle, 2.0 * np.pi)

    M_shift_clone = EmpiricalDistribution(M_shift).rvs(len(a_rep))
    M_clone = np.mod(M_shift_clone + mean_angle, 2.0 * np.pi)

    n_rep = np.repeat(n_rad_day, clone_factor)
    tp_out = (t_obs.mjd - (M_clone / n_rep)).astype(np.float32)

    # conditional K | M sampling
    M_parent_deg = np.rad2deg(M_obs) % 360.0
    M_clone_deg = np.rad2deg(M_clone) % 360.0
    K_parent_deg = (raan.astype(np.float64) + argp.astype(np.float64)) % 360.0

    K_parent_unwrapped = K_parent_deg.copy()
    K_parent_unwrapped[K_parent_unwrapped < unwrap_threshold_deg] += 360.0

    edges = np.arange(0.0, 360.0 + bin_width_deg, bin_width_deg)
    nbins = len(edges) - 1

    parent_bin = np.clip(np.digitize(M_parent_deg, edges) - 1, 0, nbins - 1)
    clone_bin = np.clip(np.digitize(M_clone_deg, edges) - 1, 0, nbins - 1)

    K_by_bin = [K_parent_unwrapped[parent_bin == i] for i in range(nbins)]
    nonempty = np.array([len(v) > 0 for v in K_by_bin])
    nonempty_idx = np.where(nonempty)[0]

    K_clone_unwrapped = np.empty_like(M_clone_deg, dtype=np.float64)
    for j, b in enumerate(clone_bin):
        if len(K_by_bin[b]) == 0:
            nearest = nonempty_idx[np.argmin(np.abs(nonempty_idx - b))]
            draw_from = K_by_bin[nearest]
        else:
            draw_from = K_by_bin[b]
        K_clone_unwrapped[j] = rng.choice(draw_from)

    K_clone_deg = K_clone_unwrapped % 360.0
    argp_new = ((K_clone_deg - raan_rep.astype(np.float64)) % 360.0).astype(np.float32)

    if return_diagnostics:
        diagnostics = {
            "M_parent_deg": M_parent_deg,
            "M_clone_deg": M_clone_deg,
            "K_parent_deg": K_parent_deg,
            "K_parent_unwrapped": K_parent_unwrapped,
            "K_clone_deg": K_clone_deg,
            "K_clone_unwrapped": K_clone_unwrapped,
            "parent_bin": parent_bin,
            "clone_bin": clone_bin,
            "mean_angle_rad": mean_angle,
        }
        return (a_rep.astype(np.float32), e_rep, inc_rep, raan_rep, argp_new,
                tp_out, diagnostics)

    return (a_rep.astype(np.float32), e_rep, inc_rep, raan_rep, argp_new, tp_out)


def clone_population_conditional_K_from_M_with_skycut(
    df,
    clone_factor,
    obstime_str,
    scorer,
    max_sep_deg=DEFAULT_MAX_SEP_DEG,
    center_mode="custom_ecliptic",
    center_lon_deg=DEFAULT_CENTER_LON_DEG,
    center_lat_deg=DEFAULT_CENTER_LAT_DEG,
    chunk=100_000,
    show_progress=True,
    bin_width_deg=10.0,
    unwrap_threshold_deg=180.0,
    rng=None,
    return_raw_clone_df=False,
    return_diagnostics=False,
):
    """Clone via conditional K|M sampling, then reapply the same sky cut."""
    if rng is None:
        rng = np.random.default_rng(0)

    if return_diagnostics:
        (a_rep, e_rep, inc_rep, raan_rep, argp_rep, tp_out,
         diagnostics) = clone_population_conditional_K_from_M(
            df=df, clone_factor=clone_factor, obstime_str=obstime_str,
            bin_width_deg=bin_width_deg,
            unwrap_threshold_deg=unwrap_threshold_deg,
            rng=rng, return_diagnostics=True,
        )
    else:
        (a_rep, e_rep, inc_rep, raan_rep, argp_rep, tp_out
         ) = clone_population_conditional_K_from_M(
            df=df, clone_factor=clone_factor, obstime_str=obstime_str,
            bin_width_deg=bin_width_deg,
            unwrap_threshold_deg=unwrap_threshold_deg,
            rng=rng, return_diagnostics=False,
        )

    raw_clone_df = pd.DataFrame({
        "a": a_rep,
        "e": e_rep,
        "i": inc_rep,
        "node": raan_rep,
        "argperi": argp_rep,
        "t_p": tp_out,
    })

    clone_visible_df = build_visible_subset_dataframe(
        raw_clone_df,
        obstime_str=obstime_str,
        scorer=scorer,
        max_sep_deg=max_sep_deg,
        chunk=chunk, show_progress=show_progress,
        center_mode=center_mode,
        center_lon_deg=center_lon_deg,
        center_lat_deg=center_lat_deg,
    )

    outputs = [clone_visible_df]
    if return_raw_clone_df:
        outputs.append(raw_clone_df)
    if return_diagnostics:
        outputs.append(diagnostics)
    return outputs[0] if len(outputs) == 1 else tuple(outputs)


# =============================================================================
# Legacy uniform-mean-anomaly path (kept for backward compatibility,
# requires nsc.clone_population_uniform_mean_anomaly to be importable;
# only used when use_conditional_cloner=False)
# =============================================================================
def elements_to_vlam_vbeta_d_h(
    a_AU, e, inc, raan, argp, tp_mjd,
    obstime_str,
    scorer,
    H_vals,
    max_sep_deg=DEFAULT_MAX_SEP_DEG,
    chunk=100_000,
    show_progress=True,
    center_mode="custom_ecliptic",
    center_lon_deg=DEFAULT_CENTER_LON_DEG,
    center_lat_deg=DEFAULT_CENTER_LAT_DEG,
):
    """Convert orbital elements -> visible (vlam, vbeta, d_au, H_visible).

    Used on the uniform-cloner code path. Retains the legacy
    `nsc.radec_rates_to_ecliptic_rates_at_obstime` call.
    """
    t_obs = Time(obstime_str, scale="tdb")

    r_obj_ecl, v_obj_ecl = nsc.elements_to_helio_ecliptic_state(
        a_AU=a_AU, e=e, inc_deg=inc,
        raan_deg=raan, argp_deg=argp, tp_mjd=tp_mjd,
        obstime_str=obstime_str,
        method="newton", n_iter=10,
        chunk=chunk, show_progress=show_progress,
    )

    eps = np.deg2rad(ECLIPTIC_OBLIQUITY_DEG)
    c, s = np.cos(eps), np.sin(eps)

    r_obj_helio = np.column_stack([
        r_obj_ecl[:, 0],
        c * r_obj_ecl[:, 1] - s * r_obj_ecl[:, 2],
        s * r_obj_ecl[:, 1] + c * r_obj_ecl[:, 2],
    ])
    v_obj_helio = np.column_stack([
        v_obj_ecl[:, 0],
        c * v_obj_ecl[:, 1] - s * v_obj_ecl[:, 2],
        s * v_obj_ecl[:, 1] + c * v_obj_ecl[:, 2],
    ])

    _, rE, vE, r_tele = scorer._get_earth_and_observer(obstime_str)
    r_rel = r_obj_helio - (rE + r_tele)
    v_rel = v_obj_helio - vE

    d_raw = np.linalg.norm(r_rel, axis=1)
    if np.nanmedian(d_raw) > 1000:
        d_au_all = d_raw / AU_KM
    else:
        d_au_all = d_raw.copy()

    x, y, z = r_rel[:, 0], r_rel[:, 1], r_rel[:, 2]
    vx, vy, vz = v_rel[:, 0], v_rel[:, 1], v_rel[:, 2]

    rho2 = x*x + y*y
    r2 = rho2 + z*z
    rho = np.sqrt(rho2)

    ra = np.arctan2(y, x)
    dec = np.arctan2(z, rho)
    dra = (x*vy - y*vx) / rho2
    ddec = (vz*rho2 - z*(x*vx + y*vy)) / (r2 * rho)

    good = (
        np.isfinite(ra) & np.isfinite(dec) &
        np.isfinite(dra) & np.isfinite(ddec) &
        (rho2 > 0) & np.isfinite(r2) & (r2 > 0)
    )

    ra_deg = (np.rad2deg(ra[good]) % 360.0)
    dec_deg = np.rad2deg(dec[good])
    dra_deg_day = np.rad2deg(dra[good]) * 86400.0
    ddec_deg_day = np.rad2deg(ddec[good]) * 86400.0

    d_au_good = d_au_all[good]
    H_good = np.asarray(H_vals, dtype=float)[good]

    sc = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame=GCRS(obstime=t_obs))

    if center_mode == "antisun":
        sun_gcrs = get_sun(t_obs).transform_to(GCRS(obstime=t_obs))
        center_sc = SkyCoord(
            x=-sun_gcrs.cartesian.x,
            y=-sun_gcrs.cartesian.y,
            z=-sun_gcrs.cartesian.z,
            representation_type="cartesian",
            frame=GCRS(obstime=t_obs),
        )
    elif center_mode == "custom_ecliptic":
        center_ecl = SkyCoord(
            lon=center_lon_deg*u.deg, lat=center_lat_deg*u.deg,
            distance=1.0*u.AU,
            frame=GeocentricTrueEcliptic(obstime=t_obs),
        )
        center_sc = center_ecl.transform_to(GCRS(obstime=t_obs))
    else:
        raise ValueError("center_mode must be 'antisun' or 'custom_ecliptic'")

    sep = sc.separation(center_sc).deg
    m = sep <= max_sep_deg

    ra_deg = ra_deg[m]; dec_deg = dec_deg[m]
    dra_deg_day = dra_deg_day[m]; ddec_deg_day = ddec_deg_day[m]
    d_au_vis = d_au_good[m]; H_vis = H_good[m]

    vlam, vbeta = nsc.radec_rates_to_ecliptic_rates_at_obstime(
        ra_deg, dec_deg, dra_deg_day, ddec_deg_day, obstime_str,
    )

    finite = np.isfinite(vlam) & np.isfinite(vbeta) & np.isfinite(d_au_vis) & np.isfinite(H_vis)
    return (
        vlam[finite].astype(np.float32),
        vbeta[finite].astype(np.float32),
        d_au_vis[finite].astype(np.float32),
        H_vis[finite].astype(np.float32),
    )


# =============================================================================
# Bayesian kNN density estimation (Appendix B 2D)
# =============================================================================
def log_posterior_d0_2d(d0_grid, d_knn):
    """log p(d0 | d_1,...,d_N) for the 2D Appendix B nearest-neighbor
    density estimator. Unnormalized."""
    d0_grid = np.asarray(d0_grid, dtype=float)
    d_knn = np.asarray(d_knn, dtype=float)
    if np.any(d0_grid <= 0):
        raise ValueError("All d0_grid values must be > 0.")
    if np.any(d_knn <= 0):
        raise ValueError("All nearest-neighbor distances must be > 0.")

    N = len(d_knn)
    k_vals = np.arange(1, N + 1)
    logp = np.zeros_like(d0_grid, dtype=float)

    for i, d0 in enumerate(d0_grid):
        term_sum = 0.0
        for k, dk in zip(k_vals, d_knn):
            term = (
                np.log(2.0)
                - np.log(d0)
                - lgamma(k)
                - (dk / d0) ** 2
                + (2 * k - 1) * np.log(dk / d0)
            )
            term_sum += term
        logp[i] = term_sum
    return logp


def normalize_posterior_from_loggrid(d0_grid, logp):
    """Trapezoid-rule normalization of an unnormalized log-posterior."""
    d0_grid = np.asarray(d0_grid, dtype=float)
    logp = np.asarray(logp, dtype=float)

    logp_shift = logp - np.max(logp)
    p_unnorm = np.exp(logp_shift)
    norm = np.trapezoid(p_unnorm, d0_grid)
    if norm <= 0 or not np.isfinite(norm):
        raise ValueError("Posterior normalization failed.")
    return p_unnorm / norm


def get_knn_distances(tree, x0, y0, k=10):
    """k-NN distances from a single evaluation point (x0, y0) into tree."""
    q = np.array([[x0, y0]])
    d, idx = tree.query(q, k=k)
    return d[0], idx[0]


def estimate_density_full_posterior_2d(
    tree,
    x0, y0,
    k=10,
    d0_min_factor=5.0,
    d0_max_factor=10.0,
    n_d0_grid=800,
    return_diagnostics=False,
):
    """Posterior-mean local 2D density at one evaluation point."""
    d_knn, idx_knn = get_knn_distances(tree, x0, y0, k=k)

    d0_min = d_knn[0] / d0_min_factor
    d0_max = d_knn[-1] * d0_max_factor
    d0_grid = np.logspace(np.log10(d0_min), np.log10(d0_max), n_d0_grid)

    logp_d0 = log_posterior_d0_2d(d0_grid, d_knn)
    p_d0 = normalize_posterior_from_loggrid(d0_grid, logp_d0)

    n0_grid = 1.0 / (np.pi * d0_grid**2)
    d0_map = d0_grid[np.argmax(p_d0)]
    d0_mean = np.trapezoid(d0_grid * p_d0, d0_grid)
    n0_map = 1.0 / (np.pi * d0_map**2)
    n0_mean = np.trapezoid(n0_grid * p_d0, d0_grid)

    if return_diagnostics:
        diag = {
            "d_knn": d_knn, "d0_grid": d0_grid, "p_d0": p_d0,
            "d0_map": d0_map, "d0_mean": d0_mean,
            "n0_map": n0_map, "n0_mean": n0_mean,
        }
        return n0_mean, diag
    return n0_mean


def evaluate_density_map_full_posterior_2d(
    tree,
    grid_points,
    k=10,
    d0_min_factor=5.0,
    d0_max_factor=10.0,
    n_d0_grid=400,
    show_progress=True,
):
    """Loop estimate_density_full_posterior_2d over all evaluation points."""
    density_vals = np.empty(len(grid_points), dtype=float)

    iterator = enumerate(grid_points)
    if show_progress:
        iterator = tqdm(iterator, total=len(grid_points), desc="Evaluating density map")

    for i, (x0, y0) in iterator:
        density_vals[i] = estimate_density_full_posterior_2d(
            tree, x0, y0, k=k,
            d0_min_factor=d0_min_factor,
            d0_max_factor=d0_max_factor,
            n_d0_grid=n_d0_grid,
            return_diagnostics=False,
        )
    return density_vals


# =============================================================================
# Build cloned density maps for one (sky-center, mag-bin) cell
# =============================================================================
def build_cloned_maps_for_center_magbin(
    center_lon_deg,
    center_lat_deg,
    center_label,
    clone_sources,
    obstime_str,
    max_sep_deg,
    rng,
    overlay_max,
    grid_points,
    X0,
    k_map,
    n_d0_grid_map,
    mag_min=None,
    mag_max=None,
):
    """Build per-population cloned density maps for one mag bin.

    Mirrors the notebook's `build_cloned_maps_for_center_magbin`.
    Each population goes: full df -> mag cut -> sky cut -> conditional
    cloner -> kNN density estimate -> downweighted map. Also stores the
    nearest-object distance map for later masking.
    """
    clone_overlay = {}
    clone_density_maps = {}
    clone_density_maps_downweighted = {}
    density_maps_downweighted_raw = {}
    nearest_dist_maps = {}
    magcut_counts = {}

    for pop_name, info in clone_sources.items():
        f = info["clone_factor"]
        df = info["df"]

        print(f"\n[{center_label}] {pop_name} with mag cut {mag_min} to {mag_max}...")

        mag_app = compute_apparent_magnitude_for_population(
            df=df, obstime_str=obstime_str, scorer=info["scorer"],
            chunk=100_000, show_progress=True,
        )

        df_sel, mag_app_cut = select_df_by_mag_bin(
            df=df, mag_app=mag_app, mag_min=mag_min, mag_max=mag_max,
        )

        magcut_counts[pop_name] = len(df_sel)
        print("  objects passing mag cut:", len(df_sel))

        if len(df_sel) == 0:
            empty_map = np.zeros_like(X0, dtype=float)
            inf_map = np.full_like(X0, np.inf, dtype=float)
            clone_overlay[pop_name] = {
                "vlam": np.array([]), "vbeta": np.array([]),
                "scatter_size": info["scatter_size"],
                "scatter_alpha": info["scatter_alpha"],
            }
            clone_density_maps[pop_name] = empty_map
            clone_density_maps_downweighted[pop_name] = empty_map
            density_maps_downweighted_raw[pop_name] = empty_map.copy()
            nearest_dist_maps[pop_name] = inf_map
            continue

        if info.get("use_conditional_cloner", True):
            df_cloner_input = build_visible_subset_dataframe(
                df_sel, obstime_str=obstime_str, scorer=info["scorer"],
                max_sep_deg=max_sep_deg,
                chunk=100_000, show_progress=True,
                center_mode="custom_ecliptic",
                center_lon_deg=center_lon_deg, center_lat_deg=center_lat_deg,
            )
            print("  visible for conditional cloner:", len(df_cloner_input))

            if len(df_cloner_input) == 0:
                empty_map = np.zeros_like(X0, dtype=float)
                inf_map = np.full_like(X0, np.inf, dtype=float)
                clone_overlay[pop_name] = {
                    "vlam": np.array([]), "vbeta": np.array([]),
                    "scatter_size": info["scatter_size"],
                    "scatter_alpha": info["scatter_alpha"],
                }
                clone_density_maps[pop_name] = empty_map
                clone_density_maps_downweighted[pop_name] = empty_map
                density_maps_downweighted_raw[pop_name] = empty_map.copy()
                nearest_dist_maps[pop_name] = inf_map
                continue

            clone_visible_df = clone_population_conditional_K_from_M_with_skycut(
                df=df_cloner_input, clone_factor=f, obstime_str=obstime_str,
                scorer=info["scorer"], max_sep_deg=max_sep_deg,
                center_mode="custom_ecliptic",
                center_lon_deg=center_lon_deg, center_lat_deg=center_lat_deg,
                chunk=100_000, show_progress=True, rng=rng,
            )
            print("  final cloned visible points:", len(clone_visible_df))

            vlam_clone = clone_visible_df["vlam"].to_numpy(dtype=float)
            vbeta_clone = clone_visible_df["vbeta"].to_numpy(dtype=float)

        else:
            # uniform-cloner fallback path; expects nsc to expose
            # clone_population_uniform_mean_anomaly
            a_rep, e_rep, inc_rep, raan_rep, argp_rep, tp_new = (
                nsc.clone_population_uniform_mean_anomaly(
                    df_sel, clone_factor=f, obstime_str=obstime_str, rng=rng,
                )
            )
            H_rep = np.repeat(df_sel["H"].to_numpy(dtype=np.float32), f)
            print("  cloned orbit count:", len(a_rep))

            vlam_clone, vbeta_clone, d_au_clone, H_clone = elements_to_vlam_vbeta_d_h(
                a_rep, e_rep, inc_rep, raan_rep, argp_rep, tp_new,
                obstime_str=obstime_str, scorer=info["scorer"],
                H_vals=H_rep, max_sep_deg=max_sep_deg,
                chunk=100_000, show_progress=True,
                center_mode="custom_ecliptic",
                center_lon_deg=center_lon_deg, center_lat_deg=center_lat_deg,
            )
            print("  final cloned visible points:", len(vlam_clone))

            del a_rep, e_rep, inc_rep, raan_rep, argp_rep, tp_new
            del H_rep, d_au_clone, H_clone
            gc.collect()

        # overlay subsample
        n_overlay = min(overlay_max, len(vlam_clone))
        if n_overlay > 0:
            idx_overlay = rng.choice(len(vlam_clone), size=n_overlay, replace=False)
            clone_overlay[pop_name] = {
                "vlam": vlam_clone[idx_overlay].copy(),
                "vbeta": vbeta_clone[idx_overlay].copy(),
                "scatter_size": info["scatter_size"],
                "scatter_alpha": info["scatter_alpha"],
            }
        else:
            clone_overlay[pop_name] = {
                "vlam": np.array([]), "vbeta": np.array([]),
                "scatter_size": info["scatter_size"],
                "scatter_alpha": info["scatter_alpha"],
            }

        pts_clone = np.column_stack([vlam_clone, vbeta_clone])
        n_visible = len(pts_clone)

        if n_visible < 2:
            print(f"  too few visible points ({n_visible}) to build density map; returning zeros.")
            density_clone_map = np.zeros_like(X0, dtype=float)
            density_downweighted_map = np.zeros_like(X0, dtype=float)
            nearest_dist = np.full_like(X0, np.inf, dtype=float)
        else:
            k_eff = min(info.get("k_map", k_map), n_visible - 1)
            if k_eff < 2:
                print(f"  too few visible points ({n_visible}) for stable kNN density; returning zeros.")
                density_clone_map = np.zeros_like(X0, dtype=float)
                density_downweighted_map = np.zeros_like(X0, dtype=float)
                nearest_dist = np.full_like(X0, np.inf, dtype=float)
            else:
                print(f"  building density map with k_eff = {k_eff} (n_visible = {n_visible})")
                tree_clone = cKDTree(pts_clone)

                density_clone_flat = evaluate_density_map_full_posterior_2d(
                    tree_clone, grid_points,
                    k=k_eff, n_d0_grid=n_d0_grid_map, show_progress=True,
                )
                density_clone_map = density_clone_flat.reshape(X0.shape)
                density_downweighted_map = density_clone_map / f

                nearest_dist_flat, _ = tree_clone.query(grid_points, k=1)
                nearest_dist = nearest_dist_flat.reshape(X0.shape)

                del tree_clone, density_clone_flat, nearest_dist_flat
                gc.collect()

        clone_density_maps[pop_name] = density_clone_map
        clone_density_maps_downweighted[pop_name] = density_downweighted_map
        density_maps_downweighted_raw[pop_name] = density_downweighted_map.copy()
        nearest_dist_maps[pop_name] = nearest_dist.copy()

        print(f"{pop_name} done.")
        print("  raw cloned max density   :", np.nanmax(density_clone_map))
        print("  downweighted max density :", np.nanmax(density_downweighted_map))

        del vlam_clone, vbeta_clone, pts_clone
        gc.collect()

    density_all_clone_downweighted = sum(
        clone_density_maps_downweighted.values()
    )

    return {
        "label": center_label,
        "center_lon_deg": center_lon_deg,
        "center_lat_deg": center_lat_deg,
        "mag_min": mag_min,
        "mag_max": mag_max,
        "magcut_counts": magcut_counts,
        "overlay": clone_overlay,
        "density_maps": clone_density_maps,
        "density_maps_downweighted": clone_density_maps_downweighted,
        "density_maps_downweighted_raw": density_maps_downweighted_raw,
        "nearest_dist_maps": nearest_dist_maps,
        "density_all_downweighted": density_all_clone_downweighted,
    }


# =============================================================================
# High-level pipeline: generate maps for all mag bins at one obstime
# =============================================================================
def generate_probability_maps(
    obstime_str,
    output_path,
    clone_sources=None,
    population_settings=None,
    mag_bins=None,
    center_lon_deg=DEFAULT_CENTER_LON_DEG,
    center_lat_deg=DEFAULT_CENTER_LAT_DEG,
    center_label=None,
    max_sep_deg=DEFAULT_MAX_SEP_DEG,
    grid_lim=DEFAULT_GRID_LIM,
    grid_step=DEFAULT_GRID_STEP,
    k_map=DEFAULT_K_MAP,
    n_d0_grid_map=DEFAULT_N_D0_GRID_MAP,
    overlay_max=200_000,
    seed=42,
    save_overlays=False,
    verbose=True,
):
    """Generate density maps for every magnitude bin at one obstime, save as .npz.

    Parameters
    ----------
    obstime_str : str
        Observation time, ISO format. e.g. "2025-03-21T00:00:00".
    output_path : str
        Destination .npz file. Will be overwritten if it exists.
    clone_sources : dict or None
        Pre-loaded clone_sources dict (as returned by load_all_populations).
        If None, populations are loaded via population_settings.
    population_settings : dict or None
        Settings used by load_all_populations if clone_sources is None.
    mag_bins : list of dict or None
        Each entry: {"label": ..., "mag_min": ..., "mag_max": ...}.
        Defaults to DEFAULT_MAG_BINS.
    center_lon_deg, center_lat_deg : float
        Sky-patch center in geocentric ecliptic coordinates.
    center_label : str or None
        Human-readable label; default is f"lon{lon}_lat{lat}".
    max_sep_deg : float
        Sky-patch radius.
    grid_lim, grid_step : tuple, float
        Velocity grid limits and spacing.
    k_map, n_d0_grid_map : int
        Density estimator parameters.
    overlay_max : int
        Maximum overlay points to retain per population per bin
        (only used in memory; not saved unless save_overlays=True).
    seed : int
        RNG seed.
    save_overlays : bool
        If True, also save the cloned-object overlay subsamples
        (can produce a large file).
    verbose : bool

    Returns
    -------
    results : dict
        In-memory results, keyed by mag-bin label.
    """
    if mag_bins is None:
        mag_bins = DEFAULT_MAG_BINS

    if center_label is None:
        center_label = f"lon{center_lon_deg:g}_lat{center_lat_deg:g}"

    # populations
    if clone_sources is None:
        if verbose:
            print("[generate_probability_maps] loading populations...")
        clone_sources = load_all_populations(
            population_settings=population_settings,
            verbose=verbose,
        )

    population_names = list(clone_sources.keys())

    # grid
    x_grid, y_grid, X0, Y0, grid_points = make_default_grid(
        grid_lim=grid_lim, grid_step=grid_step,
    )

    rng = np.random.default_rng(seed)

    # iterate mag bins
    allpop_magbin_results = {}
    for mbin in mag_bins:
        if verbose:
            print("\n" + "=" * 80)
            print(f"Running all populations for mag bin: {mbin['label']}")

        allpop_magbin_results[mbin["label"]] = build_cloned_maps_for_center_magbin(
            center_lon_deg=center_lon_deg,
            center_lat_deg=center_lat_deg,
            center_label=center_label,
            clone_sources=clone_sources,
            obstime_str=obstime_str,
            max_sep_deg=max_sep_deg,
            rng=rng,
            overlay_max=overlay_max,
            grid_points=grid_points,
            X0=X0,
            k_map=k_map,
            n_d0_grid_map=n_d0_grid_map,
            mag_min=mbin["mag_min"],
            mag_max=mbin["mag_max"],
        )

    # save
    save_maps_to_npz(
        path=output_path,
        results=allpop_magbin_results,
        x_grid=x_grid, y_grid=y_grid,
        obstime_str=obstime_str,
        center_lon_deg=center_lon_deg,
        center_lat_deg=center_lat_deg,
        center_label=center_label,
        max_sep_deg=max_sep_deg,
        mag_bins=mag_bins,
        population_names=population_names,
        save_overlays=save_overlays,
    )
    if verbose:
        print(f"\n[generate_probability_maps] wrote {output_path}")

    return allpop_magbin_results


# =============================================================================
# Save / load
# =============================================================================
def save_maps_to_npz(
    path,
    results,
    x_grid,
    y_grid,
    obstime_str,
    center_lon_deg,
    center_lat_deg,
    center_label,
    max_sep_deg,
    mag_bins,
    population_names,
    save_overlays=False,
):
    """Persist a generate_probability_maps() output dict to compressed .npz.

    Stored arrays (per population pop and per mag bin label B):

        density_raw__POP__B          (ny, nx) downweighted (raw) density map
        nearest_dist__POP__B         (ny, nx) nearest-cloned-object distance
        magcut_count__POP__B         scalar int

    Plus metadata:
        x_grid, y_grid               1D coordinate axes
        obstime_str                  the obstime
        center_lon_deg, center_lat_deg, max_sep_deg
        center_label
        population_names             1D string array
        mag_bin_labels               1D string array
        mag_bin_mins, mag_bin_maxs   1D float arrays

    Optional (only if save_overlays):
        overlay_vlam__POP__B         1D arrays (variable length)
        overlay_vbeta__POP__B        1D arrays (variable length)
    """
    out = {
        "x_grid": np.asarray(x_grid, dtype=np.float64),
        "y_grid": np.asarray(y_grid, dtype=np.float64),
        "obstime_str": np.asarray(obstime_str),
        "center_lon_deg": np.asarray(center_lon_deg, dtype=np.float64),
        "center_lat_deg": np.asarray(center_lat_deg, dtype=np.float64),
        "center_label": np.asarray(center_label),
        "max_sep_deg": np.asarray(max_sep_deg, dtype=np.float64),
        "population_names": np.asarray(population_names),
        "mag_bin_labels": np.asarray([mb["label"] for mb in mag_bins]),
        "mag_bin_mins": np.asarray([mb["mag_min"] for mb in mag_bins], dtype=np.float64),
        "mag_bin_maxs": np.asarray([mb["mag_max"] for mb in mag_bins], dtype=np.float64),
    }

    for mb in mag_bins:
        label = mb["label"]
        bin_res = results[label]
        for pop in population_names:
            density = bin_res["density_maps_downweighted_raw"][pop]
            nearest = bin_res["nearest_dist_maps"][pop]
            count = bin_res["magcut_counts"][pop]
            out[f"density_raw__{pop}__{label}"] = np.asarray(density, dtype=np.float32)
            out[f"nearest_dist__{pop}__{label}"] = np.asarray(nearest, dtype=np.float32)
            out[f"magcut_count__{pop}__{label}"] = np.asarray(count, dtype=np.int64)

            if save_overlays:
                overlay = bin_res["overlay"][pop]
                out[f"overlay_vlam__{pop}__{label}"] = np.asarray(
                    overlay["vlam"], dtype=np.float32,
                )
                out[f"overlay_vbeta__{pop}__{label}"] = np.asarray(
                    overlay["vbeta"], dtype=np.float32,
                )

    np.savez_compressed(path, **out)


def load_maps_from_npz(path):
    """Load .npz back into a Python dict mirroring generate_probability_maps output.

    Convenience wrapper. For most use, prefer ProbMapSet.from_npz().
    """
    with np.load(path, allow_pickle=False) as z:
        x_grid = z["x_grid"]
        y_grid = z["y_grid"]
        obstime_str = str(z["obstime_str"])
        center_lon_deg = float(z["center_lon_deg"])
        center_lat_deg = float(z["center_lat_deg"])
        center_label = str(z["center_label"])
        max_sep_deg = float(z["max_sep_deg"])
        population_names = [str(s) for s in z["population_names"]]
        mag_bin_labels = [str(s) for s in z["mag_bin_labels"]]
        mag_bin_mins = z["mag_bin_mins"]
        mag_bin_maxs = z["mag_bin_maxs"]

        results = {}
        for label, mn, mx in zip(mag_bin_labels, mag_bin_mins, mag_bin_maxs):
            density_raw = {}
            nearest = {}
            magcut_count = {}
            for pop in population_names:
                density_raw[pop] = np.asarray(
                    z[f"density_raw__{pop}__{label}"], dtype=np.float64,
                )
                nearest[pop] = np.asarray(
                    z[f"nearest_dist__{pop}__{label}"], dtype=np.float64,
                )
                magcut_count[pop] = int(z[f"magcut_count__{pop}__{label}"])

            results[label] = {
                "label": center_label,
                "mag_min": float(mn),
                "mag_max": float(mx),
                "magcut_counts": magcut_count,
                "density_maps_downweighted_raw": density_raw,
                "nearest_dist_maps": nearest,
            }

    return {
        "x_grid": x_grid, "y_grid": y_grid,
        "obstime_str": obstime_str,
        "center_lon_deg": center_lon_deg, "center_lat_deg": center_lat_deg,
        "center_label": center_label,
        "max_sep_deg": max_sep_deg,
        "population_names": population_names,
        "mag_bins": [{"label": l, "mag_min": float(mn), "mag_max": float(mx)}
                     for l, mn, mx in zip(mag_bin_labels, mag_bin_mins, mag_bin_maxs)],
        "results": results,
    }


# =============================================================================
# ProbMapSet: the user-facing lookup object
# =============================================================================
class ProbMapSet:
    """Loaded probability maps for one obstime + sky center.

    Provides three lookup APIs:

      * score_visible(vlam, vbeta, mag_app)
            For arrays you've already propagated. Returns
            {pop_name: prob_array}.

      * score_observation(ra_deg, dec_deg, dra_deg_day, ddec_deg_day,
                          mag_app, scorer, obstime_str=None)
            For one or more raw RA/Dec/rates observations. Converts to
            (vlam, vbeta) via manual ecliptic rotation, then scores.

      * score_orbital_df(df, scorer, obstime_str=None, return_visible=True)
            For a full S3M-style DataFrame of orbital elements. Propagates,
            applies the sky cut, computes apparent mag, then scores.
            Returns (probs_dict, visible_df).

    Probabilities are computed as:
        P(pop | vlam, vbeta) = density_pop / sum_pops density_pop
    masked to zero wherever the per-population nearest-cloned-object
    distance exceeds `mask_radius_deg_per_day` (default 0.2 deg/day).
    """

    def __init__(
        self,
        x_grid,
        y_grid,
        obstime_str,
        center_lon_deg,
        center_lat_deg,
        center_label,
        max_sep_deg,
        population_names,
        mag_bins,
        results,
        mask_radius_deg_per_day=DEFAULT_MASK_RADIUS_DEG_PER_DAY,
    ):
        self.x_grid = np.asarray(x_grid, dtype=np.float64)
        self.y_grid = np.asarray(y_grid, dtype=np.float64)
        self.obstime_str = obstime_str
        self.center_lon_deg = float(center_lon_deg)
        self.center_lat_deg = float(center_lat_deg)
        self.center_label = center_label
        self.max_sep_deg = float(max_sep_deg)
        self.population_names = list(population_names)
        self.mag_bins = list(mag_bins)
        self.results = results
        self.mask_radius_deg_per_day = float(mask_radius_deg_per_day)

        # precompute per-bin probability maps {label -> {pop -> prob_map}}
        # _prob_maps_raw: unmasked (pure density ratio) — use for ROC scoring
        # _prob_maps:     masked with mask_radius      — use for visualization
        self._prob_maps_raw = {}
        self._prob_maps = {}
        for mb in self.mag_bins:
            label = mb["label"]
            density = self.results[label]["density_maps_downweighted_raw"]
            nearest = self.results[label]["nearest_dist_maps"]
            density_total = sum(density.values())

            raw_for_label = {}
            masked_for_label = {}
            for pop in self.population_names:
                frac = np.zeros_like(density[pop], dtype=np.float64)
                good = density_total > 0
                frac[good] = density[pop][good] / density_total[good]
                raw_for_label[pop] = frac.copy()
                masked = frac.copy()
                too_far = nearest[pop] > self.mask_radius_deg_per_day
                masked[too_far] = 0.0
                masked_for_label[pop] = masked
            self._prob_maps_raw[label] = raw_for_label
            self._prob_maps[label] = masked_for_label

        # mag-bin lookup table
        self._mag_mins = np.asarray([mb["mag_min"] for mb in self.mag_bins], dtype=float)
        self._mag_maxs = np.asarray([mb["mag_max"] for mb in self.mag_bins], dtype=float)
        self._labels = [mb["label"] for mb in self.mag_bins]

    # ------------------------------------------------------------------
    # constructors
    # ------------------------------------------------------------------
    @classmethod
    def from_npz(cls, path, mask_radius_deg_per_day=DEFAULT_MASK_RADIUS_DEG_PER_DAY):
        d = load_maps_from_npz(path)
        return cls(
            x_grid=d["x_grid"], y_grid=d["y_grid"],
            obstime_str=d["obstime_str"],
            center_lon_deg=d["center_lon_deg"],
            center_lat_deg=d["center_lat_deg"],
            center_label=d["center_label"],
            max_sep_deg=d["max_sep_deg"],
            population_names=d["population_names"],
            mag_bins=d["mag_bins"],
            results=d["results"],
            mask_radius_deg_per_day=mask_radius_deg_per_day,
        )

    # ------------------------------------------------------------------
    # mag-bin selection
    # ------------------------------------------------------------------
    def _pick_bin_label(self, mag_app):
        """For each scalar mag_app return the matching bin label, or None.

        Vectorized: returns a list of length len(mag_app)."""
        mag_app = np.atleast_1d(np.asarray(mag_app, dtype=float))
        out = [None] * len(mag_app)
        for i, m in enumerate(mag_app):
            if not np.isfinite(m):
                continue
            sel = (m >= self._mag_mins) & (m < self._mag_maxs)
            if np.any(sel):
                out[i] = self._labels[np.argmax(sel)]
        return out

    def get_probability_map(self, mag_bin_label, pop_name):
        """Direct accessor for the masked probability map at one bin/pop."""
        return self._prob_maps[mag_bin_label][pop_name]

    def get_xy_meshgrid(self):
        """Return (X0, Y0) meshgrids matching the stored prob maps."""
        X0, Y0 = np.meshgrid(self.x_grid, self.y_grid, indexing="xy")
        return X0, Y0

    # ------------------------------------------------------------------
    # core lookup
    # ------------------------------------------------------------------
    def _lookup_in_map(self, prob_map, vlam, vbeta):
        """Bilinear interpolation of prob_map at (vlam, vbeta).

        Out-of-grid queries return 0.
        """
        vlam = np.atleast_1d(np.asarray(vlam, dtype=float))
        vbeta = np.atleast_1d(np.asarray(vbeta, dtype=float))

        # grid spacing (uniform)
        dx = self.x_grid[1] - self.x_grid[0]
        dy = self.y_grid[1] - self.y_grid[0]
        x0 = self.x_grid[0]
        y0 = self.y_grid[0]

        # fractional indices; replace non-finite inputs with a sentinel
        # so the float->int cast doesn't warn (these get masked out below)
        finite = np.isfinite(vlam) & np.isfinite(vbeta)
        fx = np.where(finite, (vlam - x0) / dx, -1e9)
        fy = np.where(finite, (vbeta - y0) / dy, -1e9)

        # int part + fractional part
        ix = np.floor(fx).astype(int)
        iy = np.floor(fy).astype(int)
        rx = fx - ix
        ry = fy - iy

        # in-bounds mask: need ix and ix+1 both in [0, nx-1]
        nx = len(self.x_grid)
        ny = len(self.y_grid)
        inb = (ix >= 0) & (ix + 1 < nx) & (iy >= 0) & (iy + 1 < ny) & finite

        out = np.zeros_like(vlam, dtype=float)
        if not np.any(inb):
            return out

        # bilinear interp at in-bounds points
        ix_c = np.clip(ix, 0, nx - 2)
        iy_c = np.clip(iy, 0, ny - 2)

        v00 = prob_map[iy_c, ix_c]
        v10 = prob_map[iy_c, ix_c + 1]
        v01 = prob_map[iy_c + 1, ix_c]
        v11 = prob_map[iy_c + 1, ix_c + 1]

        interp = (
            v00 * (1 - rx) * (1 - ry)
            + v10 * rx * (1 - ry)
            + v01 * (1 - rx) * ry
            + v11 * rx * ry
        )
        out[inb] = interp[inb]
        return out

    def score_visible(
        self,
        vlam,
        vbeta,
        mag_app,
        return_bin_labels=False,
        apply_mask=True,
    ):
        """Score already-propagated visible objects against the maps.

        Parameters
        ----------
        vlam, vbeta : array-like, deg/day
        mag_app : array-like, apparent magnitude
        return_bin_labels : bool
            If True, also return the assigned mag-bin label per object.
        apply_mask : bool, default True
            If True, uses the masked probability maps (0.2 deg/day nearest-
            object mask applied — correct for visualization and single-object
            scoring where you want P=0 in empty velocity-space regions).
            If False, uses the raw unmasked density ratios — use this for ROC
            curve analysis so that sparse but non-zero NEO regions still
            contribute positive P(NEO) values rather than being zeroed out.

        Returns
        -------
        probs : dict
            {pop_name: ndarray of probabilities, length N}.
            Objects whose apparent mag falls outside all bins get 0.
        bin_labels : list (optional)
            One bin label per object (None if outside all bins).
        """
        vlam = np.atleast_1d(np.asarray(vlam, dtype=float))
        vbeta = np.atleast_1d(np.asarray(vbeta, dtype=float))
        mag_app = np.atleast_1d(np.asarray(mag_app, dtype=float))
        if not (len(vlam) == len(vbeta) == len(mag_app)):
            raise ValueError("vlam, vbeta, mag_app must be the same length.")

        prob_source = self._prob_maps if apply_mask else self._prob_maps_raw

        probs = {pop: np.zeros(len(vlam), dtype=float) for pop in self.population_names}
        bin_labels = self._pick_bin_label(mag_app)

        # group by bin label so each prob_map lookup is one vectorized call
        for label in self._labels:
            idx = np.array([i for i, lbl in enumerate(bin_labels) if lbl == label],
                           dtype=int)
            if len(idx) == 0:
                continue
            for pop in self.population_names:
                pmap = prob_source[label][pop]
                probs[pop][idx] = self._lookup_in_map(pmap, vlam[idx], vbeta[idx])

        if return_bin_labels:
            return probs, bin_labels
        return probs

    # ------------------------------------------------------------------
    # convenience: score from raw RA/Dec/rates
    # ------------------------------------------------------------------
    def score_observation(
        self,
        ra_deg,
        dec_deg,
        dra_deg_day,
        ddec_deg_day,
        mag_app,
        obstime_str=None,
        return_intermediates=False,
    ):
        """Convert RA/Dec/rates to (vlam, vbeta), then look up probability.

        Uses the manual equatorial->ecliptic rotation — no scorer or
        external package needed. Just the measured sky position, rates,
        and apparent magnitude.

        Parameters
        ----------
        ra_deg, dec_deg : float or array-like, degrees
        dra_deg_day, ddec_deg_day : float or array-like, degrees/day
        mag_app : float or array-like, apparent V magnitude
        obstime_str : str or None
            Stored in return_intermediates for bookkeeping only.
            Defaults to self.obstime_str.
        return_intermediates : bool
            If True, returns dict with probs, vlam, vbeta, bin_labels.

        Returns
        -------
        probs : dict  {pop_name: ndarray shape (N,)} with P in [0,1].
        """
        if obstime_str is None:
            obstime_str = self.obstime_str

        ra_deg = np.atleast_1d(np.asarray(ra_deg, dtype=float))
        dec_deg = np.atleast_1d(np.asarray(dec_deg, dtype=float))
        dra_deg_day = np.atleast_1d(np.asarray(dra_deg_day, dtype=float))
        ddec_deg_day = np.atleast_1d(np.asarray(ddec_deg_day, dtype=float))
        mag_app = np.atleast_1d(np.asarray(mag_app, dtype=float))

        vlam, vbeta = radec_rates_to_ecliptic_rates_manual(
            ra_deg, dec_deg, dra_deg_day, ddec_deg_day,
        )

        probs, bin_labels = self.score_visible(
            vlam=vlam, vbeta=vbeta, mag_app=mag_app, return_bin_labels=True,
        )

        if return_intermediates:
            return {
                "probs": probs,
                "vlam": vlam, "vbeta": vbeta,
                "bin_labels": bin_labels,
                "obstime_str": obstime_str,
            }
        return probs

    # ------------------------------------------------------------------
    # convenience: score a full orbital-element DataFrame
    # ------------------------------------------------------------------
    def score_orbital_df(
        self,
        df,
        scorer,
        obstime_str=None,
        max_sep_deg=None,
        chunk=100_000,
        show_progress=True,
        return_visible=True,
        apply_mask=False,
    ):
        """Propagate orbital elements at obstime, sky-cut, score visible objects.

        Single propagation pass: computes apparent magnitude, sky position,
        vlam, and vbeta all in one call to avoid propagating 13M+ objects twice.

        Parameters
        ----------
        df : pandas.DataFrame
            Must have columns: a, e, i, node, argperi, t_p, H.
        scorer : nsc.NEOMODScorer
            Used for observer state (Earth position/velocity).
        obstime_str : str or None
            Defaults to self.obstime_str.
        max_sep_deg : float or None
            Defaults to self.max_sep_deg.
        chunk : int
            Propagation chunk size. Increase (e.g. 500_000) to speed up
            large populations like MBA.
        return_visible : bool
            If True, returns (probs_dict, visible_df).
        apply_mask : bool, default False
            Whether to apply the nearest-object mask. False (default) gives
            raw density-ratio probabilities suitable for ROC analysis.
            Set True for visualization-consistent scoring.
        """
        if obstime_str is None:
            obstime_str = self.obstime_str
        if max_sep_deg is None:
            max_sep_deg = self.max_sep_deg

        t_obs = Time(obstime_str, scale="tdb")

        a_AU   = df["a"].to_numpy(dtype=float)
        e      = df["e"].to_numpy(dtype=float)
        inc    = df["i"].to_numpy(dtype=float)
        raan   = df["node"].to_numpy(dtype=float)
        argp   = df["argperi"].to_numpy(dtype=float)
        tp_mjd = df["t_p"].to_numpy(dtype=float)
        H_full = df["H"].to_numpy(dtype=float)

        # ── single propagation pass ────────────────────────────────────────
        r_ecl, v_ecl = nsc.elements_to_helio_ecliptic_state(
            a_AU=a_AU, e=e, inc_deg=inc,
            raan_deg=raan, argp_deg=argp, tp_mjd=tp_mjd,
            obstime_str=obstime_str,
            method="newton", n_iter=10,
            chunk=chunk, show_progress=show_progress,
        )

        eps = np.deg2rad(ECLIPTIC_OBLIQUITY_DEG)
        c, s = np.cos(eps), np.sin(eps)

        # ecliptic -> equatorial position + velocity
        r_eq = np.column_stack([
            r_ecl[:, 0],
            c * r_ecl[:, 1] - s * r_ecl[:, 2],
            s * r_ecl[:, 1] + c * r_ecl[:, 2],
        ])
        v_eq = np.column_stack([
            v_ecl[:, 0],
            c * v_ecl[:, 1] - s * v_ecl[:, 2],
            s * v_ecl[:, 1] + c * v_ecl[:, 2],
        ])

        _, rE, vE, r_tele = scorer._get_earth_and_observer(obstime_str)

        r_rel = r_eq - (rE + r_tele)   # topocentric equatorial position
        v_rel = v_eq - vE               # topocentric equatorial velocity

        # ── apparent magnitude (full population, before sky cut) ───────────
        r_helio = r_eq
        r_sun_vec = np.linalg.norm(r_helio, axis=1)
        r_obs_vec = np.linalg.norm(r_rel, axis=1)

        r_sun_au = r_sun_vec / AU_KM if np.nanmedian(r_sun_vec) > 1000 else r_sun_vec
        Delta_au = r_obs_vec / AU_KM if np.nanmedian(r_obs_vec) > 1000 else r_obs_vec

        u_sun = -r_helio / np.linalg.norm(r_helio, axis=1, keepdims=True)
        u_obs = -r_rel   / np.linalg.norm(r_rel,   axis=1, keepdims=True)

        cos_alpha = np.clip(np.sum(u_sun * u_obs, axis=1), -1.0, 1.0)
        alpha = np.arccos(cos_alpha)

        G = 0.15
        A1, A2, B1, B2 = 3.33, 1.87, 0.63, 1.22
        phi1 = np.exp(-A1 * np.tan(0.5 * alpha) ** B1)
        phi2 = np.exp(-A2 * np.tan(0.5 * alpha) ** B2)
        Phi  = np.maximum((1 - G) * phi1 + G * phi2, 1e-30)

        mag_app_full = H_full + 5.0 * np.log10(r_sun_au * Delta_au) - 2.5 * np.log10(Phi)

        # ── sky cut ────────────────────────────────────────────────────────
        x, y, z   = r_rel[:, 0], r_rel[:, 1], r_rel[:, 2]
        vx, vy, vz = v_rel[:, 0], v_rel[:, 1], v_rel[:, 2]

        rho2 = x*x + y*y
        r2   = rho2 + z*z
        rho  = np.sqrt(rho2)

        ra  = np.arctan2(y, x)
        dec = np.arctan2(z, rho)
        dra  = (x*vy - y*vx) / rho2
        ddec = (vz*rho2 - z*(x*vx + y*vy)) / (r2 * rho)

        good = (
            np.isfinite(ra) & np.isfinite(dec) &
            np.isfinite(dra) & np.isfinite(ddec) &
            np.isfinite(rho2) & (rho2 > 0) &
            np.isfinite(r2) & (r2 > 0) &
            np.isfinite(rho) & (rho > 0)
        )

        ra_deg_g   = np.rad2deg(ra[good]) % 360.0
        dec_deg_g  = np.rad2deg(dec[good])

        sc = SkyCoord(ra=ra_deg_g * u.deg, dec=dec_deg_g * u.deg,
                      frame=GCRS(obstime=t_obs))

        center_ecl = SkyCoord(
            lon=self.center_lon_deg * u.deg,
            lat=self.center_lat_deg * u.deg,
            distance=1.0 * u.AU,
            frame=GeocentricTrueEcliptic(obstime=t_obs),
        )
        center_sc = center_ecl.transform_to(GCRS(obstime=t_obs))
        sky_sep = sc.separation(center_sc).deg
        m = sky_sep <= max_sep_deg

        # original indices of survivors
        good_idx = np.flatnonzero(good)
        sel_idx  = good_idx[m]

        # ── ecliptic rates for survivors (manual rotation) ─────────────────
        r_rel_sel = r_rel[good][m]
        v_rel_sel = v_rel[good][m]

        x_ecl  = r_rel_sel[:, 0]
        y_ecl  =  c * r_rel_sel[:, 1] + s * r_rel_sel[:, 2]
        z_ecl  = -s * r_rel_sel[:, 1] + c * r_rel_sel[:, 2]
        vx_ecl =  v_rel_sel[:, 0]
        vy_ecl =  c * v_rel_sel[:, 1] + s * v_rel_sel[:, 2]
        vz_ecl = -s * v_rel_sel[:, 1] + c * v_rel_sel[:, 2]

        rho2_ecl = x_ecl*x_ecl + y_ecl*y_ecl
        r2_ecl   = rho2_ecl + z_ecl*z_ecl
        rho_ecl  = np.sqrt(rho2_ecl)

        dlam  = (x_ecl*vy_ecl - y_ecl*vx_ecl) / rho2_ecl
        dbeta = (vz_ecl*rho2_ecl - z_ecl*(x_ecl*vx_ecl + y_ecl*vy_ecl)) / (r2_ecl * rho_ecl)

        vlam_vis  = np.rad2deg(dlam)  * 86400.0
        vbeta_vis = np.rad2deg(dbeta) * 86400.0

        mag_app_vis = mag_app_full[sel_idx]

        # ── build visible_df ───────────────────────────────────────────────
        visible_df = df.iloc[sel_idx].copy().reset_index(drop=True)
        visible_df["ra_deg"]      = ra_deg_g[m]
        visible_df["dec_deg"]     = dec_deg_g[m]
        visible_df["sky_sep_deg"] = sky_sep[m]
        visible_df["vlam"]        = vlam_vis
        visible_df["vbeta"]       = vbeta_vis
        visible_df["mag_app"]     = mag_app_vis

        # ── score ──────────────────────────────────────────────────────────
        probs, bin_labels = self.score_visible(
            vlam=vlam_vis, vbeta=vbeta_vis, mag_app=mag_app_vis,
            return_bin_labels=True,
            apply_mask=apply_mask,
        )
        visible_df["mag_bin_label"] = bin_labels
        for pop, p in probs.items():
            visible_df[f"P_{pop}"] = p

        if return_visible:
            return probs, visible_df
        return probs


# =============================================================================
# Manual RA/Dec rates -> ecliptic rates (avoids astropy differential bug)
# =============================================================================
def radec_rates_to_ecliptic_rates_manual(ra_deg, dec_deg, dra_deg_day, ddec_deg_day):
    """Compute (vlam, vbeta) from (ra, dec, dra/dt, ddec/dt) by direct
    rotation of the topocentric position+velocity unit vectors via the
    obliquity matrix.

    Inputs are arrays in degrees and deg/day. Returns (vlam, vbeta) in
    deg/day. The unit-distance assumption falls out: only the angular
    rates appear in the result.
    """
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    dra = np.deg2rad(dra_deg_day)
    ddec = np.deg2rad(ddec_deg_day)

    # unit topocentric position vector in equatorial frame
    cosd = np.cos(dec); sind = np.sin(dec)
    cosa = np.cos(ra); sina = np.sin(ra)

    x = cosd * cosa
    y = cosd * sina
    z = sind

    # velocity = d/dt of unit vector (with d/dt = 1/day units)
    vx = -sind * cosa * ddec - cosd * sina * dra
    vy = -sind * sina * ddec + cosd * cosa * dra
    vz = cosd * ddec

    # rotate equatorial -> ecliptic
    eps = np.deg2rad(ECLIPTIC_OBLIQUITY_DEG)
    c, s = np.cos(eps), np.sin(eps)

    x_ecl = x
    y_ecl = c * y + s * z
    z_ecl = -s * y + c * z

    vx_ecl = vx
    vy_ecl = c * vy + s * vz
    vz_ecl = -s * vy + c * vz

    rho2_ecl = x_ecl*x_ecl + y_ecl*y_ecl
    r2_ecl = rho2_ecl + z_ecl*z_ecl
    rho_ecl = np.sqrt(rho2_ecl)

    dlam = (x_ecl*vy_ecl - y_ecl*vx_ecl) / rho2_ecl
    dbeta = (vz_ecl*rho2_ecl - z_ecl*(x_ecl*vx_ecl + y_ecl*vy_ecl)) / (r2_ecl*rho_ecl)

    return np.rad2deg(dlam), np.rad2deg(dbeta)


# =============================================================================
# Plotting (the 0.2 deg/day mask plot from the notebook)
# =============================================================================
def plot_probability_maps(
    prob_map_set,
    mag_bin_label,
    mask_radius=DEFAULT_MASK_RADIUS_DEG_PER_DAY,
    fig=None,
    axes=None,
    cmap=None,
    show=True,
):
    """Reproduce the notebook's per-population masked probability plot.

    For each population, divides the raw downweighted density by the
    raw downweighted density summed over all populations, then sets
    pixels with `nearest_dist > mask_radius` to zero. Displayed on a
    [0, 1] color scale.

    Parameters
    ----------
    prob_map_set : ProbMapSet
    mag_bin_label : str
        Must match one of the labels in prob_map_set.mag_bins.
    mask_radius : float
        deg/day cutoff. Defaults to 0.2.
    fig, axes : matplotlib figure and ndarray of 4 axes, optional.
        If None, a 2x2 figure is created.
    cmap : str or Colormap, optional
    show : bool
        If True, calls plt.show() at the end.

    Returns
    -------
    fig, axes
    """
    import matplotlib.pyplot as plt

    if mag_bin_label not in prob_map_set.results:
        raise KeyError(
            f"{mag_bin_label!r} not in available labels "
            f"{list(prob_map_set.results.keys())}"
        )

    bin_res = prob_map_set.results[mag_bin_label]
    maps_raw = bin_res["density_maps_downweighted_raw"]
    nearest_dist_maps = bin_res["nearest_dist_maps"]
    pop_names = list(prob_map_set.population_names)

    maps_all_raw = sum(maps_raw[p] for p in pop_names)

    X0, Y0 = prob_map_set.get_xy_meshgrid()
    xlim = (prob_map_set.x_grid[0], prob_map_set.x_grid[-1])
    ylim = (prob_map_set.y_grid[0], prob_map_set.y_grid[-1])

    if fig is None or axes is None:
        n = len(pop_names)
        ncols = 2 if n > 1 else 1
        nrows = int(np.ceil(n / ncols))
        fig, axes_arr = plt.subplots(
            nrows, ncols,
            figsize=(7 * ncols, 6 * nrows),
            constrained_layout=True,
        )
        axes = np.atleast_1d(axes_arr).ravel()
    else:
        axes = np.atleast_1d(axes).ravel()

    for ax, pop_name in zip(axes, pop_names):
        # exact computation from the notebook's mask_radius=0.2 cell
        frac_map = np.zeros_like(maps_raw[pop_name], dtype=float)
        good = maps_all_raw > 0
        frac_map[good] = maps_raw[pop_name][good] / maps_all_raw[good]

        too_far = nearest_dist_maps[pop_name] > mask_radius
        frac_map[too_far] = 0.0

        im = ax.pcolormesh(
            X0, Y0, frac_map, shading="auto", vmin=0, vmax=1, cmap=cmap,
        )

        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel(r"$v_\lambda$ (deg/day)")
        ax.set_ylabel(r"$v_\beta$ (deg/day)")
        ax.set_title(
            f"{pop_name}: cloned/downweighted probability\n"
            f"{mag_bin_label}, center=({prob_map_set.center_lon_deg:g},"
            f"{prob_map_set.center_lat_deg:g}), mask>{mask_radius} deg/day"
        )
        fig.colorbar(im, ax=ax, label="probability")

    # hide unused axes if pop count doesn't fill the grid
    for ax in axes[len(pop_names):]:
        ax.set_visible(False)

    if show:
        plt.show()

    return fig, axes


# =============================================================================
# Public API
# =============================================================================
__all__ = [
    # constants
    "AU_KM", "MU_SUN", "ECLIPTIC_OBLIQUITY_DEG",
    "DEFAULT_GRID_LIM", "DEFAULT_GRID_STEP",
    "DEFAULT_K_MAP", "DEFAULT_N_D0_GRID_MAP",
    "DEFAULT_MAX_SEP_DEG",
    "DEFAULT_CENTER_LON_DEG", "DEFAULT_CENTER_LAT_DEG",
    "DEFAULT_MAG_BINS", "DEFAULT_POPULATION_SETTINGS",
    "DEFAULT_MASK_RADIUS_DEG_PER_DAY",
    # grid
    "make_default_grid",
    # S3M loading
    "load_s3m_population", "load_all_populations",
    # geometry
    "wrap360_deg",
    "mean_anomaly_deg_at_obstime", "mean_longitude_deg_at_obstime",
    # mag
    "compute_apparent_magnitude_for_population", "select_df_by_mag_bin",
    # sky cut + visible subset
    "build_visible_subset_dataframe",
    # cloning
    "clone_population_conditional_K_from_M",
    "clone_population_conditional_K_from_M_with_skycut",
    # legacy
    "elements_to_vlam_vbeta_d_h",
    # density
    "log_posterior_d0_2d", "normalize_posterior_from_loggrid",
    "get_knn_distances",
    "estimate_density_full_posterior_2d",
    "evaluate_density_map_full_posterior_2d",
    # main pipeline
    "build_cloned_maps_for_center_magbin",
    "generate_probability_maps",
    # save/load
    "save_maps_to_npz", "load_maps_from_npz",
    # lookup
    "ProbMapSet", "radec_rates_to_ecliptic_rates_manual",
    # plotting
    "plot_probability_maps",
]
