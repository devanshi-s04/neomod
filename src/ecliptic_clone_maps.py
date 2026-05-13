from __future__ import annotations

import gc
from math import lgamma
from typing import Dict, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from tqdm import tqdm

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, GCRS, GeocentricTrueEcliptic, get_sun

import neoscore as nsc
import s3m_loader as load_s3m


AU_KM = 149_597_870.7
MU_SUN = 1.32712440018e11  # km^3/s^2
EPSILON_DEG = 23.439291


# -----------------------------------------------------------------------------
# Loading populations and scorers
# -----------------------------------------------------------------------------

def load_population_and_scorer(pop: str, verbose: bool = False):
    """
    Load one S3M population and build its NEOMODScorer from its own 4D array.

    Parameters
    ----------
    pop : str
        One of: "mba", "neo", "tno", "trojan".
    verbose : bool
        Passed to load_s3m.define_s3m.

    Returns
    -------
    df : pandas.DataFrame
        Population dataframe.
    scorer : nsc.NEOMODScorer
        Scorer built from that population's 4D histogram.
    """
    df = load_s3m.define_s3m(pop=pop, verbose=verbose)
    df, array4d, Hc, ac, ec, ic = load_s3m.s3m_array(df)
    scorer = nsc.NEOMODScorer(array4d, Hc, ac, ec, ic)
    return df, scorer


def load_default_populations(verbose: bool = False) -> Dict[str, dict]:
    """
    Load MBA, NEO, TNO, and Trojan populations with their scorers.

    Returns
    -------
    populations : dict
        Dictionary keyed by display name: "MBA", "NEO", "TNO", "Trojans".
        Each value contains {"df": ..., "scorer": ...}.
    """
    mba_df, mba_scorer = load_population_and_scorer("mba", verbose=verbose)
    neo_df, neo_scorer = load_population_and_scorer("neo", verbose=verbose)
    tno_df, tno_scorer = load_population_and_scorer("tno", verbose=verbose)
    trojan_df, trojan_scorer = load_population_and_scorer("trojan", verbose=verbose)

    return {
        "MBA": {"df": mba_df, "scorer": mba_scorer},
        "NEO": {"df": neo_df, "scorer": neo_scorer},
        "TNO": {"df": tno_df, "scorer": tno_scorer},
        "Trojans": {"df": trojan_df, "scorer": trojan_scorer},
    }


# -----------------------------------------------------------------------------
# Grid helpers
# -----------------------------------------------------------------------------

def make_rate_grid(
    xlim: Tuple[float, float] = (-0.8, 0.8),
    ylim: Tuple[float, float] = (-0.8, 0.8),
    grid_step: float = 0.01,
):
    """
    Build the regular (v_lambda, v_beta) evaluation grid.

    Returns
    -------
    grid : dict
        Contains xlim, ylim, grid_step, x_grid, y_grid, X0, Y0, grid_points.
    """
    x_grid = np.arange(xlim[0], xlim[1] + grid_step, grid_step)
    y_grid = np.arange(ylim[0], ylim[1] + grid_step, grid_step)
    X0, Y0 = np.meshgrid(x_grid, y_grid, indexing="xy")
    grid_points = np.column_stack([X0.ravel(), Y0.ravel()])

    return {
        "xlim": xlim,
        "ylim": ylim,
        "grid_step": grid_step,
        "x_grid": x_grid,
        "y_grid": y_grid,
        "X0": X0,
        "Y0": Y0,
        "grid_points": grid_points,
    }


# -----------------------------------------------------------------------------
# Bayesian nearest-neighbor density estimator
# -----------------------------------------------------------------------------

def log_posterior_d0_2d(d0_grid, d_knn):
    """
    Unnormalized log posterior for the 2D Appendix-B nearest-neighbor model.
    """
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
    """
    Normalize a posterior sampled on a d0 grid.
    """
    d0_grid = np.asarray(d0_grid, dtype=float)
    logp = np.asarray(logp, dtype=float)

    logp_shift = logp - np.max(logp)
    p_unnorm = np.exp(logp_shift)
    norm = np.trapz(p_unnorm, d0_grid)

    if norm <= 0 or not np.isfinite(norm):
        raise ValueError("Posterior normalization failed.")

    return p_unnorm / norm


def get_knn_distances(tree, x0, y0, k=10):
    """
    Return distances to the k nearest neighbors of one evaluation point.
    """
    q = np.array([[x0, y0]])
    d, idx = tree.query(q, k=k)
    return d[0], idx[0]


def estimate_density_full_posterior_2d(
    tree,
    x0,
    y0,
    k=10,
    d0_min_factor=5.0,
    d0_max_factor=10.0,
    n_d0_grid=800,
    return_diagnostics=False,
):
    """
    Estimate the posterior-mean local 2D density at one point.
    """
    d_knn, _ = get_knn_distances(tree, x0, y0, k=k)

    d0_min = d_knn[0] / d0_min_factor
    d0_max = d_knn[-1] * d0_max_factor
    d0_grid = np.logspace(np.log10(d0_min), np.log10(d0_max), n_d0_grid)

    logp_d0 = log_posterior_d0_2d(d0_grid, d_knn)
    p_d0 = normalize_posterior_from_loggrid(d0_grid, logp_d0)

    n0_grid = 1.0 / (np.pi * d0_grid**2)
    d0_map = d0_grid[np.argmax(p_d0)]
    d0_mean = np.trapz(d0_grid * p_d0, d0_grid)
    n0_map = 1.0 / (np.pi * d0_map**2)
    n0_mean = np.trapz(n0_grid * p_d0, d0_grid)

    if return_diagnostics:
        return n0_mean, {
            "d_knn": d_knn,
            "d0_grid": d0_grid,
            "p_d0": p_d0,
            "d0_map": d0_map,
            "d0_mean": d0_mean,
            "n0_map": n0_map,
            "n0_mean": n0_mean,
        }

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
    """
    Evaluate the posterior-mean density on every point of a grid.
    """
    density_vals = np.empty(len(grid_points), dtype=float)

    iterator = enumerate(grid_points)
    if show_progress:
        iterator = tqdm(iterator, total=len(grid_points), desc="Evaluating density map")

    for i, (x0, y0) in iterator:
        density_vals[i] = estimate_density_full_posterior_2d(
            tree,
            x0,
            y0,
            k=k,
            d0_min_factor=d0_min_factor,
            d0_max_factor=d0_max_factor,
            n_d0_grid=n_d0_grid,
            return_diagnostics=False,
        )

    return density_vals


# -----------------------------------------------------------------------------
# Cloning and rate conversion
# -----------------------------------------------------------------------------

def clone_population_random_mean_anomaly(df, clone_factor, obstime_str, rng=None):
    """
    Clone each orbit by randomizing mean anomaly at the observation epoch.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    t_obs = Time(obstime_str, scale="tdb")

    a_AU = df["a"].to_numpy(dtype=np.float32)
    e = df["e"].to_numpy(dtype=np.float32)
    inc = df["i"].to_numpy(dtype=np.float32)
    raan = df["node"].to_numpy(dtype=np.float32)
    argp = df["argperi"].to_numpy(dtype=np.float32)

    a_rep = np.repeat(a_AU, clone_factor).astype(np.float32)
    e_rep = np.repeat(e, clone_factor).astype(np.float32)
    inc_rep = np.repeat(inc, clone_factor).astype(np.float32)
    raan_rep = np.repeat(raan, clone_factor).astype(np.float32)
    argp_rep = np.repeat(argp, clone_factor).astype(np.float32)

    a_km = a_rep.astype(np.float64) * AU_KM
    n_rad_s = np.sqrt(MU_SUN / (a_km**3))
    n_rad_day = n_rad_s * 86400.0

    M_rand = rng.uniform(0.0, 2.0 * np.pi, size=len(a_rep)).astype(np.float64)
    tp_new = (t_obs.mjd - (M_rand / n_rad_day)).astype(np.float32)

    return a_rep, e_rep, inc_rep, raan_rep, argp_rep, tp_new


def elements_to_vlam_vbeta(
    a_AU,
    e,
    inc,
    raan,
    argp,
    tp_mjd,
    obstime_str,
    scorer,
    max_sep_deg=30.0,
    chunk=100_000,
    show_progress=True,
    center_mode="custom_ecliptic",
    center_lon_deg=180.0,
    center_lat_deg=0.0,
):
    """
    Convert orbital elements to observed ecliptic rates in a selected sky patch.
    """
    t_obs = Time(obstime_str, scale="tdb")

    r_obj_ecl, v_obj_ecl = nsc.elements_to_helio_ecliptic_state(
        a_AU=a_AU,
        e=e,
        inc_deg=inc,
        raan_deg=raan,
        argp_deg=argp,
        tp_mjd=tp_mjd,
        obstime_str=obstime_str,
        method="newton",
        n_iter=10,
        chunk=chunk,
        show_progress=show_progress,
    )

    eps = np.deg2rad(EPSILON_DEG)
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

    rho2 = x * x + y * y
    r2 = rho2 + z * z
    rho = np.sqrt(rho2)

    ra = np.arctan2(y, x)
    dec = np.arctan2(z, rho)

    dra = (x * vy - y * vx) / rho2
    ddec = (vz * rho2 - z * (x * vx + y * vy)) / (r2 * rho)

    good = (
        np.isfinite(ra)
        & np.isfinite(dec)
        & np.isfinite(dra)
        & np.isfinite(ddec)
        & (rho2 > 0)
        & np.isfinite(r2)
        & (r2 > 0)
    )

    ra_deg = np.rad2deg(ra[good]) % 360.0
    dec_deg = np.rad2deg(dec[good])
    dra_deg_day = np.rad2deg(dra[good]) * 86400.0
    ddec_deg_day = np.rad2deg(ddec[good]) * 86400.0

    sc = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame=GCRS(obstime=t_obs))

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
            lon=center_lon_deg * u.deg,
            lat=center_lat_deg * u.deg,
            distance=1.0 * u.AU,
            frame=GeocentricTrueEcliptic(obstime=t_obs),
        )
        center_sc = center_ecl.transform_to(GCRS(obstime=t_obs))
    else:
        raise ValueError("center_mode must be 'antisun' or 'custom_ecliptic'")

    sep = sc.separation(center_sc).deg
    m = sep <= max_sep_deg

    ra_deg = ra_deg[m]
    dec_deg = dec_deg[m]
    dra_deg_day = dra_deg_day[m]
    ddec_deg_day = ddec_deg_day[m]

    vlam, vbeta = nsc.radec_rates_to_ecliptic_rates_at_obstime(
        ra_deg,
        dec_deg,
        dra_deg_day,
        ddec_deg_day,
        obstime_str,
    )

    finite = np.isfinite(vlam) & np.isfinite(vbeta)
    return vlam[finite].astype(np.float32), vbeta[finite].astype(np.float32)


# -----------------------------------------------------------------------------
# High-level map builder
# -----------------------------------------------------------------------------

def make_clone_sources(populations: Dict[str, dict], clone_factors: Optional[Dict[str, int]] = None):
    """
    Attach plotting defaults and clone factors to loaded populations.
    """
    if clone_factors is None:
        clone_factors = {
            "MBA": 1,
            "NEO": 10,
            "TNO": 10,
            "Trojans": 10,
        }

    style_defaults = {
        "MBA": {"scatter_size": 1, "scatter_alpha": 0.05},
        "NEO": {"scatter_size": 4, "scatter_alpha": 0.08},
        "TNO": {"scatter_size": 6, "scatter_alpha": 0.10},
        "Trojans": {"scatter_size": 6, "scatter_alpha": 0.10},
    }

    clone_sources = {}
    for pop_name, info in populations.items():
        clone_sources[pop_name] = {
            "df": info["df"],
            "scorer": info["scorer"],
            "clone_factor": clone_factors[pop_name],
            "scatter_size": style_defaults[pop_name]["scatter_size"],
            "scatter_alpha": style_defaults[pop_name]["scatter_alpha"],
        }

    return clone_sources


def build_cloned_maps_for_center(
    center_lon_deg,
    center_lat_deg,
    center_label,
    clone_sources,
    obstime_str,
    max_sep_deg=30.0,
    rng=None,
    overlay_max=100_000,
    xlim=(-0.8, 0.8),
    ylim=(-0.8, 0.8),
    grid_step=0.01,
    k_map=10,
    n_d0_grid_map=400,
    chunk=100_000,
    show_progress=True,
):
    """
    Build cloned/downweighted density maps for one ecliptic sky center.

    Returns
    -------
    result : dict
        Contains overlays, individual density maps, combined density map,
        and the grid used to evaluate them.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    grid = make_rate_grid(xlim=xlim, ylim=ylim, grid_step=grid_step)
    X0 = grid["X0"]
    Y0 = grid["Y0"]
    grid_points = grid["grid_points"]

    clone_overlay = {}
    clone_density_maps = {}
    clone_density_maps_downweighted = {}

    for pop_name, info in clone_sources.items():
        f = info["clone_factor"]
        df = info["df"]

        print(f"\n[{center_label}] Cloning {pop_name} with factor {f}...")

        a_rep, e_rep, inc_rep, raan_rep, argp_rep, tp_new = clone_population_random_mean_anomaly(
            df=df,
            clone_factor=f,
            obstime_str=obstime_str,
            rng=rng,
        )
        print("  cloned orbit count:", len(a_rep))

        vlam_clone, vbeta_clone = elements_to_vlam_vbeta(
            a_rep,
            e_rep,
            inc_rep,
            raan_rep,
            argp_rep,
            tp_new,
            obstime_str=obstime_str,
            scorer=info["scorer"],
            max_sep_deg=max_sep_deg,
            chunk=chunk,
            show_progress=show_progress,
            center_mode="custom_ecliptic",
            center_lon_deg=center_lon_deg,
            center_lat_deg=center_lat_deg,
        )
        print("  final cloned visible points:", len(vlam_clone))

        del a_rep, e_rep, inc_rep, raan_rep, argp_rep, tp_new
        gc.collect()

        n_overlay = min(overlay_max, len(vlam_clone))
        idx_overlay = rng.choice(len(vlam_clone), size=n_overlay, replace=False)
        clone_overlay[pop_name] = {
            "vlam": vlam_clone[idx_overlay].copy(),
            "vbeta": vbeta_clone[idx_overlay].copy(),
            "scatter_size": info["scatter_size"],
            "scatter_alpha": info["scatter_alpha"],
        }

        pts_clone = np.column_stack([vlam_clone, vbeta_clone])
        tree_clone = cKDTree(pts_clone)

        density_clone_flat = evaluate_density_map_full_posterior_2d(
            tree=tree_clone,
            grid_points=grid_points,
            k=k_map,
            n_d0_grid=n_d0_grid_map,
            show_progress=show_progress,
        )

        density_clone_map = density_clone_flat.reshape(X0.shape)
        density_downweighted_map = density_clone_map / f

        clone_density_maps[pop_name] = density_clone_map
        clone_density_maps_downweighted[pop_name] = density_downweighted_map

        print(f"{pop_name} done.")
        print("  raw cloned max density   :", np.nanmax(density_clone_map))
        print("  downweighted max density :", np.nanmax(density_downweighted_map))

        del vlam_clone, vbeta_clone, pts_clone, tree_clone, density_clone_flat
        gc.collect()

    density_all_downweighted = sum(clone_density_maps_downweighted.values())

    return {
        "label": center_label,
        "center_lon_deg": center_lon_deg,
        "center_lat_deg": center_lat_deg,
        "grid": grid,
        "overlay": clone_overlay,
        "density_maps": clone_density_maps,
        "density_maps_downweighted": clone_density_maps_downweighted,
        "density_all_downweighted": density_all_downweighted,
    }


# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------

def _get_grid_from_result(result):
    grid = result["grid"]
    return grid["X0"], grid["Y0"], grid["xlim"], grid["ylim"]


def plot_combined_density(result, floor=1e-30, figsize=(8, 8)):
    """
    Plot the summed downweighted density map from all populations.
    """
    X0, Y0, xlim, ylim = _get_grid_from_result(result)
    log_map = np.log10(np.maximum(result["density_all_downweighted"], floor))

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.pcolormesh(X0, Y0, log_map, shading="auto")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$v_\lambda$ (deg/day)")
    ax.set_ylabel(r"$v_\beta$ (deg/day)")
    ax.set_title(f"{result['label']}: combined cloned/downweighted log10 density")
    fig.colorbar(im, ax=ax, label=r"$\log_{10}\langle n_0 \rangle$")
    plt.show()



def plot_population_density_maps(result, floor=1e-30, figsize=(14, 12), with_overlay=False):
    """
    Plot one panel per population.

    Parameters
    ----------
    result : dict
        Output from build_cloned_maps_for_center.
    with_overlay : bool
        If True, overlay a random subset of cloned points.
    """
    X0, Y0, xlim, ylim = _get_grid_from_result(result)
    maps = result["density_maps_downweighted"]
    overlay_dict = result["overlay"]

    fig, axes = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)
    axes = axes.ravel()

    for ax, pop_name in zip(axes, maps.keys()):
        log_map = np.log10(np.maximum(maps[pop_name], floor))
        im = ax.pcolormesh(X0, Y0, log_map, shading="auto")

        if with_overlay:
            overlay = overlay_dict[pop_name]
            ax.scatter(
                overlay["vlam"],
                overlay["vbeta"],
                s=overlay["scatter_size"],
                alpha=overlay["scatter_alpha"],
                color="white",
                rasterized=True,
            )
            title = f"{pop_name}: cloned log-density + overlay sample"
        else:
            title = f"{pop_name}: cloned/downweighted log10 density"

        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel(r"$v_\lambda$ (deg/day)")
        ax.set_ylabel(r"$v_\beta$ (deg/day)")
        ax.set_title(title)

        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(r"$\log_{10}\langle n_0 \rangle$")

    plt.show()


# -----------------------------------------------------------------------------
# Convenience wrapper for the full workflow
# -----------------------------------------------------------------------------

def run_default_clone_workflow(
    obstime_str,
    center_lon_deg=180.0,
    center_lat_deg=0.0,
    center_label=None,
    max_sep_deg=30.0,
    clone_factors=None,
    xlim=(-0.8, 0.8),
    ylim=(-0.8, 0.8),
    grid_step=0.01,
    k_map=10,
    n_d0_grid_map=400,
    overlay_max=100_000,
    rng_seed=42,
    verbose=False,
    show_progress=True,
):
    """
    Full one-call workflow:
    1. load populations
    2. create clone source config
    3. build maps for one center

    Returns
    -------
    result : dict
        Output from build_cloned_maps_for_center.
    """
    if center_label is None:
        center_label = f"lon={center_lon_deg}, lat={center_lat_deg}"

    populations = load_default_populations(verbose=verbose)
    clone_sources = make_clone_sources(populations, clone_factors=clone_factors)
    rng = np.random.default_rng(rng_seed)

    return build_cloned_maps_for_center(
        center_lon_deg=center_lon_deg,
        center_lat_deg=center_lat_deg,
        center_label=center_label,
        clone_sources=clone_sources,
        obstime_str=obstime_str,
        max_sep_deg=max_sep_deg,
        rng=rng,
        overlay_max=overlay_max,
        xlim=xlim,
        ylim=ylim,
        grid_step=grid_step,
        k_map=k_map,
        n_d0_grid_map=n_d0_grid_map,
        show_progress=show_progress,
    )
