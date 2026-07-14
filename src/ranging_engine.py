#!/usr/bin/env python3
"""
ranging_engine.py — systematic-ranging NEO/non-NEO classifier (item 1A).

For each 2-detection tracklet the attributable A=(alpha,delta,alpha_dot,delta_dot)
is EXACT (4 data = 4 params, chi^2==0), so there is no per-node least-squares fit
(unlike Farnocchia 2015 / Spoto 2018, who fit >=3 noisy obs). The whole computation
is closed-form and vectorised over (tracklet x grid-node). At each admissible
(rho, rho_dot) node we convert to heliocentric (r,v) -> orbital elements (a,e,i,q)
and an absolute magnitude H(rho), look up NEO and non-NEO population densities, and
form a weighted class ratio:

    P_NEO = sum_nodes w * f_NEO / sum_nodes w * (f_NEO + f_nonNEO)

Design + paper sources: docs/D2_detail.md. Every knob is an explicit argument.

Frames/units: geometry done in km, km/s. Inputs RA/Dec are ICRS (equatorial);
we rotate r,v to the ECLIPTIC frame before computing inclination so that i matches
NEOMOD3/S3M (ecliptic) — the old NEO_H.py used equatorial i, a ~23 deg error.
"""
import os
import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import (EarthLocation, get_body_barycentric_posvel,
                                 solar_system_ephemeris)

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEFAULT_CACHE = os.path.join(_ROOT, "outputs", "pop_cache_wide.npz")

# ---- constants ----
AU_KM = 149_597_870.7
MU_SUN = 1.32712440018e11          # km^3/s^2
C_KM_S = 299_792.458
OBLIQ = np.deg2rad(23.439291)      # J2000 mean obliquity
A_MAX_AU = 100.0                   # admissible-region outer bound (Spoto a_max)
RUBIN = EarthLocation.from_geodetic(lon=-70.7366*u.deg, lat=-30.2407*u.deg, height=2647*u.m)

# equatorial(ICRS) -> ecliptic rotation matrix R_x(+obliquity)
_R_EQ2ECL = np.array([[1, 0, 0],
                      [0, np.cos(OBLIQ),  np.sin(OBLIQ)],
                      [0, -np.sin(OBLIQ), np.cos(OBLIQ)]])


# ------------------------------------------------------------------ observer/tracklet
def earth_observer_state(obstime, ephem='de432s'):
    """SUN-centred (heliocentric) observer position AND velocity at obstime.
    Two corrections the old NEO_H.py missed, each worth ~1-3% in recovered elements:
      (1) subtract the Sun's barycentric state (Sun is up to ~1.5e6 km / ~0.4 km/s off
          the SS barycentre) — elements are heliocentric, not barycentric;
      (2) include the observer's diurnal position+velocity (Earth rotation ~0.46 km/s;
          since a ~ v^2, omitting it is a ~3% error in a).
    Accepts scalar or array-like obstime. Returns r_obs_helio [km], v_obs_helio [km/s],
    each (...,3) in EQUATORIAL (ICRS/GCRS) frame."""
    if isinstance(obstime, Time):
        t = obstime
    else:
        ot = np.asarray(obstime)
        t = Time(ot.astype(float), format="mjd") if np.issubdtype(ot.dtype, np.number) else Time(ot)
    p_geo, v_geo = RUBIN.get_gcrs_posvel(t)                             # earth->observer pos+vel
    r_obs_geo = p_geo.xyz.to_value(u.km).T
    v_obs_geo = v_geo.xyz.to_value(u.km/u.s).T
    with solar_system_ephemeris.set(ephem):
        rE, vE = get_body_barycentric_posvel('earth', t)               # bary->earth
        rS, vS = get_body_barycentric_posvel('sun', t)                 # bary->sun
    rE = rE.xyz.to_value(u.km).T; vE = vE.xyz.to_value(u.km/u.s).T
    rS = rS.xyz.to_value(u.km).T; vS = vS.xyz.to_value(u.km/u.s).T
    return (rE - rS) + r_obs_geo, (vE - vS) + v_obs_geo                 # sun->observer pos, vel


def sky_basis(ra_rad, dec_rad):
    """Right-handed (l_hat, e_ra, e_dec) at each sky position. Arrays broadcast."""
    cd, sd = np.cos(dec_rad), np.sin(dec_rad)
    ca, sa = np.cos(ra_rad), np.sin(ra_rad)
    l_hat = np.stack([cd*ca, cd*sa, sd], axis=-1)
    e_ra  = np.stack([-sa, ca, np.zeros_like(sa)], axis=-1)
    e_dec = np.stack([-sd*ca, -sd*sa, cd], axis=-1)
    return l_hat, e_ra, e_dec


def tracklet_vectors(ra_deg, dec_deg, dra_deg_day, ddec_deg_day, dra_cosdec=True):
    """Unit line-of-sight l_hat and its angular-motion vector v_ang [1/s] per tracklet.
    v_ang has units 1/s so that d/dt(rho*l_hat)_tangential = rho[km]*v_ang.

    dra_cosdec: RA-RATE CONVENTION — declare it, do not guess. e_ra is a UNIT east vector,
      so the tangential velocity needs alpha_dot*cos(dec) in that direction.
        True  -> dra_deg_day is already alpha_dot*cos(dec)  (proper-motion convention;
                 e.g. Kurlander's `RARateCosDec_deg_day`).
        False -> dra_deg_day is the RAW alpha_dot; we multiply by cos(dec) here
                 (e.g. the v5 parquet's `mean_dra`, verified 2026-07-08: corr 0.99999 with
                 the wrap-safe finite difference of ra0,ra1 with NO cos(dec)).
      Getting this wrong inflates the RA tangential velocity by 1/cos(dec) (~3% at |dec|=13
      deg, ~15% at 30 deg) and corrupts the recovered elements at the few-% level.
    """
    ra, dec = np.deg2rad(ra_deg), np.deg2rad(dec_deg)
    l_hat, e_ra, e_dec = sky_basis(ra, dec)
    dra_pm = np.asarray(dra_deg_day, float)
    if not dra_cosdec:
        dra_pm = dra_pm * np.cos(dec)
    dra = np.deg2rad(dra_pm) / 86400.0                                   # rad/s
    ddec = np.deg2rad(ddec_deg_day) / 86400.0
    v_ang = dra[..., None]*e_ra + ddec[..., None]*e_dec                  # (...,3) [1/s]
    return l_hat, v_ang


# ------------------------------------------------------------------ admissible grid
def build_grid(l_hat, v_ang, r_obs, v_earth, n_rho=64, n_rhodot=32,
               rho_min=0.01, rho_max=100.0, a_max_au=A_MAX_AU):
    """Per-tracklet log-rho columns x AR-restricted rho_dot rows.

    For each rho, the bounded-orbit energy condition E < -mu/(2 a_max) is a quadratic
    in rho_dot (r_helio is fixed by rho): rho_dot^2 + 2 b rho_dot + c < 0. We sample
    n_rhodot points uniformly inside the admissible interval and carry its length as
    the node weight factor dRhoDot. Empty-interval columns contribute nothing.

    All inputs per-tracklet: l_hat (N,3), v_ang (N,3) [1/s], r_obs (3,), v_earth (3,).
    Returns dict of (N, n_rho, n_rhodot) arrays: rho [au], rho_dot [km/s], dRhoDot [km/s],
    admissible (bool), plus rho_edges for the log measure.
    """
    N = l_hat.shape[0]
    r_obs = np.atleast_2d(r_obs); v_earth = np.atleast_2d(v_earth)       # (N,3) per-tracklet
    if r_obs.shape[0] == 1:
        r_obs = np.broadcast_to(r_obs, (N, 3)); v_earth = np.broadcast_to(v_earth, (N, 3))
    log_edges = np.linspace(np.log10(rho_min), np.log10(rho_max), n_rho + 1)
    rho_au = 10 ** (0.5 * (log_edges[1:] + log_edges[:-1]))               # (n_rho,) centres
    dlog = log_edges[1] - log_edges[0]
    rho_km = rho_au * AU_KM

    # heliocentric object position r_helio = r_obs + rho*l_hat  (per tracklet, per rho)
    r_helio = r_obs[:, None, :] + rho_km[None, :, None] * l_hat[:, None, :]      # (N,n_rho,3)
    r = np.linalg.norm(r_helio, axis=-1)                                 # (N,n_rho)

    # rho_dot-independent velocity part: v0 = v_earth + rho*v_ang  (tangential)
    v0 = v_earth[:, None, :] + rho_km[None, :, None] * v_ang[:, None, :]         # (N,n_rho,3)
    b = np.einsum('nrk,nk->nr', v0, l_hat)                               # v0 . l_hat
    v0sq = np.einsum('nrk,nrk->nr', v0, v0)
    a_max_km = a_max_au * AU_KM
    c = v0sq - 2*MU_SUN/r + MU_SUN/a_max_km                              # energy quadratic const
    disc = b*b - c
    admissible_col = disc > 0
    sq = np.sqrt(np.where(admissible_col, disc, 0.0))
    lo, hi = -b - sq, -b + sq                                            # rho_dot interval [km/s]

    # sample n_rhodot uniformly inside [lo,hi] per (tracklet,rho)
    frac = (np.arange(n_rhodot) + 0.5) / n_rhodot                        # (n_rhodot,)
    rho_dot = lo[..., None] + (hi - lo)[..., None] * frac                # (N,n_rho,n_rhodot)
    dRhoDot = np.where(admissible_col, (hi - lo) / n_rhodot, 0.0)[..., None] * np.ones(n_rhodot)
    admissible = np.broadcast_to(admissible_col[..., None], rho_dot.shape)

    return dict(rho_au=np.broadcast_to(rho_au[None, :, None], rho_dot.shape).copy(),
                rho_dot=rho_dot, dRhoDot=dRhoDot, admissible=admissible,
                r_helio=r_helio, dlog=dlog, rho_km=rho_km,
                r_obs=r_obs, v_earth=v_earth)


# ------------------------------------------------------------------ elements from nodes
def elements_from_nodes(grid, l_hat, v_ang):
    """Closed-form (a,e,i,q) at every node. Vectorised, ecliptic inclination.

    Returns dict of (N,n_rho,n_rhodot): a_au, e, inc_deg, q_au, bound (bool).
    """
    rho_km = grid["rho_km"]                                             # (n_rho,)
    r_helio = grid["r_helio"]                                           # (N,n_rho,3)
    rho_dot = grid["rho_dot"]                                           # (N,n_rho,n_rhodot)
    v_earth = grid["v_earth"]                                           # (N,3)

    # heliocentric velocity: v = v_obs + rho_dot*l_hat + rho*v_ang
    v = (v_earth[:, None, None, :]
         + rho_dot[..., None] * l_hat[:, None, None, :]
         + (rho_km[None, :, None, None] * v_ang[:, None, None, :]))     # (N,n_rho,n_rhodot,3)
    r_vec = np.broadcast_to(r_helio[:, :, None, :], v.shape)            # (N,n_rho,n_rhodot,3)

    # rotate to ecliptic so inclination matches NEOMOD3/S3M
    r_ecl = r_vec @ _R_EQ2ECL.T
    v_ecl = v @ _R_EQ2ECL.T

    r = np.linalg.norm(r_ecl, axis=-1)
    vmag2 = np.einsum('...k,...k->...', v_ecl, v_ecl)
    energy = 0.5 * vmag2 - MU_SUN / r
    a_km = -MU_SUN / (2 * energy)                                       # >0 bound, <0/inf unbound

    h_vec = np.cross(r_ecl, v_ecl)
    h = np.linalg.norm(h_vec, axis=-1)
    # eccentricity vector e_vec = (v x h)/mu - r_hat
    e_vec = np.cross(v_ecl, h_vec) / MU_SUN - r_ecl / r[..., None]
    e = np.linalg.norm(e_vec, axis=-1)

    with np.errstate(invalid='ignore', divide='ignore'):
        inc = np.degrees(np.arccos(np.clip(h_vec[..., 2] / h, -1, 1)))
        a_au = a_km / AU_KM
        q_au = a_au * (1 - e)

    bound = (energy < 0) & (e < 1.0) & (a_au > 0) & np.isfinite(a_au)
    return dict(a_au=a_au, e=e, inc_deg=inc, q_au=q_au, bound=bound)


# ------------------------------------------------------------------ H from nodes
def H_from_nodes(grid, V_mag, G=0.15):
    """Absolute magnitude H at every node via HG phase function (G default 0.15).
    V_mag per tracklet (already colour-corrected to V). Returns (N,n_rho,n_rhodot)."""
    rho_km = grid["rho_km"]                                             # (n_rho,)
    r_helio = grid["r_helio"]                                           # (N,n_rho,3)
    r_obs = grid["r_obs"]                                               # (N,3)
    r_sun_obj_au = np.linalg.norm(r_helio, axis=-1) / AU_KM             # (N,n_rho)
    delta_au = rho_km / AU_KM                                           # (n_rho,) observer-obj

    # phase angle alpha at the object: between (obj->sun) and (obj->observer)
    obj_to_sun = -r_helio                                               # (N,n_rho,3)
    obj_to_obs = (r_obs[:, None, :] - r_helio)                          # (N,n_rho,3)
    cos_alpha = (np.einsum('nrk,nrk->nr', obj_to_sun, obj_to_obs)
                 / (np.linalg.norm(obj_to_sun, axis=-1) * np.linalg.norm(obj_to_obs, axis=-1)))
    alpha = np.arccos(np.clip(cos_alpha, -1, 1))                        # (N,n_rho)

    ta2 = np.tan(np.clip(alpha, 0, np.pi - 1e-6) / 2)
    phi1 = np.exp(-3.33 * ta2 ** 0.63)
    phi2 = np.exp(-1.87 * ta2 ** 1.22)
    phase = -2.5 * np.log10((1 - G) * phi1 + G * phi2)                  # magnitudes (>=0)

    # H = V - 5 log10(r_sun_obj * delta) - phase_correction ; broadcast delta over rho
    dist_term = 5 * np.log10(r_sun_obj_au * delta_au[None, :])          # (N,n_rho)
    H = V_mag[:, None] - dist_term - phase                             # (N,n_rho)
    return np.broadcast_to(H[:, :, None], grid["rho_dot"].shape)        # (N,n_rho,n_rhodot)


# ------------------------------------------------------------------ node weights (prior)
def node_weights(grid, gamma=2.0):
    """Geometric prior weight per node: log-measure (rho*ln10*dlog) x rho^gamma
    x dRhoDot (uniform-in-rho_dot). Zero where inadmissible. (N,n_rho,n_rhodot).
    Jeffreys prior deliberately NOT used (Farnocchia 2015 Table 3)."""
    rho = grid["rho_au"]
    w = (rho * np.log(10) * grid["dlog"]) * (rho ** gamma) * grid["dRhoDot"]
    return np.where(grid["admissible"], w, 0.0)


# ================================================================== scoring
def neomod3_to_wide(edges):
    """Resample the NEOMOD3 native (H,a,e,i) count table onto the wide grid
    (H, log10 a, e, i). Returns count array shaped like the wide-grid hists."""
    import neomod3_sampler as nm
    arr, cen, _ = nm.load_neomod3_array()                              # (52,42,25,22)
    Hc, ac, ec, ic = cen["H"], cen["a"], cen["e"], cen["i"]
    # cell-centre coordinates + weights, histogram onto wide edges
    HH, AA, EE, II = np.meshgrid(Hc, ac, ec, ic, indexing="ij")
    w = arr.ravel()
    m = w > 0
    pts = np.vstack([HH.ravel()[m], np.log10(AA.ravel()[m]),
                     EE.ravel()[m], II.ravel()[m]]).T
    wide, _ = np.histogramdd(
        pts, bins=[edges["H_edges"], edges["loga_edges"], edges["e_edges"], edges["i_edges"]],
        weights=w[m])
    return wide


def load_tables(level, cache_path=DEFAULT_CACHE, wagg_mba=1.0, neomod3_norm="per_H_match"):
    """Return (neo_counts, nonneo_counts, edges) on the wide grid for a given level.
      L0 -> tables unused (geometric); returns (None, None, edges).
      L1 -> neo = S3M-NEO, nonneo = S3M (MBA*wagg + Trojan + TNO). digest2 replication.
      L2 -> neo = NEOMOD3(wide), normalised per `neomod3_norm`; nonneo as L1.

    neomod3_norm (the 1A-vs-2B lever; see QA0c):
      'per_H_match' (DEFAULT, = item 1A) — rescale NEOMOD3 per H-bin so its total matches
          S3M-NEO's in that bin. L2 then differs from L1 ONLY in the conditional orbital
          shape (a,e,i | H). Clean isolated test of the debiased orbital distribution.
          Consequence: inherits S3M-NEO's H~25 cutoff, so it cannot show the faint-end gain.
      'absolute' (= item 2B territory) — use NEOMOD3's own absolute N(H), which extends to
          H=28. Exposes the faint-end debiased-count effect.
          *** CAVEAT: S3M's non-NEO census also has an H cutoff (MBA H<~24-25). At faint H the
          denominator empties while the NEOMOD3 numerator does not, driving P_NEO->1. Part of
          that is REAL (Wagg 2025: faint candidates are preferentially NEO because H-degeneracy
          removes faint MBAs) and part is a model-limit artifact. Do not report 'absolute'
          numbers without an explicit faint-H control. ***
    """
    z = np.load(cache_path)
    edges = {k: z[k] for k in ("H_edges", "loga_edges", "e_edges", "i_edges")}
    nonneo = wagg_mba * z["hist_mba"] + z["hist_trojan"] + z["hist_tno"]
    if level == "L0":
        return None, None, edges
    if level == "L1":
        return z["hist_neo"].copy(), nonneo, edges
    if level == "L2":
        nm_wide = neomod3_to_wide(edges)
        s3m_neo = z["hist_neo"]
        if neomod3_norm == "per_H_match":
            N_s3m = s3m_neo.sum(axis=(1, 2, 3))
            N_nm = nm_wide.sum(axis=(1, 2, 3))
            scale = np.where(N_nm > 0, N_s3m / np.maximum(N_nm, 1e-30), 0.0)
            return nm_wide * scale[:, None, None, None], nonneo, edges
        if neomod3_norm == "absolute":
            return nm_wide.copy(), nonneo, edges
        raise ValueError(f"neomod3_norm={neomod3_norm}")
    raise ValueError(level)


def _bin_idx(x, edges):
    idx = np.searchsorted(edges, x) - 1
    valid = (idx >= 0) & (idx < len(edges) - 1)
    return np.clip(idx, 0, len(edges) - 2), valid


def class_score(elems, H, w, neo_tab, nonneo_tab, edges, level):
    """P_NEO per tracklet from node elements/H/weights. Arrays (N,nr,nrd).
    L0 = geometric (q<1.3 indicator); L1/L2 = population count ratio in wide bins."""
    bound = elems["bound"]
    q = elems["q_au"]
    if level == "L0":
        wn = w * bound
        num = (wn * (q < 1.3)).sum(axis=(1, 2))
        den = wn.sum(axis=(1, 2))
        return np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)

    iH, vH = _bin_idx(H, edges["H_edges"])
    il, vl = _bin_idx(np.log10(np.maximum(elems["a_au"], 1e-6)), edges["loga_edges"])
    ie, ve = _bin_idx(elems["e"], edges["e_edges"])
    ii, vi = _bin_idx(elems["inc_deg"], edges["i_edges"])
    valid = bound & vH & vl & ve & vi
    fneo = np.where(valid, neo_tab[iH, il, ie, ii], 0.0)
    fnon = np.where(valid, nonneo_tab[iH, il, ie, ii], 0.0)
    num = (w * fneo).sum(axis=(1, 2))
    den = (w * (fneo + fnon)).sum(axis=(1, 2))
    return np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)


def score_tracklets(ra, dec, dra, ddec, V, obstime, level="L1",
                    n_rho=64, n_rhodot=32, gamma=2.0, wagg_mba=1.0,
                    cache_path=DEFAULT_CACHE, chunk=2000, G=0.15,
                    dra_cosdec=True, neomod3_norm="per_H_match"):
    """Top-level driver. ra/dec [deg], dra/ddec [deg/day, dra=alpha_dot*cos(dec)],
    V [mag, colour-corrected], obstime [array of MJD or ISO]. Returns P_NEO (N,)."""
    ra = np.asarray(ra, float); dec = np.asarray(dec, float)
    dra = np.asarray(dra, float); ddec = np.asarray(ddec, float); V = np.asarray(V, float)
    N = len(ra)
    neo_tab, nonneo_tab, edges = load_tables(level, cache_path, wagg_mba, neomod3_norm)
    if isinstance(obstime, Time):
        t = obstime
    else:
        ot = np.asarray(obstime)
        t = Time(ot.astype(float), format="mjd") if np.issubdtype(ot.dtype, np.number) else Time(ot)
    out = np.full(N, np.nan)
    for s in range(0, N, chunk):
        e = min(s + chunk, N)
        r_obs, v_obs = earth_observer_state(t[s:e])
        l_hat, v_ang = tracklet_vectors(ra[s:e], dec[s:e], dra[s:e], ddec[s:e], dra_cosdec)
        grid = build_grid(l_hat, v_ang, r_obs, v_obs, n_rho, n_rhodot)
        elems = elements_from_nodes(grid, l_hat, v_ang)
        Hn = H_from_nodes(grid, V[s:e], G=G)
        w = node_weights(grid, gamma=gamma)
        out[s:e] = class_score(elems, Hn, w, neo_tab, nonneo_tab, edges, level)
    return out
