import numpy as np
from numpy import sin, cos

## automatically reload any modules read below that might have changed (e.g. plots)
import sys
sys.path.append('../src')
sys.path.append('src')
import pyarrow as pa
from scipy.interpolate import RegularGridInterpolator
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation, get_body_barycentric_posvel, solar_system_ephemeris 
from adam_core.time import Timestamp
from adam_core.coordinates import CartesianCoordinates, Origin
from adam_core.constants import KM_P_AU, S_P_DAY
from adam_core.orbits import Orbits
from adam_core.utils import get_perturber_state
from adam_core.coordinates.origin import OriginCodes
import NEOMOD3 as nm3  

 


# ---------- angular-rate converter: (dra, ddec) -> (v_lambda, v_beta) ----------
def radec_rates_to_ecliptic_rates(ra_deg, dec_deg, dra_deg_day, ddec_deg_day):
    """
    Convert equatorial sky motion to ecliptic sky motion at the same direction.

    Parameters
    ----------
    ra_deg, dec_deg : float
        Right ascension and declination of the direction (deg).
    dra_deg_day, ddec_deg_day : float
        Angular rates on the sky (deg/day). dra is d(alpha)*cos(delta)/day.

    Returns
    -------
    v_lambda, v_beta : float
        Ecliptic longitude and latitude rates (deg/day).
    """
    ra  = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)

    
    e_ra  = np.array([-np.sin(ra),  np.cos(ra), 0.0])
    e_dec = np.array([-np.sin(dec)*np.cos(ra), -np.sin(dec)*np.sin(ra), np.cos(dec)])

    # Angular velocity vector in equatorial Cartesian 
    dra  = np.deg2rad(dra_deg_day)  / 86400.0
    ddec = np.deg2rad(ddec_deg_day) / 86400.0
    w_eq = dra*e_ra + ddec*e_dec

    # Rotation from equatorial to ecliptic 
    eps = np.deg2rad(23.439291)
    R = np.array([[ 1,           0,            0],
                  [ 0,  np.cos(eps),  np.sin(eps)],
                  [ 0, -np.sin(eps),  np.cos(eps)]])
    w_ecl = R @ w_eq


    
    l_hat_eq = np.array([np.cos(dec)*np.cos(ra),
                         np.cos(dec)*np.sin(ra),
                         np.sin(dec)])
    l_hat_ecl = R @ l_hat_eq
    lam = np.arctan2(l_hat_ecl[1], l_hat_ecl[0])
    bet = np.arcsin(l_hat_ecl[2])

    e_lam = np.array([-np.sin(lam),  np.cos(lam),  0.0])
    e_bet = np.array([-np.sin(bet)*np.cos(lam),
                      -np.sin(bet)*np.sin(lam),
                       np.cos(bet)])

    vlam  = np.dot(w_ecl, e_lam) * 86400.0 * 180.0/np.pi
    vbeta = np.dot(w_ecl, e_bet) * 86400.0 * 180.0/np.pi
    return vlam, vbeta





def Zmakeddotgrids():
    d_min, d_max, Nd = 0.05, 5.0, 10          # AU  
    v_min, v_max, Nv = -60.0, 60.0, 10         # km/s

    d_grid   = np.linspace(d_min, d_max, Nd)
    ddot_grid= np.linspace(v_min, v_max, Nv)
    return d_grid, ddot_grid


def weight_marg_over_d_and_ddot(array4D, H_center, a_center, e_center, i_center,
                                        d_grid, ddot_grid, ra_deg, dec_deg, dra_deg_day, ddec_deg_day,
                                        H0=19.0, i0=25.9, obstime_str="2025-09-23T00:00:00"):
    w_sum = 0.0
    # loop 
    for d in d_grid:
        for vlos in ddot_grid:
            out = compute_nd3_Weight(array4D, H_center, a_center, e_center, i_center,
                            d_au=d, ddot_kms=vlos,
                            ra_deg=ra_deg, dec_deg=dec_deg,
                            dra_deg_per_day=dra_deg_day, ddec_deg_per_day=ddec_deg_day,
                            H0=H0, i0=i0, obstime_str=obstime_str)
            w = out[0] if isinstance(out, tuple) else float(out)
            w_sum += w
    return w_sum



#converts geocentric vectors to heliocentric vectors 
def Geo_to_topo(r_geo_km, v_geo_kms, rE, vE, r_obs, obstime):
        # convert r_geo, v_geo for obj to heliocentric r_topo, v_topo for an observer at rubin
        # from center of earth to an observer now at rubin
        obstime = Time(obstime)
        # earth's barycentric (heliocentric) state was done in previous step
        # so was observer's offset from center of earth
        # adding the offsets
        r_topo = rE + r_obs + r_geo_km
        v_topo = vE + v_geo_kms
        return r_topo, v_topo # asteroid position, vec going from sun to asteroid assuming youre at rubin


# converts vectors to geocentric state vectors 
def observables_to_geocentric_state(ra_deg, dec_deg, dra_deg_per_day, ddec_deg_per_day,
                                    d_au, ddot_kms, ra_rate_is_plain_dalpha=False):
        AU_km = 149_597_870.7
        mu_sun = 1.32712440018e11  # [km^3/s^2] gm of the sun


        # converting from degrees to radians
        ra  = np.deg2rad(ra_deg) 
        dec = np.deg2rad(dec_deg)
        l_hat, e_ra, e_dec = Sky_basis(ra, dec)
         # angular rate conversion from deg/day to rad/s
        dra  = np.deg2rad(dra_deg_per_day)  / 86400.0
        ddec = np.deg2rad(ddec_deg_per_day) / 86400.0
        # if the RA rate is plain d alpha/dt, convert to u alpha⋅cos delta sky-plane rate
        if ra_rate_is_plain_dalpha:
             dra *= np.cos(dec)

        d_km = d_au * AU_km
        # multiplying l_hat by the rand dist. value gives actual geocentric pos vector of obj
        r_geo = d_km * l_hat  # vector from earth center to obj

        # speed
        v_rad = ddot_kms * l_hat
        v_tan = d_km * (dra * e_ra + ddec * e_dec)
        v_geo = v_rad + v_tan
        
        return r_geo, v_geo


    
def get_Earth_and_observer(obstime_str):
    obstime = Time(obstime_str)
    
    rubin_location = EarthLocation.from_geodetic( 
        lon=-70.7366*u.deg,   # longitude (west is negative)
        lat=-30.2407*u.deg,   # latitude
        height=2647*u.m       # elevation
        
    )
    obs_geo = rubin_location.get_gcrs(obstime)
    r_obs = np.array([obs_geo.cartesian.x.to(u.km).value,
                      obs_geo.cartesian.y.to(u.km).value,
                      obs_geo.cartesian.z.to(u.km).value])
    
    # earth's barycentric vector rel to sun aka center of solsys to earth
    with solar_system_ephemeris.set('de432s'): 
        rE_bary, vE_bary = get_body_barycentric_posvel('earth', obstime)

    # convert to km and km/s arrays (in Cartesian coordinate system) 
    rE = np.array([rE_bary.x.to(u.km).value, #rE goes from sun to earth center
                   rE_bary.y.to(u.km).value,
                   rE_bary.z.to(u.km).value])
    
    vE = np.array([vE_bary.x.to(u.km/u.s).value, #vE goes from sun to earth center
                   vE_bary.y.to(u.km/u.s).value,
                   vE_bary.z.to(u.km/u.s).value])
    return obstime, rE, vE, r_obs


# creates a right-handed vector basis that points to object from observer
def Sky_basis(ra_rad, dec_rad):
        l_hat = np.array([cos(dec_rad)*cos(ra_rad), # vector that points to object's pos on sky from observer
                          cos(dec_rad)*sin(ra_rad),
                          sin(dec_rad)])
        e_ra  = np.array([-sin(ra_rad),  cos(ra_rad), 0.0]) # vector towards increasing RA - points east, 90 deg to l_hat
        e_dec = np.array([-sin(dec_rad)*cos(ra_rad), # vector towards increasing dec - points north
                          -sin(dec_rad)*sin(ra_rad),
                          cos(dec_rad)])
        return l_hat, e_ra, e_dec # turning angular rates into linear velocities


def to_keplerian_adam(r_helio, v_helio, obstime):
    r_vec = np.array(r_helio, dtype=float)
    v_vec = np.array(v_helio, dtype=float)

    if r_vec.ndim == 1:
        r_vec = r_vec[np.newaxis, :]
        v_vec = v_vec[np.newaxis, :]

    cartesian_coordinates = CartesianCoordinates.from_kwargs(
    x=r_vec[:, 0] / KM_P_AU,
    y=r_vec[:, 1] / KM_P_AU,
    z=r_vec[:, 2] / KM_P_AU,
    vx=v_vec[:, 0] / KM_P_AU * S_P_DAY,
    vy=v_vec[:, 1] / KM_P_AU * S_P_DAY,
    vz=v_vec[:, 2] / KM_P_AU * S_P_DAY,
    time=Timestamp.from_astropy(obstime),
    origin=Origin.from_kwargs(code=pa.repeat("SUN", len(r_vec))),
    frame="equatorial" # ecliptic is the other choice
    )

    keplerian_coordinates = cartesian_coordinates.to_keplerian()


    # Might be useful to store them into an orbits object:
    orbits = Orbits.from_kwargs(
        orbit_id=[f"orbit_{i:08d}" for i in range(len(cartesian_coordinates))],
        coordinates=cartesian_coordinates,
    )

    a_val = keplerian_coordinates.a[0].as_py()
    e_val = keplerian_coordinates.e[0].as_py()
    i_val = keplerian_coordinates.i[0].as_py()
    return a_val, e_val, i_val


def score_from_elements(a_AU, e, H0, i0, array4D): #without interpolation
    # checks if orbit is bound right
    if (a_AU <= 0) or (e >= 1.0) or (e < 0):
        
        return 0.0

    # clip to model box
    a_cl = np.clip(a_AU, 0.0, 4.2 - 1e-9)
    e_cl = np.clip(e,    0.0, 1.0 - 1e-9)
    i_cl = np.clip(i0,  0.0, 88.0 - 1e-9)

    log_w = float(array4D([[H0, a_cl, e_cl, i_cl]])[0])
    return 0.0 if not np.isfinite(log_w) else float(np.exp(log_w))

def score_from_elements_interp(a_AU, e, H0, i0, interp4D):
    # checks if orbit is bound right
    if (a_AU <= 0) or (e >= 1.0) or (e < 0):
        
        return 0.0

    # clip gently to model box
    a_cl = np.clip(a_AU, 0.0, 4.2 - 1e-9)
    e_cl = np.clip(e,      0.0, 1.0 - 1e-9)
    i_cl = np.clip(i0,  0.0, 88.0 - 1e-9)

    log_w = float(interp4D([[H0, a_cl, e_cl, i_cl]])[0])
    return 0.0 if not np.isfinite(log_w) else float(np.exp(log_w))


def compute_nd3_Weight(array4D, H_center, a_center, e_center, i_center,
                           d_au, ddot_kms, 
                           ra_deg=10.0, dec_deg=20.0, 
                           dra_deg_per_day=0.12, ddec_deg_per_day=0.22,
                           H0=19.0, i0=25.9, obstime_str="2025-09-16T00:00:00", what="w"):

    # dummy   
    score, a_val, e_val, i_val = 0.1, 1.0, 0.01, 10.0 

    ## do coordinate transformations for the observer and the object
    obstime, rE, vE, r_obs = get_Earth_and_observer(obstime_str)  
    r_geo, v_geo = observables_to_geocentric_state(ra_deg, dec_deg, dra_deg_per_day, ddec_deg_per_day, d_au, ddot_kms)
    # heliocentric position and velocity for the object
    r_helio, v_helio = Geo_to_topo(r_geo, v_geo, rE, vE, r_obs, obstime) 

    ## now use ADAM to get Keplerian orbital elements from observed object's position and velocity vectors
    obstime = Time([obstime_str], scale="tdb")  # Note converted this to a list 
    a_val, e_val, i_val = to_keplerian_adam(r_helio, v_helio, obstime)


    # WILL THIS WORK ??
    #to calc center of bins
    n_H, H_min, H_max = 52, 15.0, 28.0
    n_a, a_min, a_max = 42, 0.0, 4.2
    n_e, e_min, e_max = 25, 0.0, 1.0
    n_i, i_min, i_max = 22, 0.0, 88.0

    Na, Ne = 200, 200
    a_grid = np.linspace(a_min, a_max, Na)
    e_grid = np.linspace(e_min, e_max, Ne)
    A, E = np.meshgrid(a_grid, e_grid, indexing="xy")

    interpolate = True
    if (interpolate):
        # so that log(0) doesn't happen
        _eps = 1e-300
        log_grid = np.log(array4D + _eps)
        interp4D = RegularGridInterpolator(
            (H_center, a_center, e_center, i_center),
            log_grid, bounds_error=False, fill_value=-np.inf
        )
        score = score_from_elements_interp(a_val, e_val, 19.0, i_val, interp4D)
    else:
        score = score_from_elements(a_val, e_val, 19.0, i_val, array4D)
 
    # print('score:', score) 
    return score, a_val, e_val, i_val


    
