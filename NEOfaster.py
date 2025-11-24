import numpy as np
from numpy import sin, cos
import sys
import pyarrow as pa
from scipy.interpolate import RegularGridInterpolator
from astropy.time import Time
import astropy.units as units
from astropy.coordinates import EarthLocation, get_body_barycentric_posvel, solar_system_ephemeris
from astropy.coordinates import SkyCoord, get_sun, GCRS, GeocentricTrueEcliptic
from adam_core.time import Timestamp
from adam_core.coordinates import CartesianCoordinates, Origin
from adam_core.constants import KM_P_AU, S_P_DAY
from adam_core.orbits import Orbits
from adam_core.utils import get_perturber_state
from adam_core.coordinates.origin import OriginCodes
import NEOMOD3 as nm3  


#constants
AU_km = 149_597_870.7
mu_sun = 1.32712440018e11  # [km^3/s^2] gm of the sun



# creates a right-handed vector basis that points to object from observer
def sky_basis(ra_rad, dec_rad): #its already fast w 5.32 us
        l_hat = np.array([np.cos(dec_rad)*np.cos(ra_rad), # vector that points to object's pos on sky from observer
                          np.cos(dec_rad)*np.sin(ra_rad),
                          np.sin(dec_rad)])
        e_ra  = np.array([-np.sin(ra_rad), np.cos(ra_rad), 0.0]) # vector towards increasing RA - points east, 90 deg to l_hat
        e_dec = np.array([-np.sin(dec_rad)*np.cos(ra_rad), # vector towards increasing dec - points north
                          -np.sin(dec_rad)*np.sin(ra_rad),
                          np.cos(dec_rad)])
        return l_hat, e_ra, e_dec # turning angular rates into linear velocities


    
# ---------- angular-rate converter: (dra, ddec) -> (v_lambda, v_beta) ----------
def radec_rates_to_ecliptic_rates1(ra_deg, dec_deg, dra_deg_day, ddec_deg_day):

    
    sc_eq = SkyCoord(
        ra=ra_deg * units.deg,
        dec=dec_deg * units.deg,
        pm_ra_cosdec=dra_deg_day * units.deg / units.day,
        pm_dec=ddec_deg_day * units.deg / units.day,
        frame="icrs"
    )

    
    sc_ecl = sc_eq.transform_to(GeocentricTrueEcliptic())

    
    vlam  = sc_ecl.pm_lon_coslat.to(units.deg / units.day).value
    vbeta = sc_ecl.pm_lat.to(units.deg / units.day).value

    return vlam, vbeta



def ecliptic_rates_to_radec_rates(ra_deg, dec_deg, vlam_deg_day, vbeta_deg_day, obstime_str=None):

    
    if obstime_str is not None:
        sc_ecl = SkyCoord(
            lon=ra_deg * units.deg,   
            lat=dec_deg * units.deg,
            pm_lon_coslat=vlam_deg_day * units.deg/units.day,
            pm_lat=vbeta_deg_day * units.deg/units.day,
            frame=GeocentricTrueEcliptic(obstime=obstime_str)
        )
    else:
        sc_ecl = SkyCoord(
            lon=ra_deg * units.deg, lat=dec_deg * units.deg,
            pm_lon_coslat=vlam_deg_day * units.deg/units.day,
            pm_lat=vbeta_deg_day * units.deg/units.day,
            frame="geocentrictrueecliptic"
        )

    
    sc_eq = sc_ecl.transform_to("icrs")

    dra  = sc_eq.pm_ra_cosdec.to(units.deg/units.day).value
    ddec = sc_eq.pm_dec.to(units.deg/units.day).value

    return dra, ddec





def makeddotgrids():
    d_min, d_max, Nd = 0.05, 5.0, 10          # AU  
    v_min, v_max, Nv = -60.0, 60.0, 10         # km/s

    d_grid   = np.linspace(d_min, d_max, Nd)
    ddot_grid= np.linspace(v_min, v_max, Nv)
    return d_grid, ddot_grid


def weight_marginalized_over_d_and_ddot(array4D, H_center, a_center, e_center, i_center,
                                        ra_deg, dec_deg, dra_deg_per_day, ddec_deg_per_day, obstime_str
                                        ,mag):
    w_sum = 0.0
    d_grid, ddot_grid = makeddotgrids()
    # loop 
    obstime, rE, vE, r_obs = get_earth_and_observer(obstime_str) 
    l_hat, v_hat = compute_unit_vectors(ra_deg, dec_deg, dra_deg_per_day, ddec_deg_per_day)
    for d_au in d_grid:
        for ddot_kms in ddot_grid:
            out = compute_neomod3_weight(d_au, ddot_kms, l_hat, v_hat, mag, obstime,
                           rE, vE, r_obs,
                           array4D, H_center, a_center, e_center, i_center)
            w = out[0] if isinstance(out, tuple) else float(out)
            w_sum += w
    #print(ra_deg, dec_deg, dra_deg_per_day, ddec_deg_per_day, l_hat, v_hat,w_sum)
    return w_sum
    



#converts geocentric vectors to heliocentric vectors 
def geo_to_topo(r_geo_km, v_geo_kms, rE, vE, r_obs, obstime):
        # convert r_geo, v_geo for obj to heliocentric r_topo, v_topo for an observer at rubin
        # from center of earth to an observer now at rubin
        
        # earth's barycentric (heliocentric) state was done in previous step
        # so was observer's offset from center of earth
        # adding the offsets
        r_topo = rE + r_obs + r_geo_km
        v_topo = vE + v_geo_kms
        return r_topo, v_topo # asteroid position, vec going from sun to asteroid assuming youre at rubin


# converts vectors to geocentric state vectors 
def compute_unit_vectors(ra_deg, dec_deg, dra_deg_per_day, ddec_deg_per_day): 
        
        # converting from degrees to radians
        ra  = np.deg2rad(ra_deg) 
        dec = np.deg2rad(dec_deg)
        l_hat, e_ra, e_dec = sky_basis(ra, dec)
         # angular rate conversion from deg/day to rad/s
        dra  = np.deg2rad(dra_deg_per_day)  / 86400.0
        ddec = np.deg2rad(ddec_deg_per_day) / 86400.0

        v_hat = (dra * e_ra + ddec * e_dec)

        return l_hat, v_hat

# converts vectors to geocentric state vectors 
def observables_to_geocentric_state(l_hat, v_hat, d_au, ddot_kms): 
    
        

        d_km = d_au * AU_km
        # multiplying l_hat by the rand dist. value gives actual geocentric pos vector of obj
        r_geo = np.outer(d_km, l_hat)  # vector from earth center to obj
        #r_geo = (d_au * AU_km) * l_hat
        #v_geo = ddot_kms * l_hat + (d_au * AU_km) * v_hat

        # speed
        v_rad = np.outer(ddot_kms, l_hat)
        v_tan = np.outer(d_km, v_hat)
        #v_tan = d_km * v_hat
        v_geo = v_rad + v_tan
        
        return r_geo, v_geo

   
rubin_location = EarthLocation.from_geodetic( 
        lon=-70.7366*units.deg,   # longitude (west is negative)
        lat=-30.2407*units.deg,   # latitude
        height=2647*units.m       # elevation
        
    )

def get_earth_and_observer(obstime_str): # 1.47 ms -> 1.27 ms
    obstime = Time(obstime_str)
    
    obs_geo = rubin_location.get_gcrs(obstime)
    r_obs =  obs_geo.cartesian.xyz.to_value(units.km).T
    
    # earth's barycentric vector rel to sun aka center of solsys to earth
    with solar_system_ephemeris.set('de432s'): 
        rE_bary, vE_bary = get_body_barycentric_posvel('earth', obstime)

    # convert to km and km/s arrays in cartesian coordinate system
    rE  = rE_bary.xyz.to_value(units.km).T #rE goes from sun to earth center
    vE  = vE_bary.xyz.to_value(units.km/units.s).T

      #vE goes from sun to earth center

    return obstime, rE, vE, r_obs

def equatorial_to_ecliptic(r_vec, v_vec):
    eps = np.deg2rad(23.439291)  # mean obliquity
    # rotation matrix about x-axis
    R = np.array([
        [1,          0,           0],
        [0,  np.cos(eps), np.sin(eps)],
        [0, -np.sin(eps), np.cos(eps)],
    ])
    r_ecl = r_vec @ R.T
    v_ecl = v_vec @ R.T
    return r_ecl, v_ecl


def to_keplerian_adam(r_helio, v_helio, obstime):
    obstime=obstime.tdb
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
    frame="ecliptic" # ecliptic is the other choice
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
# The Orbits class expects cartesian coordinates but its trivial to convert


def score_from_elements(a_AU, e, H0, i0, grid4D, H_center, a_center, e_center, i_center): #without interpolation
    # checks if orbit is bound right
    if (a_AU <= 0) or (e >= 1.0) or (e < 0):
        
        return 0.0

    # clip to model box
    a_cl = np.clip(a_AU, 0.0, 4.2 - 1e-9)
    e_cl = np.clip(e,    0.0, 1.0 - 1e-9)
    i_cl = np.clip(i0,  0.0, 88.0 - 1e-9)  
    h_cl = np.clip(H0, 15.0, 28.0 - 1e-9)  # added h_cl
    
    #iH = int(np.argmin(np.abs(H_center - H0))) # swapped h0 with h_cl
    iH = int(np.argmin(np.abs(H_center - h_cl)))
    ia = int(np.argmin(np.abs(a_center - a_cl)))
    ie = int(np.argmin(np.abs(e_center - e_cl)))
    ii = int(np.argmin(np.abs(i_center - i_cl)))
    #print("from H0,a,e,i:",H0,a_AU,e,i0)
    #print("from score_from_elements:",iH, ia, ie,ii)

    w=grid4D[iH, ia, ie, ii]
    
    if not np.isfinite(w) or w <= 0:
        return 0.0
    return float(w)
     
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

def compute_h0_from_distance(r_helio, r_geo, r_obs, d_au, mag):
    # Use the formulas from https://adsabs.harvard.edu/full/2007JBAA..117..342D 
    # using G = 0.015 
    # calculating reduced magnitude H(alpha) = mag - 5log(r*delta)
    #dist in AU
    r_AU = np.linalg.norm(r_helio, axis = -1) / AU_km
    delta_AU = np.linalg.norm(r_geo + r_obs, axis = -1) / AU_km
    # vectors from object to Sun and Earth
    u_hat = r_helio / np.linalg.norm(r_helio, axis=-1, keepdims=True)
    v_hat = (r_geo + r_obs) / np.linalg.norm(r_geo + r_obs, axis=-1, keepdims=True)
    
    # phase angle 
    cos_alpha = np.sum(u_hat * v_hat, axis=-1)
    cos_alpha = np.clip(cos_alpha, -1.0, 1.0)
    alpha_rad = np.arccos(cos_alpha)
    alpha_deg = np.degrees(alpha_rad)
    
    # reduced magnitude
    H_alpha = mag - 5.0 * np.log10(r_AU * delta_AU)
    # HG phase correction 
    A1, A2, B1, B2 = 3.33, 1.87, 0.63, 1.22
    phi1 = np.exp(-A1 * np.tan(0.5 * alpha_rad)**B1)
    phi2 = np.exp(-A2 * np.tan(0.5 * alpha_rad)**B2)
    G = 0.15
    phi = (1 - G) * phi1 + G * phi2
    
    # absolute magnitude
    H0 = H_alpha + 2.5 * np.log10(phi)
    
    return H0

def compute_neomod3_weight(d_au, ddot_kms, l_hat, v_hat, mag, obstime,
                           rE, vE, r_obs,
                           array4D, H_center, a_center, e_center,i_center, 
                           ):
    

    r_geo, v_geo = observables_to_geocentric_state(l_hat, v_hat, d_au, ddot_kms)

    r_helio, v_helio = geo_to_topo(r_geo, v_geo, rE, vE, r_obs, obstime) 
    r_ecl, v_ecl = equatorial_to_ecliptic(r_helio, v_helio)



    a_val, e_val, i_val = to_keplerian_adam(r_ecl, v_ecl, obstime)
    
    H0 = compute_h0_from_distance(r_helio, r_geo, r_obs, d_au, mag)
    #H0 = 19.0

    #score = score_from_elements_interp(a_val, e_val, 19.0, i_val, interp4D)
    score = score_from_elements(
        a_val, e_val, H0, i_val,
        array4D, H_center, a_center, e_center, i_center
    )
    #print("score is:",score)
    return score,a_val,e_val,i_val, H0


    
