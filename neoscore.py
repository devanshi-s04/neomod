import numpy as np
from numpy import sin, cos
import pyarrow as pa
from astropy.time import Time
import astropy.units as units
from astropy.coordinates import EarthLocation, get_body_barycentric_posvel, solar_system_ephemeris
from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic
from adam_core.time import Timestamp
from adam_core.coordinates import CartesianCoordinates, Origin
from adam_core.constants import KM_P_AU, S_P_DAY
 

# Constants
AU_km = 149_597_870.7
mu_sun = 1.32712440018e11  # [km^3/s^2] gm of the sun

# Rubin Observatory Location
rubin_location = EarthLocation.from_geodetic( 
        lon=-70.7366*units.deg,   # longitude (west is negative)
        lat=-30.2407*units.deg,   # latitude
        height=2647*units.m       # elevation
        
    )



class NEOMODScorer:
    """
    For computing NEOMOD3 scores from observables.
    """
    
    def __init__(self, array4D, H_center, a_center, e_center, i_center):
        """
        Initialize the scorer with a 4D population model grid.
        
        Parameters
        ----------
        array4D : numpy.ndarray
            4D array containing the population model
        H_center : numpy.ndarray
            Grid centers for H 
        a_center : numpy.ndarray
            Grid centers for semi-major axis 
        e_center : numpy.ndarray
            Grid centers for eccentricity 
        i_center : numpy.ndarray
            Grid centers for inclination 
        """
        self.array4D = array4D
        self.H_center = H_center
        self.a_center = a_center
        self.e_center = e_center
        self.i_center = i_center
    
    def compute_score(self, obstime_str, ra_deg, dec_deg, 
                     dra_deg_per_day, ddec_deg_per_day, mag,
                     d_au, ddot_kms):
        """
        Compute NEOMOD3 score for a single observation with given distance and range rate.
        
        Parameters
        ----------
        obstime_str : str
            Observation time in ISO format
        ra_deg : float
            Right ascension in degrees
        dec_deg : float
            Declination in degrees
        dra_deg_per_day : float
            Change in RA (dRA/dt) in degrees per day
        ddec_deg_per_day : float
            Change in Dec (dDec/dt) in degrees per day
        mag : float
            Apparent visual magnitude
        d_au : float
            Geocentric distance in AU
        ddot_kms : float
            Rate in km/s
        
        Returns
        -------
        dict
            Dictionary containing:
            - 'score': NEOMOD3 weight/score
            - 'a': Semi major axis (AU)
            - 'e': Eccentricity
            - 'i': Inclination (degrees)
            - 'H': Absolute magnitude
        """
        # Compute unit vectors
        l_hat, v_hat = self._compute_unit_vectors(ra_deg, dec_deg, 
                                                   dra_deg_per_day, ddec_deg_per_day)
        
        # Get Earth and observer states
        obstime, rE, vE, r_tele = self._get_earth_and_observer(obstime_str)
        
        # Compute geocentric state
        r_geo, v_geo = self._observables_to_geocentric_state(l_hat, v_hat, d_au, ddot_kms)
        
        # Convert to heliocentric
        r_helio, v_helio = self._geo_to_helio(r_geo, v_geo, rE, vE, r_tele)
        
        # Convert to ecliptic coordinates
        r_ecl, v_ecl = self._equatorial_to_ecliptic(r_helio, v_helio)
        
        # Compute Keplerian elements
        a_val, e_val, i_val = self._to_keplerian_adam(r_ecl, v_ecl, obstime)
        
        # Compute absolute magnitude
        H0 = self._compute_h0_from_distance(r_helio, r_geo, r_tele, d_au, mag)
        
        # Compute score
        score = self._score_from_elements(a_val, e_val, H0, i_val)
        
        return {
            'score': score,
            'a': a_val,
            'e': e_val,
            'i': i_val,
            'H': H0
        }
    #original Nd = 60, original Nv = 120
    def marginalize_over_distance(self, obstime_str, ra_deg, dec_deg,
                                  dra_deg_per_day, ddec_deg_per_day, mag,
                                  d_min=0.0, d_max=0.4, Nd=50,
                                  v_min=-40.0, v_max=40.0, Nv=100):
        """
        Compute marginalized score by integrating over distance and range rate.
        """
        d_grid = np.linspace(d_min, d_max, Nd)
        ddot_grid = np.linspace(v_min, v_max, Nv)
        
        # Pre-compute states that don't depend on d, ddot (ONLY ONCE!)
        l_hat, v_hat = self._compute_unit_vectors(ra_deg, dec_deg,
                                                   dra_deg_per_day, ddec_deg_per_day)
        obstime, rE, vE, r_tele = self._get_earth_and_observer(obstime_str)
        
        w_sum = 0.0
        n = 0
        for d_au in d_grid:
            for ddot_kms in ddot_grid:
                # without calling compute_score()
                r_geo, v_geo = self._observables_to_geocentric_state(l_hat, v_hat, d_au, ddot_kms)
                r_helio, v_helio = self._geo_to_helio(r_geo, v_geo, rE, vE, r_tele)
                r_ecl, v_ecl = self._equatorial_to_ecliptic(r_helio, v_helio)
                a_val, e_val, i_val = self._to_keplerian_adam(r_ecl, v_ecl, obstime)
                H0 = self._compute_h0_from_distance(r_helio, r_geo, r_tele, d_au, mag)
                score = self._score_from_elements(a_val, e_val, H0, i_val)
                w_sum += score
                n += 1
        
        return w_sum / max(n,1)
    
    
    
    def _compute_unit_vectors(self, ra_deg, dec_deg, dra_deg_per_day, ddec_deg_per_day):
        """Compute line of sight and velocity unit vectors."""
        ra = np.deg2rad(ra_deg)
        dec = np.deg2rad(dec_deg)
        
        l_hat, e_ra, e_dec = self._sky_basis(ra, dec)
        
        dra = np.deg2rad(dra_deg_per_day) / 86400.0
        ddec = np.deg2rad(ddec_deg_per_day) / 86400.0
        
        v_hat = (dra * e_ra + ddec * e_dec)
        
        return l_hat, v_hat
    
    def _sky_basis(self, ra_rad, dec_rad):
        """Create right handed basis pointing to object from observer."""
        l_hat = np.array([np.cos(dec_rad)*np.cos(ra_rad),
                         np.cos(dec_rad)*np.sin(ra_rad),
                         np.sin(dec_rad)])
        e_ra = np.array([-np.sin(ra_rad), np.cos(ra_rad), 0.0])
        e_dec = np.array([-np.sin(dec_rad)*np.cos(ra_rad),
                         -np.sin(dec_rad)*np.sin(ra_rad),
                         np.cos(dec_rad)])
        return l_hat, e_ra, e_dec
    
    def _get_earth_and_observer(self, obstime_str):
        """Get Earth's barycentric state and observer topocentric offset."""
        obstime = Time(obstime_str)
        
        obs_geo = rubin_location.get_gcrs(obstime)
        r_tele = obs_geo.cartesian.xyz.to_value(units.km).T
        
        with solar_system_ephemeris.set('de432s'):
            rE_bary, vE_bary = get_body_barycentric_posvel('earth', obstime)
        
        rE = rE_bary.xyz.to_value(units.km).T
        vE = vE_bary.xyz.to_value(units.km/units.s).T
        
        return obstime, rE, vE, r_tele
    
    def _observables_to_geocentric_state(self, l_hat, v_hat, d_au, ddot_kms):
        """Convert observables to geocentric position and velocity."""
        d_km = d_au * AU_km
        
        r_geo = np.outer(d_km, l_hat)
        v_rad = np.outer(ddot_kms, l_hat)
        v_tan = np.outer(d_km, v_hat)
        v_geo = v_rad + v_tan
        
        return r_geo, v_geo
    
    def _geo_to_helio(self, r_geo_km, v_geo_kms, rE, vE, r_tele):
        """Convert geocentric to heliocentric vectors."""
        r_helio = rE + r_tele + r_geo_km
        v_helio = vE + v_geo_kms
        return r_helio, v_helio
    
    def _equatorial_to_ecliptic(self, r_vec, v_vec):
        """Rotate from equatorial to ecliptic coordinates."""
        eps = np.deg2rad(23.439291)
        R = np.array([
            [1, 0, 0],
            [0, np.cos(eps), np.sin(eps)],
            [0, -np.sin(eps), np.cos(eps)],
        ])
        r_ecl = r_vec @ R.T
        v_ecl = v_vec @ R.T
        return r_ecl, v_ecl
    
    def _to_keplerian_adam(self, r_helio, v_helio, obstime):
        """Convert Cartesian to Keplerian elements using ADAM."""
        obstime = obstime.tdb
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
            frame="ecliptic"
        )
        
        keplerian_coordinates = cartesian_coordinates.to_keplerian()
        
        a_val = keplerian_coordinates.a[0].as_py()
        e_val = keplerian_coordinates.e[0].as_py()
        i_val = keplerian_coordinates.i[0].as_py()
        
        return a_val, e_val, i_val
    
    def _compute_h0_from_distance(self, r_helio, r_geo, r_tele, d_au, mag):
        """Compute absolute magnitude H from apparent magnitude."""
        # Use the formulas from https://adsabs.harvard.edu/full/2007JBAA..117..342D
        # using G = 0.015 
        # calculating reduced magnitude H(alpha) = mag - 5log(r*delta)
        
        r_AU = np.linalg.norm(r_helio, axis=-1) / AU_km
        delta_AU = np.linalg.norm(r_geo + r_tele, axis=-1) / AU_km
        
        u_hat = r_helio / np.linalg.norm(r_helio, axis=-1, keepdims=True)
        v_hat = (r_geo + r_tele) / np.linalg.norm(r_geo + r_tele, axis=-1, keepdims=True)
        
        cos_alpha = np.sum(u_hat * v_hat, axis=-1)
        cos_alpha = np.clip(cos_alpha, -1.0, 1.0)
        alpha_rad = np.arccos(cos_alpha)
        
        H_alpha = mag - 5.0 * np.log10(r_AU * delta_AU)
        
        # HG phase correction
        A1, A2, B1, B2 = 3.33, 1.87, 0.63, 1.22
        phi1 = np.exp(-A1 * np.tan(0.5 * alpha_rad)**B1)
        phi2 = np.exp(-A2 * np.tan(0.5 * alpha_rad)**B2)
        G = 0.15
        phi = (1 - G) * phi1 + G * phi2
        
        H0 = H_alpha + 2.5 * np.log10(phi)

        
        if np.ndim(H0) == 0:
            return float(H0)
        elif H0.size == 1:
            return float(H0.item())
             
        return H0
    
    def _score_from_elements(self, a_AU, e, H0, i0):
        """Look up score from orbital elements in the 4D grid."""
        if (a_AU <= 0) or (e >= 1.0) or (e < 0):
            return 0.0
        
        a_cl = np.clip(a_AU, 0.0, 4.2 - 1e-9)
        e_cl = np.clip(e, 0.0, 1.0 - 1e-9)
        i_cl = np.clip(i0, 0.0, 88.0 - 1e-9)
        h_cl = np.clip(H0, 15.0, 28.0 - 1e-9)
        
        iH = int(np.argmin(np.abs(self.H_center - h_cl)))
        ia = int(np.argmin(np.abs(self.a_center - a_cl)))
        ie = int(np.argmin(np.abs(self.e_center - e_cl)))
        ii = int(np.argmin(np.abs(self.i_center - i_cl)))
        
        w = self.array4D[iH, ia, ie, ii]
        
        if not np.isfinite(w) or w <= 0:
            return 0.0
        return float(w)



def radec_rates_to_ecliptic_rates(ra_deg, dec_deg, dra_deg_day, ddec_deg_day):
    """Convert RA/Dec rates to ecliptic longitude/latitude rates."""
    sc_eq = SkyCoord(
        ra=ra_deg * units.deg,
        dec=dec_deg * units.deg,
        pm_ra_cosdec=dra_deg_day * units.deg / units.day,
        pm_dec=ddec_deg_day * units.deg / units.day,
        frame="icrs"
    )
    
    sc_ecl = sc_eq.transform_to(GeocentricTrueEcliptic())
    
    vlam = sc_ecl.pm_lon_coslat.to(units.deg / units.day).value
    vbeta = sc_ecl.pm_lat.to(units.deg / units.day).value
    
    return vlam, vbeta


def ecliptic_rates_to_radec_rates(ra_deg, dec_deg, vlam_deg_day, vbeta_deg_day, obstime_str=None):
    """Convert ecliptic longitude/latitude rates to RA/Dec rates."""
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
            lon=ra_deg * units.deg,
            lat=dec_deg * units.deg,
            pm_lon_coslat=vlam_deg_day * units.deg/units.day,
            pm_lat=vbeta_deg_day * units.deg/units.day,
            frame="geocentrictrueecliptic"
        )
    
    sc_eq = sc_ecl.transform_to("icrs")
    
    dra = sc_eq.pm_ra_cosdec.to(units.deg/units.day).value
    ddec = sc_eq.pm_dec.to(units.deg/units.day).value
    
    return dra, ddec
# added after the error that messed up sides
def ecliptic_rates_to_radec_rates1(ra_deg, dec_deg, vlam_deg_day, vbeta_deg_day, obstime_str):
    """
    Convert ecliptic rates (pm_lon_coslat, pm_lat) at the sky position (ra,dec)
    into ICRS rates (pm_ra_cosdec, pm_dec).

    IMPORTANT: vlam_deg_day here is pm_lon_coslat (dλ/dt * cosβ), consistent with astropy.
    """
    t = Time(obstime_str)

    # 1) True sky position in ICRS
    sc_icrs = SkyCoord(ra=ra_deg*units.deg, dec=dec_deg*units.deg, frame="icrs", obstime=t)

    # 2) Transform that position to geocentric true ecliptic to get lon/lat
    sc_ecl_pos = sc_icrs.transform_to(GeocentricTrueEcliptic(obstime=t))

    # 3) Rebuild ecliptic coord at same lon/lat, now with ecliptic proper motions
    sc_ecl = SkyCoord(
        lon=sc_ecl_pos.lon,
        lat=sc_ecl_pos.lat,
        pm_lon_coslat=vlam_deg_day * units.deg/units.day,
        pm_lat=vbeta_deg_day * units.deg/units.day,
        frame=GeocentricTrueEcliptic(obstime=t)
    )

    # 4) Transform back to ICRS and extract RA/Dec rates
    sc_back = sc_ecl.transform_to("icrs")

    dra  = sc_back.pm_ra_cosdec.to(units.deg/units.day).value
    ddec = sc_back.pm_dec.to(units.deg/units.day).value
    return dra, ddec