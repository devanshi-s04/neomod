#!/usr/bin/env python3
"""QA1a — geometry truth round-trip. Feed TRUE (rho, rho_dot) from Kurlander's
per-detection truth into ranging_engine's element math; recover (a,e,i,q) and
compare to the stored true elements. Validates earth state, ecliptic rotation,
unit vectors, element formulas. Pass = <~1% scatter (residual = light-time)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "src")
import ranging_engine as re
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris

D = "/astro/users/jkurla/public_html/LSST_Sorcha_predictions"
df = pd.read_hdf(f"{D}/one_day_neo.h5").reset_index(drop=True)
# truth columns
ra, dec = df.RATrue_deg.values, df.DecTrue_deg.values
dra = df.RARateCosDec_deg_day.values          # already alpha_dot*cos(dec)
ddec = df.DecRate_deg_day.values
rho_km = df.Range_LTC_km.values
rhodot = df.RangeRate_LTC_km_s.values
mjd = df.fieldMJD_TAI.values
a_true, e_true, i_true = df.a.values, df.e.values, df.inc.values

# per-detection sun-centred observer state (vectorised over times)
t = Time(mjd, format="mjd", scale="tai")
r_obs, vE = re.earth_observer_state(t)                                    # sun->observer, sun->earth vel (N,3)

l_hat, v_ang = re.tracklet_vectors(ra, dec, dra, ddec)                    # (N,3),(N,3)[1/s]

# heliocentric state at true (rho, rho_dot)
r_helio = r_obs + rho_km[:, None] * l_hat
v_helio = vE + rhodot[:, None] * l_hat + rho_km[:, None] * v_ang
r_ecl = r_helio @ re._R_EQ2ECL.T
v_ecl = v_helio @ re._R_EQ2ECL.T

r = np.linalg.norm(r_ecl, axis=-1)
v2 = np.einsum('nk,nk->n', v_ecl, v_ecl)
a_km = -re.MU_SUN / (2 * (0.5*v2 - re.MU_SUN/r))
h_vec = np.cross(r_ecl, v_ecl); h = np.linalg.norm(h_vec, axis=-1)
e_vec = np.cross(v_ecl, h_vec)/re.MU_SUN - r_ecl/r[:, None]
e_rec = np.linalg.norm(e_vec, axis=-1)
i_rec = np.degrees(np.arccos(np.clip(h_vec[:, 2]/h, -1, 1)))
a_rec = a_km / re.AU_KM
q_rec = a_rec*(1-e_rec); q_true = a_true*(1-e_true)

def stats(name, rec, tru):
    d = rec - tru
    rel = np.abs(d)/np.maximum(np.abs(tru), 1e-6)
    print(f"  {name:4s} median|rel|={np.nanmedian(rel)*100:.3f}%  "
          f"max|rel|={np.nanpercentile(rel,99)*100:.3f}%(p99)  "
          f"median|abs|={np.nanmedian(np.abs(d)):.4g}")

print(f"QA1a round-trip on {len(df)} Kurlander NEO detections:")
stats("a", a_rec, a_true); stats("e", e_rec, e_true)
stats("i", i_rec, i_true); stats("q", q_rec, q_true)
# bound-orbit + NEO recovery
print(f"  recovered q<1.3 fraction: {np.mean(q_rec<1.3):.3f} (truth: all NEO)")
print(f"  worst a outliers (rec,true):",
      *[f"({a_rec[k]:.2f},{a_true[k]:.2f})" for k in np.argsort(-np.abs(a_rec-a_true))[:3]])
