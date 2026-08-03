#!/usr/bin/env python3
"""(1) Where does the uniform ~13% density overshoot come from? (docs §9.9)
   (2) Is the faint-end NEO excess entirely the H>25 tail S3M does not have?

(1) The map density is a Bayesian posterior MEAN of n0 = 1/(pi*d0^2). E[1/d0^2] != 1/E[d0]^2
    (Jensen), so a posterior mean of an inverse square is biased HIGH by construction. Measured here
    against a field of KNOWN uniform density -- if the ratio reproduces ~1.13 the bias is the
    estimator, is analytic, and is fixable by one constant.
(2) Restrict the NEOMOD3 clones to H<25 (S3M's range) and re-compare to the S3M count. If the ratio
    collapses to ~1, the faint excess is entirely the H>25 population and nothing else is wrong.
"""
import sys, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, "/mmfs1/gscratch/dirac/ds2004/sorcha/neomod/src")
from scipy.spatial import cKDTree
import velocity_density_pipeline_neomod_clone_only as v
W = "/mmfs1/gscratch/dirac/ds2004/sorcha"
pd.set_option("display.width", 200)

print("="*80); print("(1) ESTIMATOR BIAS on a field of KNOWN density")
rng = np.random.default_rng(0)
rows = []
for k in [5, 10, 20, 40]:
    for rho_true in [200.0, 1000.0]:
        L = 10.0; N = int(rho_true*L*L)
        pts = rng.uniform(0, L, size=(N, 2))
        tree = cKDTree(pts)
        ev = rng.uniform(3, 7, size=(300, 2))       # far from edges
        est = np.array([v.estimate_density_full_posterior_2d(tree, x, y, k=k, n_d0_grid=400)
                        for x, y in ev])
        rows.append(dict(k=k, rho_true=rho_true, mean_est=est.mean(),
                         ratio=est.mean()/rho_true, k_over_km1=k/(k-1)))
b = pd.DataFrame(rows)
print(b.to_string(index=False, float_format=lambda x: f"{x:,.4f}"))
print("\n  ratio = measured/true. Compare with k/(k-1), the classic kNN normalisation bias.")

print("\n" + "="*80); print("(2) IS THE FAINT NEO EXCESS JUST H>25?")
import astropy.units as u
from astropy.time import Time
from astropy.utils import iers; iers.conf.auto_max_age = None
from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic, get_sun
import healpy as hp, pyarrow.dataset as pads, pyarrow.compute as pc
EPOCH = "2027-08-25T00:00:00"; MAXSEP = 30.0
BINS = {"mag21": (21, 22), "mag22": (22, 23), "mag23": (23, 24), "mag24+": (24, 25)}
w_abs = json.load(open(f"{W}/outputs/neomod3_projection_cache/cache_metadata.json"))["w_abs_objects_per_clone"]
t = Time(EPOCH, scale="utc")
sun = get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
out = []
for label, dlon, lat in [("dlon+000_lat+00", 0.0, 0.0), ("dlon+020_lat-12", 20.0, -12.0)]:
    clon = (sun.lon.deg + 180.0 + dlon) % 360.0
    ctr_ecl = SkyCoord(lon=clon*u.deg, lat=lat*u.deg, frame=GeocentricTrueEcliptic(obstime=t))
    ctr_icrs = ctr_ecl.transform_to("icrs")
    px = hp.query_disc(8, hp.ang2vec(ctr_icrs.ra.deg, ctr_icrs.dec.deg, lonlat=True),
                       np.radians(MAXSEP + np.degrees(hp.max_pixrad(8))), inclusive=True)
    ds = pads.dataset(f"{W}/outputs/neomod3_projection_cache/by_pixel", format="parquet", partitioning="hive")
    nm = ds.to_table(filter=pc.field("pix").isin(px.tolist()),
                     columns=["ra_deg", "dec_deg", "mag_app", "H"]).to_pandas()
    c = SkyCoord(ra=nm.ra_deg.to_numpy()*u.deg, dec=nm.dec_deg.to_numpy()*u.deg, frame="icrs")
    nm = nm[c.separation(ctr_icrs).deg <= MAXSEP]
    cache = pd.read_parquet(f"{W}/outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet",
                            columns=["population", "lam_deg", "beta_deg", "mag_app"])
    cache = cache[cache.population == "NEO"]
    co = SkyCoord(lon=cache.lam_deg.to_numpy()*u.deg, lat=cache.beta_deg.to_numpy()*u.deg,
                  frame=GeocentricTrueEcliptic(obstime=t))
    cache = cache[co.separation(ctr_ecl).deg <= MAXSEP]
    for bl, (m0, m1) in BINS.items():
        n_all = ((nm.mag_app >= m0) & (nm.mag_app < m1)).sum()*w_abs
        n_h25 = ((nm.mag_app >= m0) & (nm.mag_app < m1) & (nm.H < 25)).sum()*w_abs
        s3m = int(((cache.mag_app >= m0) & (cache.mag_app < m1)).sum())
        out.append(dict(center=label[:15], magbin=bl, NEOMOD3_all=n_all, NEOMOD3_Hlt25=n_h25,
                        S3M=s3m, ratio_all=n_all/max(s3m, 1), ratio_Hlt25=n_h25/max(s3m, 1)))
o = pd.DataFrame(out)
print(o.to_string(index=False, float_format=lambda x: f"{x:,.3f}"))
print("\n  ratio_all    = NEOMOD3 (H 15-28) / S3M   -- the 2-4x from §9.9")
print("  ratio_Hlt25  = NEOMOD3 restricted to S3M's own H range / S3M")
print("                 ~1 => the excess is ENTIRELY the H>25 population S3M lacks.")
