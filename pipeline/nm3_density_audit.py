#!/usr/bin/env python3
"""Where does the 24% NEO under-count (Sum p / true = 0.763) come from? (docs §9.8)

density_raw is downweighted, so integral(rho) dA over the whole velocity grid = the expected number
of objects of that population visible in the patch, in that magnitude bin. That is checkable against
the actual objects, with NO truth labels and no simulation-transfer assumption.

Three numbers per (center, magbin, population):
  A  int(rho) dA              -- what the MAP says
  B  clones x weight          -- what the clones IN the map imply  (A/B tests the estimator+mask)
  C  true objects in patch    -- counted from the Stage-0 n-body epoch cache (A/C tests reality)

MBA IS THE CONTROL: its density is built from real S3M objects, so A/C must be ~1 for MBA. If it is
not, the method is wrong rather than the map. Only if MBA passes does a NEO deviation mean anything.
"""
import json, sys
import numpy as np, pandas as pd
import astropy.units as u
from astropy.time import Time
from astropy.utils import iers; iers.conf.auto_max_age = None
from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic, get_sun
W = "/mmfs1/gscratch/dirac/ds2004/sorcha"
pd.set_option("display.width", 250)
EPOCH = "2027-08-25T00:00:00"; MAXSEP = 30.0
BINS = {"14_16": (14, 16), "16_18": (16, 18), "18_20": (18, 20), "mag20": (20, 21),
        "mag21": (21, 22), "mag22": (22, 23), "mag23": (23, 24), "mag24+": (24, 25)}
CENTERS = [("dlon+000_lat+00", 0.0, 0.0), ("dlon+020_lat-12", 20.0, -12.0),
           ("dlon+050_lat-01", 50.0, -1.0), ("dlon-070_lat-25", -70.0, -25.0)]
w_abs = json.load(open(f"{W}/outputs/neomod3_projection_cache/cache_metadata.json"))["w_abs_objects_per_clone"]
print(f"w_abs = {w_abs:.5f} objects/clone")

t = Time(EPOCH, scale="utc")
sun = get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
cache = pd.read_parquet(f"{W}/outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet",
                        columns=["population", "lam_deg", "beta_deg", "mag_app", "vlam", "vbeta"])
print(f"epoch cache: {len(cache):,} objects")
obj = SkyCoord(lon=cache.lam_deg.to_numpy()*u.deg, lat=cache.beta_deg.to_numpy()*u.deg,
               frame=GeocentricTrueEcliptic(obstime=t))

rows = []
for label, dlon, lat in CENTERS:
    clon = (sun.lon.deg + 180.0 + dlon) % 360.0
    ctr = SkyCoord(lon=clon*u.deg, lat=lat*u.deg, frame=GeocentricTrueEcliptic(obstime=t))
    sep = obj.separation(ctr).deg
    inpatch = sep <= MAXSEP
    z = np.load(f"{W}/prob_maps_grid_neomod3_full/prob_maps_grid_{label}.npz", allow_pickle=True)
    x, y = z["x_grid"], z["y_grid"]; dA = float((x[1]-x[0])*(y[1]-y[0]))
    for b, (m0, m1) in BINS.items():
        for pop, w in (("NEO", w_abs), ("MBA", 1.0)):
            A = np.nan_to_num(np.asarray(z[f"density_raw__{pop}__{b}"], float)).sum()*dA
            B = np.nan_to_num(np.asarray(z[f"support_count__{pop}__{b}"], float)).sum()*w
            sel = inpatch & (cache.population.to_numpy() == pop) & \
                  (cache.mag_app.to_numpy() >= m0) & (cache.mag_app.to_numpy() < m1)
            C = int(sel.sum())
            v = np.maximum(np.abs(cache.vlam.to_numpy()[sel]), np.abs(cache.vbeta.to_numpy()[sel]))
            rows.append(dict(center=label, magbin=b, pop=pop, int_rho=A, clones_x_w=B, true_count=C,
                             frac_v_gt5=float((v > 5).mean()) if C else np.nan,
                             A_over_B=A/B if B > 0 else np.nan, A_over_C=A/C if C > 0 else np.nan))
r = pd.DataFrame(rows); r.to_csv(f"{W}/outputs/neomod3_fullgrid/density_audit.csv", index=False)

for pop in ["MBA", "NEO"]:
    s = r[r["pop"] == pop]
    print(f"\n{'='*118}\n{pop}" + ("   <-- CONTROL: A/C must be ~1" if pop == "MBA" else
          "   <-- the question: does the map predict the right number of NEOs?"))
    piv = s.pivot_table(index="magbin", values=["int_rho", "clones_x_w", "true_count", "A_over_B", "A_over_C"],
                        aggfunc="mean").reindex(list(BINS))
    print(piv[["int_rho", "clones_x_w", "true_count", "A_over_B", "A_over_C"]].to_string(
        float_format=lambda v: f"{v:,.4g}"))
    tot_A = s.int_rho.sum(); tot_C = s.true_count.sum()
    print(f"  TOTAL over all bins/centers: int(rho)={tot_A:,.0f}  true={tot_C:,}  ratio={tot_A/tot_C:.3f}")
print("\nfrac_v_gt5 = fraction of true objects moving faster than the ±5 grid edge (density lost off-grid)")
print(r.groupby("pop", observed=True).frac_v_gt5.mean().to_string())
