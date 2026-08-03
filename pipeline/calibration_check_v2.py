#!/usr/bin/env python3
"""Calibration check v2 — three-way comparison, with the ground truth computed directly.

v1 was wrong twice and this fixes both mistakes:
  * `magcut_count` in the npz is the magnitude cut applied to the WHOLE-SKY population; the 30 deg sky
    cut happens afterwards. Comparing int(rho) (a patch quantity) to it (an all-sky quantity) is
    meaningless, and the apparent "NEO recovers 3x less than MBA" was just MBAs being far more
    concentrated toward the ecliptic/opposition than the more isotropic NEOs.
  * the instantaneous-vs-cumulative ("NEOs cycle through the patch") hypothesis was tested by
    restricting to a single night and did NOT explain the gap.

So compute the ground truth directly from the Stage 0 n-body cache: how many objects of each
population are ACTUALLY visible in this patch, in this magnitude bin, at the map epoch.

Three quantities per (center, mag bin):
  1. TRUE VISIBLE  -- count per population from the Stage 0 cache, same 30 deg cut, same mag bin.
                      This is what the map SHOULD encode.
  2. MAP IMPLIED   -- int rho_pop dv from the density maps. What the map DOES encode.
  3. DETECTED      -- Sorcha tracklets in that map + V-band mag bin.

  (1) vs (2) tests the cross-population NORMALISATION -- the actual open question.
  (1) vs (3) tests DETECTION EFFICIENCY (trailing losses, cadence) and is a physical result, not a bug.
"""
from __future__ import annotations
import argparse, re
from pathlib import Path
import numpy as np, pandas as pd
from astropy.time import Time
from astropy.coordinates import SkyCoord, GCRS, GeocentricTrueEcliptic, get_sun
import astropy.units as u

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
CACHE = W / "outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet"
SORCHA = W / "outputs/s3m_linking/case1/sorcha_comparison_case1_nbody_Vband.parquet"
OUT = W / "outputs" / "calibration_check"
EPOCH = "2027-08-25T00:00:00"
MAX_SEP = 30.0
POP_MAP = {"NEO": "NEO", "MBA": "MBA", "TNO": "TNO", "Trojan": "Trojans"}


def center_from_name(name, t_obs):
    dlon = float(re.search(r"dlon([+-]\d+)", name).group(1))
    lat = float(re.search(r"lat([+-]\d+)", name).group(1))
    sun = get_sun(t_obs).transform_to(GeocentricTrueEcliptic(obstime=t_obs))
    anti = (sun.lon.deg + 180.0) % 360.0
    return (anti + dlon) % 360.0, lat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maps-dir", default="prob_maps_grid_s3m_nbody")
    ap.add_argument("--center", required=True)
    ap.add_argument("--tag", default="production")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    t_obs = Time(EPOCH, scale="tdb")
    clon, clat = center_from_name(a.center, t_obs)
    print(f"center {a.center} -> ecliptic ({clon:.3f}, {clat:.1f}), {MAX_SEP} deg cut", flush=True)

    # ---------- 1. TRUE VISIBLE, from the Stage 0 cache ----------
    cache = pd.read_parquet(CACHE, columns=["population", "ra_deg", "dec_deg", "mag_app"])
    print(f"cache: {len(cache):,} objects", flush=True)
    center_sc = SkyCoord(lon=clon*u.deg, lat=clat*u.deg, distance=1.0*u.AU,
                         frame=GeocentricTrueEcliptic(obstime=t_obs)).transform_to(GCRS(obstime=t_obs))
    sc = SkyCoord(ra=cache.ra_deg.to_numpy()*u.deg, dec=cache.dec_deg.to_numpy()*u.deg,
                  frame=GCRS(obstime=t_obs))
    inpatch = sc.separation(center_sc).deg <= MAX_SEP
    vis = cache[inpatch]
    print(f"visible in patch: {len(vis):,}  ({100*len(vis)/len(cache):.2f}% of sky population)", flush=True)

    # ---------- 2. MAP IMPLIED ----------
    z = np.load(W / a.maps_dir / a.center, allow_pickle=True)
    xg, yg = z["x_grid"], z["y_grid"]
    cell = float((xg[1]-xg[0])*(yg[1]-yg[0]))
    labels = [str(l) for l in z["mag_bin_labels"]]
    mins = [float(v) for v in z["mag_bin_mins"]]
    maxs = [float(v) for v in z["mag_bin_maxs"]]
    pops = [str(p) for p in z["population_names"]]

    # ---------- 3. DETECTED ----------
    trk = pd.read_parquet(SORCHA, columns=["population", "prob_map_file", "night", "mag_bin_label_Vband"])
    trk = trk[trk.prob_map_file == a.center]

    rows = []
    for lab, mlo, mhi in zip(labels, mins, maxs):
        vb = vis[(vis.mag_app >= mlo) & (vis.mag_app < mhi)]
        true_n = {p: int((vb.population == p).sum()) for p in pops}
        map_n = {p: float(z[f"density_raw__{p}__{lab}"].sum())*cell for p in pops}
        det = trk[trk.mag_bin_label_Vband == lab]
        det_n = {p: int((det.population.map(POP_MAP) == p).sum()) for p in pops}

        def frac(d):
            t = sum(d.values())
            return (d.get("NEO", 0)/t) if t > 0 else np.nan
        rows.append(dict(
            mag_bin=lab,
            true_NEO=true_n.get("NEO", 0), true_MBA=true_n.get("MBA", 0), true_NEO_frac=frac(true_n),
            map_NEO=map_n.get("NEO", 0.0), map_MBA=map_n.get("MBA", 0.0), map_NEO_frac=frac(map_n),
            det_NEO=det_n.get("NEO", 0), det_MBA=det_n.get("MBA", 0), det_NEO_frac=frac(det_n),
        ))
    df = pd.DataFrame(rows)
    df["map_over_true"] = df.map_NEO_frac / df.true_NEO_frac
    df["det_over_true"] = df.det_NEO_frac / df.true_NEO_frac
    df.insert(0, "center", a.center)
    df.to_csv(OUT / f"calibration_v2_{a.tag}.csv", index=False)

    pd.set_option("display.width", 250)
    print(f"\n{'='*110}")
    print("NEO fraction:  TRUE visible (cache)  vs  MAP implied (int rho)  vs  DETECTED (Sorcha 2yr)")
    print(df[["mag_bin", "true_NEO", "true_MBA", "true_NEO_frac",
              "map_NEO_frac", "map_over_true", "det_NEO", "det_MBA", "det_NEO_frac", "det_over_true"]]
          .to_string(index=False, float_format=lambda v: f"{v:,.4g}"))
    ok = df[np.isfinite(df.map_over_true)]
    if len(ok):
        print(f"\n[NORMALISATION]  map/true NEO fraction: median {ok.map_over_true.median():.3f}x   "
              f"range {ok.map_over_true.min():.3f}-{ok.map_over_true.max():.3f}")
        print("   1.0 = correctly normalised.  <1 = NEO under-weighted in the maps.")
    ok2 = df[np.isfinite(df.det_over_true)]
    if len(ok2):
        print(f"[DETECTION]      det/true NEO fraction: median {ok2.det_over_true.median():.3f}x   "
              f"(physical: cadence/trailing/linking, NOT a normalisation error)")
    print(f"\nwrote {OUT/f'calibration_v2_{a.tag}.csv'}")


if __name__ == "__main__":
    main()
