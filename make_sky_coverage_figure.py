#!/usr/bin/env python
"""
Mollweide sky-density map of Sorcha tracklets ("base LSST survey" figure),
in the style of Wagg et al. (2025) and Kurlander et al. (LSST yield).

Preview version uses the Arnor comparison subsample
(outputs/phase2/sorcha_comparison.parquet, v3.3 window, antisun-footprint,
non-NEO tracklets for a faithful footprint). For the publication figure,
point --input at the full Hyak v5.0 tracklet set (outputs/tracklets_v5/...)
once Hyak is available, and drop --drop-neo (the full set is unbiased).

Env:  /astro/users/ds2004/.conda/envs/neofast_py310/bin/python
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import healpy as hp
from astropy.coordinates import SkyCoord, BarycentricMeanEcliptic, ICRS
import astropy.units as u


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="outputs/phase2/sorcha_comparison.parquet")
    p.add_argument("--output", default="Figures/sorcha_sky_coverage.png")
    p.add_argument("--ra-col", default="mean_ra")
    p.add_argument("--dec-col", default="mean_dec")
    p.add_argument("--nside", type=int, default=64)
    p.add_argument("--drop-neo", action="store_true",
                   help="Exclude NEOs (oversampled in the comparison subsample).")
    p.add_argument("--smooth-deg", type=float, default=0.0,
                   help="Optional Gaussian smoothing FWHM in degrees.")
    p.add_argument("--title", default="Sorcha simulated tracklet density (2-yr, Rubin cadence)")
    args = p.parse_args()

    cols = [args.ra_col, args.dec_col]
    has_pop = False
    try:
        df = pd.read_parquet(args.input, columns=cols + ["population"])
        has_pop = True
    except Exception:
        df = pd.read_parquet(args.input, columns=cols)
    if args.drop_neo and has_pop:
        df = df[df.population != "NEO"]

    ra = df[args.ra_col].to_numpy()
    dec = df[args.dec_col].to_numpy()

    npix = hp.nside2npix(args.nside)
    pix_area = hp.nside2pixarea(args.nside, degrees=True)
    pix = hp.ang2pix(args.nside, ra, dec, lonlat=True)
    counts = np.bincount(pix, minlength=npix).astype(float)
    dens = counts / pix_area
    if args.smooth_deg > 0:
        dens = hp.smoothing(np.nan_to_num(dens), fwhm=np.radians(args.smooth_deg))
    dens[counts == 0] = np.nan

    # ecliptic plane (lat = 0) in equatorial RA/Dec, sorted by RA and broken at the
    # map edge (RA=180) so the line does not draw a spurious connector across the map.
    elon = np.linspace(0, 360, 2000)
    ecl = SkyCoord(lon=elon * u.deg, lat=np.zeros_like(elon) * u.deg,
                   frame=BarycentricMeanEcliptic).transform_to(ICRS())
    e_ra, e_dec = ecl.ra.deg, ecl.dec.deg
    o = np.argsort(e_ra)
    e_ra, e_dec = e_ra[o], e_dec[o]
    ecl_lon180 = ((e_ra + 180) % 360) - 180          # longitude in [-180, 180]
    brk = np.where(np.abs(np.diff(ecl_lon180)) > 180)[0]
    ecl_lon180 = np.insert(ecl_lon180, brk + 1, np.nan)
    ecl_dec = np.insert(e_dec, brk + 1, np.nan)

    vmin = max(np.nanpercentile(dens, 2), 1e-1)
    vmax = np.nanpercentile(dens, 99.5)
    plt.figure(figsize=(11, 6))
    hp.projview(dens, coord=["C"], projection_type="mollweide", flip="astro",
                cmap="magma", norm="log", min=vmin, max=vmax,
                graticule=True, graticule_labels=True,
                xlabel="Right Ascension", ylabel="Declination",
                unit="tracklets per square degree", cb_orientation="horizontal",
                title=args.title, badcolor="white", bgcolor="white",
                xtick_label_color="0.3", ytick_label_color="0.3")
    hp.newprojplot(theta=np.radians(90 - ecl_dec), phi=np.radians(ecl_lon180),
                   color="deepskyblue", lw=1.8)
    plt.savefig(args.output, dpi=200, bbox_inches="tight", facecolor="white")
    print("saved", args.output, "| tracklets:", len(ra),
          "| density/sqdeg: %.2f-%.1f" % (np.nanmin(dens), np.nanmax(dens)))


if __name__ == "__main__":
    main()
