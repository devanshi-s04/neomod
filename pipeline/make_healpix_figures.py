#!/usr/bin/env python3
"""Figures explaining the HEALPix-partitioned NEOMOD3 cache, for neomod3_vs_s3m_findings.ipynb.

Uses the SAME query_disc call as velocity_density_pipeline_neomod_clone_only._load_neomod3_cache_healpix,
so the pictures show what the pipeline actually reads, not an idealisation.
"""
import json, glob, os
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import healpy as hp
import pyarrow.parquet as pq
import astropy.units as u
from astropy.time import Time
from astropy.utils import iers
iers.conf.auto_max_age = None
from astropy.coordinates import SkyCoord, GCRS, GeocentricTrueEcliptic, get_sun

W = Path("/mmfs1/gscratch/dirac/ds2004/sorcha")
OUT = W/"outputs/neomod3_vs_s3m_comparison"; OUT.mkdir(parents=True, exist_ok=True)
BYPIX = W/"outputs/neomod3_projection_cache/by_pixel"
META = json.load(open(BYPIX/"_healpix_meta.json"))
NSIDE, NPIX = META["nside"], META["npix"]
EPOCH = "2027-08-25T00:00:00"; MAXSEP = 30.0
MARGIN = np.degrees(hp.max_pixrad(NSIDE))

# ---- per-pixel row counts, straight from parquet footers (no data read) ----
counts = np.zeros(NPIX, dtype=np.int64)
for d in sorted(glob.glob(str(BYPIX/"pix=*"))):
    p = int(os.path.basename(d).split("=")[1])
    counts[p] = sum(pq.ParquetFile(f).metadata.num_rows for f in glob.glob(f"{d}/*.parquet"))
sizes = np.array([sum(os.path.getsize(f) for f in glob.glob(f"{BYPIX}/pix={p}/*.parquet"))
                  for p in range(NPIX)], dtype=float)
print(f"rows {counts.sum():,}  (meta says {META['n_rows']:,})   total {sizes.sum()/2**30:.2f} GiB")

t = Time(EPOCH, scale="tdb")
_SUNLON = get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t)).lon.deg   # ONCE, not per center

def centers_radec(dlon, lat):
    """Vectorised: arrays of delta_lon/lat -> arrays of (ra, dec). One frame transform for all."""
    dlon = np.atleast_1d(np.asarray(dlon, float)); lat = np.atleast_1d(np.asarray(lat, float))
    lon = (_SUNLON + 180.0 + dlon) % 360.0
    c = SkyCoord(lon=lon*u.deg, lat=lat*u.deg, distance=1.0*u.AU,
                 frame=GeocentricTrueEcliptic(obstime=t)).transform_to(GCRS(obstime=t))
    return np.atleast_1d(c.ra.deg), np.atleast_1d(c.dec.deg)

def pixels_for(dlon, lat, _t=None):
    ra, dec = centers_radec(dlon, lat); ra, dec = float(ra[0]), float(dec[0])
    vec = hp.ang2vec(ra, dec, lonlat=True)
    return hp.query_disc(NSIDE, vec, np.radians(MAXSEP + MARGIN), inclusive=True), ra, dec

# =============== FIGURE 1: what the partition looks like, and what one task reads ===============
fig = plt.figure(figsize=(14, 9))
dens = np.where(counts > 0, counts, np.nan)
plt.axes([0.02, 0.55, 0.96, 0.40])
hp.mollview(np.log10(dens), fig=fig.number, hold=True, cmap="viridis",
            title=f"(a)  NEOMOD3 cache partitioned into {NPIX} HEALPix pixels (nside={NSIDE})\n"
                  f"colour = log$_{{10}}$(clones in pixel);  {counts.sum():,} clones total",
            unit="log10(clones per pixel)", cbar=True)
hp.graticule(dpar=30, dmer=30, color="white", alpha=0.35, verbose=False)

sel, ra, dec = pixels_for(20.0, -12.0, t)
mask = np.full(NPIX, np.nan)
mask[counts > 0] = 0.0
mask[sel] = 1.0
plt.axes([0.02, 0.06, 0.96, 0.40])
hp.mollview(mask, fig=fig.number, hold=True, cmap=ListedColormap(["0.82", "crimson"]),
            min=0, max=1, cbar=False,
            title=f"(b)  what ONE map task actually reads — center dlon+020 lat−12\n"
                  f"red = the {len(sel)} pixels touching its {MAXSEP:.0f}° patch "
                  f"({100*len(sel)/NPIX:.0f}% of the sky tiles);  grey = the other {NPIX-len(sel)} pixels, never opened")
hp.graticule(dpar=30, dmer=30, color="0.4", alpha=0.35, verbose=False)
for r, lab, c in [(MAXSEP, f"{MAXSEP:.0f}° patch", "black"),
                  (MAXSEP+MARGIN, f"+{MARGIN:.2f}° max_pixrad margin", "lime")]:
    th = np.linspace(0, 2*np.pi, 400)
    v = hp.ang2vec(ra, dec, lonlat=True)
    e1 = np.cross(v, [0, 0, 1.0]); e1 /= np.linalg.norm(e1); e2 = np.cross(v, e1)
    pts = (np.cos(np.radians(r))*v[None, :] + np.sin(np.radians(r))*(np.cos(th)[:, None]*e1 + np.sin(th)[:, None]*e2))
    lo, la = hp.vec2ang(pts, lonlat=True)
    hp.projplot(lo, la, lonlat=True, color=c, lw=2.0, ls="-", label=lab)
plt.legend(loc="lower right", fontsize=9, framealpha=0.9)
fig.savefig(OUT/"healpix_partition_overview.png", dpi=140, bbox_inches="tight")
plt.close(fig)

# =============== FIGURE 2: cost, across all 667 real map centers ===============
man = pd.read_csv(W/"neomod/pipeline/slurm/grid_map_manifest.csv")
ras, decs = centers_radec(man.delta_lon_from_antisun_deg.to_numpy(), man.lat_deg.to_numpy())
rad = np.radians(MAXSEP + MARGIN)
rows = []
for (_, r), ra_i, dec_i in zip(man.iterrows(), ras, decs):
    px = hp.query_disc(NSIDE, hp.ang2vec(float(ra_i), float(dec_i), lonlat=True), rad, inclusive=True)
    rows.append(dict(dlon=r.delta_lon_from_antisun_deg, lat=r.lat_deg, npix=len(px),
                     clones=int(counts[px].sum()), gib=sizes[px].sum()/2**30))
cen = pd.DataFrame(rows); cen.to_csv(OUT/"healpix_per_center_cost.csv", index=False)
mono_gib = os.path.getsize(W/"outputs/neomod3_projection_cache/neomod3_projection_20270825T000000.parquet")/2**30

fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
ax[0].hist(cen.npix, bins=np.arange(cen.npix.min()-0.5, cen.npix.max()+1.5), color="tab:blue")
ax[0].set_xlabel("HEALPix pixels read per map task"); ax[0].set_ylabel("map centers")
ax[0].set_title(f"(a) pixels touched per center\nmedian {cen.npix.median():.0f} of {NPIX} "
                f"({100*cen.npix.median()/NPIX:.0f}%)")
sc = ax[1].scatter(cen.dlon, cen.lat, c=cen.gib, cmap="magma", s=26)
ax[1].set_xlabel("Δlon from antisun (deg)"); ax[1].set_ylabel("ecliptic lat (deg)")
ax[1].set_title("(b) GiB read per center\n(bright = NEO-dense ecliptic tiles)")
plt.colorbar(sc, ax=ax[1], label="GiB")
ax[2].bar(["monolithic\n(read whole file)", "HEALPix\n(median center)"],
          [mono_gib, cen.gib.median()], color=["tab:red", "tab:green"])
ax[2].set_ylabel("GiB read per map task")
ax[2].set_title(f"(c) per-task read\n{mono_gib/cen.gib.median():.1f}× less I/O")
for i, v in enumerate([mono_gib, cen.gib.median()]):
    ax[2].text(i, v, f"{v:.2f}", ha="center", va="bottom", fontweight="bold")
fig.tight_layout(); fig.savefig(OUT/"healpix_cost.png", dpi=140, bbox_inches="tight"); plt.close(fig)

summary = dict(nside=NSIDE, npix=NPIX, clones=int(counts.sum()),
               partitioned_gib=float(sizes.sum()/2**30), monolithic_gib=float(mono_gib),
               median_pixels_per_center=int(cen.npix.median()),
               min_pixels=int(cen.npix.min()), max_pixels=int(cen.npix.max()),
               median_gib_per_center=float(cen.gib.median()),
               speedup=float(mono_gib/cen.gib.median()), max_pixrad_deg=float(MARGIN),
               emptiest_pixel=int(counts.min()), densest_pixel=int(counts.max()))
json.dump(summary, open(OUT/"healpix_summary.json", "w"), indent=2)
print(json.dumps(summary, indent=2))
