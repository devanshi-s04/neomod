#!/usr/bin/env python3
"""Build nightly two-detection tracklets from one Sorcha HDF5 output file.

This is the Phase 1 Hyak worker for the Sorcha/VDP/digest2 comparison. Each
invocation reads one Sorcha ``sorcha_results`` HDF5 file, groups detections by
object and night, keeps object-nights whose first/last detections span 3 to
90 minutes, and writes one compact parquet row per retained tracklet.

The output intentionally stores structured columns rather than digest2 MPC text
or VDP scores. Phase 2 can concatenate these parquet files, score VDP in bulk,
format MPC 80-column observations, and run digest2 in chunks.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
Path("/tmp/astropy").mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import GCRS, GeocentricTrueEcliptic, SkyCoord
from astropy.time import Time
from astropy.utils import iers


iers.conf.auto_download = False


SORCHA_KEY = "sorcha_results"
MIN_PAIR_DT_MIN = 3.0
MAX_PAIR_DT_MIN = 90.0
MU_SUN_AU3_DAY2 = 0.01720209895**2

REQUIRED_COLUMNS = [
    "ObjID",
    "fieldMJD_TAI",
    "RA_deg",
    "Dec_deg",
    "RARateCosDec_deg_day",
    "DecRate_deg_day",
    "PSFMag",
]

# Cartesian state vectors (hybrid/Cartesian Sorcha runs) — used to derive a, e, q
_CARTESIAN_COLUMNS = ["x", "y", "z", "xdot", "ydot", "zdot"]
# COM orbital elements (S3M/COM Sorcha runs) — used as fallback to derive a, e, q
_COM_COLUMNS = ["q", "e"]

OPTIONAL_COLUMNS = [
    "optFilter",
    "SNR",
    "H_r",
    "H_filter",
    "Range_LTC_km",
    "Obj_Sun_LTC_km",
] + _CARTESIAN_COLUMNS + _COM_COLUMNS

OUTPUT_COLUMNS = [
    "tracklet_id",
    "source_file",
    "source_file_index",
    "ObjID",
    "night",
    "n_det",
    "night_span_min",
    "dt_min",
    "mean_ra",
    "mean_dec",
    "mean_dra",
    "mean_ddec",
    "mean_mag",
    "ra0",
    "dec0",
    "mjd0_tai",
    "mjd0_utc",
    "mag0",
    "filter0",
    "snr0",
    "ra1",
    "dec1",
    "mjd1_tai",
    "mjd1_utc",
    "mag1",
    "filter1",
    "snr1",
    "H_r",
    "H_filter",
    "range_km",
    "obj_sun_km",
    "ecl_lon",
    "ecl_lat",
    "prob_map",
    "prob_map_file",
    "prob_map_obstime_str",
    "prob_map_center_lon_deg",
    "prob_map_center_lat_deg",
    "prob_map_dist_deg",
    "population",
    "a_au",
    "e",
    "q_au",
    "n_det_per_night",
]


@dataclass(frozen=True)
class ProbMapMeta:
    path: Path
    label: str
    obstime_str: str
    obstime_mjd: float        # epoch as MJD UTC for fast epoch-proximity lookup
    center_lon_deg: float
    center_lat_deg: float
    max_sep_deg: float
    # Present only for antisun-relative grid maps (sorcha_gen_maps_grid.py).
    # When set, the map is assigned by nearest (delta_lon, lat) offset rather
    # than by epoch proximity to the current antisun.
    delta_lon_from_antisun_deg: float | None = None
    grid_lat_deg: float | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build two-detection nightly tracklets from one Sorcha HDF5 file."
    )
    parser.add_argument("--infile", type=Path, help="Sorcha HDF5 file to process.")
    parser.add_argument(
        "--outfile",
        type=Path,
        help="Output parquet path. Defaults to outdir/tracklets_<input-stem>.parquet.",
    )
    parser.add_argument(
        "--file_index",
        type=int,
        help="Index into sorted inst*_part*.h5 files under --indir.",
    )
    parser.add_argument(
        "--indir",
        type=Path,
        default=Path("outputs/production_2yr"),
        help="Directory containing production Sorcha HDF5 files.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("outputs/tracklets"),
        help="Directory for default parquet outputs.",
    )
    parser.add_argument(
        "--prob_maps_dir",
        type=Path,
        default=Path("prob_maps"),
        help="Directory containing VDP prob_maps_*.npz files for footprint assignment.",
    )
    parser.add_argument(
        "--min-dt-min",
        type=float,
        default=MIN_PAIR_DT_MIN,
        help="Minimum first/last detection separation in minutes.",
    )
    parser.add_argument(
        "--max-dt-min",
        type=float,
        default=MAX_PAIR_DT_MIN,
        help="Maximum first/last detection separation in minutes.",
    )
    parser.add_argument(
        "--keep-unmapped",
        action="store_true",
        help="Keep tracklets outside all probability-map footprints.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing output parquet.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip this task if the output parquet already exists and is readable.",
    )
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path, int | None]:
    if args.infile is None:
        if args.file_index is None:
            raise ValueError("Provide either --infile or --file_index.")
        files = sorted(args.indir.glob("inst*_part*.h5"))
        if not files:
            raise FileNotFoundError(f"No inst*_part*.h5 files found in {args.indir}")
        if args.file_index < 0 or args.file_index >= len(files):
            raise IndexError(
                f"--file_index {args.file_index} is outside 0..{len(files) - 1}"
            )
        infile = files[args.file_index]
        source_file_index = args.file_index
    else:
        infile = args.infile
        source_file_index = args.file_index

    if not infile.exists():
        raise FileNotFoundError(infile)

    if args.outfile is None:
        outfile = args.outdir / f"tracklets_{infile.stem}.parquet"
    else:
        outfile = args.outfile

    if outfile.exists() and not args.overwrite and not args.skip_existing:
        raise FileExistsError(f"{outfile} exists; pass --overwrite to replace it.")

    return infile, outfile, source_file_index


def output_is_readable(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        pd.read_parquet(path, columns=["tracklet_id"])
    except Exception:
        return False
    return True


def read_sorcha_results(path: Path) -> pd.DataFrame:
    df = pd.read_hdf(path, key=SORCHA_KEY)
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise KeyError(f"{path} is missing required Sorcha columns: {missing}")

    keep = REQUIRED_COLUMNS + [col for col in OPTIONAL_COLUMNS if col in df.columns]
    df = df[keep].copy()

    finite = np.ones(len(df), dtype=bool)
    for col in [
        "fieldMJD_TAI",
        "RA_deg",
        "Dec_deg",
        "RARateCosDec_deg_day",
        "DecRate_deg_day",
        "PSFMag",
    ]:
        finite &= np.isfinite(df[col].to_numpy(float))
    df = df[finite].copy()

    dec_rad = np.deg2rad(df["Dec_deg"].to_numpy(float))
    cos_dec = np.cos(dec_rad)
    df["dra_deg_day"] = np.where(
        np.abs(cos_dec) > 1e-12,
        df["RARateCosDec_deg_day"].to_numpy(float) / cos_dec,
        np.nan,
    )
    df["ddec_deg_day"] = df["DecRate_deg_day"].to_numpy(float)
    df["night"] = np.floor(df["fieldMJD_TAI"].to_numpy(float)).astype(np.int64)
    return df[np.isfinite(df["dra_deg_day"].to_numpy(float))].copy()


def circular_mean_deg(values: pd.Series | np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    ang = np.deg2rad(arr)
    return float(np.rad2deg(np.arctan2(np.mean(np.sin(ang)), np.mean(np.cos(ang)))) % 360.0)


def _pair_col(df_subset: pd.DataFrame, col: str) -> np.ndarray:
    """Column values from a first/last subset DataFrame, or NaN array if absent."""
    if col not in df_subset.columns:
        return np.full(len(df_subset), np.nan)
    return df_subset[col].to_numpy()


def _pair_nanmean(
    first_df: pd.DataFrame, last_df: pd.DataFrame, col: str
) -> np.ndarray:
    """Per-row mean of col across first and last detection; NaN when both absent."""
    if col not in first_df.columns:
        return np.full(len(first_df), np.nan)
    a = first_df[col].to_numpy(float)
    b = last_df[col].to_numpy(float)
    with np.errstate(all="ignore"):
        return np.nanmean(np.stack([a, b], axis=1), axis=1)


def classify_population(states: pd.DataFrame) -> pd.DataFrame:
    states = states.drop_duplicates("ObjID").copy()

    if "x" in states.columns and "xdot" in states.columns:
        # Cartesian path (hybrid/Cartesian Sorcha output)
        r_vec = states[["x", "y", "z"]].to_numpy(float)
        v_vec = states[["xdot", "ydot", "zdot"]].to_numpy(float)
        r = np.linalg.norm(r_vec, axis=1)
        v2 = np.einsum("ij,ij->i", v_vec, v_vec)
        energy = 0.5 * v2 - MU_SUN_AU3_DAY2 / r
        a = np.where(energy != 0.0, -MU_SUN_AU3_DAY2 / (2.0 * energy), np.nan)
        h_vec = np.cross(r_vec, v_vec)
        e_vec = np.cross(v_vec, h_vec) / MU_SUN_AU3_DAY2 - r_vec / r[:, None]
        ecc = np.linalg.norm(e_vec, axis=1)
        q = np.where(np.isfinite(a), a * (1.0 - ecc), np.nan)
    else:
        # COM orbital elements path (S3M/COM Sorcha output — q and e are stored directly)
        q = states["q"].to_numpy(float)
        ecc = states["e"].to_numpy(float)
        a = np.where(ecc != 1.0, q / (1.0 - ecc), np.nan)

    pop = np.full(len(states), "other", dtype=object)
    pop[np.isfinite(q) & (q < 1.3)] = "NEO"
    pop[np.isfinite(a) & (a > 30.0)] = "TNO"
    pop[np.isfinite(a) & (a >= 4.8) & (a <= 5.6)] = "Trojan"
    pop[np.isfinite(a) & (a > 2.0) & (a < 3.3) & np.isfinite(ecc) & (ecc < 0.3)] = "MBA"

    return pd.DataFrame(
        {
            "ObjID": states["ObjID"].to_numpy(),
            "population": pop,
            "a_au": a,
            "e": ecc,
            "q_au": q,
        }
    )


def build_tracklets(
    df: pd.DataFrame,
    infile: Path,
    source_file_index: int | None,
    min_dt_min: float,
    max_dt_min: float,
) -> pd.DataFrame:
    df = df.sort_values(["ObjID", "night", "fieldMJD_TAI"], kind="mergesort")
    g = df.groupby(["ObjID", "night"], sort=False)

    # Identify the first and last row of every (ObjID, night) group without a
    # Python loop.  cumcount(ascending=True)==0 marks the first row;
    # cumcount(ascending=False)==0 marks the last.
    cum_fwd = g.cumcount()
    cum_rev = g.cumcount(ascending=False)
    first_df = df[cum_fwd == 0].set_index(["ObjID", "night"])
    last_df  = df[cum_rev == 0].set_index(["ObjID", "night"])
    n_det    = g.size().rename("n_det")

    # Time-separation filter: first/last must span 3–90 min, group must have ≥2.
    dt_min = (last_df["fieldMJD_TAI"] - first_df["fieldMJD_TAI"]) * 1440.0
    valid  = (n_det >= 2) & (dt_min >= min_dt_min) & (dt_min <= max_dt_min)
    if not valid.any():
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    first_v = first_df[valid]
    last_v  = last_df[valid]
    n_v     = n_det[valid]
    dt_v    = dt_min[valid]
    idx     = dt_v.index
    obj_ids = idx.get_level_values("ObjID")
    nights  = idx.get_level_values("night")

    # Circular mean of RA across the pair, handling the 0°/360° wrap.
    ra0 = first_v["RA_deg"].to_numpy(float)
    ra1 = last_v["RA_deg"].to_numpy(float)
    mean_ra = (
        np.rad2deg(
            np.arctan2(
                np.sin(np.deg2rad(ra0)) + np.sin(np.deg2rad(ra1)),
                np.cos(np.deg2rad(ra0)) + np.cos(np.deg2rad(ra1)),
            )
        ) % 360.0
    )

    mjd0 = first_v["fieldMJD_TAI"].to_numpy(float)
    mjd1 = last_v["fieldMJD_TAI"].to_numpy(float)

    out = pd.DataFrame(
        {
            "tracklet_id": [
                f"{infile.stem}:{oid}:{n}" for oid, n in zip(obj_ids, nights)
            ],
            "source_file": infile.name,
            "source_file_index": -1 if source_file_index is None else source_file_index,
            "ObjID": obj_ids,
            "night": nights.astype(int),
            "n_det": n_v.values,
            "n_det_per_night": n_v.values,
            "night_span_min": dt_v.values,
            "dt_min": dt_v.values,
            "mean_ra": mean_ra,
            "mean_dec": (first_v["Dec_deg"].to_numpy(float) + last_v["Dec_deg"].to_numpy(float)) / 2.0,
            "mean_dra": (first_v["dra_deg_day"].to_numpy(float) + last_v["dra_deg_day"].to_numpy(float)) / 2.0,
            "mean_ddec": (first_v["ddec_deg_day"].to_numpy(float) + last_v["ddec_deg_day"].to_numpy(float)) / 2.0,
            "mean_mag": (first_v["PSFMag"].to_numpy(float) + last_v["PSFMag"].to_numpy(float)) / 2.0,
            "ra0": ra0,
            "dec0": first_v["Dec_deg"].to_numpy(float),
            "mjd0_tai": mjd0,
            "mag0": first_v["PSFMag"].to_numpy(float),
            "filter0": _pair_col(first_v, "optFilter"),
            "snr0":    _pair_col(first_v, "SNR"),
            "ra1": ra1,
            "dec1": last_v["Dec_deg"].to_numpy(float),
            "mjd1_tai": mjd1,
            "mag1": last_v["PSFMag"].to_numpy(float),
            "filter1": _pair_col(last_v, "optFilter"),
            "snr1":    _pair_col(last_v, "SNR"),
            "H_r":       _pair_nanmean(first_v, last_v, "H_r"),
            "H_filter":  _pair_nanmean(first_v, last_v, "H_filter"),
            "range_km":  _pair_nanmean(first_v, last_v, "Range_LTC_km"),
            "obj_sun_km": _pair_nanmean(first_v, last_v, "Obj_Sun_LTC_km"),
        }
    )

    out["mjd0_utc"] = Time(mjd0, format="mjd", scale="tai").utc.mjd
    out["mjd1_utc"] = Time(mjd1, format="mjd", scale="tai").utc.mjd
    if "x" in df.columns and "xdot" in df.columns:
        pop = classify_population(df[["ObjID", "x", "y", "z", "xdot", "ydot", "zdot"]])
    else:
        pop = classify_population(df[["ObjID", "q", "e"]])
    out = out.merge(pop, on="ObjID", how="left")
    return out


def load_prob_map_metadata(prob_maps_dir: Path) -> list[ProbMapMeta]:
    if not prob_maps_dir.exists():
        raise FileNotFoundError(prob_maps_dir)

    metas = []
    for path in sorted(prob_maps_dir.glob("*.npz")):
        with np.load(path, allow_pickle=True) as z:
            label = str(z["center_label"]) if "center_label" in z.files else path.stem
            obstime_str = str(z["obstime_str"])
            metas.append(
                ProbMapMeta(
                    path=path,
                    label=label,
                    obstime_str=obstime_str,
                    obstime_mjd=Time(obstime_str).utc.mjd,
                    center_lon_deg=float(z["center_lon_deg"]),
                    center_lat_deg=float(z["center_lat_deg"]),
                    max_sep_deg=float(z["max_sep_deg"]),
                    delta_lon_from_antisun_deg=(
                        float(z["delta_lon_from_antisun_deg"])
                        if "delta_lon_from_antisun_deg" in z.files else None
                    ),
                    grid_lat_deg=(
                        float(z["grid_lat_deg"])
                        if "grid_lat_deg" in z.files else None
                    ),
                )
            )

    if not metas:
        raise FileNotFoundError(f"No .npz probability maps found in {prob_maps_dir}")
    return metas


_ECLIPTIC_OBL_RAD = np.deg2rad(23.439291)
_J2000_MJD = 51544.5


def _antisun_ecl_lon(mjd: np.ndarray) -> np.ndarray:
    """Approximate geocentric ecliptic longitude of the antisun (deg) for each MJD.

    Uses the low-precision sun-longitude formula (accurate to ~1°), fast and
    vectorised — no astropy call needed.  The antisun is opposite the Sun in
    ecliptic longitude.
    """
    sun_lon = (280.46 + 0.9856474 * (mjd - _J2000_MJD)) % 360.0
    return (sun_lon + 180.0) % 360.0


def _ecl_lon_from_radec(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    """Convert equatorial (RA, Dec) to geocentric ecliptic longitude (deg).

    Uses a fixed obliquity — sufficient for a footprint proximity test.
    """
    ra  = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    c, s = np.cos(_ECLIPTIC_OBL_RAD), np.sin(_ECLIPTIC_OBL_RAD)
    x    = np.cos(dec) * np.cos(ra)
    y    = np.cos(dec) * np.sin(ra)
    z    = np.sin(dec)
    y_ecl = c * y + s * z
    return np.rad2deg(np.arctan2(y_ecl, x)) % 360.0


def _ecl_lat_from_radec(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    """Convert equatorial (RA, Dec) to geocentric ecliptic latitude (deg).

    Uses a fixed obliquity — sufficient for grid-cell proximity in latitude.
    """
    ra  = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    c, s = np.cos(_ECLIPTIC_OBL_RAD), np.sin(_ECLIPTIC_OBL_RAD)
    y     = np.cos(dec) * np.sin(ra)
    z     = np.sin(dec)
    z_ecl = -s * y + c * z
    return np.rad2deg(np.arcsin(np.clip(z_ecl, -1.0, 1.0)))


def _angular_sep_ecl_lon(lon_a: np.ndarray, lon_b: np.ndarray) -> np.ndarray:
    """Great-circle-approximate angular separation between two ecliptic longitudes
    at lat ≈ 0 (deg). Fast scalar-or-array operation."""
    diff = np.abs(((lon_a - lon_b) + 180.0) % 360.0 - 180.0)
    return diff


def assign_probability_maps(tracklets: pd.DataFrame, metas: list[ProbMapMeta]) -> pd.DataFrame:
    """Assign each tracklet to its best-matching probability map.

    For 'antisun'-labelled maps the footprint criterion is whether the tracklet
    is within max_sep_deg of the **current antisun direction at observation
    time** — not the map's past sky position.  This prevents off-axis tracklets
    (those whose elongation from the Sun has changed since the map epoch) from
    being assigned to antisun maps with wrong velocity calibrations.

    For all other maps the original criterion (sky separation from the fixed map
    center) is retained.

    Among maps whose footprint contains a tracklet we pick the one whose epoch
    is closest to the tracklet's observation time.

    Grid mode (antisun-relative sky grid, sorcha_gen_maps_grid.py): if the maps
    carry ``delta_lon_from_antisun_deg`` metadata, each tracklet is placed in the
    antisun-relative frame (Δlon = ecl_lon − antisun_lon at observation time,
    plus ecliptic lat) and assigned to the nearest grid center in that 2-D space.
    Tracklets within ``SUN_EXCL_DEG`` of the Sun are left unassigned.
    """
    if len(tracklets) == 0:
        tracklets["ecl_lon"] = pd.Series(dtype=float)
        tracklets["ecl_lat"] = pd.Series(dtype=float)
        tracklets["prob_map"] = pd.Series(dtype=object)
        tracklets["prob_map_dist_deg"] = pd.Series(dtype=float)
        return tracklets

    ra  = tracklets["mean_ra"].to_numpy(float)
    dec = tracklets["mean_dec"].to_numpy(float)
    mjd_obs = tracklets["mjd0_utc"].to_numpy(float)
    N = len(tracklets)
    M = len(metas)

    # Precompute per-tracklet quantities used in footprint checks
    ecl_lon_obs    = _ecl_lon_from_radec(ra, dec)           # ecliptic lon at obs time
    ecl_lat_obs    = _ecl_lat_from_radec(ra, dec)           # ecliptic lat at obs time
    antisun_lon    = _antisun_ecl_lon(mjd_obs)               # current antisun lon

    is_grid = any(m.delta_lon_from_antisun_deg is not None for m in metas)

    # (N, M) footprint membership, plus a per-(tracklet,map) "score" we minimise
    # to choose the assigned map:
    #   grid maps          -> 2-D (Δlon, lat) distance to the grid center
    #   antisun monthly    -> |epoch − obs_time| (days)
    in_footprint = np.zeros((N, M), dtype=bool)
    score        = np.full((N, M), np.inf, dtype=float)

    if is_grid:
        SUN_EXCL_DEG = 40.0
        # antisun-relative longitude offset of each tracklet, wrapped to [-180,180]
        delta_lon_obs = ((ecl_lon_obs - antisun_lon) + 180.0) % 360.0 - 180.0
        in_sun_zone   = np.abs(delta_lon_obs) > (180.0 - SUN_EXCL_DEG)
        for j, meta in enumerate(metas):
            d_lon = np.abs(
                ((delta_lon_obs - meta.delta_lon_from_antisun_deg) + 180.0) % 360.0 - 180.0
            )
            d_lat = np.abs(ecl_lat_obs - meta.grid_lat_deg)
            score[:, j] = np.hypot(d_lon, d_lat)
        # every non-sun-zone tracklet is assigned to its nearest grid center
        in_footprint[~in_sun_zone, :] = True
    else:
        for j, meta in enumerate(metas):
            score[:, j] = np.abs(mjd_obs - meta.obstime_mjd)
            # All antisun maps are calibrated for the antisun geometry, so the
            # footprint criterion for every map is: is the tracklet within
            # max_sep of the CURRENT antisun at observation time?  This prevents
            # off-axis tracklets from being scored with a near-opposition map.
            dist = _angular_sep_ecl_lon(ecl_lon_obs, antisun_lon)
            in_footprint[:, j] = dist <= meta.max_sep_deg

    # Among in-footprint maps, pick the one with the smallest score.
    score_masked = np.where(in_footprint, score, np.inf)
    best_j = np.argmin(score_masked, axis=1)                     # (N,)
    any_in_footprint = in_footprint.any(axis=1)                  # (N,)

    # Build per-map attribute arrays for fast indexing
    labels      = np.array([m.label          for m in metas], dtype=object)
    filenames   = np.array([m.path.name      for m in metas], dtype=object)
    obstimes    = np.array([m.obstime_str    for m in metas], dtype=object)
    center_lons = np.array([m.center_lon_deg for m in metas], dtype=float)
    center_lats = np.array([m.center_lat_deg for m in metas], dtype=float)

    best_label      = np.where(any_in_footprint, labels[best_j],      None)
    best_file       = np.where(any_in_footprint, filenames[best_j],   None)
    best_obstime    = np.where(any_in_footprint, obstimes[best_j],    None)
    best_center_lon = np.where(any_in_footprint, center_lons[best_j], np.nan)
    best_center_lat = np.where(any_in_footprint, center_lats[best_j], np.nan)
    # Distance-to-assigned-center diagnostic:
    #   grid maps -> the 2-D (Δlon, lat) distance to the chosen grid center
    #   antisun   -> approximate distance from the current antisun
    if is_grid:
        best_dist = np.where(
            any_in_footprint, score_masked[np.arange(N), best_j], np.nan,
        )
    else:
        best_dist = np.where(
            any_in_footprint,
            _angular_sep_ecl_lon(ecl_lon_obs, antisun_lon),   # approx dist from antisun
            np.nan,
        )

    # Ecliptic lon already computed above from obs-time ra/dec (ecl_lon_obs).
    # For ecl_lat use the first map's frame (small error, diagnostic only).
    t0  = Time(metas[0].obstime_str)
    sc0 = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame=GCRS(obstime=t0))
    ecl0 = sc0.transform_to(GeocentricTrueEcliptic(obstime=t0))

    tracklets = tracklets.copy()
    tracklets["ecl_lon"]               = ecl0.lon.deg
    tracklets["ecl_lat"]               = ecl0.lat.deg
    tracklets["prob_map"]              = best_label
    tracklets["prob_map_file"]         = best_file
    tracklets["prob_map_obstime_str"]  = best_obstime
    tracklets["prob_map_center_lon_deg"] = best_center_lon
    tracklets["prob_map_center_lat_deg"] = best_center_lat
    tracklets["prob_map_dist_deg"]     = best_dist
    return tracklets


def main() -> None:
    args = parse_args()
    if args.min_dt_min < 0 or args.max_dt_min <= args.min_dt_min:
        raise ValueError("--max-dt-min must be greater than --min-dt-min.")

    infile, outfile, source_file_index = resolve_paths(args)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    if args.skip_existing and output_is_readable(outfile):
        print(f"Skipping existing readable output: {outfile}", flush=True)
        return

    print(f"Reading {infile}", flush=True)
    detections = read_sorcha_results(infile)
    print(f"  detections after finite-rate filter: {len(detections):,}", flush=True)

    print("Building two-detection nightly tracklets", flush=True)
    tracklets = build_tracklets(
        detections,
        infile=infile,
        source_file_index=source_file_index,
        min_dt_min=args.min_dt_min,
        max_dt_min=args.max_dt_min,
    )
    n_before_maps = len(tracklets)
    print(f"  tracklets after pair cuts: {n_before_maps:,}", flush=True)

    metas = load_prob_map_metadata(args.prob_maps_dir)
    print(f"Assigning epoch-aware probability map from {len(metas)} map files", flush=True)
    tracklets = assign_probability_maps(tracklets, metas)

    if not args.keep_unmapped:
        tracklets = tracklets[tracklets["prob_map"].notna()].reset_index(drop=True)
        print(
            f"  kept {len(tracklets):,}/{n_before_maps:,} tracklets inside map footprints",
            flush=True,
        )

    for col in OUTPUT_COLUMNS:
        if col not in tracklets.columns:
            tracklets[col] = np.nan
    tracklets = tracklets[OUTPUT_COLUMNS]

    tmp_outfile = outfile.with_name(f".{outfile.name}.tmp_{os.getpid()}")
    try:
        tracklets.to_parquet(tmp_outfile, index=False)
        tmp_outfile.replace(outfile)
    finally:
        if tmp_outfile.exists():
            tmp_outfile.unlink()
    print(f"Wrote {len(tracklets):,} tracklets to {outfile}", flush=True)


if __name__ == "__main__":
    main()
