"""
Generate a comprehensive antisun-relative sky grid of VDP probability maps.

Unlike the 24 monthly antisun maps (all at ecliptic lat=0, one per month), this
builds a full sky grid of ~500-700 maps defined in **antisun-relative ecliptic
coordinates** (Delta-lon-from-antisun, ecliptic-lat).  The grid is
time-independent: because NEO/MBA velocity structure is set by geometry relative
to the Sun, a map built at offset (dlon, lat) from the antisun is reusable at any
epoch.  At scoring time each LSST observation is placed in this frame by computing
the antisun position at the observation time and taking the offset.

Latitude sampling is intentionally NON-UNIFORM: fine (1-2 deg) steps near the
ecliptic where NEOs/MBAs concentrate, coarser steps at high latitude.

Grid geometry (all configurable via CLI):
    longitude : --lon-step deg, antisun-relative, from -180 to +180,
                excluding any center within --sun-exclusion deg of the Sun
                (the Sun sits at |dlon| = 180).
    latitude  : --lat-base is a comma list of NON-NEGATIVE values, expanded
                symmetrically about the ecliptic (0 once, +/- for the rest).
                Default packs fine steps near 0 and coarse steps at high lat.

Each map is generated exactly like the monthly antisun maps
(velocity_density_pipeline_gmm.generate_probability_maps) but centered at the
antisun-relative (dlon, lat) offset rather than always at (antisun, 0).  After
generation the .npz is augmented with two exact metadata keys so that
sorcha_postprocess.py can assign tracklets to the nearest grid center without
re-deriving the antisun position (which the low-precision formula would do with
~1 deg error -- too coarse for the 1 deg ecliptic spacing).

Usage
-----
List the grid (and write the CSV manifest) without generating anything:
    python sorcha_gen_maps_grid.py --list-only --map-grid-file map_grid.csv

Generate the Nth map (driven by a Slurm array task id):
    python sorcha_gen_maps_grid.py --task-id 0 --n-jobs 16 \
        --prob-maps-dir prob_maps_grid

Generate a single explicit map (manual test):
    python sorcha_gen_maps_grid.py --delta-lon 0 --lat 0 --n-jobs 16 \
        --output prob_maps_grid/prob_maps_grid_dlon+000_lat+00.npz
"""
import argparse
import csv
import os
import sys

import numpy as np

WORKDIR = "/mmfs1/gscratch/dirac/ds2004/sorcha"
NEOMOD  = os.path.join(WORKDIR, "neomod")

sys.path.insert(0, os.path.join(NEOMOD, "adam_core_stub"))
sys.path.insert(0, os.path.join(NEOMOD, "src"))

# Default reference epoch for ALL grid maps.  Mid-point of the v5.0 window
# (2025-11 .. 2027-11).  The grid is antisun-relative so the absolute epoch only
# fixes where the patches land on the real sky for this one generation; the
# velocity statistics are epoch-independent in the antisun frame.
DEFAULT_REF_OBSTIME = "2026-01-01T00:00:00"

# Non-negative latitude knots; expanded symmetrically about the ecliptic.
# Fine 1 deg steps to 5 deg, then progressively coarser out to 50 deg.
DEFAULT_LAT_BASE = "0,1,2,3,4,5,8,12,18,25,35,50"


def make_symmetric_lats(base_values):
    """Expand non-negative latitude knots into a sorted symmetric list."""
    vals = set()
    for v in base_values:
        vals.add(round(float(v), 6))
        vals.add(round(-float(v), 6))
    return sorted(vals)


def build_grid(lon_step, lat_values, sun_exclusion):
    """Return the ordered list of (delta_lon, lat) grid centers.

    delta_lon runs over antisun-relative longitudes -180..+180 (exclusive of the
    +180 duplicate) in lon_step increments, dropping any center within
    sun_exclusion degrees of the Sun (Sun at |delta_lon| = 180).
    Order: latitude outer (low to high), longitude inner -- deterministic so a
    Slurm --task-id maps to a stable grid point.
    """
    dlon_limit = 180.0 - sun_exclusion          # keep |delta_lon| <= this
    dlons = []
    d = -180.0
    while d < 180.0 - 1e-9:
        if abs(d) <= dlon_limit + 1e-9:
            dlons.append(round(d, 6))
        d += lon_step
    grid = []
    for lat in lat_values:
        for dlon in dlons:
            grid.append((dlon, lat))
    return grid, dlons


def map_filename(delta_lon, lat):
    """Stable, human-readable filename for a grid center."""
    return f"prob_maps_grid_dlon{int(round(delta_lon)):+04d}_lat{int(round(lat)):+03d}.npz"


def map_label(delta_lon, lat):
    return f"grid_dlon{int(round(delta_lon)):+04d}_lat{int(round(lat)):+03d}"


def antisun_lon_at(obstime_str):
    """High-precision antisun ecliptic longitude (deg) via astropy get_sun."""
    from astropy.time import Time
    from astropy.coordinates import get_sun, GeocentricTrueEcliptic
    t = Time(obstime_str, scale="utc")
    sun = get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
    return (sun.lon.deg + 180.0) % 360.0


def augment_npz(out_path, delta_lon, lat, ref_obstime):
    """Re-save the generated .npz with exact antisun-relative metadata keys."""
    with np.load(out_path, allow_pickle=True) as z:
        data = {k: z[k] for k in z.files}
    data["delta_lon_from_antisun_deg"] = np.asarray(float(delta_lon), dtype=np.float64)
    data["grid_lat_deg"]               = np.asarray(float(lat),       dtype=np.float64)
    data["ref_obstime_str"]            = np.asarray(ref_obstime)
    np.savez_compressed(out_path, **data)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # --- grid geometry (configurable) ---
    p.add_argument("--lon-step",      type=float, default=10.0,
                   help="Longitude step in degrees (advisor: 5 or 10). Default 10.")
    p.add_argument("--lat-base",      type=str, default=DEFAULT_LAT_BASE,
                   help="Comma list of NON-NEGATIVE latitude knots, expanded "
                        "symmetrically. Default packs fine steps near the ecliptic.")
    p.add_argument("--lat-points",    type=str, default=None,
                   help="Explicit comma list of latitudes (overrides --lat-base; "
                        "used verbatim, not symmetrised).")
    p.add_argument("--sun-exclusion", type=float, default=40.0,
                   help="Exclude grid centers within this many deg of the Sun. Default 40.")
    p.add_argument("--ref-obstime",   type=str, default=DEFAULT_REF_OBSTIME,
                   help="Reference epoch for all grid maps (ISO).")
    # --- selection / output ---
    p.add_argument("--task-id",       type=int, default=None,
                   help="Generate only the Nth grid map (Slurm array index).")
    p.add_argument("--delta-lon",     type=float, default=None,
                   help="Generate a single explicit map at this antisun-relative lon.")
    p.add_argument("--lat",           type=float, default=None,
                   help="Latitude for the single explicit map (with --delta-lon).")
    p.add_argument("--prob-maps-dir", type=str, default="prob_maps_grid",
                   help="Output directory for grid maps.")
    p.add_argument("--output",        type=str, default=None,
                   help="Explicit output .npz path (overrides --prob-maps-dir naming).")
    p.add_argument("--map-grid-file", type=str, default=None,
                   help="Write the (index, delta_lon, lat, filename) manifest CSV here.")
    p.add_argument("--list-only",     action="store_true",
                   help="Print the grid summary, write manifest if requested, exit.")
    p.add_argument("--n-jobs",        type=int, default=1,
                   help="Parallel workers for the VDP density estimator.")
    p.add_argument("--save-overlays", action="store_true",
                   help="Persist cloned-object overlays (vlam/vbeta/heliocentric) "
                        "for diagnostic plots. Inflates the .npz — use for test "
                        "maps only, not the full production grid.")
    p.add_argument("--overwrite",     action="store_true")
    args = p.parse_args()

    # Build the grid
    if args.lat_points is not None:
        lat_values = [round(float(x), 6) for x in args.lat_points.split(",") if x.strip() != ""]
        lat_values = sorted(set(lat_values))
    else:
        base = [float(x) for x in args.lat_base.split(",") if x.strip() != ""]
        lat_values = make_symmetric_lats(base)

    grid, dlons = build_grid(args.lon_step, lat_values, args.sun_exclusion)

    print(f"grid geometry:", flush=True)
    print(f"  lon-step       = {args.lon_step} deg  -> {len(dlons)} usable longitudes", flush=True)
    print(f"  sun-exclusion  = {args.sun_exclusion} deg", flush=True)
    print(f"  latitudes ({len(lat_values)}) = {lat_values}", flush=True)
    print(f"  TOTAL MAPS     = {len(grid)}", flush=True)
    print(f"  ref-obstime    = {args.ref_obstime}", flush=True)

    # Manifest CSV
    if args.map_grid_file:
        os.makedirs(os.path.dirname(os.path.abspath(args.map_grid_file)) or ".", exist_ok=True)
        with open(args.map_grid_file, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["index", "delta_lon_from_antisun_deg", "lat_deg", "filename"])
            for i, (dlon, lat) in enumerate(grid):
                w.writerow([i, dlon, lat, map_filename(dlon, lat)])
        print(f"  manifest -> {args.map_grid_file}", flush=True)

    if args.list_only:
        sys.exit(0)

    # Resolve which single map to generate
    if args.delta_lon is not None:
        if args.lat is None:
            p.error("--lat is required with --delta-lon")
        delta_lon, lat = float(args.delta_lon), float(args.lat)
    elif args.task_id is not None:
        if not (0 <= args.task_id < len(grid)):
            p.error(f"--task-id {args.task_id} out of range 0..{len(grid)-1}")
        delta_lon, lat = grid[args.task_id]
    else:
        p.error("provide one of --task-id, (--delta-lon and --lat), or --list-only")

    # Output path
    if args.output:
        out = args.output if os.path.isabs(args.output) else os.path.join(WORKDIR, args.output)
    else:
        out = os.path.join(WORKDIR, args.prob_maps_dir, map_filename(delta_lon, lat))

    if not args.overwrite and os.path.exists(out):
        print(f"[skip] already exists: {out}", flush=True)
        sys.exit(0)

    os.makedirs(os.path.dirname(out), exist_ok=True)

    # Antisun-relative -> absolute ecliptic center for this reference epoch
    antisun_lon = antisun_lon_at(args.ref_obstime)
    center_lon  = (antisun_lon + delta_lon) % 360.0
    center_lat  = lat
    label       = map_label(delta_lon, lat)

    print(f"antisun_lon(ref) = {antisun_lon:.4f} deg", flush=True)
    print(f"center           = ({center_lon:.4f}, {center_lat:.4f})  "
          f"[dlon={delta_lon:+.1f}, lat={lat:+.1f}]", flush=True)
    print(f"label            = {label}", flush=True)
    print(f"output           = {out}", flush=True)
    sys.stdout.flush()

    # s3m_loader searches [".", "S3Mdata"] relative to CWD -- must chdir to neomod/
    os.chdir(NEOMOD)
    import velocity_density_pipeline_gmm as vdp  # noqa: E402

    vdp.generate_probability_maps(
        obstime_str=args.ref_obstime,
        output_path=out,
        center_lon_deg=center_lon,
        center_lat_deg=center_lat,
        center_label=label,
        n_jobs=args.n_jobs,
        save_overlays=args.save_overlays,
    )

    augment_npz(out, delta_lon, lat, args.ref_obstime)
    print(f"[done] {out}  (augmented with delta_lon_from_antisun_deg, grid_lat_deg)", flush=True)


if __name__ == "__main__":
    main()
