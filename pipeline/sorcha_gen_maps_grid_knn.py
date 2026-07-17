"""
Generate the 667-map antisun-relative sky grid using the OLD kNN/K|M pipeline.

Identical grid geometry to sorcha_gen_maps_grid.py (same lon/lat sampling,
same reference epoch, same MBA clone_factor override). The only difference is
the density estimator and NEO cloner:

    sorcha_gen_maps_grid.py     -> velocity_density_pipeline_gmm   (GMM NEO cloner, wider grid)
    sorcha_gen_maps_grid_knn.py -> velocity_density_pipeline_fast  (kNN density, K|M NEO cloner)

Pipeline differences to be aware of:
  - velocity_density_pipeline_fast uses s3m_loader directly (no VDP_LOADER env needed)
  - Default velocity grid: (-1.5, 1.5) deg/day, step=0.01 (301x301)
    vs GMM grid: (-2.0, 2.0), step=0.01 (401x401)
  - NEO clone_factor default: 300 (vs GMM 80 + NEOMD3 augmentation)
  - NEO density estimator: kNN (vs GMM)
  - MBA/TNO/Trojans: same conditional K|M cloner in both pipelines

Usage
-----
    python sorcha_gen_maps_grid_knn.py --task-id 0 --n-jobs 16 \
        --prob-maps-dir prob_maps_grid_knn_s3m

    python sorcha_gen_maps_grid_knn.py --list-only --map-grid-file map_grid_knn.csv
"""
import argparse
import copy
import csv
import os
import sys

import numpy as np

WORKDIR = "/mmfs1/gscratch/dirac/ds2004/sorcha"
NEOMOD  = os.path.join(WORKDIR, "neomod")

sys.path.insert(0, os.path.join(NEOMOD, "adam_core_stub"))
sys.path.insert(0, os.path.join(NEOMOD, "src"))

DEFAULT_REF_OBSTIME = "2026-01-01T00:00:00"

DEFAULT_LAT_BASE = "0,1,2,3,4,5,8,12,18,25,35,50"


def make_symmetric_lats(base_values):
    vals = set()
    for v in base_values:
        vals.add(round(float(v), 6))
        vals.add(round(-float(v), 6))
    return sorted(vals)


def build_grid(lon_step, lat_values, sun_exclusion):
    dlon_limit = 180.0 - sun_exclusion
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
    return f"prob_maps_grid_dlon{int(round(delta_lon)):+04d}_lat{int(round(lat)):+03d}.npz"


def map_label(delta_lon, lat):
    return f"grid_dlon{int(round(delta_lon)):+04d}_lat{int(round(lat)):+03d}"


def antisun_lon_at(obstime_str):
    from astropy.time import Time
    from astropy.coordinates import get_sun, GeocentricTrueEcliptic
    t = Time(obstime_str, scale="utc")
    sun = get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
    return (sun.lon.deg + 180.0) % 360.0


def augment_npz(out_path, delta_lon, lat, ref_obstime):
    with np.load(out_path, allow_pickle=True) as z:
        data = {k: z[k] for k in z.files}
    data["delta_lon_from_antisun_deg"] = np.asarray(float(delta_lon), dtype=np.float64)
    data["grid_lat_deg"]               = np.asarray(float(lat),       dtype=np.float64)
    data["ref_obstime_str"]            = np.asarray(ref_obstime)
    np.savez_compressed(out_path, **data)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lon-step",      type=float, default=10.0)
    p.add_argument("--lat-base",      type=str, default=DEFAULT_LAT_BASE)
    p.add_argument("--lat-points",    type=str, default=None)
    p.add_argument("--sun-exclusion", type=float, default=40.0)
    p.add_argument("--ref-obstime",   type=str, default=DEFAULT_REF_OBSTIME)
    p.add_argument("--task-id",       type=int, default=None)
    p.add_argument("--delta-lon",     type=float, default=None)
    p.add_argument("--lat",           type=float, default=None)
    p.add_argument("--prob-maps-dir", type=str, default="prob_maps_grid_knn_s3m")
    p.add_argument("--output",        type=str, default=None)
    p.add_argument("--map-grid-file", type=str, default=None)
    p.add_argument("--list-only",     action="store_true")
    p.add_argument("--n-jobs",        type=int, default=1)
    p.add_argument("--mba-clone-factor", type=int, default=5,
                   help="MBA clone_factor. Default 5 matches the GMM grid for a fair comparison.")
    p.add_argument("--save-overlays", action="store_true")
    p.add_argument("--overwrite",     action="store_true")
    args = p.parse_args()

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
    print(f"  pipeline       = velocity_density_pipeline_fast (kNN density, K|M NEO cloner)", flush=True)
    print(f"  velocity grid  = (-1.5, 1.5) deg/day, step=0.01  [301x301]", flush=True)

    if args.map_grid_file:
        os.makedirs(os.path.dirname(os.path.abspath(args.map_grid_file)) or ".", exist_ok=True)
        with open(args.map_grid_file, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["index", "delta_lon_from_antisun_deg", "lat_deg", "filename", "map_branch"])
            for i, (dlon, lat) in enumerate(grid):
                w.writerow([i, dlon, lat, map_filename(dlon, lat), "knn_s3m"])
        print(f"  manifest -> {args.map_grid_file}", flush=True)

    if args.list_only:
        sys.exit(0)

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

    if args.output:
        out = args.output if os.path.isabs(args.output) else os.path.join(WORKDIR, args.output)
    else:
        out = os.path.join(WORKDIR, args.prob_maps_dir, map_filename(delta_lon, lat))

    if not args.overwrite and os.path.exists(out):
        print(f"[skip] already exists: {out}", flush=True)
        sys.exit(0)

    os.makedirs(os.path.dirname(out), exist_ok=True)

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

    # kNN pipeline uses s3m_loader directly — chdir to neomod/ so it can find S3Mdata/
    os.chdir(NEOMOD)
    import velocity_density_pipeline_fast as vdp  # noqa: E402

    pop_settings = copy.deepcopy(vdp.DEFAULT_POPULATION_SETTINGS)
    pop_settings["MBA"]["clone_factor"] = args.mba_clone_factor
    print(f"MBA clone_factor = {args.mba_clone_factor} "
          f"(NEO={pop_settings['NEO']['clone_factor']}, "
          f"TNO={pop_settings['TNO']['clone_factor']}, "
          f"Trojans={pop_settings['Trojans']['clone_factor']})", flush=True)

    vdp.generate_probability_maps(
        obstime_str=args.ref_obstime,
        output_path=out,
        population_settings=pop_settings,
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
