"""
Generate one VDP probability map for a given epoch and sky center.

Usage:
    python sorcha_gen_map.py \
        --obstime 2025-05-01T00:00:00 \
        --center-lon 21.5 --center-lat 0.0 \
        --center-label antisun \
        --output /path/to/prob_maps/prob_maps_2025-05-01_antisun.npz
"""
import argparse
import os
import sys

WORKDIR = "/mmfs1/gscratch/dirac/ds2004/sorcha"
NEOMOD  = os.path.join(WORKDIR, "neomod")

sys.path.insert(0, os.path.join(NEOMOD, "adam_core_stub"))
sys.path.insert(0, os.path.join(NEOMOD, "src"))
# s3m_loader searches [".", "S3Mdata"] relative to CWD — must chdir to neomod/
os.chdir(NEOMOD)

import velocity_density_pipeline_hybrid as vdp  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--obstime",      required=True,
                   help="ISO obstime, e.g. 2025-05-01T00:00:00")
    p.add_argument("--center-lon",   type=float, required=True,
                   help="Sky-patch center ecliptic longitude (deg)")
    p.add_argument("--center-lat",   type=float, default=0.0,
                   help="Sky-patch center ecliptic latitude (deg)")
    p.add_argument("--center-label", required=True,
                   help="Short label stored in .npz metadata")
    p.add_argument("--output",       required=True,
                   help="Output .npz path")
    p.add_argument("--overwrite",    action="store_true")
    p.add_argument("--n-jobs",       type=int, default=1,
                   help="Parallel workers for density estimation (default 1)")
    args = p.parse_args()

    # make output path absolute before the chdir takes effect
    out = args.output if os.path.isabs(args.output) \
          else os.path.join(WORKDIR, args.output)

    if not args.overwrite and os.path.exists(out):
        print(f"[skip] already exists: {out}", flush=True)
        sys.exit(0)

    os.makedirs(os.path.dirname(out), exist_ok=True)

    print(f"obstime    = {args.obstime}",                              flush=True)
    print(f"center     = ({args.center_lon:.4f}, {args.center_lat})", flush=True)
    print(f"label      = {args.center_label}",                         flush=True)
    print(f"output     = {out}",                                        flush=True)
    sys.stdout.flush()

    vdp.generate_probability_maps(
        obstime_str=args.obstime,
        output_path=out,
        center_lon_deg=args.center_lon,
        center_lat_deg=args.center_lat,
        center_label=args.center_label,
        n_jobs=args.n_jobs,
    )
    print(f"[done] {out}", flush=True)


if __name__ == "__main__":
    main()
