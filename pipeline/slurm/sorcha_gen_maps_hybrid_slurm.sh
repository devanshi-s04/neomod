#!/bin/bash
#SBATCH --job-name=sorcha_g_hybrid
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-23%24
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Generate 24 VDP antisun probability maps: one per month covering the
# full 2-year Sorcha simulation (2025-05 through 2027-04).
#
# Why antisun only (not the off-axis patches):
#   The antisun direction rotates 360° over the year, so 24 monthly maps
#   cover every ecliptic longitude at the correct epoch.  Off-axis patches
#   at ±45/90/120° from antisun performed poorly (epoch-matched F1 ≈ 0.39–0.53)
#   and are not worth generating.
#
# Speed-ups vs original:
#   1. velocity_density_pipeline_fast.py — vectorized Bayesian kNN integral
#      (no Python double-loop; numpy broadcasting + scipy.special.gammaln)
#   2. Apparent magnitudes pre-computed once per population per map,
#      not redundantly for each of the 8 mag bins
#   3. --cpus-per-task=8 passed as --n-jobs to joblib density evaluator
#
# Expected wall time: ~20–40 min per map  (was ~2 h with original code)

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
N_JOBS=16

# 24 monthly epochs spanning the simulation window
MONTHS=(
  "2025-05-01T00:00:00"
  "2025-06-01T00:00:00"
  "2025-07-01T00:00:00"
  "2025-08-01T00:00:00"
  "2025-09-01T00:00:00"
  "2025-10-01T00:00:00"
  "2025-11-01T00:00:00"
  "2025-12-01T00:00:00"
  "2026-01-01T00:00:00"
  "2026-02-01T00:00:00"
  "2026-03-01T00:00:00"
  "2026-04-01T00:00:00"
  "2026-05-01T00:00:00"
  "2026-06-01T00:00:00"
  "2026-07-01T00:00:00"
  "2026-08-01T00:00:00"
  "2026-09-01T00:00:00"
  "2026-10-01T00:00:00"
  "2026-11-01T00:00:00"
  "2026-12-01T00:00:00"
  "2027-01-01T00:00:00"
  "2027-02-01T00:00:00"
  "2027-03-01T00:00:00"
  "2027-04-01T00:00:00"
)

OBSTIME=${MONTHS[$SLURM_ARRAY_TASK_ID]}
DATE=${OBSTIME%%T*}

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}"
echo "obstime=${OBSTIME}"
echo "host=$(hostname)  cpus=${SLURM_CPUS_PER_TASK}"
echo "start=$(date -Is)"

mkdir -p "$WORKDIR/logs" "$WORKDIR/prob_maps_hybrid"

# Compute antisun ecliptic longitude at this epoch
ANTISUN_LON=$("$PY" -c "
from astropy.time import Time
from astropy.coordinates import get_sun, GeocentricTrueEcliptic
t = Time('${OBSTIME}', scale='utc')
sun = get_sun(t).transform_to(GeocentricTrueEcliptic(obstime=t))
print(f'{(sun.lon.deg + 180.0) % 360.0:.6f}')
")

echo "antisun_lon=${ANTISUN_LON}"

OUTPUT="${WORKDIR}/prob_maps_hybrid/prob_maps_${DATE}_antisun.npz"
echo "output=${OUTPUT}"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_gen_map_hybrid.py" \
    --obstime    "${OBSTIME}"    \
    --center-lon "${ANTISUN_LON}" \
    --center-lat 0.0             \
    --center-label "antisun"     \
    --output     "${OUTPUT}"     \
    --n-jobs     "${N_JOBS}"     \
    --overwrite

echo "end=$(date -Is)"
