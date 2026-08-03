#!/bin/bash
#SBATCH --job-name=nm3_grid
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --array=0-666%200
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/nm3grid/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/nm3grid/%x_%A_%a.err
#
# FULL 667-map grid, NEOMOD3 clone-only, velocity domain +-5 deg/day.
# Gated on T1-T4 (docs/new_neomod_cloning.md §8) -- all PASS.
#
# NEO clones: NEOMOD3 projection cache, HEALPix-partitioned (nside=8) -> each task reads only the
#   ~100 pixels its 30-deg patch touches (~1.4 GB) instead of the whole 7.9 GB monolithic cache.
# MBA/TNO/Trojans: Stage-0 n-body epoch-state cache, NOT cloned (clone_factor=1).
#
# Writes to a NEW directory. prob_maps_grid_s3m_nbody/ (production) is NEVER touched.
# Idempotent: no --overwrite, so a task whose .npz already exists exits immediately.
# ckpt-all is preemptible -> just resubmit this same script to mop up killed tasks.
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_nm3grid_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
export NEOMOD3_CACHE_DIR="$W/outputs/neomod3_projection_cache/by_pixel"

echo "task=${SLURM_ARRAY_TASK_ID} host=$(hostname) start=$(date -Is)"
"$PY" "$W/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --task-id       "${SLURM_ARRAY_TASK_ID}" \
    --ref-obstime   "2027-08-25T00:00:00" \
    --cache         "$W/outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet" \
    --prob-maps-dir prob_maps_grid_neomod3_full \
    --vdp-module    velocity_density_pipeline_neomod_clone_only \
    --velocity-grid-limit 5.0 --velocity-grid-step 0.01 \
    --mba-clone-factor 1 \
    --n-jobs        16
echo "end=$(date -Is)"
