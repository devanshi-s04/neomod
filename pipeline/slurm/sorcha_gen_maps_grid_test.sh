#!/bin/bash
#SBATCH --job-name=sorcha_g_gridtest
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=333,652
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# TEST maps before launching the full 667-map grid.
#   index 333 -> (Δlon=0, lat=0)   the antisun center; should reproduce the
#                already-validated monthly map prob_maps_2026-01-01_antisun.npz
#                (center_lon=100.2). Strong correctness check.
#   index 652 -> (Δlon=0, lat=50)  a sparse high-latitude patch; confirms the
#                GMM/K|M fallback does not crash where few objects are visible.
# --save-overlays so the diagnostic plots (heliocentric x-y, velocity coverage)
# can be made from these test maps.

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
N_JOBS=16

mkdir -p "$WORKDIR/logs" "$WORKDIR/prob_maps_grid"

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --task-id       "${SLURM_ARRAY_TASK_ID}" \
    --prob-maps-dir prob_maps_grid \
    --n-jobs        "${N_JOBS}" \
    --save-overlays \
    --overwrite

echo "end=$(date -Is)"
