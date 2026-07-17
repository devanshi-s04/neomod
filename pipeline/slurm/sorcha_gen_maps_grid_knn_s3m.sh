#!/bin/bash
#SBATCH --job-name=sorcha_g_knn
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-666%48
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Generate the 667-map antisun-relative sky grid using the OLD kNN/K|M pipeline
# with PURE S3M data (s3m_loader, hardcoded in velocity_density_pipeline_fast.py).
# Output: prob_maps_grid_knn_s3m/   (prob_maps_grid_s3m/ GMM maps are UNTOUCHED)
#
# Key differences from sorcha_gen_maps_grid_s3m.sh (GMM branch):
#   - Uses velocity_density_pipeline_fast (kNN density, K|M NEO cloner)
#   - Velocity grid: (-1.5, 1.5) deg/day 301x301 (vs GMM (-2.0, 2.0) 401x401)
#   - NEO clone_factor: 300 (fast pipeline default, vs GMM 80+NEOMD3)
#   - MBA clone_factor: 5 (same as GMM branch for fair comparison)
#   - No VDP_LOADER env needed (fast pipeline always uses s3m_loader)
#
# Resubmit same array to mop up preempted tasks (skip-existing logic).
# DO NOT add --requeue.

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
N_JOBS=16

mkdir -p "$WORKDIR/logs" "$WORKDIR/prob_maps_grid_knn_s3m"

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}"
echo "host=$(hostname)  cpus=${SLURM_CPUS_PER_TASK}"
echo "pipeline=velocity_density_pipeline_fast (kNN/K|M, pure S3M)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_gen_maps_grid_knn.py" \
    --task-id       "${SLURM_ARRAY_TASK_ID}" \
    --prob-maps-dir prob_maps_grid_knn_s3m \
    --n-jobs        "${N_JOBS}" \
    --mba-clone-factor 5

echo "end=$(date -Is)"
