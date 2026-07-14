#!/bin/bash
#SBATCH --job-name=sorcha_g_s3m
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

# Generate the antisun-relative sky GRID of VDP probability maps using PURE S3M
# population data (raw .s3m files via s3m_loader, NOT hybrid_elements.parquet).
# Output: prob_maps_grid_s3m/   (prob_maps_grid/ is UNTOUCHED)
#
# Identical grid to sorcha_gen_maps_grid_slurm.sh (667 maps, 29 lon x 23 lat,
# MBA clone_factor=5). Only difference: VDP_LOADER=s3m forces the pipeline to
# load orbital distributions from S3Mdata/*.s3m instead of hybrid_elements.parquet.
#
# Resubmit same array to mop up any preempted tasks (skip-existing: tasks skip
# when output .npz already exists). DO NOT add --requeue.

set -euo pipefail

export VDP_LOADER=s3m          # use pure .s3m files, not hybrid_elements.parquet

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
N_JOBS=16

mkdir -p "$WORKDIR/logs" "$WORKDIR/prob_maps_grid_s3m"

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}"
echo "host=$(hostname)  cpus=${SLURM_CPUS_PER_TASK}"
echo "VDP_LOADER=${VDP_LOADER}"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --task-id       "${SLURM_ARRAY_TASK_ID}" \
    --prob-maps-dir prob_maps_grid_s3m \
    --n-jobs        "${N_JOBS}"

echo "end=$(date -Is)"
