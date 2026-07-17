#!/bin/bash
#SBATCH --job-name=sorcha_vdp_s3m
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-112%32
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Phase 2a — VDP scoring for the pure-S3M run.
# Scores every tracklet in outputs/tracklets_s3m/ against its assigned
# antisun-relative grid map in prob_maps_grid_s3m/ (pure-S3M maps).
#
# Shards: 14,381 files / batch-size 128 = 113 shards (array 0-112).
# prob_maps_grid_s3m/ uses VDP_LOADER=s3m (pure .s3m files, not hybrid_elements.parquet).
# outputs/phase2_v5/ and outputs/tracklets_v5_grid/ are UNTOUCHED.
#
# Output: outputs/phase2_s3m/vdp_shards/vdp_00000.parquet ... vdp_00112.parquet

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha

cd "$WORKDIR"
mkdir -p logs outputs/phase2_s3m

export XDG_CACHE_HOME=/tmp/${USER}/sorcha_phase2_xdg_cache
mkdir -p "$XDG_CACHE_HOME"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" score-vdp \
  --tracklet-dir outputs/tracklets_s3m \
  --prob-maps-dir prob_maps_grid_s3m \
  --work-dir outputs/phase2_s3m \
  --batch-size 128 \
  --shard-index "${SLURM_ARRAY_TASK_ID}" \
  --overwrite \
  --no-nearest-dist-mask \
  --support-mask-min 1

echo "end=$(date -Is)"
