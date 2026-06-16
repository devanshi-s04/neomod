#!/bin/bash
#SBATCH --job-name=sorcha_post_v5g
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-999%64
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Phase 1 re-run for v5.0 using the antisun-relative GRID maps (prob_maps_grid/).
# Identical to sorcha_postprocess_v5.sh except:
#   --prob_maps_dir prob_maps_grid   (grid maps instead of 24 monthly antisun maps)
#   --outdir outputs/tracklets_v5_grid
# sorcha_postprocess.py auto-detects grid maps (delta_lon_from_antisun_deg key)
# and assigns each tracklet to the nearest (Δlon, lat) grid center.
#
# Submit in batches of 1000 under the 2000-task QOS limit (14,445 files total):
#   sbatch --export=ALL,INDEX_OFFSET=0     --array=0-999%64 sorcha_postprocess_v5_grid.sh
#   sbatch --export=ALL,INDEX_OFFSET=1000  --array=0-999%64 sorcha_postprocess_v5_grid.sh
#   ... through INDEX_OFFSET=14000 --array=0-444%64
#
# DO NOT add --requeue: causes "user env retrieval failed requeued held" on ckpt.

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha

cd "$WORKDIR"

mkdir -p logs outputs/tracklets_v5_grid

export XDG_CACHE_HOME=/tmp/${USER}/sorcha_postprocess_xdg_cache
mkdir -p "$XDG_CACHE_HOME"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "job_id=${SLURM_JOB_ID}"
echo "array_job_id=${SLURM_ARRAY_JOB_ID}"
echo "array_task_id=${SLURM_ARRAY_TASK_ID}"
INDEX_OFFSET=${INDEX_OFFSET:-0}
FILE_INDEX=$((SLURM_ARRAY_TASK_ID + INDEX_OFFSET))
echo "index_offset=${INDEX_OFFSET}"
echo "file_index=${FILE_INDEX}"
echo "host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_postprocess.py" \
  --file_index "${FILE_INDEX}" \
  --indir outputs/production_2yr_v5 \
  --outdir outputs/tracklets_v5_grid \
  --prob_maps_dir prob_maps_grid \
  --overwrite

echo "end=$(date -Is)"
