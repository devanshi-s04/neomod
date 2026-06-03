#!/bin/bash
#SBATCH --job-name=sorcha_vdp
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-112%32
#SBATCH --chdir=/mmfs1/gscratch/astro/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

PY=/mmfs1/gscratch/astro/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/astro/ds2004/sorcha

cd "$WORKDIR"
mkdir -p logs outputs/phase2

export XDG_CACHE_HOME=/tmp/${USER}/sorcha_phase2_xdg_cache
mkdir -p "$XDG_CACHE_HOME"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "job_id=${SLURM_JOB_ID}"
echo "array_job_id=${SLURM_ARRAY_JOB_ID}"
echo "array_task_id=${SLURM_ARRAY_TASK_ID}"
echo "host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" score-vdp \
  --tracklet-dir outputs/tracklets \
  --prob-maps-dir prob_maps \
  --work-dir outputs/phase2 \
  --batch-size 128 \
  --shard-index "${SLURM_ARRAY_TASK_ID}" \
  --overwrite

echo "end=$(date -Is)"
