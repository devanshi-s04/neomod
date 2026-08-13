#!/bin/bash
#SBATCH --job-name=knnmaps
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=08:00:00
#SBATCH --array=0-7
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/knnmaps/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/knnmaps/%x_%A_%a.err
#
# The eight experimental cases: {BASE,HIGH} x k in {10,25,50,100}.
#   task 0..3 -> BASE k=10,25,50,100
#   task 4..7 -> HIGH k=10,25,50,100
# Every k within a row reads the SAME frozen source parquet, so k is the only difference in a row.
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_knnmap_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
export NEOMOD3_CACHE_DIR="$W/outputs/neomod3_projection_cache/by_pixel"

KS=(10 25 50 100 10 25 50 100)
SRCS=(BASE BASE BASE BASE HIGH HIGH HIGH HIGH)
K=${KS[$SLURM_ARRAY_TASK_ID]}
SRC=${SRCS[$SLURM_ARRAY_TASK_ID]}

echo "task=${SLURM_ARRAY_TASK_ID} source=${SRC} k=${K} cpus=${SLURM_CPUS_PER_TASK} host=$(hostname) start=$(date -Is)"
"$PY" "$W/neomod/pipeline/knn_sample_experiment.py" map \
    --source "${SRC}" --k "${K}" --n-jobs "${SLURM_CPUS_PER_TASK}"
echo "end=$(date -Is)"
