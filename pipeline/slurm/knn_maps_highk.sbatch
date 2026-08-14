#!/bin/bash
#SBATCH --job-name=knnhighk
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=12:00:00
#SBATCH --array=0-5
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/knnhighk/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/knnhighk/%x_%A_%a.err
#
# INCREMENTAL extension of the more_neomod_samples_knn experiment to k = 150, 200, 250.
#   task 0..2 -> BASE k=150,200,250
#   task 3..5 -> HIGH k=150,200,250
# Reuses the frozen source_base.parquet / source_high.parquet unchanged: no new NEOMOD3 draws.
# The original eight tasks (k = 10, 25, 50, 100) are NOT rerun and NOT modified.
# N_D0_GRID stays 8000; verified accurate to 1.8e-14 at k=250 by `integration-check`.
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_knnhk_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
export NEOMOD3_CACHE_DIR="$W/outputs/neomod3_projection_cache/by_pixel"

KS=(150 200 250 150 200 250)
SRCS=(BASE BASE BASE HIGH HIGH HIGH)
K=${KS[$SLURM_ARRAY_TASK_ID]}
SRC=${SRCS[$SLURM_ARRAY_TASK_ID]}

echo "task=${SLURM_ARRAY_TASK_ID} source=${SRC} k=${K} cpus=${SLURM_CPUS_PER_TASK} host=$(hostname) start=$(date -Is)"
"$PY" "$W/neomod/pipeline/knn_sample_experiment.py" map \
    --source "${SRC}" --k "${K}" --n-jobs "${SLURM_CPUS_PER_TASK}"
echo "end=$(date -Is)"
