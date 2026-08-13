#!/bin/bash
#SBATCH --job-name=knnhigh
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=02:00:00
#SBATCH --array=0-319%160
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/knnhigh/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/knnhigh/%x_%A_%a.err
#
# Additional independent NEOMOD3 draws for the HIGH source case.
# Seeds are HIGH_SEED_BASE + shard = 1,000,000 + 0..319, disjoint from the original cache's 42..141.
# Each shard draws, projects all-sky exactly as neomod3_projection_cache.build does, casts to
# float32 as the cache does, then applies the production magnitude + 30-deg sky selection and keeps
# ONLY the survivors -- so 6.4e8 draws never materialise as rows on disk.
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_knnhigh_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
export NEOMOD3_CACHE_DIR="$W/outputs/neomod3_projection_cache/by_pixel"
echo "shard=${SLURM_ARRAY_TASK_ID} host=$(hostname) start=$(date -Is)"
"$PY" "$W/neomod/pipeline/knn_sample_experiment.py" draw-high --shard "${SLURM_ARRAY_TASK_ID}"
echo "end=$(date -Is)"
