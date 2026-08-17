#!/bin/bash
#SBATCH --job-name=nm3allsky
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=02:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/nm3allsky/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/nm3allsky/%x_%A_%a.err
#
# One global all-sky NEOMOD3 GEN realization. STAGE is "draw" (HIGH seeds) or "rebase" (BASE shards).
#   sbatch --array=0-319%160 --export=ALL,STAGE=draw   ...
#   sbatch --array=0-99%100  --export=ALL,STAGE=rebase ...
# Each task writes its OWN shard; nothing appends concurrently to a shared Parquet file.
set -euo pipefail
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_allsky_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
STAGE="${STAGE:-draw}"
echo "stage=${STAGE} shard=${SLURM_ARRAY_TASK_ID} host=$(hostname) start=$(date -Is)"
"$W/conda_prep/bin/python" "$W/neomod/pipeline/build_neomod3_high_allsky_cache.py" \
    "${STAGE}" --shard "${SLURM_ARRAY_TASK_ID}"
echo "end=$(date -Is)"
