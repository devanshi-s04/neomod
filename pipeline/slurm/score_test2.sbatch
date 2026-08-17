#!/bin/bash
#SBATCH --job-name=t2score
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/t2score/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/t2score/%x_%A_%a.err
set -euo pipefail
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_t2s_${SLURM_ARRAY_TASK_ID:-0}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
MODE="${MODE:?set MODE}"
NSH="${NSH:?set NSH}"
echo "mode=${MODE} shard=${SLURM_ARRAY_TASK_ID} nshards=${NSH} host=$(hostname) start=$(date -Is)"
"$W/conda_prep/bin/python" "$W/neomod/pipeline/score_test2.py" \
    --mode "${MODE}" --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${NSH}"
echo "end=$(date -Is)"
