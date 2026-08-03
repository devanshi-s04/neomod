#!/bin/bash
#SBATCH --job-name=bm_nm3
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --array=0-9
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%A_%a.err
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_bmnm3_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
"$PY" "$W/neomod/pipeline/gen_benchmark_tracklets_neomod3.py" neo-shard \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards 10 --n-orbits-total 10000000
