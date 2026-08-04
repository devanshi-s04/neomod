#!/bin/bash
#SBATCH --job-name=testneo
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --array=0-29
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/testneo/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/testneo/%x_%A_%a.err
# TEST NEO draw per TEST_DATA_SEAL.json: seed 31415926, 30 shards x 1M orbits.
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_tn_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 VDP_LOADER=s3m
export BM_NEO_SEED=31415926
export BM_OUT_DIR="$W/outputs/test_tracklets_neomod3"
"$PY" "$W/neomod/pipeline/gen_benchmark_tracklets_neomod3.py" neo-shard \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards 30 --n-orbits-total 30000000
