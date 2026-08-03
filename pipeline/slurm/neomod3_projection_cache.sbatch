#!/bin/bash
#SBATCH --job-name=nm3_cache
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:30:00
#SBATCH --array=0-39
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%A_%a.err
#
# NEOMOD3 projection cache: sample + project ONCE, all-sky, at the map epoch.
# Each of the 667 map centers later slices this same cache via the existing
# build_visible_subset_dataframe cache branch (fixing_integrator.md §12.2).
#
# Default: 40 shards x 1,000,000 orbits = 40,000,000 total.
# Override with --export=ALL,N_TOTAL=...,NSHARDS=... (keep --array matching NSHARDS-1).
#
#   sbatch neomod/pipeline/slurm/neomod3_projection_cache.sbatch
#   # then, once all shards are done:
#   conda_prep/bin/python neomod/pipeline/neomod3_projection_cache.py merge

set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
N_TOTAL="${N_TOTAL:-40000000}"
NSHARDS="${NSHARDS:-40}"
EPOCH="${EPOCH:-2027-08-25T00:00:00}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="/tmp/${USER}/mpl_nm3cache"; mkdir -p "$MPLCONFIGDIR"
export VDP_LOADER=s3m   # scorer source; pure-S3M (t_0=2008), never the hybrid catalog

echo "shard=${SLURM_ARRAY_TASK_ID}/${NSHARDS} n_total=${N_TOTAL} epoch=${EPOCH} start=$(date -Is)"
"$PY" neomod/pipeline/neomod3_projection_cache.py build \
    --shard "${SLURM_ARRAY_TASK_ID}" \
    --nshards "${NSHARDS}" \
    --n-orbits-total "${N_TOTAL}" \
    --epoch "${EPOCH}"
echo "end=$(date -Is)"
