#!/bin/bash
#SBATCH --job-name=t2neo
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=02:00:00
#SBATCH --array=0-149%148
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/t2neo/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/t2neo/%x_%A_%a.err
#
# TEST2 NEO parents: an INDEPENDENT NEOMOD3 realization, 150 shards x 2e6 = 3e8 draws.
# Seed = BM_NEO_SEED + shard = 777,000,000 + 0..149.
# Disjoint from: map GEN cache 42..141, all-sky HIGH 1,000,000..1,000,319,
#                prior (inspected) TEST 31,415,926+.
# 10x the prior TEST draw so the oversampled NEO rows are drawn from ~650k clones, not ~65k.
set -euo pipefail
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_t2neo_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
export BM_NEO_SEED=777000000
export BM_SPLIT_ROLE=TEST
export BM_OUT_DIR="$W/outputs/test2_geometric"
export BM_NEO_SHARDS="$W/outputs/test2_geometric/neo_shards"
echo "shard=${SLURM_ARRAY_TASK_ID} seed=$((777000000+SLURM_ARRAY_TASK_ID)) host=$(hostname) start=$(date -Is)"
"$W/conda_prep/bin/python" "$W/neomod/pipeline/gen_benchmark_tracklets_neomod3.py" \
    neo-shard --shard "${SLURM_ARRAY_TASK_ID}" --nshards 150 --n-orbits-total 300000000
echo "end=$(date -Is)"
