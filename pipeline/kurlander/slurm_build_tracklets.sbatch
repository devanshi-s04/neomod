#!/bin/bash
#SBATCH --job-name=kurl_tracklets
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --array=0-39
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/logs/kurl_tracklets_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/logs/kurl_tracklets_%A_%a.err
#
# FULL Kurlander referee tracklet build (all 2778 s3m chunks + NEO + Trojan + TNO).
# 40-way array; each task takes every 40th file. Vectorised builder ~1.9 s/file
# => ~70 files/task * 2 s ~= 3 min/task (plus HDF read). Generous 4 h wall for safety.
#
# Launch:
#   mkdir -p /mmfs1/gscratch/dirac/ds2004/logs
#   sbatch pipeline/kurlander/slurm_build_tracklets.sbatch
# Then merge:
#   python pipeline/kurlander/merge_shards.py \
#       --glob '/mmfs1/gscratch/dirac/ds2004/kurlander2025/tracklets/shard_*.parquet' \
#       --out  /mmfs1/gscratch/dirac/ds2004/kurlander2025/tracklets/all_tracklets.parquet

set -e
ROOT=/mmfs1/gscratch/dirac/ds2004/kurlander2025
OUTDIR=$ROOT/tracklets
mkdir -p "$OUTDIR"

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python

srun $PY neomod/pipeline/kurlander/build_tracklets.py \
    --root "$ROOT" \
    --pops neomod s3m trojanmod cfeps \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards 40 \
    --out "$OUTDIR/shard_${SLURM_ARRAY_TASK_ID}.parquet"

echo "shard $SLURM_ARRAY_TASK_ID done"
