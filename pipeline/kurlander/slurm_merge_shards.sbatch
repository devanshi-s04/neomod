#!/bin/bash
#SBATCH --job-name=kurl_merge
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=00:30:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/logs/kurl_merge_%j.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/logs/kurl_merge_%j.err
#
# Merge the 40 Kurlander tracklet shards (~19 GB on disk, well over login-node memory
# once loaded/concatenated in pandas -- the direct login-node run OOM-killed at 137).

set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python

echo "start=$(date -Is)  host=$(hostname)"

"$PY" neomod/pipeline/kurlander/merge_shards.py \
    --glob '/mmfs1/gscratch/dirac/ds2004/kurlander2025/tracklets/shard_*.parquet' \
    --out  /mmfs1/gscratch/dirac/ds2004/kurlander2025/tracklets/referee.parquet

echo "end=$(date -Is)"
