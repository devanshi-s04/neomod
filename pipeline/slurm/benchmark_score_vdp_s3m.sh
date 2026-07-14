#!/bin/bash
#SBATCH --job-name=bm_vdp_s3m
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Score benchmark pure-S3M tracklets with VDP (GMM 667-map grid).
#
# Input:   outputs/benchmark_tracklets_s3m/tracklets_benchmark.parquet
# Maps:    prob_maps_grid_s3m/  (667 GMM maps)
# Output:  outputs/phase2_benchmark_s3m/vdp_shards/vdp_00000.parquet
#
# Runs as a single job (no array): the one tracklet file is one shard.
# Loads all required maps sequentially from prob_maps_grid_s3m/ —
# expect ~30-60 min of map I/O for ~200-300 unique cells in the parquet.
#
# Flags:
#   --no-nearest-dist-mask  GMM maps do not need kNN nearest-clone mask
#   --support-mask-min 1    zero P(NEO) where in-cell clone count < 1
#
# Run AFTER benchmark_gen_tracklets_s3m.sh completes.
# Run BEFORE benchmark_run_digest2_s3m.sh.

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha

mkdir -p logs outputs/phase2_benchmark_s3m

export XDG_CACHE_HOME=/tmp/${USER}/bm_vdp_s3m_xdg
mkdir -p "$XDG_CACHE_HOME"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "job_id=${SLURM_JOB_ID}  host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" score-vdp \
  --tracklet-dir outputs/benchmark_tracklets_s3m \
  --prob-maps-dir prob_maps_grid_s3m \
  --work-dir outputs/phase2_benchmark_s3m \
  --batch-size 1 \
  --shard-index 0 \
  --overwrite \
  --no-nearest-dist-mask \
  --support-mask-min 1

echo ""
echo "=== VDP shard summary ==="
"$PY" -c "
import pandas as pd, numpy as np
df = pd.read_parquet('outputs/phase2_benchmark_s3m/vdp_shards/vdp_00000.parquet',
                      columns=['population','P_NEO_vdp'])
print(f'rows: {len(df):,}')
print(f'NaN P_NEO_vdp: {df.P_NEO_vdp.isna().sum():,}')
is_neo = df.population == 'NEO'
print(f'NEO: {is_neo.sum():,}  non-NEO: {(~is_neo).sum():,}')
import math; chunk=5000
n = math.ceil(len(df)/chunk)
print(f'digest2 array needs: 0-{n-1}  (chunk={chunk}, set in benchmark_run_digest2_s3m.sh)')
"

echo "end=$(date -Is)"
