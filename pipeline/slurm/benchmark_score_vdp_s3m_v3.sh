#!/bin/bash
#SBATCH --job-name=bm_vdp_s3m_v3
#SBATCH --partition=ckpt-all
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

# Score v3 (proportional-cap + busiest-night epoch) benchmark pure-S3M tracklets
# with VDP (GMM 667-map grid).
#
# Input:   outputs/benchmark_tracklets_s3m_v3/tracklets_benchmark_v3.parquet
# Maps:    prob_maps_grid_s3m/  (667 GMM maps)
# Output:  outputs/phase2_benchmark_s3m_v3/vdp_shards/vdp_00000.parquet
#
# Run AFTER benchmark_gen_tracklets_s3m_v3.sh + --combine-only completes.
# Run BEFORE benchmark_run_digest2_s3m_v3.sh.

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha

mkdir -p logs outputs/phase2_benchmark_s3m_v3

export XDG_CACHE_HOME=/tmp/${USER}/bm_vdp_s3m_v3_xdg
mkdir -p "$XDG_CACHE_HOME"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "job_id=${SLURM_JOB_ID}  host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" score-vdp \
  --tracklet-dir outputs/benchmark_tracklets_s3m_v3 \
  --prob-maps-dir prob_maps_grid_s3m \
  --work-dir outputs/phase2_benchmark_s3m_v3 \
  --batch-size 1 \
  --shard-index 0 \
  --overwrite \
  --no-nearest-dist-mask \
  --support-mask-min 1

echo ""
echo "=== VDP shard summary ==="
"$PY" -c "
import pandas as pd, numpy as np, math
df = pd.read_parquet('outputs/phase2_benchmark_s3m_v3/vdp_shards/vdp_00000.parquet',
                      columns=['population','P_NEO_vdp'])
print(f'rows: {len(df):,}')
print(f'NaN P_NEO_vdp: {df.P_NEO_vdp.isna().sum():,}')
print(df['population'].value_counts().to_string())
chunk=4000
n = math.ceil(len(df)/chunk)
print(f'digest2 array needs: 0-{n-1}  (chunk={chunk}, set in benchmark_run_digest2_s3m_v3.sh)')
"

echo "end=$(date -Is)"
