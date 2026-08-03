#!/bin/bash
#SBATCH --job-name=bm_vdp_nm3
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# VDP-score the NEOMOD3 benchmark against the NEOMOD3 clone-only 667-map grid (docs §11.3-§11.5).
# Production scoring settings, identical to the S3M benchmark: support_mask_min=1, nearest-dist
# mask OFF.
#
# NOTE: sorcha_phase2.py globs tracklets_*.parquet, so the per-population files live in
# per_population/ -- leaving them alongside the combined file would DOUBLE-COUNT every row.
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
mkdir -p logs outputs/phase2_benchmark_neomod3
export XDG_CACHE_HOME=/tmp/${USER}/bm_vdp_nm3_xdg; mkdir -p "$XDG_CACHE_HOME"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK} NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export VDP_LOADER=s3m
echo "start=$(date -Is)"
"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" score-vdp \
  --tracklet-dir outputs/benchmark_tracklets_neomod3 \
  --prob-maps-dir prob_maps_grid_neomod3_full \
  --work-dir outputs/phase2_benchmark_neomod3 \
  --batch-size 1 --shard-index 0 --overwrite \
  --no-nearest-dist-mask --support-mask-min 1
echo "=== summary ==="
"$PY" -c "
import pandas as pd
d=pd.read_parquet('outputs/phase2_benchmark_neomod3/vdp_shards/vdp_00000.parquet',columns=['population','P_NEO_vdp'])
print(f'rows {len(d):,}   NaN {d.P_NEO_vdp.isna().sum():,}')
print(d.population.value_counts().to_string())
print('NEO prior: %.3f%%' % (100*(d.population==\"NEO\").mean()))"
echo "end=$(date -Is)"
