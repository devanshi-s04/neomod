#!/bin/bash
#SBATCH --job-name=bisect820v5
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --array=0-15
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
SORCHA_RUN=$WORKDIR/conda_prep/bin/sorcha-run
BISECT=$WORKDIR/work/production_2yr_v5/instance_00820/bisect
B=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

echo "start=$(date -Is)  host=$(hostname)  batch=${B}"

"$SORCHA_RUN" \
  -c neomod/pipeline/config/Rubin_full_footprint_detections.ini \
  --pd baseline_v5.0.0_2yrs.db \
  -o outputs/production_2yr_v5/bisect \
  -t "bisect820_b${B}" \
  --ob "$BISECT/orbits_b${B}.csv" \
  -p  "$BISECT/physical_b${B}.csv" \
  --ar sorcha_cache_2025-07-06 \
  -f

echo "end=$(date -Is)"
