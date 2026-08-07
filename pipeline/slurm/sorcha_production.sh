#!/bin/bash
#SBATCH --job-name=sorchaprod
#SBATCH --partition=ckpt-all
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=04:00:00
#SBATCH --array=0-319%320
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
W=/mmfs1/gscratch/dirac/ds2004/sorcha
I=$(printf "%04d" "$SLURM_ARRAY_TASK_ID")
# FROZEN per SORCHA_TEST_PREPRODUCTION_SEAL.json: SORCHA_SEED = 20260806 + task_id
export SORCHA_SEED=$((20260806 + SLURM_ARRAY_TASK_ID))
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}
OUT=$W/outputs/sorcha_production/production
# idempotent: a completed chunk is never redone, so the array can be resubmitted to fill gaps
if [ -s "$OUT/prod_$I.h5" ]; then echo "chunk $I already complete"; exit 0; fi
echo "chunk $I  SORCHA_SEED=$SORCHA_SEED  cpus=$SLURM_CPUS_PER_TASK"
$W/conda_prep/bin/sorcha-run \
  -c $W/neomod/pipeline/config/neomod3_test_detections.ini \
  --pd $W/baseline_v5.0.0_2yrs.db \
  -o $OUT -t prod_$I \
  --ob $W/outputs/sorcha_production/work/chunk_$I/orbits.csv \
  -p  $W/outputs/sorcha_production/work/chunk_$I/physical.csv \
  --ar $W/sorcha_cache_2025-07-06 -f
echo "chunk $I DONE"
