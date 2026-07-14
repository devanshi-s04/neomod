#!/bin/bash
#SBATCH --job-name=sorcha_d2_s3m
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --array=0-129%64
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Phase 2b — digest2 scoring for the pure-S3M run.
# Runs AFTER score-vdp + sample steps.
# 648,908 subsample rows / 5000 per chunk = 130 tasks (array 0-129).
# outputs/phase2_v5/ is UNTOUCHED.
#
# Output: outputs/phase2_s3m/digest2_shards/

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
CHUNK_SIZE=5000

cd "$WORKDIR"
mkdir -p logs outputs/phase2_s3m/digest2_shards

ROW_START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
ROW_STOP=$(( (SLURM_ARRAY_TASK_ID + 1) * CHUNK_SIZE ))

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  rows=${ROW_START}-${ROW_STOP}  host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" run-digest2 \
  --work-dir outputs/phase2_s3m \
  --subsample-file outputs/phase2_s3m/sorcha_subsample.parquet \
  --digest2-dir "$WORKDIR/digest2" \
  --digest2-chunk-tracklets 5000 \
  --digest2-timeout-sec 3600 \
  --d2-row-start "$ROW_START" \
  --d2-row-stop  "$ROW_STOP" \
  --overwrite

echo "end=$(date -Is)"
