#!/bin/bash
#SBATCH --job-name=sorcha_d2_s3mr
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# RETRY for S3M digest2 tasks that timed out or failed.
# Same row mapping as sorcha_digest2_s3m.sh (5000 rows per task index), but
# digest2-chunk-tracklets=1000 (5 subprocess calls instead of 1) to avoid
# TimeoutExpired on slow ckpt nodes. Longer wall (03:00:00) too.
#
# Submit ONLY the failed indices:
#   sbatch --array=6,7,8,9,10,11,12,13,18,19,20,21,22,23,24,25,33,34,35,36,37,38,39,40,55,56,57,58,59,60,61,62,63,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,113,114,115,120,121,122,123,124%64 sorcha_digest2_s3m_retry.sh

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
  --digest2-chunk-tracklets 1000 \
  --digest2-timeout-sec 3600 \
  --d2-row-start "$ROW_START" \
  --d2-row-stop  "$ROW_STOP" \
  --overwrite

echo "end=$(date -Is)"
