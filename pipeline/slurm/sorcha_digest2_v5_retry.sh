#!/bin/bash
#SBATCH --job-name=sorcha_d2_v5r
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

# RETRY for digest2 tasks that hit TimeoutExpired on slow ckpt nodes.
# Same row mapping as sorcha_digest2_v5.sh (5000 rows per task index), but runs
# digest2 in 5 x 1000-tracklet subprocess calls instead of 1 x 5000, so each
# binary call finishes well under the 3600s timeout. Longer wall (03:00:00) too.
#
# Submit ONLY the failed indices via an explicit --array list, e.g.:
#   sbatch --array=3,4,5,6,32,...%64 sorcha_digest2_v5_retry.sh

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
CHUNK_SIZE=5000

cd "$WORKDIR"
mkdir -p logs outputs/phase2_v5/digest2_shards

ROW_START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
ROW_STOP=$(( (SLURM_ARRAY_TASK_ID + 1) * CHUNK_SIZE ))

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  rows=${ROW_START}-${ROW_STOP}  host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" run-digest2 \
  --work-dir outputs/phase2_v5 \
  --subsample-file outputs/phase2_v5/sorcha_subsample.parquet \
  --digest2-dir "$WORKDIR/digest2" \
  --digest2-chunk-tracklets 1000 \
  --digest2-timeout-sec 3600 \
  --d2-row-start "$ROW_START" \
  --d2-row-stop  "$ROW_STOP" \
  --overwrite

echo "end=$(date -Is)"
