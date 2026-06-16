#!/bin/bash
#SBATCH --job-name=sorcha_d2_retry
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-35%36
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Retry for the 36 chunks that timed out in the first run.
# Uses 1000-tracklet sub-chunks (5 per task) so no single digest2 call
# can hit the timeout, even on the slowest ckpt nodes.

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
CHUNK_SIZE=5000

# Ordered list of the 36 missing original-task IDs
MISSING_TASKS=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 33 41 51 59 60 61 62 63 65 93 98 114)

cd "$WORKDIR"
mkdir -p logs outputs/phase2/digest2_shards

ORIG_TASK=${MISSING_TASKS[$SLURM_ARRAY_TASK_ID]}
ROW_START=$(( ORIG_TASK * CHUNK_SIZE ))
ROW_STOP=$(( (ORIG_TASK + 1) * CHUNK_SIZE ))

echo "job_id=${SLURM_JOB_ID}"
echo "retry_array_task=${SLURM_ARRAY_TASK_ID}  orig_task=${ORIG_TASK}"
echo "rows=${ROW_START}-${ROW_STOP}"
echo "host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" run-digest2 \
  --work-dir outputs/phase2 \
  --digest2-dir "$WORKDIR/digest2" \
  --digest2-chunk-tracklets 1000 \
  --digest2-timeout-sec 1800 \
  --d2-row-start "$ROW_START" \
  --d2-row-stop  "$ROW_STOP" \
  --overwrite

echo "end=$(date -Is)"
