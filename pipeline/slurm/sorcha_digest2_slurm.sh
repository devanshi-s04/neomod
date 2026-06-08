#!/bin/bash
#SBATCH --job-name=sorcha_d2
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --array=0-141%64
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# 598,670 rows / 5,000 per chunk = 120 chunks (array tasks 0-119)
# Each task processes rows [TASK*5000, (TASK+1)*5000) from sorcha_subsample.parquet
# ~28 min per chunk; 64-way concurrency → ~30 min total wall time

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
CHUNK_SIZE=5000

cd "$WORKDIR"
mkdir -p logs outputs/phase2/digest2_shards

ROW_START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
ROW_STOP=$(( (SLURM_ARRAY_TASK_ID + 1) * CHUNK_SIZE ))

echo "job_id=${SLURM_JOB_ID}"
echo "array_task=${SLURM_ARRAY_TASK_ID}"
echo "rows=${ROW_START}-${ROW_STOP}"
echo "host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" run-digest2 \
  --work-dir outputs/phase2 \
  --digest2-dir "$WORKDIR/digest2" \
  --digest2-chunk-tracklets 5000 \
  --digest2-timeout-sec 3600 \
  --d2-row-start "$ROW_START" \
  --d2-row-stop  "$ROW_STOP" \
  --overwrite

echo "end=$(date -Is)"
