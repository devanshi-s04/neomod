#!/bin/bash
#SBATCH --job-name=s3mlink_d2
#SBATCH --partition=ckpt-all
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --array=0-129%130
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# STAGE 5 — digest2 scoring of the subsample.
#
# Speed: each 5000-row task is split into 4 parallel workers of 1250 rows via
# GNU parallel (--d2-row-start/stop), 4 cpu/task, and %130 array concurrency.
# Original single-worker 5000-row tasks took ~56 min and 56/130 TIMED OUT; the
# 4x1250 split cuts per-worker load ~4x -> ~15 min wall and far fewer timeouts.
#
# IMPORTANT: --subsample-file is passed explicitly (points at THIS case's
# subsample). This is the guard against the contamination bug where an absent
# flag silently fell back to outputs/phase2/sorcha_subsample.parquet.
#
# The array upper bound must cover ceil(subsample_rows / 5000). The default
# 0-129 matches ~648,908 rows (case1/case3). case2 (tk=3) has far fewer rows —
# set a smaller --array to avoid empty tasks. Stage 4 prints the exact size.
#
# Submit per case (after stage 4 completes), adjusting --array as needed:
#   sbatch --export=ALL,CASE=case1 --array=0-129%130 5_digest2.sh

module load gnu/parallel/20210422
set -euo pipefail
source /mmfs1/gscratch/dirac/ds2004/sorcha/neomod/pipeline/slurm/s3m_linking/_case_env.sh

OUTER_CHUNK=5000   # rows per array task
INNER_CHUNK=1250   # rows per parallel worker (4 workers per task)
BASE_ROW=$(( SLURM_ARRAY_TASK_ID * OUTER_CHUNK ))
echo "base_row=${BASE_ROW}  workers=4x${INNER_CHUNK}"

run_worker() {
  local sub_id=$1
  local row_start=$(( BASE_ROW + sub_id * INNER_CHUNK ))
  local row_stop=$(( row_start + INNER_CHUNK ))
  "$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" run-digest2 \
    --work-dir "$PHASE2_DIR" \
    --subsample-file "$SUBSAMPLE" \
    --digest2-dir "$DIGEST2_DIR" \
    --digest2-chunk-tracklets "$INNER_CHUNK" \
    --digest2-timeout-sec 3600 \
    --d2-row-start "$row_start" \
    --d2-row-stop  "$row_stop" \
    --overwrite
}
export -f run_worker
export PY WORKDIR PHASE2_DIR SUBSAMPLE DIGEST2_DIR BASE_ROW INNER_CHUNK

parallel --jobs 4 run_worker ::: 0 1 2 3

echo "end=$(date -Is)"
