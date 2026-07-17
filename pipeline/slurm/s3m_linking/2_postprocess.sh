#!/bin/bash
#SBATCH --job-name=s3mlink_post
#SBATCH --partition=ckpt-all
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=0-899%450
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# STAGE 2 — build nightly two-detection tracklets from the linked Sorcha .h5.
# Same worker (sorcha_postprocess.py) as the original. Each array task processes a
# contiguous block of FILES_PER_TASK files from the global sorted inst*_part*.h5
# list, preserving the exact integer --file_index (and source_file_index).
#
# Speed: 8 cpu + GNU parallel --jobs 8 fan-out over the block (files are
# independent), and %450 array concurrency. Original per-file ~2 min (max 9) ->
# ~10-15 min wall.
#
# Footprint/map assignment uses prob_maps_grid (geometry only) exactly as the
# original sorcha_postprocess_s3m.sh did.
#
# Submit per case (after stage 1 completes):
#   sbatch --export=ALL,CASE=case1 2_postprocess.sh

module load gnu/parallel/20210422
set -euo pipefail
source /mmfs1/gscratch/dirac/ds2004/sorcha/neomod/pipeline/slurm/s3m_linking/_case_env.sh

export XDG_CACHE_HOME=/tmp/${USER}/s3mlink_post_xdg
mkdir -p "$XDG_CACHE_HOME"

FILES_PER_TASK=16
NFILES=$(ls "$PROD_DIR"/inst*_part*.h5 2>/dev/null | wc -l)
echo "NFILES=$NFILES  FILES_PER_TASK=$FILES_PER_TASK  jobs=8"

START=$(( SLURM_ARRAY_TASK_ID * FILES_PER_TASK ))

# One worker = one file_index. Independent outputs, safe to run in parallel.
run_one() {
  local IDX=$1
  echo "--- file_index $IDX ---"
  "$PY" "$WORKDIR/neomod/pipeline/sorcha_postprocess.py" \
    --file_index "$IDX" \
    --indir "$PROD_DIR" \
    --outdir "$TRACK_DIR" \
    --prob_maps_dir "$PROB_MAPS_FOOTPRINT" \
    --skip-existing
}
export -f run_one
export PY WORKDIR PROD_DIR TRACK_DIR PROB_MAPS_FOOTPRINT

# Build the list of in-range indices for this task, then fan out 8-wide.
INDICES=()
for (( k=0; k<FILES_PER_TASK; k++ )); do
  IDX=$(( START + k ))
  [ "$IDX" -ge "$NFILES" ] && break
  INDICES+=( "$IDX" )
done

if [ "${#INDICES[@]}" -gt 0 ]; then
  printf '%s\n' "${INDICES[@]}" | parallel --jobs 8 run_one
fi

echo "end=$(date -Is)"
