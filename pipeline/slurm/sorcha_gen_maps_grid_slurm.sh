#!/bin/bash
#SBATCH --job-name=sorcha_g_grid
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-666%32
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Generate the antisun-relative sky GRID of VDP probability maps (667 maps).
#
# Grid (see sorcha_gen_maps_grid.py for the configurable defaults):
#   longitude : 10 deg steps, antisun-relative, 40 deg sun exclusion -> 29 usable
#   latitude  : NON-UNIFORM, fine 1 deg near the ecliptic, coarse at high lat:
#               0, +/-1, +/-2, +/-3, +/-4, +/-5, +/-8, +/-12, +/-18, +/-25, +/-35, +/-50
#               (23 values)  ->  29 x 23 = 667 maps
#   ref epoch : 2026-01-01T00:00:00 (maps are antisun-relative / time-independent)
#
# Each task generates ONE map by its grid index (--task-id $SLURM_ARRAY_TASK_ID).
# The grid is rebuilt deterministically inside the script, so no manifest file is
# needed at runtime (one was written to slurm/grid_map_manifest.csv for reference).
#
# NO --overwrite: on ckpt preemption a re-submit skips already-finished maps
# (the script skips when the output .npz already exists).  Resubmit the same
# array to mop up any preempted tasks.
#
# Output: prob_maps_grid/prob_maps_grid_dlon{+DDD}_lat{+DD}.npz
# Expected wall time: ~13-40 min per map; ~667/32 batches -> ~7-8 h total.

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
N_JOBS=16

mkdir -p "$WORKDIR/logs" "$WORKDIR/prob_maps_grid"

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}"
echo "host=$(hostname)  cpus=${SLURM_CPUS_PER_TASK}"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --task-id       "${SLURM_ARRAY_TASK_ID}" \
    --prob-maps-dir prob_maps_grid \
    --n-jobs        "${N_JOBS}"

echo "end=$(date -Is)"
