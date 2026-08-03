#!/bin/bash
#SBATCH --job-name=vlim5_pilot
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --array=0-28
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# STAGE D: rebuild the 29 map centers touched by the out-of-±2 night-61642 NEOs at velocity
# limit +/-5.0 deg/day (0.01 step -> 1001x1001), into a NEW dir. Does NOT touch the production
# +/-2.0 grid (prob_maps_grid_s3m_nbody/). Same cache/epoch/clone-factor as production; only the
# velocity domain changes. One center per array task, read from vlim5_pilot_centers.txt.

set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/codex_integrator:$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_vlim5"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$W/neomod/pipeline/slurm/vlim5_pilot_centers.txt")
DLON=$(echo "$LINE" | awk '{print $1}')
LAT=$(echo "$LINE" | awk '{print $2}')
echo "task=${SLURM_ARRAY_TASK_ID} dlon=${DLON} lat=${LAT} start=$(date -Is) host=$(hostname)"

"$PY" "$W/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --delta-lon "${DLON}" --lat "${LAT}" \
    --ref-obstime "2027-08-25T00:00:00" \
    --cache "$W/outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet" \
    --prob-maps-dir prob_maps_grid_s3m_nbody_vlim5_pilot \
    --mba-clone-factor 5 \
    --velocity-grid-limit 5.0 \
    --velocity-grid-step 0.01 \
    --n-jobs 8 \
    --overwrite
echo "end=$(date -Is)"
