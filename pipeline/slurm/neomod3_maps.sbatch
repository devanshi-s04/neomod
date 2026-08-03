#!/bin/bash
#SBATCH --job-name=nm3_map
#SBATCH --account=astro
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=03:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%A_%a.err
#SBATCH --array=0-1
#
# NEOMOD3 clone-only maps at +-5 deg/day, for the two test centers.
#   task 0 -> dlon +000 lat +00  (antisun; the geometry of interest)
#   task 1 -> dlon +020 lat -12  (busiest center; statistical power for the ROC)
# NEO clones come from the NEOMOD3 projection cache with absolute normalisation; MBA/TNO/Trojans
# come from the Stage-0 n-body cache and are NOT cloned. Writes to a NEW directory.
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_nm3map"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m

case "$SLURM_ARRAY_TASK_ID" in
  0) DLON=0;  LAT=0  ;;
  1) DLON=20; LAT=-12 ;;
esac
echo "task=$SLURM_ARRAY_TASK_ID dlon=$DLON lat=$LAT start=$(date -Is)"
"$PY" "$W/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --delta-lon "$DLON" --lat "$LAT" \
    --ref-obstime "2027-08-25T00:00:00" \
    --cache "$W/outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet" \
    --prob-maps-dir prob_maps_grid_neomod3_vlim5 \
    --vdp-module velocity_density_pipeline_neomod_clone_only \
    --velocity-grid-limit 5.0 --velocity-grid-step 0.01 \
    --mba-clone-factor 1 \
    --n-jobs 16 --overwrite
echo "end=$(date -Is)"
