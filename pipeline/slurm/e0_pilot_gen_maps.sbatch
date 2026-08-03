#!/bin/bash
#SBATCH --job-name=e0pilot
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --array=0-15
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/e0pilot/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/e0pilot/%x_%A_%a.err
#
# E0 PILOT -- 16 preregistered centers (docs/E0_PILOT_PREREGISTRATION.md §1, sha 2608c472c04f2d81).
# Maps built from the GEN split ONLY; non-NEO density corrected by the exact per-(population,
# magbin) retained fraction so it represents the FULL population (protocol v1.2 §0.2).
# Gates the full 667-map rebuild. TEST is not touched.
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_e0_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
export NEOMOD3_CACHE_DIR="$W/outputs/neomod3_projection_cache/by_pixel"

DLONS=(  0  20  10 -10   0  50  50 -50  90 140 -140   0   0   0   0 -50)
LATS=(   0 -12   2  -2   8  -1   1   1   0   0    0  25 -25  50 -50 -50)
D=${DLONS[$SLURM_ARRAY_TASK_ID]}
L=${LATS[$SLURM_ARRAY_TASK_ID]}
echo "task=$SLURM_ARRAY_TASK_ID dlon=$D lat=$L start=$(date -Is)"
"$PY" "$W/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --delta-lon "$D" --lat "$L" \
    --ref-obstime "2027-08-25T00:00:00" \
    --cache "$W/outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet" \
    --split-manifest "$W/outputs/splits/nonneo_split_manifest.parquet" \
    --split-role GEN \
    --split-provenance "$W/outputs/splits/split_provenance.json" \
    --prob-maps-dir prob_maps_grid_neomod3_GEN \
    --vdp-module velocity_density_pipeline_neomod_clone_only \
    --velocity-grid-limit 5.0 --velocity-grid-step 0.01 \
    --mba-clone-factor 1 --n-jobs 16 --overwrite
echo "end=$(date -Is)"
