#!/bin/bash
#SBATCH --job-name=e0prod
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --array=0-63
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/e0prod/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/e0prod/%x_%A_%a.err
#
# E0 candidate products: 16 preregistered centers x support thresholds {2,3,5,10} raw clones.
# ONLY --smooth-support-threshold varies. Grid, sigma, truncation, passes, normalisation, split,
# k, and the GEN artifact are identical across all 64. Every task requires a PASSing GEN
# verification receipt tied to the live manifest.
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_e0p_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
export NEOMOD3_CACHE_DIR="$W/outputs/neomod3_projection_cache/by_pixel"

DLONS=(  0  20  10 -10   0  50  50 -50  90 140 -140   0   0   0   0 -50)
LATS=(   0 -12   2  -2   8  -1   1   1   0   0    0  25 -25  50 -50 -50)
THRS=(2 3 5 10)
CI=$(( SLURM_ARRAY_TASK_ID % 16 ))
TI=$(( SLURM_ARRAY_TASK_ID / 16 ))
D=${DLONS[$CI]}; L=${LATS[$CI]}; T=${THRS[$TI]}
echo "task=$SLURM_ARRAY_TASK_ID center=($D,$L) threshold=$T start=$(date -Is)"
"$PY" "$W/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --delta-lon "$D" --lat "$L" \
    --ref-obstime "2027-08-25T00:00:00" \
    --cache "$W/outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet" \
    --split-manifest "$W/outputs/splits/nonneo_split_manifest.parquet" \
    --split-role GEN \
    --split-provenance "$W/outputs/splits/split_provenance.json" \
    --require-receipt \
    --smooth-support-threshold "$T" \
    --prob-maps-dir "prob_maps_e0_thr${T}" \
    --vdp-module velocity_density_pipeline_neomod_clone_only \
    --velocity-grid-limit 5.0 --velocity-grid-step 0.01 \
    --mba-clone-factor 1 --n-jobs 16 --overwrite
echo "end=$(date -Is)"
