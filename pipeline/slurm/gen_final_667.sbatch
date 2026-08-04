#!/bin/bash
#SBATCH --job-name=genfinal
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=0-666%250
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/genfinal/%x_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/genfinal/%x_%A_%a.err
#
# FULL 667-MAP GEN BUILD, sealed by outputs/splits/MAP_BUILD_SEAL.json
#   sha256 c30cc29d24ead70da1111782092b37557f21d6b90f581609bf2490907a474960
# Frozen policy: step 0.01, +-5 deg/day, smoothing threshold 10 raw clones (scale_by_clone_factor
# False), GEN split with per-(pop,magbin) 1/f normalisation, symmetric unmasked densities at scoring
# time, technical-only abstention. Every task requires the GEN verification receipt.
# Idempotent: no --overwrite, so a resubmit only fills gaps. Preemptible partition -- verify with an
# open-and-read integrity scan afterwards, NOT just a file count (docs new_neomod_cloning.md 9.2).
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_gf_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
export NEOMOD3_CACHE_DIR="$W/outputs/neomod3_projection_cache/by_pixel"
echo "task=${SLURM_ARRAY_TASK_ID} host=$(hostname) start=$(date -Is)"
"$PY" "$W/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --task-id "${SLURM_ARRAY_TASK_ID}" \
    --ref-obstime "2027-08-25T00:00:00" \
    --cache "$W/outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet" \
    --split-manifest "$W/outputs/splits/nonneo_split_manifest.parquet" --split-role GEN \
    --split-provenance "$W/outputs/splits/split_provenance.json" --require-receipt \
    --smooth-support-threshold 10 \
    --prob-maps-dir prob_maps_grid_neomod3_GEN_final \
    --vdp-module velocity_density_pipeline_neomod_clone_only \
    --velocity-grid-limit 5.0 --velocity-grid-step 0.01 \
    --mba-clone-factor 1 --n-jobs 16
echo "end=$(date -Is)"
