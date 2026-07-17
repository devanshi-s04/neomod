#!/bin/bash
#SBATCH --job-name=s3mlink_prod
#SBATCH --partition=ckpt-all
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --array=0-898%450
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# STAGE 1 — Sorcha production run WITH native SSP linking.
# Identical resources to the original multi_sorcha_production_s3m.sh
# (16 cpu, 96G, 4h, array 0-898%128). Only --config and the output dirs change.
#
# Submit per case:
#   sbatch --export=ALL,CASE=case1 1_production.sh
#   sbatch --export=ALL,CASE=case2 1_production.sh
#   sbatch --export=ALL,CASE=case3 1_production.sh
#
# 14,380,436 S3M objects / 16000 per instance = 899 instances (array 0-898).
# Speed: %450 concurrency (was %128) on ckpt-all (~10.5k idle cpu). Original
# per-task ~13 min (max 28) -> ~2 waves -> ~1-1.5 h wall.
# DO NOT add --requeue (ckpt env-retrieval hold); h5_is_good handles resubmission.
# If preempted tasks leave gaps, just re-sbatch this script (skips finished h5).

set -euo pipefail
source /mmfs1/gscratch/dirac/ds2004/sorcha/neomod/pipeline/slurm/s3m_linking/_case_env.sh

"$PY" "$WORKDIR/production_run/multi_sorcha_production.py" \
  --input_orbits   "$INPUT_ORBITS" \
  --input_physical "$INPUT_PHYS" \
  --pointings      "$POINTINGS" \
  --config         "$CONFIG" \
  --ar_data_path   "$AR_CACHE" \
  --outdir         "$PROD_DIR" \
  --workdir        "$SORCHA_WORK" \
  --chunksize 16000 \
  --norbits 1000 \
  --cores 16 \
  --instance ${SLURM_ARRAY_TASK_ID} \
  --sorcha_run "$SORCHA_RUN"

echo "end=$(date -Is)"
