#!/bin/bash
#SBATCH --job-name=sorcha_s3m
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --array=0-898%128
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# PURE-S3M production run (original 14,380,436 S3M objects, NO MPCORB), v5.0.0 cadence.
# Writes to NEW directories outputs/production_2yr_s3m / work/production_2yr_s3m so the
# existing hybrid run (outputs/production_2yr_v5) is UNTOUCHED.
#
# 14,380,436 objects / 16000 per instance = 899 instances (array 0-898; instance 898
# reads the final ~12.4k objects, 899+ would be empty). Same config / DB / cache /
# chunking as the hybrid run.
#
# DO NOT add --requeue (ckpt env-retrieval hold). Skip-existing (h5_is_good) handles
# resubmission. NB: S3M epoch is ~2007 (MJD 54466) vs the hybrid ~2023, so propagation
# to the 2025-2027 window is longer -> watch for a few integrator-hang instances
# (recover them the same way as inst00820 in the hybrid run).

set -euo pipefail

WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
PY=$WORKDIR/conda_prep/bin/python
SORCHA_RUN=$WORKDIR/conda_prep/bin/sorcha-run

cd "$WORKDIR"
mkdir -p outputs/production_2yr_s3m work/production_2yr_s3m logs

"$PY" "$WORKDIR/production_run/multi_sorcha_production.py" \
  --input_orbits inputs/s3m_sorcha_orbits.csv \
  --input_physical inputs/s3m_sorcha_phys.csv \
  --pointings baseline_v5.0.0_2yrs.db \
  --config neomod/pipeline/config/Rubin_full_footprint_detections.ini \
  --ar_data_path sorcha_cache_2025-07-06 \
  --outdir outputs/production_2yr_s3m \
  --workdir work/production_2yr_s3m \
  --chunksize 16000 \
  --norbits 1000 \
  --cores 16 \
  --instance ${SLURM_ARRAY_TASK_ID} \
  --sorcha_run $SORCHA_RUN
