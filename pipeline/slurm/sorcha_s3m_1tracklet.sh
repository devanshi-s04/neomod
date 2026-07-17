#!/bin/bash
#SBATCH --job-name=s3m_1tr
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

# Case 1: SSP_number_tracklets=1, SSP_separation_threshold=0.5 arcsec.
# Isolates the effect of the 3-night linking requirement — only 1 nightly tracklet
# needed for detection; motion threshold unchanged.
# Compare to baseline (sorcha_comparison_s3m.parquet) to see how many more TNOs
# appear when the multi-night linking requirement is removed.

set -euo pipefail

WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
PY=$WORKDIR/conda_prep/bin/python
SORCHA_RUN=$WORKDIR/conda_prep/bin/sorcha-run

mkdir -p outputs/production_2yr_s3m_1tracklet work/production_2yr_s3m_1tracklet logs

"$PY" "$WORKDIR/production_run/multi_sorcha_production.py" \
  --input_orbits inputs/s3m_sorcha_orbits.csv \
  --input_physical inputs/s3m_sorcha_phys.csv \
  --pointings baseline_v5.0.0_2yrs.db \
  --config neomod/pipeline/config/Rubin_full_footprint_1tracklet.ini \
  --ar_data_path sorcha_cache_2025-07-06 \
  --outdir outputs/production_2yr_s3m_1tracklet \
  --workdir work/production_2yr_s3m_1tracklet \
  --chunksize 16000 \
  --norbits 1000 \
  --cores 16 \
  --instance ${SLURM_ARRAY_TASK_ID} \
  --sorcha_run $SORCHA_RUN
