#!/bin/bash
#SBATCH --job-name=s3m_1tr_nosep
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

# Case 2: SSP_number_tracklets=1, SSP_separation_threshold=0.001 arcsec (effectively 0).
# Removes BOTH the multi-night linking requirement AND the motion threshold.
# Any object with 2 detections in one night is kept, regardless of how little it moved.
# Compare to Case 1 to isolate how much the 0.5 arcsec motion threshold costs.
# Note: 0.001 arcsec used instead of 0.0 because Sorcha rejects ssp_separation_threshold <= 0.

set -euo pipefail

WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
PY=$WORKDIR/conda_prep/bin/python
SORCHA_RUN=$WORKDIR/conda_prep/bin/sorcha-run

mkdir -p outputs/production_2yr_s3m_1tracklet_nosep work/production_2yr_s3m_1tracklet_nosep logs

"$PY" "$WORKDIR/production_run/multi_sorcha_production.py" \
  --input_orbits inputs/s3m_sorcha_orbits.csv \
  --input_physical inputs/s3m_sorcha_phys.csv \
  --pointings baseline_v5.0.0_2yrs.db \
  --config neomod/pipeline/config/Rubin_full_footprint_1tracklet_nosep.ini \
  --ar_data_path sorcha_cache_2025-07-06 \
  --outdir outputs/production_2yr_s3m_1tracklet_nosep \
  --workdir work/production_2yr_s3m_1tracklet_nosep \
  --chunksize 16000 \
  --norbits 1000 \
  --cores 16 \
  --instance ${SLURM_ARRAY_TASK_ID} \
  --sorcha_run $SORCHA_RUN
