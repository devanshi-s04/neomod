#!/bin/bash
#SBATCH --job-name=sorcha_v5
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --array=0-903%16
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# DO NOT add --requeue: on ckpt, preempted jobs with --requeue fail to re-fetch
# --export=ALL env → stuck in "user env retrieval failed requeued held".
# Use --skip-existing logic (h5_is_good check in py) + manual resubmit instead.

set -euo pipefail

WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
PY=$WORKDIR/conda_prep/bin/python
SORCHA_RUN=$WORKDIR/conda_prep/bin/sorcha-run

cd "$WORKDIR"

mkdir -p outputs/production_2yr_v5 work/production_2yr_v5 logs

"$PY" "$WORKDIR/production_run/multi_sorcha_production.py" \
  --input_orbits inputs/hybrid_sorcha_orbits.csv \
  --input_physical inputs/hybrid_sorcha_phys.csv \
  --pointings baseline_v5.0.0_2yrs.db \
  --config neomod/pipeline/config/Rubin_full_footprint_detections.ini \
  --ar_data_path sorcha_cache_2025-07-06 \
  --outdir outputs/production_2yr_v5 \
  --workdir work/production_2yr_v5 \
  --chunksize 16000 \
  --norbits 1000 \
  --cores 16 \
  --instance ${SLURM_ARRAY_TASK_ID} \
  --sorcha_run $SORCHA_RUN
