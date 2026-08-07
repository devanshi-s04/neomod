#!/bin/bash
#SBATCH --job-name=prodtrk
#SBATCH --partition=ckpt-all
#SBATCH --account=astro
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=06:00:00
#SBATCH --array=0-39
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
W=/mmfs1/gscratch/dirac/ds2004/sorcha
$W/conda_prep/bin/python $W/neomod/pipeline/sorcha_test_build_tracklets.py \
  --prod-dir $W/outputs/sorcha_production/production \
  --objects  $W/outputs/sorcha_production/SORCHA_TEST_OBJECTS_production.parquet \
  --out-dir  $W/outputs/sorcha_production/tracklets \
  --tag production --shard $SLURM_ARRAY_TASK_ID --nshards 40
