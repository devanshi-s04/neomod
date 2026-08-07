#!/bin/bash
#SBATCH --job-name=prodscore
#SBATCH --partition=ckpt-all
#SBATCH --account=astro
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=0-79
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
W=/mmfs1/gscratch/dirac/ds2004/sorcha
$W/conda_prep/bin/python $W/neomod/pipeline/sorcha_test_score_shard.py \
  --first-tracklet $W/outputs/sorcha_production/SORCHA_TEST_FIRST_TRACKLET_production.parquet \
  --shard $SLURM_ARRAY_TASK_ID --nshards 80 \
  --out-dir $W/outputs/sorcha_production/scored
