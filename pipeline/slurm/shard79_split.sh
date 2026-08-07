#!/bin/bash
#SBATCH --job-name=s79sub
#SBATCH --partition=cpu-g2
#SBATCH --account=astro
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=0-15
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
W=/mmfs1/gscratch/dirac/ds2004/sorcha
N=14891; K=16
S=$(( SLURM_ARRAY_TASK_ID * N / K ))
E=$(( (SLURM_ARRAY_TASK_ID+1) * N / K ))
echo "sub-range $S..$E"
$W/conda_prep/bin/python $W/neomod/pipeline/sorcha_test_score_shard.py \
  --first-tracklet $W/outputs/sorcha_production/SORCHA_TEST_FIRST_TRACKLET_production.parquet \
  --shard 79 --nshards 80 --out-dir $W/outputs/sorcha_production/scored_sub79 \
  --row-start $S --row-stop $E --d2-batch 200
