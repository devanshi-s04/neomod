#!/bin/bash
#SBATCH --job-name=srcdig
#SBATCH --partition=cpu-g2
#SBATCH --account=astro
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --array=0-47
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python \
  /mmfs1/gscratch/dirac/ds2004/sorcha/neomod/pipeline/source_selection_digest.py \
  --shard $SLURM_ARRAY_TASK_ID --nshards 48
