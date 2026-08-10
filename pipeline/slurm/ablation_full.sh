#!/bin/bash
#SBATCH --job-name=ablfull
#SBATCH --partition=cpu-g2
#SBATCH --account=astro
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-666
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
W=/mmfs1/gscratch/dirac/ds2004/sorcha
B=$W/outputs/geometric_density_estimator_ablation
CEN=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $B/frozen_center_list.txt)
echo "mode=$ABL_MODE center=$CEN"
$W/conda_prep/bin/python $W/neomod/pipeline/run_density_ablation.py \
  --mode $ABL_MODE --center $CEN --out-dir $B/full/maps_$ABL_MODE --n-jobs 8
