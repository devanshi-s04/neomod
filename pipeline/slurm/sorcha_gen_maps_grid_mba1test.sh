#!/bin/bash
#SBATCH --job-name=sorcha_g_mba1
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=333,338,342
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# A/B test: MBA clone_factor 5 -> 1, to test whether the residual over-weighting
# of MBA density in the high-|vbeta| NEO wings (pure-NEO cells score ~0.68/0.56
# instead of ~1 under training-aligned labels) is driven by cf=5 fattening the
# kNN density tails. NEO stays on GMM (default). Maps: 333 antisun / 338 dlon+50
# / 342 dlon+90 -> prob_maps_grid_mba1test/ (production maps untouched).

set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
mkdir -p "$WORKDIR/logs" "$WORKDIR/prob_maps_grid_mba1test"
echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  host=$(hostname)"
echo "start=$(date -Is)"
"$PY" "$WORKDIR/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --task-id          "${SLURM_ARRAY_TASK_ID}" \
    --prob-maps-dir    prob_maps_grid_mba1test \
    --mba-clone-factor 1 \
    --n-jobs           16 \
    --save-overlays \
    --overwrite
echo "end=$(date -Is)"
