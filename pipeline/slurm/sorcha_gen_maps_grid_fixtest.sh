#!/bin/bash
#SBATCH --job-name=sorcha_g_fixtest
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

# Validate the cross-population normalisation fix (Bayesian mixture: rho_pop =
# N_source * pdf_pop) on 3 representative maps BEFORE the full 667-map regen:
#   333 = dlon+000 lat+00 (antisun)  -> control, must stay correct (no regression)
#   338 = dlon+050 lat+00            -> previously broken, must recover
#   342 = dlon+090 lat+00            -> previously broken, must recover
# Writes to prob_maps_grid_fixtest/ so the production maps are untouched for
# side-by-side comparison. --save-overlays for density/P(NEO) inspection.

set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
mkdir -p "$WORKDIR/logs" "$WORKDIR/prob_maps_grid_fixtest"
echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  host=$(hostname)"
echo "start=$(date -Is)"
"$PY" "$WORKDIR/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --task-id       "${SLURM_ARRAY_TASK_ID}" \
    --prob-maps-dir prob_maps_grid_fixtest \
    --n-jobs        16 \
    --save-overlays \
    --overwrite
echo "end=$(date -Is)"
