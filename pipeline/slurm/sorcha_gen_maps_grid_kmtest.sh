#!/bin/bash
#SBATCH --job-name=sorcha_g_kmtest
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

# A/B TEST: regenerate 3 maps with the K|M/kNN NEO cloner (VDP_NEO_CLONER=km)
# instead of the GMM, to test whether it recovers the NEO velocity wings that
# the GMM under-covers (advisor's pure-NEO-cell test).
#   333 = dlon+000 lat+00 (antisun)  338 = dlon+050 lat+00   342 = dlon+090 lat+00
# Writes to prob_maps_grid_kmtest/ so the GMM production maps stay intact for
# side-by-side comparison.

set -euo pipefail
export VDP_NEO_CLONER=km            # <-- force the K|M/kNN cloner for NEO
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
mkdir -p "$WORKDIR/logs" "$WORKDIR/prob_maps_grid_kmtest"
echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  host=$(hostname)  NEO_CLONER=${VDP_NEO_CLONER}"
echo "start=$(date -Is)"
"$PY" "$WORKDIR/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --task-id       "${SLURM_ARRAY_TASK_ID}" \
    --prob-maps-dir prob_maps_grid_kmtest \
    --n-jobs        16 \
    --save-overlays \
    --overwrite
echo "end=$(date -Is)"
