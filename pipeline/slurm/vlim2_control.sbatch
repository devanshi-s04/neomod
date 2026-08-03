#!/bin/bash
#SBATCH --job-name=vlim2_ctrl
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# STAGE C control: rebuild ONE center at the PRODUCTION velocity domain (+/-2.0, 0.01) via the NEW
# CLI options, into a separate dir, to prove the --velocity-grid-limit/-step change reproduces the
# existing production map (i.e. passing the defaults explicitly == not passing them).
# Center dlon=0 lat=-1 (first pilot center; exists in production).

set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/codex_integrator:$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_vlim2ctrl"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

echo "control dlon=0 lat=-1 at +/-2.0 start=$(date -Is) host=$(hostname)"
"$PY" "$W/neomod/pipeline/sorcha_gen_maps_grid.py" \
    --delta-lon 0 --lat -1 \
    --ref-obstime "2027-08-25T00:00:00" \
    --cache "$W/outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet" \
    --prob-maps-dir prob_maps_grid_s3m_nbody_control2 \
    --mba-clone-factor 5 \
    --velocity-grid-limit 2.0 \
    --velocity-grid-step 0.01 \
    --n-jobs 8 \
    --overwrite
echo "end=$(date -Is)"
