#!/bin/bash
#SBATCH --job-name=mag025gate
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=64
#SBATCH --mem=300G
#SBATCH --time=08:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.err
set -euo pipefail
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m MPLCONFIGDIR="$W/.nbexec/mpl"
export NEOMOD3_CACHE_DIR="$W/outputs/neomod3_projection_cache/by_pixel"
echo "host=$(hostname) cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
"$W/conda_prep/bin/python" "$W/neomod/pipeline/build_neomod3_mag025_k150_maps.py" \
  --center dlon+000_lat+00 --n-jobs "${SLURM_CPUS_PER_TASK}" \
  --density-engine closed_form --neo-source high_cell \
  --only-bins V024.00_024.25,V024.25_024.50,V024.50_024.75,V024.75_025.00 \
  --out-dir "$W/outputs/neomod3_mag025_k150_maps_v2_gate" --overwrite
echo "end=$(date -Is)"
