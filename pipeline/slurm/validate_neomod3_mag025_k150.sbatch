#!/bin/bash
#SBATCH --job-name=mag025val
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
"$W/conda_prep/bin/python" "$W/neomod/pipeline/validate_neomod3_mag025_k150_maps.py" \
  --n-jobs "${SLURM_CPUS_PER_TASK}" "$@"
