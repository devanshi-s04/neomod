#!/bin/bash
#SBATCH --job-name=nm3allskyfin
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.err
set -euo pipefail
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m MPLCONFIGDIR="$W/.nbexec/mpl"
P="$W/conda_prep/bin/python"
S="$W/neomod/pipeline/build_neomod3_high_allsky_cache.py"
echo "=== partition ==="; "$P" "$S" partition
echo "=== validate ===";  "$P" "$S" validate
