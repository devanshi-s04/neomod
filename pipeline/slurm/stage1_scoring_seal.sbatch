#!/bin/bash
#SBATCH --job-name=t2seal
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.err
set -euo pipefail
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m MPLCONFIGDIR="$W/.nbexec/mpl"
"$W/conda_prep/bin/python" "$W/neomod/pipeline/stage1_scoring_seal.py"
