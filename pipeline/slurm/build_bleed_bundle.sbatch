#!/bin/bash
#SBATCH --job-name=bleedbundle
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=00:40:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export BUNDLE_JOBS="${SLURM_CPUS_PER_TASK}"
export MPLCONFIGDIR="/tmp/${USER}/mpl_bundle"; mkdir -p "$MPLCONFIGDIR"
/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python neomod/pipeline/build_bleed_bundle.py
