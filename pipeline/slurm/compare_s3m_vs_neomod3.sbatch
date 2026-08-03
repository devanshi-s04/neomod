#!/bin/bash
#SBATCH --job-name=s3m_vs_nm3
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.err
set -euo pipefail
export VDP_LOADER=s3m   # pure S3M elements (t_0=2008), NOT the hybrid catalog default (t_0=2023, fixing_integrator.md §13.6)
export MPLCONFIGDIR="/tmp/${USER}/mpl_nm3"; mkdir -p "$MPLCONFIGDIR"
/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python neomod/pipeline/compare_s3m_vs_neomod3.py
