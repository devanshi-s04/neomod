#!/bin/bash
#SBATCH --job-name=oneB_falsify
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha/neomod
#SBATCH --export=all
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.err

set -euo pipefail

WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
PY="$WORKDIR/conda_prep/bin/python"

mkdir -p "$WORKDIR/logs" "$WORKDIR/outputs"

export XDG_CACHE_HOME=/tmp/${USER}/oneB_falsification_xdg_cache
mkdir -p "$XDG_CACHE_HOME"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

export ONEB_PARQUET="$WORKDIR/outputs/phase2_v5/sorcha_comparison_v5_masked.parquet"
export ONEB_MAPDIR="$WORKDIR/prob_maps_grid"
export ONEB_OUTFILE="$WORKDIR/outputs/oneB_falsification_scores.parquet"

echo "job_id=${SLURM_JOB_ID} host=$(hostname)"
echo "start=$(date -Is)"
echo "ONEB_PARQUET=$ONEB_PARQUET"
echo "ONEB_MAPDIR=$ONEB_MAPDIR"
echo "ONEB_OUTFILE=$ONEB_OUTFILE"

"$PY" -u oneB_falsification.py

echo "end=$(date -Is)"
