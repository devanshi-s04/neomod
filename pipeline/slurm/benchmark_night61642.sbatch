#!/bin/bash
#SBATCH --job-name=bm_n61642
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/logs/bm_n61642_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/logs/bm_n61642_%A_%a.err
#
# Targeted benchmark on Sorcha's night-61642 object set (Arnor ObjID-mapping request).
# One task per population (0=NEO 1=MBA 2=TNO 3=Trojans). MBA (task 1) loads the full
# 13.88M S3M population then filters to the ~10.6k targets -> the heavy one; 48G mem.
# After all 4 finish:
#   conda_prep/bin/python neomod/pipeline/gen_benchmark_night61642.py --combine-only

set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python

export XDG_CACHE_HOME=/tmp/${USER}/bm_n61642_xdg_${SLURM_ARRAY_TASK_ID}
mkdir -p "$XDG_CACHE_HOME"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "job_id=${SLURM_JOB_ID}  pop=${SLURM_ARRAY_TASK_ID}  host=$(hostname)  start=$(date -Is)"

"$PY" neomod/pipeline/gen_benchmark_night61642.py --pop-id "${SLURM_ARRAY_TASK_ID}"

echo "end=$(date -Is)"
