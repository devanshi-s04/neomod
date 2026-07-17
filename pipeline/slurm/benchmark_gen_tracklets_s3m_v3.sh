#!/bin/bash
#SBATCH --job-name=bm_gen_s3m_v3
#SBATCH --partition=ckpt-all
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# v3 benchmark: proportional population caps (same as v2, Option A) PLUS epoch
# matched to MJD 61642 (2027-08-25) -- the busiest single Sorcha night across
# case1/2/3 (4,316 / 4,490 / 4,353 tracklets, highest of any night in all three).
# See WAGG_SORCHA_HYAK_CONTEXT.md "benchmark population-cap discovery" and
# sorcha_full_pipeline.md sec 9. Writes to a NEW directory (v3), leaving v1/v2
# untouched. Maps are antisun-relative / epoch-independent, so only REF_OBSTIME
# changed in gen_benchmark_tracklets_s3m.py -- grid geometry is unaffected.
#
#   task 0 = NEO      (cap 12,900,  was uncapped)
#   task 1 = MBA      (cap 650,000, was 200,000)
#   task 2 = TNO      (cap 1,600,   was uncapped)
#   task 3 = Trojans  (cap 6,000,   was 100,000)
#
# Each task writes:
#   outputs/benchmark_tracklets_s3m_v3/tracklets_{POP}.parquet
#
# After ALL 4 tasks complete, combine on login node:
#   conda_prep/bin/python neomod/pipeline/gen_benchmark_tracklets_s3m.py --combine-only
#
# Then submit:
#   sbatch neomod/pipeline/slurm/benchmark_score_vdp_s3m_v3.sh

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha

mkdir -p logs outputs/benchmark_tracklets_s3m_v3

export XDG_CACHE_HOME=/tmp/${USER}/bm_gen_s3m_v3_xdg_${SLURM_ARRAY_TASK_ID}
mkdir -p "$XDG_CACHE_HOME"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/gen_benchmark_tracklets_s3m.py" \
    --pop-id "${SLURM_ARRAY_TASK_ID}"

echo "end=$(date -Is)"
