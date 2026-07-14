#!/bin/bash
#SBATCH --job-name=bm_d2_s3m
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --array=0-118%100
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# digest2 scoring for benchmark pure-S3M tracklets (667-map GMM grid).
#
# IMPORTANT: --subsample-file must point to a nonexistent path so run-digest2
# reads from vdp_shards/ instead of the default outputs/phase2/sorcha_subsample.parquet.
#
# Each SLURM task covers 4000 rows split across 4 parallel digest2 workers (1000 each).
# 119 tasks x ~3 min each, %100 concurrency -> ~4 min total (all run at once).
#
# Reads:  outputs/phase2_benchmark_s3m/vdp_shards/vdp_00000.parquet (475,027 rows)
# Output: outputs/phase2_benchmark_s3m/digest2_shards/d2_vdp_00000_r*.parquet
#
# Combine after all tasks finish with:
#
#   python neomod/pipeline/sorcha_phase2.py combine \
#       --work-dir outputs/phase2_benchmark_s3m \
#       --outfile outputs/phase2_benchmark_s3m/benchmark_comparison_s3m.parquet \
#       --overwrite

module load gnu/parallel/20210422

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
OUTER_CHUNK=4000   # rows per SLURM task (4 workers x 1000)
INNER_CHUNK=1000   # rows per digest2 worker

mkdir -p logs outputs/phase2_benchmark_s3m/digest2_shards

BASE_ROW=$(( SLURM_ARRAY_TASK_ID * OUTER_CHUNK ))

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  base_row=${BASE_ROW}  host=$(hostname)"
echo "start=$(date -Is)"

# Run 4 digest2 workers in parallel, one per 1000-row sub-chunk.
run_worker() {
    local sub_id=$1
    local row_start=$(( BASE_ROW + sub_id * INNER_CHUNK ))
    local row_stop=$(( row_start + INNER_CHUNK ))
    "$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" run-digest2 \
        --work-dir outputs/phase2_benchmark_s3m \
        --subsample-file outputs/phase2_benchmark_s3m/no_subsample.parquet \
        --digest2-dir "$WORKDIR/digest2" \
        --digest2-chunk-tracklets 1000 \
        --digest2-timeout-sec 900 \
        --d2-row-start "$row_start" \
        --d2-row-stop  "$row_stop" \
        --overwrite
}

export -f run_worker
export PY WORKDIR BASE_ROW INNER_CHUNK

parallel --jobs 4 run_worker ::: 0 1 2 3

echo "end=$(date -Is)"
