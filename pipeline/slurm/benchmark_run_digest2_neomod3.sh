#!/bin/bash
#SBATCH --job-name=bm_d2_nm3
#SBATCH --partition=ckpt-all
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-805%810
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/d2nm3/%x_%A_%a.out
#SBATCH --error=logs/d2nm3/%x_%A_%a.err
#
# digest2 on the FULL NEOMOD3 benchmark -- all 3,220,714 rows, NO subsampling, so the NEO prior
# (0.776%) is correct by construction and no reweighting is needed (EVALUATION_PROTOCOL.md v1.0 §4.2).
# 806 tasks x 4000 rows (4 parallel workers x 1000).
#
# digest2 runs in `repeatable` mode (sorcha_phase2.py cfg_text) -- deterministic, per protocol §1.
# ckpt-all is preemptible; resubmit the same array to mop up killed shards (shard files are skipped
# when present unless --overwrite).
module load gnu/parallel/20210422
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
OUTER_CHUNK=4000
INNER_CHUNK=1000
mkdir -p logs/d2nm3 outputs/phase2_benchmark_neomod3/digest2_shards
BASE_ROW=$(( SLURM_ARRAY_TASK_ID * OUTER_CHUNK ))
echo "task=${SLURM_ARRAY_TASK_ID} base_row=${BASE_ROW} host=$(hostname) start=$(date -Is)"
run_worker() {
    local sub_id=$1
    local row_start=$(( BASE_ROW + sub_id * INNER_CHUNK ))
    local row_stop=$(( row_start + INNER_CHUNK ))
    "$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" run-digest2 \
        --work-dir outputs/phase2_benchmark_neomod3 \
        --subsample-file outputs/phase2_benchmark_neomod3/no_subsample.parquet \
        --digest2-dir "$WORKDIR/digest2" \
        --digest2-chunk-tracklets 1000 \
        --digest2-timeout-sec 900 \
        --d2-row-start "$row_start" --d2-row-stop "$row_stop"
}
export -f run_worker
export PY WORKDIR BASE_ROW INNER_CHUNK
parallel --jobs 4 run_worker ::: 0 1 2 3
echo "end=$(date -Is)"
