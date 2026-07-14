#!/bin/bash
#SBATCH --job-name=kurl_d2
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-23
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/logs/kurl_d2_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/logs/kurl_d2_%A_%a.err
#
# digest2 scoring for the Kurlander referee (referee_eval.parquet, 120,000 tracklets).
# The original serial call (score_digest2.py, --chunk 20000, run directly) took ~80 min
# for ONE 20,000-row chunk (digest2 itself, not I/O, scales poorly with a single huge
# call) -- ~8h projected serial, and was wrongly run on the shared login node. Fixed by
# mirroring the proven pattern from case1/2/3 + benchmark v1/v2/v3 digest2 stages this
# session: OUTER_CHUNK=5000/task, 4 parallel GNU-parallel workers of INNER_CHUNK=1250
# rows each, 24 tasks x 5000 = 120,000 rows exactly.
#
# Writes: outputs/kurlander/d2_shards/d2_r<start>.parquet  (one per 1250-row worker)
# Merge after all tasks finish (see bottom of this file for the one-liner).

module load gnu/parallel/20210422
set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
IN="$WORKDIR/outputs/kurlander/referee_eval.parquet"
OUTDIR="$WORKDIR/outputs/kurlander/d2_shards"
D2DIR="$WORKDIR/digest2"
OUTER_CHUNK=5000   # rows per SLURM task (4 workers x 1250)
INNER_CHUNK=1250   # rows per digest2 worker

mkdir -p "$OUTDIR"

BASE_ROW=$(( SLURM_ARRAY_TASK_ID * OUTER_CHUNK ))
echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  base_row=${BASE_ROW}  host=$(hostname)"
echo "start=$(date -Is)"

run_worker() {
    local sub_id=$1
    local row_start=$(( BASE_ROW + sub_id * INNER_CHUNK ))
    local row_stop=$(( row_start + INNER_CHUNK ))
    "$PY" "$WORKDIR/neomod/pipeline/kurlander/score_digest2.py" \
        --in "$IN" \
        --out "$OUTDIR/d2_r$(printf %07d $row_start).parquet" \
        --d2-dir "$D2DIR" \
        --chunk "$INNER_CHUNK" \
        --row-start "$row_start" \
        --row-stop  "$row_stop"
}

export -f run_worker
export PY WORKDIR IN OUTDIR D2DIR BASE_ROW INNER_CHUNK

parallel --jobs 4 run_worker ::: 0 1 2 3

echo "end=$(date -Is)"

# ---- after all 24 tasks finish, merge on the login node ----
#   python -c "
#   import pandas as pd, glob
#   parts = [pd.read_parquet(f) for f in sorted(glob.glob('outputs/kurlander/d2_shards/d2_r*.parquet'))]
#   out = pd.concat(parts, ignore_index=True)
#   out.to_parquet('outputs/kurlander/referee_eval_d2.parquet', index=False)
#   print(len(out), 'rows; nan frac', out.P_NEO_d2.isna().mean())
#   "
