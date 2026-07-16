#!/bin/bash
#SBATCH --job-name=n61642_score
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/logs/n61642_score_%j.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/logs/n61642_score_%j.err
#
# Add P_NEO_vdp + P_NEO_d2 to the identity-matched benchmark_night61642 tracklets
# (Arnor request, docs/NOTE_FOR_HYAK_benchmark_night61642_scores.md). One node, chains:
#   1) score-vdp   with prob_maps_grid_s3m + --no-nearest-dist-mask --support-mask-min 1
#      (IDENTICAL to Sorcha stage 3 / 3_vdp.sh, so scores are comparable to Sorcha's P_NEO_vdp)
#   2) run-digest2 (4 GNU-parallel workers over the single 8,436-row vdp shard)
#   3) combine     -> benchmark_night61642_scored.parquet  (P_NEO_vdp + P_NEO_d2 + s3m_objid)
# 8,436 rows: VDP a few min (loads only referenced cells), digest2 ~8 min at 4x.

module load gnu/parallel/20210422
set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
BASE=$WORKDIR/outputs/benchmark_night61642_scored
TRACK_DIR=$BASE/tracklets_in          # contains ONLY tracklets_night61642.parquet
FINAL=$WORKDIR/outputs/benchmark_night61642/benchmark_night61642_scored.parquet

export XDG_CACHE_HOME=/tmp/${USER}/n61642_xdg; mkdir -p "$XDG_CACHE_HOME"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK} NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "host=$(hostname)  start=$(date -Is)"

# ---- 1. VDP scoring (S3M densities, Sorcha-stage-3 flags) ----
"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" score-vdp \
  --tracklet-dir "$TRACK_DIR" \
  --prob-maps-dir prob_maps_grid_s3m \
  --work-dir "$BASE" \
  --batch-size 1 --shard-index 0 \
  --overwrite --no-nearest-dist-mask --support-mask-min 1
echo "vdp done $(date -Is)"

# ---- 2. digest2 scoring (4 parallel workers over the one vdp shard) ----
N=$("$PY" -c "import pyarrow.parquet as pq,glob; print(pq.ParquetFile(sorted(glob.glob('$BASE/vdp_shards/vdp_*.parquet'))[0]).metadata.num_rows)")
CH=$(( (N + 3) / 4 ))                  # 4 near-equal row chunks
echo "digest2: $N rows, chunk=$CH"

d2_worker() {
    local rs=$1 re=$2
    "$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" run-digest2 \
        --work-dir "$BASE" \
        --subsample-file "$BASE/no_subsample.parquet" \
        --digest2-dir "$WORKDIR/digest2" \
        --digest2-chunk-tracklets 1000 --digest2-timeout-sec 1800 \
        --d2-row-start "$rs" --d2-row-stop "$re" --overwrite
}
export -f d2_worker; export PY WORKDIR BASE
parallel --jobs 4 --colsep ' ' d2_worker {1} {2} ::: \
  "0 $CH" "$CH $((2*CH))" "$((2*CH)) $((3*CH))" "$((3*CH)) $N"
echo "digest2 done $(date -Is)"

# ---- 3. combine ----
"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" combine \
  --work-dir "$BASE" --outfile "$FINAL" --overwrite

"$PY" - <<PYEOF
import pandas as pd
df = pd.read_parquet("$FINAL")
print("rows:", len(df))
print("has P_NEO_vdp:", "P_NEO_vdp" in df, " has P_NEO_d2:", "P_NEO_d2" in df, " has s3m_objid:", "s3m_objid" in df)
print("nan P_NEO_vdp:", int(df.P_NEO_vdp.isna().sum()), " nan P_NEO_d2:", int(df.P_NEO_d2.isna().sum()))
print(df["population"].value_counts().to_dict())
PYEOF
echo "end=$(date -Is)"
