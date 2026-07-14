# Shared environment for the S3M linking experiment.
# Sourced by every stage script AFTER its #SBATCH block.
# The ONLY thing that varies between the three cases is $CASE, which selects
# the config file. Every path and every other pipeline parameter is identical.
#
# Required: export CASE=case1|case2|case3   (via  sbatch --export=ALL,CASE=caseN)

: "${CASE:?Set CASE=case1|case2|case3 via  --export=ALL,CASE=caseN}"

WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
PY=$WORKDIR/conda_prep/bin/python
SORCHA_RUN=$WORKDIR/conda_prep/bin/sorcha-run

case "$CASE" in
  case1) CONFIG=$WORKDIR/neomod/pipeline/config/s3m_case1_tk1_sep05.ini ;;  # tk=1  sep=0.5
  case2) CONFIG=$WORKDIR/neomod/pipeline/config/s3m_case2_tk3_sep05.ini ;;  # tk=3  sep=0.5
  case3) CONFIG=$WORKDIR/neomod/pipeline/config/s3m_case3_tk1_sep0.ini  ;;  # tk=1  sep=0.001
  *) echo "ERROR: unknown CASE=$CASE (expected case1|case2|case3)"; exit 1 ;;
esac

# Per-case directory tree (identical structure for every case).
BASE=$WORKDIR/outputs/s3m_linking/$CASE
PROD_DIR=$BASE/production                       # Sorcha .h5 (deletable after stage 2)
TRACK_DIR=$BASE/tracklets                       # per-file tracklet parquets
PHASE2_DIR=$BASE/phase2                         # vdp_shards/, digest2_shards/
SUBSAMPLE=$PHASE2_DIR/subsample.parquet
COMPARISON=$BASE/sorcha_comparison_$CASE.parquet
SORCHA_WORK=$WORKDIR/work/s3m_linking/$CASE     # Sorcha scratch (deletable after stage 1)

# Fixed inputs — identical to the original production run.
INPUT_ORBITS=$WORKDIR/inputs/s3m_sorcha_orbits.csv
INPUT_PHYS=$WORKDIR/inputs/s3m_sorcha_phys.csv
POINTINGS=$WORKDIR/baseline_v5.0.0_2yrs.db
AR_CACHE=$WORKDIR/sorcha_cache_2025-07-06
PROB_MAPS_FOOTPRINT=$WORKDIR/prob_maps_grid       # geometry only, for map assignment
PROB_MAPS_SCORE=$WORKDIR/prob_maps_grid_s3m       # pure-S3M densities, for VDP scoring
DIGEST2_DIR=$WORKDIR/digest2

cd "$WORKDIR"
mkdir -p logs "$BASE" "$PROD_DIR" "$TRACK_DIR" "$PHASE2_DIR" "$SORCHA_WORK"

echo "================ S3M LINKING PIPELINE ================"
echo "CASE     = $CASE"
echo "CONFIG   = $CONFIG"
echo "BASE     = $BASE"
echo "host     = $(hostname)   start=$(date -Is)"
echo "====================================================="
