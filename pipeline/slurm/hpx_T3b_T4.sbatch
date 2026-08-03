#!/bin/bash
#SBATCH --job-name=hpxT3bT4
#SBATCH --account=astro
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=400G
#SBATCH --time=03:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/sorcha/logs/%x_%j.err
set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
W=/mmfs1/gscratch/dirac/ds2004/sorcha
export PYTHONPATH="$W/neomod/src"
export MPLCONFIGDIR="/tmp/${USER}/mpl_t3b"; mkdir -p "$MPLCONFIGDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VDP_LOADER=s3m
# THE point of this job: read the partitioned cache, not the monolithic one.
export NEOMOD3_CACHE_DIR="$W/outputs/neomod3_projection_cache/by_pixel"

# both centers concurrently, n_jobs=16 each to match how the reference maps were built
for spec in "0 0" "20 -12"; do
  set -- $spec
  "$PY" "$W/neomod/pipeline/sorcha_gen_maps_grid.py" \
      --delta-lon "$1" --lat "$2" --ref-obstime "2027-08-25T00:00:00" \
      --cache "$W/outputs/epoch_state_cache/epoch_state_2027-08-25T000000.parquet" \
      --prob-maps-dir prob_maps_grid_neomod3_vlim5_hpx \
      --vdp-module velocity_density_pipeline_neomod_clone_only \
      --velocity-grid-limit 5.0 --velocity-grid-step 0.01 \
      --mba-clone-factor 1 --n-jobs 16 --overwrite &
done
wait
echo "=== builds done $(date -Is) ==="
ls -la "$W/prob_maps_grid_neomod3_vlim5_hpx/"
echo; echo "############ T3b: map comparison ############"
"$PY" "$W/neomod/pipeline/hpx_T3b_map_compare.py"
echo; echo "############ T4: scoring regression ############"
NM3_MAPS_DIR="$W/prob_maps_grid_neomod3_vlim5_hpx" \
NM3_OUT_DIR="$W/outputs/neomod3_score_roc_hpx" \
  "$PY" "$W/neomod/pipeline/nm3_score_roc.py"
echo; echo "############ reference (monolithic-built) numbers for comparison ############"
cat "$W/outputs/neomod3_score_roc/neomod3_roc_summary.csv"
