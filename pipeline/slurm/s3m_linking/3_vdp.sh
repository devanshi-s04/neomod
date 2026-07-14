#!/bin/bash
#SBATCH --job-name=s3mlink_vdp
#SBATCH --partition=ckpt-all
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-112%113
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# STAGE 3 — VDP scoring. Scores every tracklet against its assigned
# antisun-relative grid map in prob_maps_grid_s3m (pure-S3M densities).
# Same resources / flags as the original sorcha_phase2_vdp_s3m.sh
# (2 cpu, 32G, batch-size 128, --no-nearest-dist-mask, --support-mask-min 1).
#
# Speed: %113 runs all shards at once (was %32) on ckpt-all -> ~10 min wall.
# The array upper bound (0-112) matches ~14,381 tracklet files / 128.
# case2 (tk=3) produces fewer tracklet files; extra shards simply no-op.
# For an exact bound:  NSHARDS = ceil(#tracklet_parquets / 128).
#
# Submit per case (after stage 2 completes):
#   sbatch --export=ALL,CASE=case1 3_vdp.sh

set -euo pipefail
source /mmfs1/gscratch/dirac/ds2004/sorcha/neomod/pipeline/slurm/s3m_linking/_case_env.sh

export XDG_CACHE_HOME=/tmp/${USER}/s3mlink_vdp_xdg
mkdir -p "$XDG_CACHE_HOME"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" score-vdp \
  --tracklet-dir "$TRACK_DIR" \
  --prob-maps-dir "$PROB_MAPS_SCORE" \
  --work-dir "$PHASE2_DIR" \
  --batch-size 128 \
  --shard-index "${SLURM_ARRAY_TASK_ID}" \
  --overwrite \
  --no-nearest-dist-mask \
  --support-mask-min 1

echo "end=$(date -Is)"
