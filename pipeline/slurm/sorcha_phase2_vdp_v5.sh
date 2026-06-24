#!/bin/bash
#SBATCH --job-name=sorcha_vdp_v5
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-112%32
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Phase 2a — VDP scoring for the v5.0 GRID run.
# Scores every tracklet in outputs/tracklets_v5_grid/ against its assigned
# antisun-relative grid map in prob_maps_grid/ (GMM maps -> --no-nearest-dist-mask).
#
# Shards: 14,445 files / batch-size 128 = 113 shards (array 0-112).
# Requires the sorcha_phase2.py fix that accepts "grid" map filenames
# (otherwise every grid-assigned tracklet scores NaN).
#
# Output: outputs/phase2_v5/vdp_shards/vdp_00000.parquet ... vdp_00112.parquet
#         (adds P_NEO_vdp, vlam, vbeta, mag_bin_label)

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha

cd "$WORKDIR"
mkdir -p logs outputs/phase2_v5

export XDG_CACHE_HOME=/tmp/${USER}/sorcha_phase2_xdg_cache
mkdir -p "$XDG_CACHE_HOME"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" score-vdp \
  --tracklet-dir outputs/tracklets_v5_grid \
  --prob-maps-dir prob_maps_grid \
  --work-dir outputs/phase2_v5 \
  --batch-size 128 \
  --shard-index "${SLURM_ARRAY_TASK_ID}" \
  --overwrite \
  --no-nearest-dist-mask \
  --support-mask-min 1

echo "end=$(date -Is)"

# --support-mask-min 1 (adopted 2026-06-22): zeros non-smoothed populations' density in
# velocity cells with no in-cell clone support, cutting the kNN estimator's bleed of the
# dense MBA core into the sparse NEO wings (the confirmed wing-suppression bug). NEO is
# smoothed and exempt. Makes the full vdp_shards canonical-masked.
