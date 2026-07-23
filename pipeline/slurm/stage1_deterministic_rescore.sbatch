#!/bin/bash
#SBATCH --job-name=d2s1
#SBATCH --partition=ckpt-all
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=00:40:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# STAGE 1 deterministic digest2 rescore, one row-chunk per array task.
# digest2 `repeatable` is cpu-count-independent, so 16 cpus just speeds it up (scores identical).
#
# Submit (chunk-size 20000):
#   sbatch --export=ALL,SOURCE=bench         --array=0-33 neomod/pipeline/slurm/stage1_deterministic_rescore.sbatch
#   sbatch --export=ALL,SOURCE=sorcha_rawfix --array=0-32 neomod/pipeline/slurm/stage1_deterministic_rescore.sbatch
#   sbatch --export=ALL,SOURCE=sorcha_vcorr  --array=0-32 neomod/pipeline/slurm/stage1_deterministic_rescore.sbatch
# Then per source: combine + analyze on the login node.

set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
SOURCE="${SOURCE:?set SOURCE=bench|sorcha_rawfix|sorcha_vcorr}"
CHUNK_SIZE="${CHUNK_SIZE:-20000}"
echo "start=$(date -Is) source=$SOURCE chunk=$SLURM_ARRAY_TASK_ID size=$CHUNK_SIZE cpus=$SLURM_CPUS_PER_TASK"
DIGEST2_AUDIT_CPUS="$SLURM_CPUS_PER_TASK" "$PY" \
  /mmfs1/gscratch/dirac/ds2004/sorcha/neomod/pipeline/stage1_deterministic_rescore.py run \
  --source "$SOURCE" --chunk-index "$SLURM_ARRAY_TASK_ID" --chunk-size "$CHUNK_SIZE"
echo "end=$(date -Is)"
