#!/bin/bash
#SBATCH --job-name=d2audit
#SBATCH --partition=ckpt-all
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# Auditable digest2 re-scoring, one chunk of objects per array task (all 3 variants B/S/V).
# The runner reads outputs/digest2_audit_night61642/digest2_audit_pairs.parquet (build it first
# on the login node: `audit_digest2.py build-pairs`).
#
# SELECT=neo   -> 794 matched NEOs      -> chunk-size 200 => --array=0-3
# SELECT=all   -> ~4316 matched objects -> chunk-size 200 => --array=0-21
# (empty tail chunks no-op cleanly, so an oversized --array is harmless.)
#
# Submit:
#   sbatch --export=ALL,SELECT=neo,CHUNK_SIZE=200 --array=0-3   neomod/pipeline/slurm/audit_digest2.sbatch
#   sbatch --export=ALL,SELECT=all,CHUNK_SIZE=200 --array=0-21  neomod/pipeline/slurm/audit_digest2.sbatch
# Then combine on the login node:
#   $PY neomod/pipeline/audit_digest2.py combine --select neo
#   $PY neomod/pipeline/audit_digest2.py combine --select all

set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
SELECT="${SELECT:-neo}"
CHUNK_SIZE="${CHUNK_SIZE:-200}"

echo "start=$(date -Is)  select=${SELECT}  chunk=${SLURM_ARRAY_TASK_ID}  chunk_size=${CHUNK_SIZE}"
"$PY" /mmfs1/gscratch/dirac/ds2004/sorcha/neomod/pipeline/audit_digest2.py run \
  --select "$SELECT" \
  --chunk-index "$SLURM_ARRAY_TASK_ID" \
  --chunk-size "$CHUNK_SIZE"
echo "end=$(date -Is)"
