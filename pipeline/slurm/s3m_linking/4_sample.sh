#!/bin/bash
#SBATCH --job-name=s3mlink_sample
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err
#
# STAGE 4 — stratified subsample for digest2 (keeps ALL NEOs, samples 500k
# non-NEOs). Single job. Matches the original run (n-sample-nonneo=500000,
# seed=42). Writes phase2/subsample.parquet.
#
# Submit per case (after stage 3 completes):
#   sbatch --export=ALL,CASE=case1 4_sample.sh

set -euo pipefail
source /mmfs1/gscratch/dirac/ds2004/sorcha/neomod/pipeline/slurm/s3m_linking/_case_env.sh

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" sample \
  --work-dir "$PHASE2_DIR" \
  --subsample-file "$SUBSAMPLE" \
  --n-sample-nonneo 500000 \
  --sample-seed 42 \
  --overwrite

echo "Subsample rows: $("$PY" -c "import pyarrow.parquet as pq; print(pq.ParquetFile('$SUBSAMPLE').metadata.num_rows)")"
echo "  -> digest2 array size = ceil(rows / 5000)"
echo "end=$(date -Is)"
