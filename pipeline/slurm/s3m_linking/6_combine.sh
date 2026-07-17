#!/bin/bash
#SBATCH --job-name=s3mlink_combine
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
# STAGE 6 — concatenate the digest2 shards into the final per-case comparison
# parquet:  outputs/s3m_linking/<case>/sorcha_comparison_<case>.parquet
# This is the analogue of sorcha_comparison_s3m.parquet, one per case.
#
# Submit per case (after stage 5 completes):
#   sbatch --export=ALL,CASE=case1 6_combine.sh

set -euo pipefail
source /mmfs1/gscratch/dirac/ds2004/sorcha/neomod/pipeline/slurm/s3m_linking/_case_env.sh

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" combine \
  --work-dir "$PHASE2_DIR" \
  --outfile "$COMPARISON" \
  --overwrite

echo "Wrote $COMPARISON"
"$PY" - <<PYEOF
import pandas as pd
df = pd.read_parquet("$COMPARISON", columns=["population"])
print("rows:", len(df))
print(df["population"].value_counts().to_string())
PYEOF
echo "end=$(date -Is)"
