#!/bin/bash
#SBATCH --job-name=split820v5
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Recovery for inst00820_part003: same 3 objects that hung v3.3 (A804, A854, A868).
# Strips those rows from the work CSVs and runs Sorcha on the remaining 997 objects.
# See production_run/production_2yr_excluded_objects.txt for background.

set -euo pipefail

WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
PY=$WORKDIR/conda_prep/bin/python
SORCHA_RUN=$WORKDIR/conda_prep/bin/sorcha-run

ORBS_IN=$WORKDIR/work/production_2yr_v5/instance_00820/orbits_003.csv
PHYS_IN=$WORKDIR/work/production_2yr_v5/instance_00820/physical_003.csv
ORBS_OUT=$WORKDIR/work/production_2yr_v5/instance_00820/orbits_003_skip3.csv
PHYS_OUT=$WORKDIR/work/production_2yr_v5/instance_00820/physical_003_skip3.csv

echo "start=$(date -Is)  host=$(hostname)"

# Drop rows at 0-indexed data positions 514, 541, 617 (A804, A854, A868).
# With header at line 1, these are file lines 516, 543, 619.
"$PY" - <<'EOF'
import pandas as pd
orbs = pd.read_csv("/mmfs1/gscratch/dirac/ds2004/sorcha/work/production_2yr_v5/instance_00820/orbits_003.csv")
phys = pd.read_csv("/mmfs1/gscratch/dirac/ds2004/sorcha/work/production_2yr_v5/instance_00820/physical_003.csv")
drop_idx = [514, 541, 617]
print(f"Dropping rows {drop_idx}: {list(orbs.iloc[drop_idx]['ObjID'])}")
orbs = orbs.drop(index=drop_idx).reset_index(drop=True)
phys = phys.drop(index=drop_idx).reset_index(drop=True)
print(f"Remaining: {len(orbs)} objects")
orbs.to_csv("/mmfs1/gscratch/dirac/ds2004/sorcha/work/production_2yr_v5/instance_00820/orbits_003_skip3.csv", index=False)
phys.to_csv("/mmfs1/gscratch/dirac/ds2004/sorcha/work/production_2yr_v5/instance_00820/physical_003_skip3.csv", index=False)
EOF

echo "Filtered CSVs written. Running sorcha-run..."

mkdir -p "$WORKDIR/outputs/production_2yr_v5"

"$SORCHA_RUN" \
  -c neomod/pipeline/config/Rubin_full_footprint_detections.ini \
  --pd baseline_v5.0.0_2yrs.db \
  -o outputs/production_2yr_v5 \
  -t inst00820_part003 \
  --ob "$ORBS_OUT" \
  -p "$PHYS_OUT" \
  --ar sorcha_cache_2025-07-06 \
  -f

echo "end=$(date -Is)"
echo "Output: outputs/production_2yr_v5/inst00820_part003.h5"
