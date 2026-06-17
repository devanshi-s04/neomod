#!/bin/bash
#SBATCH --job-name=split820v5f
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

# Recovery for inst00820_part003: 5 pathological objects removed total.
# From original part003: A804 RA, A854 OA, A868 TA (v3.3 known)
# From skip3:            A854 RA, A868 WA (found via bisection on v5.0)
# Remaining: 995 objects.

set -euo pipefail

WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha
SORCHA_RUN=$WORKDIR/conda_prep/bin/sorcha-run

echo "start=$(date -Is)  host=$(hostname)"

"$SORCHA_RUN" \
  -c neomod/pipeline/config/Rubin_full_footprint_detections.ini \
  --pd baseline_v5.0.0_2yrs.db \
  -o outputs/production_2yr_v5 \
  -t inst00820_part003 \
  --ob work/production_2yr_v5/instance_00820/orbits_003_skip5.csv \
  -p  work/production_2yr_v5/instance_00820/physical_003_skip5.csv \
  --ar sorcha_cache_2025-07-06 \
  -f

echo "end=$(date -Is)"
echo "Output: outputs/production_2yr_v5/inst00820_part003.h5"
