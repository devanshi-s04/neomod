#!/bin/bash
#SBATCH --job-name=vdp_Vband
#SBATCH --account=astro
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:30:00
#SBATCH --array=0-2
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --output=/mmfs1/gscratch/dirac/ds2004/logs/vdp_Vband_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/dirac/ds2004/logs/vdp_Vband_%A_%a.err
#
# Option A (SORCHA_MAG_BAND_FIX_PLAN.md): re-score Sorcha VDP in Johnson-V magnitude.
# One task per case. Adds mean_mag_V + P_NEO_vdp_Vband to a NEW *_Vband.parquet;
# the original comparison parquets are untouched. ~40-65 min/task, all 3 in parallel.

set -euo pipefail
PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
CASES=(case1 case2 case3)
C=${CASES[$SLURM_ARRAY_TASK_ID]}

IN=outputs/s3m_linking/$C/sorcha_comparison_$C.parquet
OUT=outputs/s3m_linking/$C/sorcha_comparison_${C}_Vband.parquet

echo "case=$C  host=$(hostname)  start=$(date -Is)"
"$PY" neomod/pipeline/rescore_vdp_Vband.py --comparison "$IN" --out "$OUT"
echo "end=$(date -Is)"
