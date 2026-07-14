#!/bin/bash
#SBATCH --job-name=bm_d2_troj
#SBATCH --partition=ckpt
#SBATCH --account=astro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --array=0-49%50
#SBATCH --chdir=/mmfs1/gscratch/dirac/ds2004/sorcha
#SBATCH --export=all
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Retry for the 50 Trojan digest2 shards that timed out in job 36607184.
# 1 worker per task, 45-min limit, all 50 run simultaneously.

set -euo pipefail

PY=/mmfs1/gscratch/dirac/ds2004/sorcha/conda_prep/bin/python
WORKDIR=/mmfs1/gscratch/dirac/ds2004/sorcha

ROW_STARTS=(376000 377000 378000 379000 380000 381000 382000 383000 384000 385000 386000 387000 400000 401000 402000 403000 404000 405000 406000 407000 408000 409000 410000 411000 412000 413000 414000 415000 416000 417000 418000 419000 420000 423000 444000 445000 446000 447000 448000 449000 450000 451000 452000 453000 454000 455000 456000 457000 458000 459000)

mkdir -p logs outputs/phase2_benchmark_s3m/digest2_shards

ROW_START=${ROW_STARTS[$SLURM_ARRAY_TASK_ID]}
ROW_STOP=$(( ROW_START + 1000 ))

echo "job_id=${SLURM_JOB_ID}  array_task=${SLURM_ARRAY_TASK_ID}  rows=${ROW_START}-${ROW_STOP}  host=$(hostname)"
echo "start=$(date -Is)"

"$PY" "$WORKDIR/neomod/pipeline/sorcha_phase2.py" run-digest2 \
  --work-dir outputs/phase2_benchmark_s3m \
  --subsample-file outputs/phase2_benchmark_s3m/no_subsample.parquet \
  --digest2-dir "$WORKDIR/digest2" \
  --digest2-chunk-tracklets 1000 \
  --digest2-timeout-sec 2400 \
  --d2-row-start "$ROW_START" \
  --d2-row-stop  "$ROW_STOP" \
  --overwrite

echo "end=$(date -Is)"
