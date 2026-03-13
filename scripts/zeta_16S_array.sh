#!/bin/bash
#SBATCH --job-name=zeta_16S
#SBATCH -p shared
#SBATCH --array=1-20          # orders 1 to 20; change upper limit as needed
#SBATCH --mem=8G
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=logs/zeta/zeta_16S_%a.out
#SBATCH --error=logs/zeta/zeta_16S_%a.err

# Zeta-diversity decline – 16S bacteria
# Each array task computes ONE zeta order (bagged + unbagged) using MC sampling
# Results saved as per-order RDS files, then collected by zeta_collect.R
#
# Usage:
#   mkdir -p logs/zeta
#   sbatch scripts/zeta_16S_array.sh

echo "=== Zeta 16S | order ${SLURM_ARRAY_TASK_ID} ==="
echo "Date:  $(date)"
echo "Host:  $(hostname)"
echo "JobID: ${SLURM_JOB_ID}, TaskID: ${SLURM_ARRAY_TASK_ID}"

# Load R module
# NOTE: check available versions with: module avail r
module load r/4.4.3

# Create tmp output directory if it doesn't exist
mkdir -p results/rds/zeta_tmp

Rscript scripts/zeta_array.R \
    --amplicon 16S \
    --order    ${SLURM_ARRAY_TASK_ID} \
    --sam      10000

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: order ${SLURM_ARRAY_TASK_ID} completed at $(date)"
else
    echo "ERROR: order ${SLURM_ARRAY_TASK_ID} failed"
    exit 1
fi
