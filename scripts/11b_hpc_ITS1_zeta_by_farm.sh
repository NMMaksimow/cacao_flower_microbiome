#!/bin/bash
#SBATCH --job-name=zeta_ITS1_byfarm
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 14
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --output=logs/zeta/11b_ITS1_byfarm_%j.out
#SBATCH --error=logs/zeta/11b_ITS1_byfarm_%j.err
#SBATCH --mail-type=END,FAIL

# 14 cores = one per farm × group combination (7 farms × 2 treatments).
# ITS1 OTU97 has fewer features than 16S ASVs so expected runtime is shorter.

module load r/4.5.1

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/zeta

Rscript scripts/11b_hpc_ITS1_zeta_by_farm.R
