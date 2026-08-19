#!/bin/bash
#SBATCH --job-name=zeta_16S_byfarm
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 14
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --output=logs/zeta/11a_16S_byfarm_%j.out
#SBATCH --error=logs/zeta/11a_16S_byfarm_%j.err
#SBATCH --mail-type=END,FAIL

# 14 cores = one per farm × group combination (7 farms × 2 treatments).
# Wall time 4h is generous; actual runtime expected ~30–60 min.

module load r/4.5.1

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/zeta

Rscript scripts/11a_hpc_16S_zeta_by_farm.R
