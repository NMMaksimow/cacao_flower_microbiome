#!/bin/bash
#SBATCH --job-name=net_spring
#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -t 16:00:00
#SBATCH --output=logs/network/21_spring_%j.out
#SBATCH --error=logs/network/21_spring_%j.err
#SBATCH --mail-type=END,FAIL

# SPRING networks, 16S (21a) + ITS1 (21b), nPerm = 200.
# NOTE: netCompare currently FAILS for SPRING (nearPD does not converge on
# permuted data). Do not submit until script 23's diagnostic is resolved.
# 1 core: parallel workers measured ~2x slower.

module load bioinformatics
module load r/4.6.0

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/network results/rds results/figures results/tables

Rscript scripts/21a_16S_network_by_farm_SPRING.R
Rscript scripts/21b_ITS1_network_by_farm_SPRING.R

Rscript -e "sink('logs/network/21_session_info.txt'); sessionInfo(); sink()"
