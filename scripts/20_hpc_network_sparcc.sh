#!/bin/bash
#SBATCH --job-name=net_sparcc
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -t 16:00:00
#SBATCH --output=logs/network/20_sparcc_%j.out
#SBATCH --error=logs/network/20_sparcc_%j.err
#SBATCH --mail-type=END,FAIL

# SparCC networks, 16S (20a) + ITS1 (20b), nPerm = 200.
# Measured 6.47 s per permutation on one core -> ~5 h for 14 farm-pairs.
# 16 h is ~3x headroom. 1 core: parallel workers measured ~2x slower.

module load bioinformatics
module load r/4.6.0

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/network results/rds results/figures results/tables

Rscript scripts/20a_16S_network_by_farm_SparCC.R
Rscript scripts/20b_ITS1_network_by_farm_SparCC.R

Rscript -e "sink('logs/network/20_session_info.txt'); sessionInfo(); sink()"
