#!/bin/bash
#SBATCH --job-name=net_timing
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --output=logs/network/22_timing_%j.out
#SBATCH --error=logs/network/22_timing_%j.err

# Short probe: measures netCompare permutation cost on one farm so the real
# SparCC/SPRING jobs can be given a realistic wall time and core count.
# 4 cores so the 1-core vs 4-core speedup can be measured.

module load bioinformatics
module load r/4.6.0

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/network

Rscript scripts/22_hpc_netcompare_timing.R
