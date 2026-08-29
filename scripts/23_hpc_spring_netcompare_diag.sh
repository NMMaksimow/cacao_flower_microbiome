#!/bin/bash
#SBATCH --job-name=net_diag
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH --output=logs/network/23_diag_%j.out
#SBATCH --error=logs/network/23_diag_%j.err

# Diagnostic: is the netCompare failure specific to SPRING, and to which farms?
# cores = 1 because the timing probe showed parallel workers are slower here.

module load bioinformatics
module load r/4.6.0

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/network

Rscript scripts/23_hpc_spring_netcompare_diag.R
