#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 32
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH --job-name=ITS1_alpha_rare
#SBATCH --output=logs/07b_ITS1_alpha_%j.out
#SBATCH --error=logs/07b_ITS1_alpha_%j.err

module load bioinformatics
module load r/4.5.2

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Cores allocated: $SLURM_CPUS_PER_TASK"

Rscript scripts/07b_hpc_ITS1_alpha_div_rarefaction.R

echo "Job finished: $(date)"
