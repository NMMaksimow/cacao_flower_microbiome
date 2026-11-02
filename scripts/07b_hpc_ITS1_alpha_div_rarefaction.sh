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
module load r/4.6.1    # same R minor version as the laptop; packages in ~/R/library

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs

# Data version (v1 | v2) and sample set (main = 288 flowers, with_deviations = 294),
# read by scripts/utils_paths.R. Set them at submission, e.g.
#   CFM_SAMPLE_SET=with_deviations sbatch scripts/07b_hpc_ITS1_alpha_div_rarefaction.sh
export CFM_VERSION=${CFM_VERSION:-v2}
export CFM_SAMPLE_SET=${CFM_SAMPLE_SET:-main}

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Cores allocated: $SLURM_CPUS_PER_TASK"
echo "CFM_VERSION: $CFM_VERSION, sample set: $CFM_SAMPLE_SET"

Rscript scripts/07b_hpc_ITS1_alpha_div_rarefaction.R

echo "Job finished: $(date)"
