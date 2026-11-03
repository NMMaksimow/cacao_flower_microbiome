#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 32
#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH --job-name=16S_alpha_rare
#SBATCH --output=logs/07a_16S_alpha_%j.out
#SBATCH --error=logs/07a_16S_alpha_%j.err

module load bioinformatics
module load r/4.6.1    # same R minor version as the laptop; packages in ~/R/library

# Without this the whole job runs on ONE core, silently. Hamilton sets
# OMP_PROC_BIND=spread in every job; the bioinformatics R then binds its main
# process to one core and all PSOCK workers inherit that binding (seen 2026-11-02:
# 33 R processes on one CPU). Same fix as in 03_hpc_rfc.sh.
export OMP_PROC_BIND=false

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs

# Data version (v1 | v2) and sample set (main = 288 flowers, with_deviations = 294),
# read by scripts/utils_paths.R. Set them at submission, e.g.
#   CFM_SAMPLE_SET=with_deviations sbatch scripts/07a_hpc_16S_alpha_div_rarefaction.sh
export CFM_VERSION=${CFM_VERSION:-v2}
export CFM_SAMPLE_SET=${CFM_SAMPLE_SET:-main}

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Cores allocated: $SLURM_CPUS_PER_TASK"
echo "CFM_VERSION: $CFM_VERSION, sample set: $CFM_SAMPLE_SET"

Rscript scripts/07a_hpc_16S_alpha_div_rarefaction.R

echo "Job finished: $(date)"
