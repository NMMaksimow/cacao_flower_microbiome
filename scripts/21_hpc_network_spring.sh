#!/bin/bash
#SBATCH --job-name=net_spring
#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH --output=logs/network/21_spring_%j.out
#SBATCH --error=logs/network/21_spring_%j.err
#SBATCH --mail-type=END,FAIL

# SPRING networks, 16S (21a) + ITS1 (21b), PERM_TEST = FALSE.
# netCompare's permutation test fails for SPRING on every farm (nearPD does not
# converge on permuted data), so SPRING is descriptive: topology and hubs, no
# p-values. Without permutations the cost is construction only, which measured
# ~2 h for both markers, so 6 h is ample.
# 1 core: parallel workers measured ~2x slower.

module load bioinformatics
module load r/4.6.1

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/network

# Data version (v1 | v2) and sample set (main = 288 flowers, with_deviations = 294),
# read by scripts/utils_paths.R. Set them at submission, e.g.
#   CFM_SAMPLE_SET=with_deviations sbatch scripts/21_hpc_network_spring.sh
export CFM_VERSION=${CFM_VERSION:-v2}
export CFM_SAMPLE_SET=${CFM_SAMPLE_SET:-main}

echo "CFM_VERSION: $CFM_VERSION, sample set: $CFM_SAMPLE_SET"

Rscript scripts/21a_16S_network_by_farm_SPRING.R
Rscript scripts/21b_ITS1_network_by_farm_SPRING.R

Rscript -e "sink('logs/network/21_session_info.txt'); sessionInfo(); sink()"
