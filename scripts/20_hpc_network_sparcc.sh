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
module load r/4.6.1

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/network

# Data version (v1 | v2) and sample set (main = 288 flowers, with_deviations = 294),
# read by scripts/utils_paths.R. Set them at submission, e.g.
#   CFM_SAMPLE_SET=with_deviations sbatch scripts/20_hpc_network_sparcc.sh
export CFM_VERSION=${CFM_VERSION:-v2}
export CFM_SAMPLE_SET=${CFM_SAMPLE_SET:-main}

echo "CFM_VERSION: $CFM_VERSION, sample set: $CFM_SAMPLE_SET"

Rscript scripts/20a_16S_network_by_farm_SparCC.R
Rscript scripts/20b_ITS1_network_by_farm_SparCC.R

Rscript -e "sink('logs/network/20_session_info.txt'); sessionInfo(); sink()"
