#!/bin/bash
#SBATCH --job-name=net_spieceasi
#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --output=logs/network/19_spieceasi_%j.out
#SBATCH --error=logs/network/19_spieceasi_%j.err
#SBATCH --mail-type=END,FAIL

# SpiecEasi-MB networks, 16S (19a) + ITS1 (19b), nPerm = 200.
# netCompare dominates runtime: measured ~2.5-3 h per farm at nPerm=1000,
# so ~0.5-0.6 h per farm at 200 -> ~8 h for 14 farm-pairs. 24 h is ~3x headroom.
# 1 core: the timing probe (script 22) showed netCompare's parallel workers
# are about 2x SLOWER than serial, so extra cores only cost queue time.

module load bioinformatics
module load r/4.6.0

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/network results/rds results/figures results/tables

Rscript scripts/19a_16S_network_by_farm_SpiecEasi.R
Rscript scripts/19b_ITS1_network_by_farm_SpiecEasi.R

Rscript -e "sink('logs/network/19_session_info.txt'); sessionInfo(); sink()"
