#!/bin/bash
#SBATCH --job-name=zeta_16S_byfarm
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 14
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --output=logs/zeta/11a_16S_byfarm_%j.out
#SBATCH --error=logs/zeta/11a_16S_byfarm_%j.err
#SBATCH --mail-type=END,FAIL

# 14 cores = one per farm × group combination (7 farms × 2 treatments).
# Wall time 4h is generous; actual runtime expected ~30–60 min.

# zetadiv was installed under R 4.5.1; it loads and runs under the bioinformatics
# R 4.6.1 (tested 2026-11-02), which also loads phyloseq and dplyr from ~/R/library.
# r/4.5.1 alone no longer works: dplyr cannot find rlang.
module load bioinformatics
module load r/4.6.1

# Without this the whole job runs on ONE core, silently. Hamilton sets
# OMP_PROC_BIND=spread in every job; the bioinformatics R then binds its main
# process to one core and all PSOCK workers inherit that binding (seen 2026-11-02
# in the 07a/07b alpha jobs). Same fix as in 03_hpc_rfc.sh.
export OMP_PROC_BIND=false

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/zeta

# Data version (v1 | v2) and sample set (main = 288 flowers, with_deviations = 294),
# read by scripts/utils_paths.R. Set them at submission, e.g.
#   CFM_SAMPLE_SET=with_deviations sbatch scripts/11a_hpc_16S_zeta_by_farm.sh
export CFM_VERSION=${CFM_VERSION:-v2}
export CFM_SAMPLE_SET=${CFM_SAMPLE_SET:-main}
echo "CFM_VERSION: $CFM_VERSION, sample set: $CFM_SAMPLE_SET"

Rscript scripts/11a_hpc_16S_zeta_by_farm.R
