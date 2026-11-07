#!/bin/bash
#SBATCH --job-name=ancom_16S_farm
#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --output=logs/ancombc2/19a_16S_farm_%j.out
#SBATCH --error=logs/ancombc2/19a_16S_farm_%j.err
#SBATCH --mail-type=END,FAIL

# 16 cores = N_CORES in the R script. ANCOMBC 2.14.0, as on the laptop, in ~/R/library.
module load bioinformatics
module load r/4.6.1

# Without this the job runs on one core: see 11a_hpc_16S_zeta_by_farm.sh
export OMP_PROC_BIND=false

cd ~/nobackup/cacao_flower_microbiome
mkdir -p logs/ancombc2

# Data version (v1 | v2) and sample set (main = 288 flowers, with_deviations = 294),
# read by scripts/utils_paths.R. Set them at submission, e.g.
#   CFM_SAMPLE_SET=with_deviations sbatch scripts/19a_hpc_16S_ancombc2_farm.sh
export CFM_VERSION=${CFM_VERSION:-v2}
export CFM_SAMPLE_SET=${CFM_SAMPLE_SET:-main}
echo "CFM_VERSION: $CFM_VERSION, sample set: $CFM_SAMPLE_SET"

Rscript scripts/19a_hpc_16S_ancombc2_farm.R
