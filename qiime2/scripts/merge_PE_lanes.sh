#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=4G
#SBATCH -c 2
#SBATCH -n 1
#SBATCH -t 01:00:00

# Script to merge L1 and L2 files for each sub-library
# Input: data/raw_data/X204SC25042477-Z01-F001/01.RawData/
# Output: qiime2/import/merged_files/

echo "Starting lane merging..."
echo "Job started at: $(date)"

# Create output directory
mkdir -p qiime2/import/merged_files

# Set input directory
INPUT_DIR="data/raw_data/X204SC25042477-Z01-F001/01.RawData"

# Loop through each sub-library directory
for sublibrary in $(ls $INPUT_DIR); do
    echo "Processing sublibrary: $sublibrary"
    
    # Define file paths
    L1_R1="${INPUT_DIR}/${sublibrary}/${sublibrary}_FKDL250197829-1A_HVGYYDRX5_L1_1.fq.gz"
    L1_R2="${INPUT_DIR}/${sublibrary}/${sublibrary}_FKDL250197829-1A_HVGYYDRX5_L1_2.fq.gz"
    L2_R1="${INPUT_DIR}/${sublibrary}/${sublibrary}_FKDL250197829-1A_HWYVWDRX5_L2_1.fq.gz"
    L2_R2="${INPUT_DIR}/${sublibrary}/${sublibrary}_FKDL250197829-1A_HWYVWDRX5_L2_2.fq.gz"
    
    # Output paths
    MERGED_R1="qiime2/import/merged_files/${sublibrary}_R1.fq.gz"
    MERGED_R2="qiime2/import/merged_files/${sublibrary}_R2.fq.gz"
    
    # Merge forward reads (L1_R1 + L2_R1)
    echo "  Merging forward reads..."
    cat "$L1_R1" "$L2_R1" > "$MERGED_R1"
    
    # Merge reverse reads (L1_R2 + L2_R2)
    echo "  Merging reverse reads..."
    cat "$L1_R2" "$L2_R2" > "$MERGED_R2"
    
    echo "  Completed: $sublibrary"
done

echo "Lane merging completed!"
echo "Merged files are in: qiime2/import/merged_files/"
echo "Job completed at: $(date)"
