#!/bin/bash
# Request resources
#SBATCH -p shared           # queue name
#SBATCH --mem=2G            # mem per node, units=M,G,T
#SBATCH -c 1                # CPU cores, up to 128
#SBATCH -n 1                # tasks/ranks (MPI jobs)
#SBATCH -t 01:00:00         # time limit, format dd-hh:mm:ss
#SBATCH --gres=tmp:3G       # TMPDIR space, up to 400G
#SBATCH --mail-type=ALL     # Choose from BEGIN, END, FAIL, ALL
#SBATCH --job-name=trim_adapters
#SBATCH --output=logs/trim_adapters_%j.out
#SBATCH --error=logs/trim_adapters_%j.err

# Script to trim Illumina sequencing adapters from multiplexed reads
# Input: Lane merged raw FASTQ sublibrary files from Novogene
# Output: Adapter-trimmed FASTQ sublibrary files ready for demultiplexing

# Load required modules
module load bioinformatics
module load cutadapt/4.9

# Set input and output directories (constant/configuration variables in capital letters; all directories already created)
INPUT_DIR="qiime2/import/merged_files"
OUTPUT_DIR="qiime2/import/trimmed_reads"

# Get current working directory before changing directories
CURRENT_DIR=$(pwd)

echo "Starting adapter trimming..."
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Processing $(ls $INPUT_DIR/*R1.fq.gz | wc -l) sublibrary pairs"

# Change to input directory
cd $INPUT_DIR

# Process each R1 file and find its corresponding R2 (* is a wildcard; backwards ending: do ... done)
for r1_file in *R1.fq.gz; do
    # Get base filename without _R1.fq.gz (% removes shortest match from the end)
    base_name=${r1_file%_R1.fq.gz} 
    r2_file="${base_name}_R2.fq.gz"
    
    # Check if R2 file exists (if [[ condition ]]; then #do something fi)
    if [[ -f "$r2_file" ]]; then
        echo "Processing: $base_name"
        
        # Define output files using absolute paths to avoid confusion
        output_r1="${CURRENT_DIR}/${OUTPUT_DIR}/${base_name}_R1.trimmed.fq.gz"
        output_r2="${CURRENT_DIR}/${OUTPUT_DIR}/${base_name}_R2.trimmed.fq.gz"
        
        # Run cutadapt to remove YOUR Illumina overhangs/adapters
        # -a: Forward overhang (remove from 3' end of R1 reads)
        # -A: Reverse overhang (remove from 3' end of R2 reads)
        # -o: output R1 file (stays as a separate file)
        # -p: output R2 file (stays as a separate file)
        cutadapt \
            -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
            -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
            -o "$output_r1" \
            -p "$output_r2" \
            "$r1_file" "$r2_file"
        
        # Check exit code of the last command ($?) -eq (equal to) 0 (succeeded)
        if [[ $? -eq 0 ]]; then
            echo "  ✓ Successfully trimmed $base_name"
        else
            echo "  ✗ Error trimming $base_name"
        fi
    else
        echo "Warning: R2 file not found for $r1_file"
    fi
done

echo "Adapter trimming complete!"
echo "Trimmed files saved in: $OUTPUT_DIR"

# Summary statistics
echo ""
echo "Summary:"
echo "Input files: $(ls *R1.fq.gz | wc -l) sublibrary pairs"
echo "Output files: $(ls ${CURRENT_DIR}/${OUTPUT_DIR}/*R1.trimmed.fq.gz 2>/dev/null | wc -l) trimmed pairs"
