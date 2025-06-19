#!/bin/bash
#SBATCH --job-name=demux_stacks
#SBATCH -p shared
#SBATCH --mem=8G
#SBATCH -c 4
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --gres=tmp:10G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/04_demultiplex_stacks_no_cqr.out
#SBATCH --error=logs/04_demultiplex_stacks_no_cqr.err

# Script to demultiplex trimmed sublibrary files using Stacks process_radtags
# Input: 48 adapter-trimmed sublibrary pairs + 48 mapping files
# Output: Individual sample files separated by amplicon (16S and ITS1)

echo "Starting Stacks demultiplexing..."
echo "Date: $(date)"
echo "Host: $(hostname)"

# Load required modules
module load bioinformatics
module load stacks/2.64

# Set directories
INPUT_DIR="qiime2/import/trimmed_reads"
MAPPING_DIR="qiime2/import/demux/internal_tag_mappings"
OUTPUT_16S="qiime2/import/demux/demultiplexed_sample_files_no_cqr/16S"
OUTPUT_ITS1="qiime2/import/demux/demultiplexed_sample_files_no_cqr/ITS1"

echo "Input directory: $INPUT_DIR"
echo "Mapping directory: $MAPPING_DIR"
echo "Output 16S: $OUTPUT_16S"
echo "Output ITS1: $OUTPUT_ITS1"

# Count total sublibraries to process
total_sublibraries=$(ls $INPUT_DIR/*R1.trimmed.fq.gz | wc -l)
echo "Total sublibraries to process: $total_sublibraries"

# Initialize counters
processed=0
successful=0
failed=0

# Process each sublibrary
for r1_file in $INPUT_DIR/*R1.trimmed.fq.gz; do
    # Extract sublibrary name
    base_name=$(basename "$r1_file" _R1.trimmed.fq.gz)
    r2_file="$INPUT_DIR/${base_name}_R2.trimmed.fq.gz"
    mapping_file="$MAPPING_DIR/${base_name}_tags.tsv"
    
    processed=$((processed + 1))
    echo ""
    echo "=== Processing sublibrary $processed/$total_sublibraries: $base_name ==="
    echo "R1 file: $r1_file"
    echo "R2 file: $r2_file"
    echo "Mapping: $mapping_file"
    
    # Check if all required files exist
    if [[ ! -f "$r2_file" ]]; then
        echo "ERROR: R2 file not found: $r2_file"
        failed=$((failed + 1))
        continue
    fi
    
    if [[ ! -f "$mapping_file" ]]; then
        echo "ERROR: Mapping file not found: $mapping_file"
        failed=$((failed + 1))
        continue
    fi
    
    # Show mapping file stats
    mapping_entries=$(wc -l < "$mapping_file")
    echo "Mapping entries: $mapping_entries"
    echo "Expected samples: $((mapping_entries / 2)) (16S + ITS1 per sample)"
    
    # Create temporary output directory for this sublibrary
    temp_output="temp_${base_name}"
    mkdir -p "$temp_output"
    
    echo "Running process_radtags..."
    
    # Run process_radtags
    process_radtags \
        -P \
        -1 "$r1_file" \
        -2 "$r2_file" \
        -b "$mapping_file" \
        -o "$temp_output" \
        --inline_inline \
        --disable-rad-check \
        --retain-header
    
    # Check if process_radtags succeeded
    if [[ $? -eq 0 ]]; then
        echo "SUCCESS: process_radtags completed successfully for $base_name"
        
        # Count and organize output files
        output_files=$(ls "$temp_output"/*.1.fq.gz "$temp_output"/*.2.fq.gz 2>/dev/null | wc -l)
        echo "Output files created: $output_files"
        
        # Separate files by amplicon and move to final directories
        echo "Organizing files by amplicon..."
        
        moved_16s=0
        moved_its1=0
        
        # Move 16S files
        for file in "$temp_output"/*_16S.*.fq.gz; do
            if [[ -f "$file" ]]; then
                mv "$file" "$OUTPUT_16S/"
                moved_16s=$((moved_16s + 1))
            fi
        done
        
        # Move ITS1 files
        for file in "$temp_output"/*_ITS1.*.fq.gz; do
            if [[ -f "$file" ]]; then
                mv "$file" "$OUTPUT_ITS1/"
                moved_its1=$((moved_its1 + 1))
            fi
        done
        
        echo "Moved to 16S directory: $moved_16s files"
        echo "Moved to ITS1 directory: $moved_its1 files"
        
        # Clean up temporary directory
        rm -rf "$temp_output"
        
        successful=$((successful + 1))
        
    else
        echo "ERROR: process_radtags failed for $base_name"
        echo "Check temporary directory: $temp_output"
        failed=$((failed + 1))
    fi
    
    echo "Progress: $processed/$total_sublibraries processed ($successful successful, $failed failed)"
done

echo ""
echo "=== DEMULTIPLEXING COMPLETE ==="
echo "Total sublibraries processed: $processed"
echo "Successful: $successful"
echo "Failed: $failed"

# Final file counts
echo ""
echo "=== FINAL OUTPUT SUMMARY ==="
s16_files=$(ls "$OUTPUT_16S"/*.fq.gz 2>/dev/null | wc -l)
its1_files=$(ls "$OUTPUT_ITS1"/*.fq.gz 2>/dev/null | wc -l)

echo "16S files created: $s16_files"
echo "ITS1 files created: $its1_files"
echo "Total sample files: $((s16_files + its1_files))"

# Expected: ~336 samples × 2 amplicons × 2 read directions = ~1,344 files
expected_total=$((336 * 2 * 2))
echo "Expected total files: ~$expected_total"

if [[ $((s16_files + its1_files)) -gt 1000 ]]; then
    echo "SUCCESS: Demultiplexing appears successful!"
    echo "Ready for QIIME2 import!"
else
    echo "WARNING: Lower than expected file count - check for issues"
fi

echo "Demultiplexing completed at: $(date)"
