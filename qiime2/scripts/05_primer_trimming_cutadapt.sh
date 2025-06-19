#!/bin/bash
#SBATCH --job-name=primer_trim
#SBATCH -p shared
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH --gres=tmp:20G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/05_primer_trimming.out
#SBATCH --error=logs/05_primer_trimming.err

# Script to trim locus-specific primers from demultiplexed samples
# Step 5: Remove 16S and ITS1 primer sequences using cutadapt
# Input: Demultiplexed sample files separated by amplicon
# Output: Primer-trimmed files ready for QIIME2 import

echo "Starting primer trimming with cutadapt..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"

# Load required modules
module load bioinformatics
module load cutadapt/4.9

# Define primer sequences (locus-specific parts only)
# 16S V4 region primers (515f_modified/806r_modified)
PRIMER_16S_FWD="GTGYCAGCMGCCGCGGTAA"     # 515f_modified
PRIMER_16S_REV="GGACTACNVGGGTWTCTAAT"    # 806r_modified

# ITS1 region primers (ITS1f/ITS2)
PRIMER_ITS1_FWD="CTTGGTCATTTAGAGGAAGTAA"  # ITS1f
PRIMER_ITS1_REV="GCTGCGTTCTTCATCGATGC"    # ITS2

echo "Primer sequences:"
echo "16S Forward: $PRIMER_16S_FWD"
echo "16S Reverse: $PRIMER_16S_REV"
echo "ITS1 Forward: $PRIMER_ITS1_FWD"
echo "ITS1 Reverse: $PRIMER_ITS1_REV"

# Set directories
INPUT_16S="qiime2/import/demux/demultiplexed_sample_files/16S"
INPUT_ITS1="qiime2/import/demux/demultiplexed_sample_files/ITS1"
OUTPUT_16S="qiime2/import/primer_trimmed/16S"
OUTPUT_ITS1="qiime2/import/primer_trimmed/ITS1"

# Create output directories
mkdir -p "$OUTPUT_16S"
mkdir -p "$OUTPUT_ITS1"

echo "Input 16S: $INPUT_16S"
echo "Input ITS1: $INPUT_ITS1"
echo "Output 16S: $OUTPUT_16S"
echo "Output ITS1: $OUTPUT_ITS1"

# Initialize counters
total_16s_samples=0
total_its1_samples=0
processed_16s=0
processed_its1=0
failed_16s=0
failed_its1=0

# Function to trim primers for a sample pair
trim_primers() {
    local r1_file=$1
    local r2_file=$2
    local output_dir=$3
    local fwd_primer=$4
    local rev_primer=$5
    local amplicon=$6
    
    local sample_name=$(basename "$r1_file" .1.fq.gz)
    local output_r1="${output_dir}/${sample_name}.1.trimmed.fq.gz"
    local output_r2="${output_dir}/${sample_name}.2.trimmed.fq.gz"
    
    echo "  Processing: $sample_name"
    echo "    Input R1: $r1_file"
    echo "    Input R2: $r2_file"
    echo "    Output R1: $output_r1"
    echo "    Output R2: $output_r2"
    
    # Count input reads
    input_reads_r1=$(zcat "$r1_file" | wc -l)
    input_reads_r1=$((input_reads_r1 / 4))
    
    # Run cutadapt with paired-end primer trimming
    # Note: --quality-cutoff removed to let QIIME2 handle quality filtering
    cutadapt \
        -g "$fwd_primer" \
        -G "$rev_primer" \
        -a "$rev_primer" \
        -A "$fwd_primer" \
        --match-read-wildcards \
        --minimum-length 50 \
        --maximum-length 300 \
        --cores $SLURM_CPUS_PER_TASK \
        --output "$output_r1" \
        --paired-output "$output_r2" \
        "$r1_file" "$r2_file"
    
    # Check if cutadapt succeeded
    if [[ $? -eq 0 && -f "$output_r1" && -f "$output_r2" ]]; then
        # Count output reads
        output_reads_r1=$(zcat "$output_r1" | wc -l)
        output_reads_r1=$((output_reads_r1 / 4))
        
        retention_rate=$(echo "scale=2; $output_reads_r1 * 100 / $input_reads_r1" | bc)
        
        echo "    SUCCESS: $input_reads_r1 → $output_reads_r1 reads (${retention_rate}% retained)"
        return 0
    else
        echo "    ERROR: cutadapt failed for $sample_name"
        return 1
    fi
}

echo ""
echo "=== PROCESSING 16S SAMPLES ==="

# Count 16S sample pairs (exclude .rem files)
total_16s_samples=$(ls "$INPUT_16S"/*.1.fq.gz 2>/dev/null | grep -v "\.rem\." | wc -l)
echo "Total 16S sample pairs to process: $total_16s_samples"

# Process 16S samples
for r1_file in "$INPUT_16S"/*.1.fq.gz; do
    # Skip .rem files
    if [[ "$r1_file" == *".rem."* ]]; then
        continue
    fi
    
    # Check if file exists
    if [[ ! -f "$r1_file" ]]; then
        continue
    fi
    
    # Get corresponding R2 file
    r2_file="${r1_file/.1.fq.gz/.2.fq.gz}"
    
    if [[ ! -f "$r2_file" ]]; then
        echo "ERROR: R2 file not found for $r1_file"
        failed_16s=$((failed_16s + 1))
        continue
    fi
    
    processed_16s=$((processed_16s + 1))
    echo ""
    echo "16S Sample $processed_16s/$total_16s_samples"
    
    if trim_primers "$r1_file" "$r2_file" "$OUTPUT_16S" "$PRIMER_16S_FWD" "$PRIMER_16S_REV" "16S"; then
        echo "16S sample processed successfully"
    else
        failed_16s=$((failed_16s + 1))
    fi
    
    # Progress update every 50 samples
    if [[ $((processed_16s % 50)) -eq 0 ]]; then
        echo "16S Progress: $processed_16s/$total_16s_samples processed"
    fi
done

echo ""
echo "=== PROCESSING ITS1 SAMPLES ==="

# Count ITS1 sample pairs (exclude .rem files)
total_its1_samples=$(ls "$INPUT_ITS1"/*.1.fq.gz 2>/dev/null | grep -v "\.rem\." | wc -l)
echo "Total ITS1 sample pairs to process: $total_its1_samples"

# Process ITS1 samples
for r1_file in "$INPUT_ITS1"/*.1.fq.gz; do
    # Skip .rem files
    if [[ "$r1_file" == *".rem."* ]]; then
        continue
    fi
    
    # Check if file exists
    if [[ ! -f "$r1_file" ]]; then
        continue
    fi
    
    # Get corresponding R2 file
    r2_file="${r1_file/.1.fq.gz/.2.fq.gz}"
    
    if [[ ! -f "$r2_file" ]]; then
        echo "ERROR: R2 file not found for $r1_file"
        failed_its1=$((failed_its1 + 1))
        continue
    fi
    
    processed_its1=$((processed_its1 + 1))
    echo ""
    echo "ITS1 Sample $processed_its1/$total_its1_samples"
    
    if trim_primers "$r1_file" "$r2_file" "$OUTPUT_ITS1" "$PRIMER_ITS1_FWD" "$PRIMER_ITS1_REV" "ITS1"; then
        echo "ITS1 sample processed successfully"
    else
        failed_its1=$((failed_its1 + 1))
    fi
    
    # Progress update every 50 samples
    if [[ $((processed_its1 % 50)) -eq 0 ]]; then
        echo "ITS1 Progress: $processed_its1/$total_its1_samples processed"
    fi
done

echo ""
echo "=== PRIMER TRIMMING COMPLETE ==="
echo "16S Results:"
echo "  Total samples: $total_16s_samples"
echo "  Successfully processed: $((processed_16s - failed_16s))"
echo "  Failed: $failed_16s"

echo ""
echo "ITS1 Results:"
echo "  Total samples: $total_its1_samples"
echo "  Successfully processed: $((processed_its1 - failed_its1))"
echo "  Failed: $failed_its1"

# Final file counts
echo ""
echo "=== OUTPUT SUMMARY ==="
trimmed_16s=$(ls "$OUTPUT_16S"/*.trimmed.fq.gz 2>/dev/null | wc -l)
trimmed_its1=$(ls "$OUTPUT_ITS1"/*.trimmed.fq.gz 2>/dev/null | wc -l)

echo "16S trimmed files created: $trimmed_16s"
echo "ITS1 trimmed files created: $trimmed_its1"
echo "Total trimmed files: $((trimmed_16s + trimmed_its1))"

# Expected: 336 samples × 2 amplicons × 2 read directions = 1,344 files
expected_total=$((336 * 2 * 2))
echo "Expected total files: $expected_total"

if [[ $((trimmed_16s + trimmed_its1)) -gt 1200 ]]; then
    echo "SUCCESS: Primer trimming appears successful!"
    echo "Ready for QIIME2 import!"
else
    echo "WARNING: Lower than expected file count - check for issues"
fi

echo ""
echo "Next steps:"
echo "1. Check quality of trimmed sequences"
echo "2. Import into QIIME2"
echo "3. Run DADA2 denoising"

echo "Primer trimming completed at: $(date)"
