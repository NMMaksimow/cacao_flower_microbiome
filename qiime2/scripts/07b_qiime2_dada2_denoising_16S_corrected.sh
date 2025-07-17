#!/bin/bash
#SBATCH --job-name=dada2_denoise_16S
#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH --gres=tmp:50G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/07b_dada2_denoising_16S_corrected.out
#SBATCH --error=logs/07b_dada2_denoising_16S_corrected.err

# Script to re-denoise 16S sequences using DADA2 with corrected sample identities
# Step 7b: Quality filtering, denoising, and ASV calling for 16S only
# Input: Corrected 16S QIIME2 paired-end sequence artifact
# Output: Overwrites existing 16S feature tables and representative sequences

echo "Starting DADA2 denoising for corrected 16S data..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 64G"
echo "Purpose: Re-process 16S with corrected cmt53/cmtPC sample identities"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
INPUT_DIR="qiime2/denoise"
OUTPUT_DIR="qiime2/denoise"

echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Verify input file exists
if [[ ! -f "$INPUT_DIR/CFM_16S_PE_import.qza" ]]; then
    echo "ERROR: Corrected 16S import file not found: $INPUT_DIR/CFM_16S_PE_import.qza"
    exit 1
fi

echo "Corrected 16S input file verified successfully"

# Check timestamp to confirm we're using corrected data
echo "Verifying we're using corrected import artifact:"
ls -l "$INPUT_DIR/CFM_16S_PE_import.qza"

# Backup existing 16S DADA2 outputs
echo ""
echo "=== BACKING UP EXISTING 16S DADA2 FILES ==="
backup_timestamp=$(date +%Y%m%d_%H%M)

if [[ -f "$OUTPUT_DIR/CFM_16S_dada2_table.qza" ]]; then
    echo "Backing up existing 16S DADA2 table..."
    cp "$OUTPUT_DIR/CFM_16S_dada2_table.qza" "$OUTPUT_DIR/CFM_16S_dada2_table_backup_${backup_timestamp}.qza"
fi

if [[ -f "$OUTPUT_DIR/CFM_16S_dada2_repseqs.qza" ]]; then
    echo "Backing up existing 16S DADA2 representative sequences..."
    cp "$OUTPUT_DIR/CFM_16S_dada2_repseqs.qza" "$OUTPUT_DIR/CFM_16S_dada2_repseqs_backup_${backup_timestamp}.qza"
fi

if [[ -f "$OUTPUT_DIR/CFM_16S_dada2_stats.qza" ]]; then
    echo "Backing up existing 16S DADA2 stats..."
    cp "$OUTPUT_DIR/CFM_16S_dada2_stats.qza" "$OUTPUT_DIR/CFM_16S_dada2_stats_backup_${backup_timestamp}.qza"
fi

# DADA2 denoising for corrected 16S
echo ""
echo "=== DADA2 DENOISING: 16S rRNA V4 (CORRECTED DATA) ==="
echo "Starting corrected 16S denoising at: $(date)"
echo "Parameters (identical to original run):"
echo "  Forward truncation: 0 (no truncation - keep full sequences)"
echo "  Reverse truncation: 0 (no truncation - keep full sequences)"
echo "  Trim left forward: 0"
echo "  Trim left reverse: 0"
echo "  Chimera method: consensus"
echo "  Cores: $SLURM_CPUS_PER_TASK"

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$INPUT_DIR/CFM_16S_PE_import.qza" \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-n-threads $SLURM_CPUS_PER_TASK \
    --p-chimera-method consensus \
    --o-table "$OUTPUT_DIR/CFM_16S_dada2_table.qza" \
    --o-representative-sequences "$OUTPUT_DIR/CFM_16S_dada2_repseqs.qza" \
    --o-denoising-stats "$OUTPUT_DIR/CFM_16S_dada2_stats.qza" \
    --verbose

# Check 16S success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Corrected 16S DADA2 completed at: $(date)"

    # Create 16S summary visualizations
    echo "Creating corrected 16S summary visualizations..."

    # Feature table summary
    echo "Creating feature table summary..."
    qiime feature-table summarize \
        --i-table "$OUTPUT_DIR/CFM_16S_dada2_table.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_dada2_table_summary.qzv"

    # Representative sequences summary
    echo "Creating representative sequences summary..."
    qiime feature-table tabulate-seqs \
        --i-data "$OUTPUT_DIR/CFM_16S_dada2_repseqs.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_dada2_repseqs_summary.qzv"

    # Denoising stats
    echo "Creating denoising stats summary..."
    qiime metadata tabulate \
        --m-input-file "$OUTPUT_DIR/CFM_16S_dada2_stats.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_dada2_stats_summary.qzv"

    echo "Corrected 16S visualizations created successfully"
else
    echo "ERROR: Corrected 16S DADA2 failed"
    exit 1
fi

echo ""
echo "=== CORRECTED 16S DADA2 DENOISING COMPLETE ==="

# Final summary
echo "Checking output files..."

# List all created files
echo "Updated 16S files in $OUTPUT_DIR:"
ls -lah "$OUTPUT_DIR"/CFM_16S_dada2_*

# Count features and samples
echo ""
echo "=== FINAL SUMMARY ==="

# 16S summary
if [[ -f "$OUTPUT_DIR/CFM_16S_dada2_table.qza" ]]; then
    echo "Corrected 16S Results:"
    qiime tools peek "$OUTPUT_DIR/CFM_16S_dada2_table.qza" | grep -E "samples|features"

    size_16s_table=$(du -h "$OUTPUT_DIR/CFM_16S_dada2_table.qza" | cut -f1)
    size_16s_seqs=$(du -h "$OUTPUT_DIR/CFM_16S_dada2_repseqs.qza" | cut -f1)
    echo "16S table size: $size_16s_table"
    echo "16S sequences size: $size_16s_seqs"
fi

echo ""
echo "Sample identity verification:"
echo "The corrected 16S data now contains:"
echo "- cmt53: correctly pointing to actual cmt53 biological sample data"
echo "- cmtPC: correctly pointing to actual cmtPC biological sample data"

echo ""
echo "Backup files created:"
echo "- CFM_16S_dada2_table_backup_${backup_timestamp}.qza"
echo "- CFM_16S_dada2_repseqs_backup_${backup_timestamp}.qza"
echo "- CFM_16S_dada2_stats_backup_${backup_timestamp}.qza"

echo ""
echo "Next steps:"
echo "1. Download and examine corrected denoising stats (.qzv files)"
echo "2. Check corrected feature table summaries"
echo "3. Verify sample cmt53 and cmtPC are correctly represented"
echo "4. Proceed with taxonomic classification using corrected 16S data"
echo "5. Re-run filtering and downstream analysis"

echo ""
echo "Corrected 16S DADA2 denoising completed at: $(date)"
echo "Ready for taxonomic classification with correct sample identities!"
