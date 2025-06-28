#!/bin/bash
#SBATCH --job-name=dada2_denoise
#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --gres=tmp:50G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/07_dada2_denoising.out
#SBATCH --error=logs/07_dada2_denoising.err

# Script to denoise paired-end sequences using DADA2
# Step 7: Quality filtering, denoising, and ASV calling
# Input: QIIME2 paired-end sequence artifacts
# Output: Feature tables and representative sequences

echo "Starting DADA2 denoising..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 64G"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
INPUT_DIR="qiime2/denoise"
OUTPUT_DIR="qiime2/denoise"

echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Verify input files exist
if [[ ! -f "$INPUT_DIR/CFM_16S_PE_import.qza" ]]; then
    echo "ERROR: 16S import file not found: $INPUT_DIR/CFM_16S_PE_import.qza"
    exit 1
fi

if [[ ! -f "$INPUT_DIR/CFM_ITS1_PE_import.qza" ]]; then
    echo "ERROR: ITS1 import file not found: $INPUT_DIR/CFM_ITS1_PE_import.qza"
    exit 1
fi

echo "Input files verified successfully"

# DADA2 denoising for 16S
echo ""
echo "=== DADA2 DENOISING: 16S rRNA V4 ==="
echo "Starting 16S denoising at: $(date)"
echo "Parameters:"
echo "  Forward truncation: 0 (no truncation - keep full sequences)"
echo "  Reverse truncation: 0 (no truncation - keep full sequences)"
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
    echo "SUCCESS: 16S DADA2 completed at: $(date)"
    
    # Create 16S summary visualizations
    echo "Creating 16S summary visualizations..."
    
    # Feature table summary
    qiime feature-table summarize \
        --i-table "$OUTPUT_DIR/CFM_16S_dada2_table.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_dada2_table_summary.qzv"
    
    # Representative sequences summary
    qiime feature-table tabulate-seqs \
        --i-data "$OUTPUT_DIR/CFM_16S_dada2_repseqs.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_dada2_repseqs_summary.qzv"
    
    # Denoising stats
    qiime metadata tabulate \
        --m-input-file "$OUTPUT_DIR/CFM_16S_dada2_stats.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_dada2_stats_summary.qzv"
    
    echo "16S visualizations created successfully"
else
    echo "ERROR: 16S DADA2 failed"
    exit 1
fi

echo ""

# DADA2 denoising for ITS1
echo "=== DADA2 DENOISING: ITS1 Fungal ==="
echo "Starting ITS1 denoising at: $(date)"
echo "Parameters:"
echo "  Forward truncation: 0 (no truncation - keep full sequences)"
echo "  Reverse truncation: 0 (no truncation - keep full sequences)"
echo "  Cores: $SLURM_CPUS_PER_TASK"

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$INPUT_DIR/CFM_ITS1_PE_import.qza" \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-n-threads $SLURM_CPUS_PER_TASK \
    --p-chimera-method consensus \
    --o-table "$OUTPUT_DIR/CFM_ITS1_dada2_table.qza" \
    --o-representative-sequences "$OUTPUT_DIR/CFM_ITS1_dada2_repseqs.qza" \
    --o-denoising-stats "$OUTPUT_DIR/CFM_ITS1_dada2_stats.qza" \
    --verbose

# Check ITS1 success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: ITS1 DADA2 completed at: $(date)"
    
    # Create ITS1 summary visualizations
    echo "Creating ITS1 summary visualizations..."
    
    # Feature table summary
    qiime feature-table summarize \
        --i-table "$OUTPUT_DIR/CFM_ITS1_dada2_table.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_ITS1_dada2_table_summary.qzv"
    
    # Representative sequences summary
    qiime feature-table tabulate-seqs \
        --i-data "$OUTPUT_DIR/CFM_ITS1_dada2_repseqs.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_ITS1_dada2_repseqs_summary.qzv"
    
    # Denoising stats
    qiime metadata tabulate \
        --m-input-file "$OUTPUT_DIR/CFM_ITS1_dada2_stats.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_ITS1_dada2_stats_summary.qzv"
    
    echo "ITS1 visualizations created successfully"
else
    echo "ERROR: ITS1 DADA2 failed"
    exit 1
fi

echo ""
echo "=== DADA2 DENOISING COMPLETE ==="

# Final summary
echo "Checking output files..."

# List all created files
echo "Files created in $OUTPUT_DIR:"
ls -lah "$OUTPUT_DIR"/CFM_*_dada2_*

# Count features and samples
echo ""
echo "=== FINAL SUMMARY ==="

# 16S summary
if [[ -f "$OUTPUT_DIR/CFM_16S_dada2_table.qza" ]]; then
    echo "16S Results:"
    qiime tools peek "$OUTPUT_DIR/CFM_16S_dada2_table.qza" | grep -E "samples|features"
    
    size_16s_table=$(du -h "$OUTPUT_DIR/CFM_16S_dada2_table.qza" | cut -f1)
    size_16s_seqs=$(du -h "$OUTPUT_DIR/CFM_16S_dada2_repseqs.qza" | cut -f1)
    echo "16S table size: $size_16s_table"
    echo "16S sequences size: $size_16s_seqs"
fi

echo ""

# ITS1 summary
if [[ -f "$OUTPUT_DIR/CFM_ITS1_dada2_table.qza" ]]; then
    echo "ITS1 Results:"
    qiime tools peek "$OUTPUT_DIR/CFM_ITS1_dada2_table.qza" | grep -E "samples|features"
    
    size_its1_table=$(du -h "$OUTPUT_DIR/CFM_ITS1_dada2_table.qza" | cut -f1)
    size_its1_seqs=$(du -h "$OUTPUT_DIR/CFM_ITS1_dada2_repseqs.qza" | cut -f1)
    echo "ITS1 table size: $size_its1_table"
    echo "ITS1 sequences size: $size_its1_seqs"
fi

echo ""
echo "Next steps:"
echo "1. Download and examine denoising stats (.qzv files)"
echo "2. Check feature table summaries"
echo "3. Proceed with taxonomic classification"
echo "4. Run diversity analysis"

echo ""
echo "DADA2 denoising completed at: $(date)"
echo "Ready for taxonomic classification!"
