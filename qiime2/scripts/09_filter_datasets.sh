#!/bin/bash
#SBATCH --job-name=filter_datasets
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --gres=tmp:20G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/09_filter_datasets.out
#SBATCH --error=logs/09_filter_datasets.err

# Script to filter datasets removing contaminants and non-target taxa
# Step 9: Dataset filtering (plastids/mitochondria from 16S, non-fungi from ITS1)
# Input: DADA2 feature tables and taxonomy
# Output: Filtered feature tables and representative sequences

echo "Starting dataset filtering..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 32G"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
INPUT_DIR_DENOISE="qiime2/denoise"
INPUT_DIR_TAXONOMY="qiime2/taxonomy"
OUTPUT_DIR="qiime2/filtered"

echo "Input denoise directory: $INPUT_DIR_DENOISE"
echo "Input taxonomy directory: $INPUT_DIR_TAXONOMY"
echo "Output directory: $OUTPUT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Verify input files exist
if [[ ! -f "$INPUT_DIR_DENOISE/CFM_16S_dada2_table.qza" ]]; then
    echo "ERROR: 16S feature table not found: $INPUT_DIR_DENOISE/CFM_16S_dada2_table.qza"
    exit 1
fi

if [[ ! -f "$INPUT_DIR_TAXONOMY/CFM_16S_taxonomy.qza" ]]; then
    echo "ERROR: 16S taxonomy not found: $INPUT_DIR_TAXONOMY/CFM_16S_taxonomy.qza"
    exit 1
fi

echo "Input files verified successfully"

# ============================================================================
# 16S BACTERIAL DATASET FILTERING
# ============================================================================

echo ""
echo "=== 16S BACTERIAL DATASET FILTERING ==="
echo "Starting 16S filtering at: $(date)"
echo "Removing plastids (chloroplasts) and mitochondria..."

# Remove plastids (chloroplasts) and mitochondria from 16S data
qiime taxa filter-table \
    --i-table "$INPUT_DIR_DENOISE/CFM_16S_dada2_table.qza" \
    --i-taxonomy "$INPUT_DIR_TAXONOMY/CFM_16S_taxonomy.qza" \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table "$OUTPUT_DIR/CFM_16S_filtered_table.qza" \
    --verbose

# Check 16S table filtering success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: 16S table filtering completed"

    # Filter representative sequences to match filtered table
    echo "Filtering 16S representative sequences..."
    qiime taxa filter-seqs \
        --i-sequences "$INPUT_DIR_DENOISE/CFM_16S_dada2_repseqs.qza" \
        --i-taxonomy "$INPUT_DIR_TAXONOMY/CFM_16S_taxonomy.qza" \
        --p-exclude mitochondria,chloroplast \
        --o-filtered-sequences "$OUTPUT_DIR/CFM_16S_filtered_repseqs.qza" \
        --verbose

    # Generate summary of filtered 16S table
    echo "Creating 16S filtered table summary..."
    qiime feature-table summarize \
        --i-table "$OUTPUT_DIR/CFM_16S_filtered_table.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_filtered_table_summary.qzv"

    echo "16S filtering complete. Plastids and mitochondria removed."
else
    echo "ERROR: 16S table filtering failed"
    exit 1
fi

echo ""

# ============================================================================
# ITS1 FUNGAL DATASET FILTERING
# ============================================================================

echo "=== ITS1 FUNGAL DATASET FILTERING ==="
echo "Starting ITS1 filtering at: $(date)"
echo "Keeping only Kingdom Fungi..."

# Keep only Kingdom Fungi from ITS1 data
qiime taxa filter-table \
    --i-table "$INPUT_DIR_DENOISE/CFM_ITS1_dada2_table.qza" \
    --i-taxonomy "$INPUT_DIR_TAXONOMY/CFM_ITS1_taxonomy.qza" \
    --p-include k__Fungi \
    --o-filtered-table "$OUTPUT_DIR/CFM_ITS1_filtered_table.qza" \
    --verbose

# Check ITS1 table filtering success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: ITS1 table filtering completed"

    # Filter representative sequences to match filtered table
    echo "Filtering ITS1 representative sequences..."
    qiime taxa filter-seqs \
        --i-sequences "$INPUT_DIR_DENOISE/CFM_ITS1_dada2_repseqs.qza" \
        --i-taxonomy "$INPUT_DIR_TAXONOMY/CFM_ITS1_taxonomy.qza" \
        --p-include k__Fungi \
        --o-filtered-sequences "$OUTPUT_DIR/CFM_ITS1_filtered_repseqs.qza" \
        --verbose

    # Generate summary of filtered ITS1 table
    echo "Creating ITS1 filtered table summary..."
    qiime feature-table summarize \
        --i-table "$OUTPUT_DIR/CFM_ITS1_filtered_table.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_ITS1_filtered_table_summary.qzv"

    echo "ITS1 filtering complete. Only Fungi retained."
else
    echo "ERROR: ITS1 table filtering failed"
    exit 1
fi

echo ""

# ============================================================================
# ADDITIONAL FILTERING: SINGLETONS AND LOW-PREVALENCE ASVs
# ============================================================================

echo "=== ADDITIONAL FILTERING ==="
echo "Removing singleton ASVs (frequency < 2)..."

# 16S singletons removal
qiime feature-table filter-features \
    --i-table "$OUTPUT_DIR/CFM_16S_filtered_table.qza" \
    --p-min-frequency 2 \
    --o-filtered-table "$OUTPUT_DIR/CFM_16S_filtered_no_singletons.qza" \
    --verbose

# ITS1 singletons removal
qiime feature-table filter-features \
    --i-table "$OUTPUT_DIR/CFM_ITS1_filtered_table.qza" \
    --p-min-frequency 2 \
    --o-filtered-table "$OUTPUT_DIR/CFM_ITS1_filtered_no_singletons.qza" \
    --verbose

echo "Removing low-prevalence ASVs (present in <5% of samples)..."

# Calculate 5% of 336 samples = ~17 samples
MIN_SAMPLES=17

# 16S low-prevalence removal
qiime feature-table filter-features \
    --i-table "$OUTPUT_DIR/CFM_16S_filtered_no_singletons.qza" \
    --p-min-samples $MIN_SAMPLES \
    --o-filtered-table "$OUTPUT_DIR/CFM_16S_final_filtered.qza" \
    --verbose

# ITS1 low-prevalence removal
qiime feature-table filter-features \
    --i-table "$OUTPUT_DIR/CFM_ITS1_filtered_no_singletons.qza" \
    --p-min-samples $MIN_SAMPLES \
    --o-filtered-table "$OUTPUT_DIR/CFM_ITS1_final_filtered.qza" \
    --verbose

# Generate final summaries
echo "Creating final filtered table summaries..."

qiime feature-table summarize \
    --i-table "$OUTPUT_DIR/CFM_16S_final_filtered.qza" \
    --o-visualization "$OUTPUT_DIR/CFM_16S_final_filtered_summary.qzv"

qiime feature-table summarize \
    --i-table "$OUTPUT_DIR/CFM_ITS1_final_filtered.qza" \
    --o-visualization "$OUTPUT_DIR/CFM_ITS1_final_filtered_summary.qzv"

# ============================================================================
# OPTIONAL: PROKARYOTES-ONLY 16S ANALYSIS (for comparison)
# ============================================================================

echo ""
echo "=== CREATING PROKARYOTES-ONLY 16S DATASET ==="
echo "Purpose: Standard bacterial microbiome analysis excluding eukaryotic DNA"
echo "Note: This will remove Phytophthora sequences found in the 16S data"

# Remove all Eukaryota from 16S data for standard bacterial analysis
echo "Removing Domain Eukaryota from 16S dataset..."
qiime taxa filter-table \
    --i-table "$OUTPUT_DIR/CFM_16S_final_filtered.qza" \
    --i-taxonomy "$INPUT_DIR_TAXONOMY/CFM_16S_taxonomy.qza" \
    --p-exclude "d__Eukaryota" \
    --o-filtered-table "$OUTPUT_DIR/CFM_16S_prokaryotes_only.qza" \
    --verbose

# Check prokaryotes-only filtering success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Prokaryotes-only table created"

    # Filter representative sequences to match prokaryotes-only table
    echo "Filtering representative sequences for prokaryotes-only dataset..."
    qiime taxa filter-seqs \
        --i-sequences "$OUTPUT_DIR/CFM_16S_filtered_repseqs.qza" \
        --i-taxonomy "$INPUT_DIR_TAXONOMY/CFM_16S_taxonomy.qza" \
        --p-exclude "d__Eukaryota" \
        --o-filtered-sequences "$OUTPUT_DIR/CFM_16S_prokaryotes_only_repseqs.qza" \
        --verbose

    # Generate summary for prokaryotes-only dataset
    echo "Creating prokaryotes-only summary..."
    qiime feature-table summarize \
        --i-table "$OUTPUT_DIR/CFM_16S_prokaryotes_only.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_prokaryotes_only_summary.qzv"

    echo "Prokaryotes-only dataset created successfully"
else
    echo "ERROR: Prokaryotes-only filtering failed"
    echo "Continuing with standard analysis..."
fi

echo ""
echo "=== DUAL 16S ANALYSIS APPROACH COMPLETE ==="
echo "Two 16S datasets now available:"
echo "1. CFM_16S_final_filtered.qza - includes eukaryotic DNA (Phytophthora)"
echo "2. CFM_16S_prokaryotes_only.qza - standard bacterial microbiome"

echo ""
echo "=== DATASET FILTERING COMPLETE ==="

# Final summary
echo "Checking output files..."

# List all created files
echo "Files created in $OUTPUT_DIR:"
ls -lah "$OUTPUT_DIR"/CFM_*

echo ""
echo "=== FILTERING SUMMARY ==="

# 16S summary - both datasets
if [[ -f "$OUTPUT_DIR/CFM_16S_final_filtered.qza" ]]; then
    echo "16S Full Dataset (includes eukaryotic DNA):"
    qiime tools peek "$OUTPUT_DIR/CFM_16S_final_filtered.qza" | grep -E "samples|features"
    size_16s_final=$(du -h "$OUTPUT_DIR/CFM_16S_final_filtered.qza" | cut -f1)
    echo "16S full dataset size: $size_16s_final"
fi

if [[ -f "$OUTPUT_DIR/CFM_16S_prokaryotes_only.qza" ]]; then
    echo ""
    echo "16S Prokaryotes-Only Dataset:"
    qiime tools peek "$OUTPUT_DIR/CFM_16S_prokaryotes_only.qza" | grep -E "samples|features"
    size_16s_prok=$(du -h "$OUTPUT_DIR/CFM_16S_prokaryotes_only.qza" | cut -f1)
    echo "16S prokaryotes-only size: $size_16s_prok"
fi

echo ""

# ITS1 summary
if [[ -f "$OUTPUT_DIR/CFM_ITS1_final_filtered.qza" ]]; then
    echo "ITS1 Filtering Results:"
    qiime tools peek "$OUTPUT_DIR/CFM_ITS1_final_filtered.qza" | grep -E "samples|features"
    size_its1_final=$(du -h "$OUTPUT_DIR/CFM_ITS1_final_filtered.qza" | cut -f1)
    echo "ITS1 final filtered table size: $size_its1_final"
fi

echo ""
echo "Output files generated:"
echo "16S DATASETS:"
echo "- CFM_16S_final_filtered.qza (includes eukaryotic DNA like Phytophthora)"
echo "- CFM_16S_filtered_repseqs.qza (representative sequences for full dataset)"
echo "- CFM_16S_final_filtered_summary.qzv (full dataset summary)"
echo "- CFM_16S_prokaryotes_only.qza (standard bacterial microbiome)"
echo "- CFM_16S_prokaryotes_only_repseqs.qza (prokaryotic sequences only)"
echo "- CFM_16S_prokaryotes_only_summary.qzv (prokaryotes-only summary)"
echo ""
echo "ITS1 DATASET:"
echo "- CFM_ITS1_final_filtered.qza (fungal communities only)"
echo "- CFM_ITS1_filtered_repseqs.qza (fungal representative sequences)"
echo "- CFM_ITS1_final_filtered_summary.qzv (fungal dataset summary)"

echo ""
echo "Filtering applied:"
echo "- 16S FULL: Removed plastids and mitochondria, kept eukaryotic DNA (Phytophthora)"
echo "- 16S PROKARYOTES: Additionally removed all Domain Eukaryota"
echo "- ITS1: Kept only Kingdom Fungi"
echo "- ALL: Removed singletons and low-prevalence ASVs (<5% samples)"

echo ""
echo "Analytical approach:"
echo "1. Use 16S FULL dataset for comprehensive analysis including pathogens"
echo "2. Use 16S PROKARYOTES dataset for standard bacterial microbiome comparison"
echo "3. Use ITS1 dataset for fungal community analysis"
echo "4. Compare results between 16S approaches in your publication"

echo ""
echo "Next steps:"
echo "1. Download and examine all filtering summary (.qzv) files"
echo "2. Compare feature counts between 16S full vs prokaryotes-only datasets"
echo "3. Proceed with alpha rarefaction analysis on BOTH 16S datasets"
echo "4. Document Phytophthora detection via 16S amplification in your methods"

echo ""
echo "Dataset filtering completed at: $(date)"
echo "Ready for dual-approach alpha rarefaction analysis!"
