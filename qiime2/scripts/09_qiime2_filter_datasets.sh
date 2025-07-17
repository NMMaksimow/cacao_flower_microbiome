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

# Calculate minimum samples threshold (1% of 336 samples = ~3 samples)
MIN_SAMPLES=3

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
    --o-filtered-table "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito.qza" \
    --verbose

# Check 16S table filtering success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: 16S taxonomic filtering completed"

    # Filter representative sequences to match taxonomically filtered table
    echo "Filtering 16S representative sequences..."
    qiime taxa filter-seqs \
        --i-sequences "$INPUT_DIR_DENOISE/CFM_16S_dada2_repseqs.qza" \
        --i-taxonomy "$INPUT_DIR_TAXONOMY/CFM_16S_taxonomy.qza" \
        --p-exclude mitochondria,chloroplast \
        --o-filtered-sequences "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito_repseqs.qza" \
        --verbose

    echo "16S taxonomic filtering complete. Plastids and mitochondria removed."
else
    echo "ERROR: 16S table filtering failed"
    exit 1
fi

# Combined singleton and prevalence filtering for 16S
echo "Removing singletons (frequency < 2) and low-prevalence ASVs (present in < $MIN_SAMPLES samples)..."
qiime feature-table filter-features \
    --i-table "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito.qza" \
    --p-min-frequency 2 \
    --p-min-samples $MIN_SAMPLES \
    --o-filtered-table "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito_singl_lowprev3.qza" \
    --verbose

# Check 16S combined filtering success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: 16S combined filtering completed"

    # Filter representative sequences to match filtered table
    echo "Filtering 16S representative sequences to match filtered table..."
    qiime feature-table filter-seqs \
        --i-data "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito_repseqs.qza" \
        --i-table "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito_singl_lowprev3.qza" \
        --o-filtered-data "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito_singl_lowprev3_repseqs.qza" \
        --verbose

    # Generate summary of final 16S table
    echo "Creating 16S filtered table summary..."
    qiime feature-table summarize \
        --i-table "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito_singl_lowprev3.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito_singl_lowprev3_summary.qzv"

    echo "16S filtering pipeline complete."
else
    echo "ERROR: 16S combined filtering failed"
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
    --o-filtered-table "$OUTPUT_DIR/CFM_ITS1_filtered_fungi_only.qza" \
    --verbose

# Check ITS1 table filtering success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: ITS1 taxonomic filtering completed"

    # Filter representative sequences to match taxonomically filtered table
    echo "Filtering ITS1 representative sequences..."
    qiime taxa filter-seqs \
        --i-sequences "$INPUT_DIR_DENOISE/CFM_ITS1_dada2_repseqs.qza" \
        --i-taxonomy "$INPUT_DIR_TAXONOMY/CFM_ITS1_taxonomy.qza" \
        --p-include k__Fungi \
        --o-filtered-sequences "$OUTPUT_DIR/CFM_ITS1_filtered_fungi_only_repseqs.qza" \
        --verbose

    echo "ITS1 taxonomic filtering complete. Only Fungi retained."
else
    echo "ERROR: ITS1 table filtering failed"
    exit 1
fi

# Combined singleton and prevalence filtering for ITS1
echo "Removing singletons (frequency < 2) and low-prevalence ASVs (present in < $MIN_SAMPLES samples)..."
qiime feature-table filter-features \
    --i-table "$OUTPUT_DIR/CFM_ITS1_filtered_fungi_only.qza" \
    --p-min-frequency 2 \
    --p-min-samples $MIN_SAMPLES \
    --o-filtered-table "$OUTPUT_DIR/CFM_ITS1_filtered_fungi_only_singl_lowprev3.qza" \
    --verbose

# Check ITS1 combined filtering success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: ITS1 combined filtering completed"

    # Filter representative sequences to match filtered table
    echo "Filtering ITS1 representative sequences to match filtered table..."
    qiime feature-table filter-seqs \
        --i-data "$OUTPUT_DIR/CFM_ITS1_filtered_fungi_only_repseqs.qza" \
        --i-table "$OUTPUT_DIR/CFM_ITS1_filtered_fungi_only_singl_lowprev3.qza" \
        --o-filtered-data "$OUTPUT_DIR/CFM_ITS1_filtered_fungi_only_singl_lowprev3_repseqs.qza" \
        --verbose

    # Generate summary of final ITS1 table
    echo "Creating ITS1 filtered table summary..."
    qiime feature-table summarize \
        --i-table "$OUTPUT_DIR/CFM_ITS1_filtered_fungi_only_singl_lowprev3.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_ITS1_filtered_fungi_only_singl_lowprev3_summary.qzv"

    echo "ITS1 filtering pipeline complete."
else
    echo "ERROR: ITS1 combined filtering failed"
    exit 1
fi

echo ""

# ============================================================================
# OPTIONAL: PROKARYOTES-ONLY 16S ANALYSIS (for comparison)
# ============================================================================

echo "=== CREATING PROKARYOTES-ONLY 16S DATASET ==="
echo "Purpose: Standard bacterial microbiome analysis excluding eukaryotic DNA"
echo "Note: This will remove Phytophthora sequences found in the 16S data"

# Remove all Eukaryota from 16S data for standard bacterial analysis
echo "Removing Domain Eukaryota from 16S dataset..."
qiime taxa filter-table \
    --i-table "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito_singl_lowprev3.qza" \
    --i-taxonomy "$INPUT_DIR_TAXONOMY/CFM_16S_taxonomy.qza" \
    --p-exclude "d__Eukaryota" \
    --o-filtered-table "$OUTPUT_DIR/CFM_16S_filtered_prokaryotes_only.qza" \
    --verbose

# Apply singleton and prevalence filtering to prokaryotes-only dataset
echo "Applying singleton and prevalence filtering to prokaryotes-only dataset..."
qiime feature-table filter-features \
    --i-table "$OUTPUT_DIR/CFM_16S_filtered_prokaryotes_only.qza" \
    --p-min-frequency 2 \
    --p-min-samples $MIN_SAMPLES \
    --o-filtered-table "$OUTPUT_DIR/CFM_16S_filtered_prokaryotes_only_singl_lowprev3.qza" \
    --verbose

# Check prokaryotes-only filtering success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Prokaryotes-only table created and filtered"

    # Filter representative sequences by taxonomy first (remove eukaryotes)
    echo "Filtering representative sequences to remove eukaryotes..."
    qiime taxa filter-seqs \
        --i-sequences "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito_singl_lowprev3_repseqs.qza" \
        --i-taxonomy "$INPUT_DIR_TAXONOMY/CFM_16S_taxonomy.qza" \
        --p-exclude "d__Eukaryota" \
        --o-filtered-sequences "$OUTPUT_DIR/CFM_16S_filtered_prokaryotes_only_repseqs_temp.qza" \
        --verbose

    # Then filter to match the final prokaryotes-only table
    echo "Filtering representative sequences to match prokaryotes-only table..."
    qiime feature-table filter-seqs \
        --i-data "$OUTPUT_DIR/CFM_16S_filtered_prokaryotes_only_repseqs_temp.qza" \
        --i-table "$OUTPUT_DIR/CFM_16S_filtered_prokaryotes_only_singl_lowprev3.qza" \
        --o-filtered-data "$OUTPUT_DIR/CFM_16S_filtered_prokaryotes_only_singl_lowprev3_repseqs.qza" \
        --verbose

    # Generate summary for prokaryotes-only dataset
    echo "Creating prokaryotes-only summary..."
    qiime feature-table summarize \
        --i-table "$OUTPUT_DIR/CFM_16S_filtered_prokaryotes_only_singl_lowprev3.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_filtered_prokaryotes_only_singl_lowprev3_summary.qzv"

    # Clean up temporary file
    rm -f "$OUTPUT_DIR/CFM_16S_filtered_prokaryotes_only_repseqs_temp.qza"

    echo "Prokaryotes-only dataset created successfully"
else
    echo "ERROR: Prokaryotes-only filtering failed"
fi

echo ""
echo "=== DATASET FILTERING COMPLETE ==="

echo ""
echo "Key files for downstream analysis:"
echo "- CFM_16S_filtered_chloro_mito_singl_lowprev3.qza (16S feature table - all domains)"
echo "- CFM_16S_filtered_chloro_mito_singl_lowprev3_repseqs.qza (16S sequences - all domains)"
echo "- CFM_ITS1_filtered_fungi_only_singl_lowprev3.qza (ITS1 feature table - fungi only)"
echo "- CFM_ITS1_filtered_fungi_only_singl_lowprev3_repseqs.qza (ITS1 sequences - fungi only)"
echo "- CFM_16S_filtered_prokaryotes_only_singl_lowprev3.qza (16S feature table - prokaryotes only)"
echo "- CFM_16S_filtered_prokaryotes_only_singl_lowprev3_repseqs.qza (16S sequences - prokaryotes only)"

echo ""
echo "Dataset filtering completed at: $(date)"
echo "Ready for alpha rarefaction analysis!"
