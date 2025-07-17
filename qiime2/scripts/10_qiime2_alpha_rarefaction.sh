#!/bin/bash
#SBATCH --job-name=alpha_rarefaction
#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --gres=tmp:30G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/10_alpha_rarefaction.out
#SBATCH --error=logs/10_alpha_rarefaction.err

# Script to perform alpha rarefaction analysis for determining optimal sampling depth
# Step 10: Alpha rarefaction curves for unfiltered and fully filtered datasets
# Input: Unfiltered and fully filtered feature tables
# Output: 5 alpha rarefaction curves optimized for actual data distributions

echo "Starting alpha rarefaction analysis..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 64G"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
INPUT_DIR_DENOISE="qiime2/denoise"
INPUT_DIR_FILTERED="qiime2/filtered"
OUTPUT_DIR="qiime2/rarefaction"
METADATA_FILE="data/cfm_qiime2_metadata.txt"

echo "Input denoise directory: $INPUT_DIR_DENOISE"
echo "Input filtered directory: $INPUT_DIR_FILTERED"
echo "Output directory: $OUTPUT_DIR"
echo "Metadata file: $METADATA_FILE"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if metadata file exists
if [[ -f "$METADATA_FILE" ]]; then
    echo "Metadata file found, will include in rarefaction analysis"
    METADATA_OPTION="--m-metadata-file $METADATA_FILE"
else
    echo "WARNING: Metadata file not found, proceeding without metadata"
    METADATA_OPTION=""
fi

# Verify input files exist
if [[ ! -f "$INPUT_DIR_DENOISE/CFM_16S_dada2_table.qza" ]]; then
    echo "ERROR: 16S unfiltered table not found"
    exit 1
fi

if [[ ! -f "$INPUT_DIR_DENOISE/CFM_ITS1_dada2_table.qza" ]]; then
    echo "ERROR: ITS1 unfiltered table not found"
    exit 1
fi

echo "Input files verified successfully"

# ============================================================================
# UNFILTERED DATA RAREFACTION (for comparison)
# ============================================================================

echo ""
echo "=== UNFILTERED DATA RAREFACTION ANALYSIS ==="
echo "Purpose: Quality assessment and comparison with filtered data"
echo "Starting unfiltered rarefaction analysis at: $(date)"

# 16S Unfiltered Alpha Rarefaction
echo ""
echo "--- 16S UNFILTERED RAREFACTION ---"
echo "Data summary: Min=71, Median=84,909, Max=244,992 reads"
echo "Parameters: Max depth: 200,000 reads, Steps: 75, Iterations: 10"

qiime diversity alpha-rarefaction \
    --i-table "$INPUT_DIR_DENOISE/CFM_16S_dada2_table.qza" \
    --p-max-depth 200000 \
    --p-steps 75 \
    --p-iterations 10 \
    $METADATA_OPTION \
    --o-visualization "$OUTPUT_DIR/CFM_16S_unfiltered_alpha_rarefaction.qzv" \
    --verbose

# Check 16S unfiltered rarefaction success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: 16S unfiltered rarefaction completed"
else
    echo "ERROR: 16S unfiltered rarefaction failed"
    exit 1
fi

# ITS1 Unfiltered Alpha Rarefaction
echo ""
echo "--- ITS1 UNFILTERED RAREFACTION ---"
echo "Data summary: Min=0, Median=50,423, Max=148,702 reads"
echo "Parameters: Max depth: 130,000 reads, Steps: 65, Iterations: 10"

qiime diversity alpha-rarefaction \
    --i-table "$INPUT_DIR_DENOISE/CFM_ITS1_dada2_table.qza" \
    --p-max-depth 130000 \
    --p-steps 65 \
    --p-iterations 10 \
    $METADATA_OPTION \
    --o-visualization "$OUTPUT_DIR/CFM_ITS1_unfiltered_alpha_rarefaction.qzv" \
    --verbose

# Check ITS1 unfiltered rarefaction success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: ITS1 unfiltered rarefaction completed"
else
    echo "ERROR: ITS1 unfiltered rarefaction failed"
    exit 1
fi

echo "Unfiltered rarefaction analysis completed at: $(date)"

# ============================================================================
# FILTERED DATA RAREFACTION ANALYSIS
# ============================================================================

echo ""
echo "=== FILTERED DATA RAREFACTION ANALYSIS ==="
echo "Purpose: Determine optimal sampling depth for diversity analysis"
echo "Starting filtered rarefaction analysis at: $(date)"

# 16S Fully Filtered (chloroplast/mitochondria + singleton/low-prevalence)
if [[ -f "$INPUT_DIR_FILTERED/CFM_16S_filtered_chloro_mito_singl_lowprev3.qza" ]]; then
    echo ""
    echo "--- 16S FULLY FILTERED RAREFACTION ---"
    echo "Data summary: Min=31, Median=17,424, Max=173,188 reads"
    echo "Parameters: Max depth: 150,000 reads, Steps: 75, Iterations: 10"
    
    qiime diversity alpha-rarefaction \
        --i-table "$INPUT_DIR_FILTERED/CFM_16S_filtered_chloro_mito_singl_lowprev3.qza" \
        --p-max-depth 150000 \
        --p-steps 75 \
        --p-iterations 10 \
        $METADATA_OPTION \
        --o-visualization "$OUTPUT_DIR/CFM_16S_filtered_chloro_mito_singl_lowprev3_alpha_rarefaction.qzv" \
        --verbose

    if [[ $? -eq 0 ]]; then
        echo "SUCCESS: 16S fully filtered rarefaction completed"
    else
        echo "ERROR: 16S fully filtered rarefaction failed"
    fi
else
    echo "WARNING: CFM_16S_filtered_chloro_mito_singl_lowprev3.qza not found, skipping rarefaction"
fi

# 16S Prokaryotes-Only Filtered
if [[ -f "$INPUT_DIR_FILTERED/CFM_16S_filtered_prokaryotes_only_singl_lowprev3.qza" ]]; then
    echo ""
    echo "--- 16S PROKARYOTES-ONLY FILTERED RAREFACTION ---"
    echo "Data summary: Min=31, Median=16,678, Max=173,188 reads"
    echo "Parameters: Max depth: 150,000 reads, Steps: 75, Iterations: 10"
    
    qiime diversity alpha-rarefaction \
        --i-table "$INPUT_DIR_FILTERED/CFM_16S_filtered_prokaryotes_only_singl_lowprev3.qza" \
        --p-max-depth 150000 \
        --p-steps 75 \
        --p-iterations 10 \
        $METADATA_OPTION \
        --o-visualization "$OUTPUT_DIR/CFM_16S_filtered_prokaryotes_only_singl_lowprev3_alpha_rarefaction.qzv" \
        --verbose

    if [[ $? -eq 0 ]]; then
        echo "SUCCESS: 16S prokaryotes-only rarefaction completed"
    else
        echo "ERROR: 16S prokaryotes-only rarefaction failed"
    fi
else
    echo "WARNING: CFM_16S_filtered_prokaryotes_only_singl_lowprev3.qza not found, skipping rarefaction"
fi

# ITS1 Fully Filtered (fungi-only + singleton/low-prevalence)
if [[ -f "$INPUT_DIR_FILTERED/CFM_ITS1_filtered_fungi_only_singl_lowprev3.qza" ]]; then
    echo ""
    echo "--- ITS1 FULLY FILTERED RAREFACTION ---"
    echo "Data summary: Min=16, Median=48,052, Max=146,897 reads"
    echo "Parameters: Max depth: 130,000 reads, Steps: 65, Iterations: 10"
    
    qiime diversity alpha-rarefaction \
        --i-table "$INPUT_DIR_FILTERED/CFM_ITS1_filtered_fungi_only_singl_lowprev3.qza" \
        --p-max-depth 130000 \
        --p-steps 65 \
        --p-iterations 10 \
        $METADATA_OPTION \
        --o-visualization "$OUTPUT_DIR/CFM_ITS1_filtered_fungi_only_singl_lowprev3_alpha_rarefaction.qzv" \
        --verbose

    if [[ $? -eq 0 ]]; then
        echo "SUCCESS: ITS1 fully filtered rarefaction completed"
    else
        echo "ERROR: ITS1 fully filtered rarefaction failed"
    fi
else
    echo "WARNING: CFM_ITS1_filtered_fungi_only_singl_lowprev3.qza not found, skipping rarefaction"
fi

echo ""
echo "Filtered rarefaction analysis completed at: $(date)"

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "=== ALPHA RAREFACTION ANALYSIS COMPLETE ==="

echo ""
echo "Checking output files..."
ls -lah "$OUTPUT_DIR"/*.qzv 2>/dev/null || echo "No .qzv files found"

echo ""
echo "Output files generated (5 rarefaction curves):"
echo ""
echo "UNFILTERED DATA (for comparison):"
echo "- CFM_16S_unfiltered_alpha_rarefaction.qzv"
echo "- CFM_ITS1_unfiltered_alpha_rarefaction.qzv"
echo ""
echo "FILTERED DATA (for analysis):"
echo "- CFM_16S_filtered_chloro_mito_singl_lowprev3_alpha_rarefaction.qzv"
echo "- CFM_16S_filtered_prokaryotes_only_singl_lowprev3_alpha_rarefaction.qzv"
echo "- CFM_ITS1_filtered_fungi_only_singl_lowprev3_alpha_rarefaction.qzv"

echo ""
echo "Next steps:"
echo "1. Download and view ALL .qzv files at https://view.qiime2.org"
echo "2. Compare unfiltered vs filtered rarefaction curves"
echo "3. Compare 16S full vs prokaryotes-only datasets"
echo "4. Choose optimal rarefaction depths for downstream diversity analysis"
echo "5. Document your chosen depths and rationale for methods section"

echo ""
echo "Alpha rarefaction analysis completed at: $(date)"
echo "Ready for phylogenetic tree construction with optimal depths!"
