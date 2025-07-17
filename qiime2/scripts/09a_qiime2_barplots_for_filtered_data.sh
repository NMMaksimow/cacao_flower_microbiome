#!/bin/bash
#SBATCH --job-name=filtered_barplots
#SBATCH -p shared
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --gres=tmp:10G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/09a_barplots_for_filtered_data.out
#SBATCH --error=logs/09a_barplots_for_filtered_data.err

# Script to create taxonomy barplots visualizations for filtered datasets
# Step 9a: Barplot visualizations for all filtered feature tables
# Input: Filtered feature tables from qiime2/filtered/
# Output: Taxonomy barplots

echo "Starting filtered data visualization creation..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 16G"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
FILTERED_DIR="qiime2/filtered"
TAXONOMY_DIR="qiime2/taxonomy"
METADATA_FILE="data/cfm_qiime2_metadata.txt"

echo "Filtered data directory: $FILTERED_DIR"
echo "Taxonomy directory: $TAXONOMY_DIR"
echo "Metadata file: $METADATA_FILE"

# Verify required files exist
if [[ ! -f "$TAXONOMY_DIR/CFM_16S_taxonomy.qza" ]]; then
    echo "ERROR: 16S taxonomy not found: $TAXONOMY_DIR/CFM_16S_taxonomy.qza"
    exit 1
fi

if [[ ! -f "$TAXONOMY_DIR/CFM_ITS1_taxonomy.qza" ]]; then
    echo "ERROR: ITS1 taxonomy not found: $TAXONOMY_DIR/CFM_ITS1_taxonomy.qza"
    exit 1
fi

if [[ ! -f "$METADATA_FILE" ]]; then
    echo "WARNING: Metadata file not found: $METADATA_FILE"
    echo "Proceeding without metadata..."
    METADATA_OPTION=""
else
    echo "Metadata file found, will include in visualizations"
    METADATA_OPTION="--m-metadata-file $METADATA_FILE"
fi

echo "Required files verified"

# ============================================================================
# 16S BARPLOTS - ALL FILTERED DATASETS
# ============================================================================

echo ""
echo "=== CREATING 16S TAXONOMY BARPLOTS ==="

# 1. 16S Chloroplast/Mitochondria filtered (taxonomic filtering only)
if [[ -f "$FILTERED_DIR/CFM_16S_filtered_chloro_mito.qza" ]]; then
    echo "Creating 16S chloroplast/mitochondria filtered barplot..."
    qiime taxa barplot \
        --i-table "$FILTERED_DIR/CFM_16S_filtered_chloro_mito.qza" \
        --i-taxonomy "$TAXONOMY_DIR/CFM_16S_taxonomy.qza" \
        $METADATA_OPTION \
        --o-visualization "$FILTERED_DIR/CFM_16S_filtered_chloro_mito_barplot.qzv" \
        --verbose
    echo "SUCCESS: 16S chloroplast/mitochondria barplot created"
else
    echo "WARNING: CFM_16S_filtered_chloro_mito.qza not found, skipping barplot"
fi

# 2. 16S Chloroplast/Mitochondria + Singleton/Low-prevalence filtered
if [[ -f "$FILTERED_DIR/CFM_16S_filtered_chloro_mito_singl_lowprev3.qza" ]]; then
    echo "Creating 16S fully filtered barplot..."
    qiime taxa barplot \
        --i-table "$FILTERED_DIR/CFM_16S_filtered_chloro_mito_singl_lowprev3.qza" \
        --i-taxonomy "$TAXONOMY_DIR/CFM_16S_taxonomy.qza" \
        $METADATA_OPTION \
        --o-visualization "$FILTERED_DIR/CFM_16S_filtered_chloro_mito_singl_lowprev3_barplot.qzv" \
        --verbose
    echo "SUCCESS: 16S fully filtered barplot created"
else
    echo "WARNING: CFM_16S_filtered_chloro_mito_singl_lowprev3.qza not found, skipping barplot"
fi

# 3. 16S Prokaryotes-only + Singleton/Low-prevalence filtered
if [[ -f "$FILTERED_DIR/CFM_16S_filtered_prokaryotes_only_singl_lowprev3.qza" ]]; then
    echo "Creating 16S prokaryotes-only barplot..."
    qiime taxa barplot \
        --i-table "$FILTERED_DIR/CFM_16S_filtered_prokaryotes_only_singl_lowprev3.qza" \
        --i-taxonomy "$TAXONOMY_DIR/CFM_16S_taxonomy.qza" \
        $METADATA_OPTION \
        --o-visualization "$FILTERED_DIR/CFM_16S_filtered_prokaryotes_only_singl_lowprev3_barplot.qzv" \
        --verbose
    echo "SUCCESS: 16S prokaryotes-only barplot created"
else
    echo "WARNING: CFM_16S_filtered_prokaryotes_only_singl_lowprev3.qza not found, skipping barplot"
fi

# ============================================================================
# ITS1 BARPLOTS - ALL FILTERED DATASETS
# ============================================================================

echo ""
echo "=== CREATING ITS1 TAXONOMY BARPLOTS ==="

# 1. ITS1 Fungi-only filtered (taxonomic filtering only)
if [[ -f "$FILTERED_DIR/CFM_ITS1_filtered_fungi_only.qza" ]]; then
    echo "Creating ITS1 fungi-only filtered barplot..."
    qiime taxa barplot \
        --i-table "$FILTERED_DIR/CFM_ITS1_filtered_fungi_only.qza" \
        --i-taxonomy "$TAXONOMY_DIR/CFM_ITS1_taxonomy.qza" \
        $METADATA_OPTION \
        --o-visualization "$FILTERED_DIR/CFM_ITS1_filtered_fungi_only_barplot.qzv" \
        --verbose
    echo "SUCCESS: ITS1 fungi-only barplot created"
else
    echo "WARNING: CFM_ITS1_filtered_fungi_only.qza not found, skipping barplot"
fi

# 2. ITS1 Fungi-only + Singleton/Low-prevalence filtered
if [[ -f "$FILTERED_DIR/CFM_ITS1_filtered_fungi_only_singl_lowprev3.qza" ]]; then
    echo "Creating ITS1 fully filtered barplot..."
    qiime taxa barplot \
        --i-table "$FILTERED_DIR/CFM_ITS1_filtered_fungi_only_singl_lowprev3.qza" \
        --i-taxonomy "$TAXONOMY_DIR/CFM_ITS1_taxonomy.qza" \
        $METADATA_OPTION \
        --o-visualization "$FILTERED_DIR/CFM_ITS1_filtered_fungi_only_singl_lowprev3_barplot.qzv" \
        --verbose
    echo "SUCCESS: ITS1 fully filtered barplot created"
else
    echo "WARNING: CFM_ITS1_filtered_fungi_only_singl_lowprev3.qza not found, skipping barplot"
fi

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "=== VISUALIZATION CREATION COMPLETE ==="
echo ""
echo "Created barplot visualizations:"
echo "- CFM_16S_filtered_chloro_mito_barplot.qzv"
echo "- CFM_16S_filtered_chloro_mito_singl_lowprev3_barplot.qzv"
echo "- CFM_16S_filtered_prokaryotes_only_singl_lowprev3_barplot.qzv"
echo "- CFM_ITS1_filtered_fungi_only_barplot.qzv"
echo "- CFM_ITS1_filtered_fungi_only_singl_lowprev3_barplot.qzv"
echo ""
echo "All visualizations saved in: $FILTERED_DIR"
echo ""
echo "Visualization creation completed at: $(date)"
