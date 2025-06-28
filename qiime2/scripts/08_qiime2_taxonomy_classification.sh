#!/bin/bash
#SBATCH --job-name=taxonomy_classify
#SBATCH -p bigmem
#SBATCH --mem=500G
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH --gres=tmp:100G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/08_taxonomy_classification.out
#SBATCH --error=logs/08_taxonomy_classification.err

# Script to assign taxonomy to ASVs using SILVA (16S) and UNITE (ITS1) databases
# Step 8: Taxonomic classification of representative sequences
# Input: DADA2 representative sequences
# Output: Taxonomic assignments and visualizations

echo "Starting taxonomic classification..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 500G"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
INPUT_DIR="qiime2/denoise"
OUTPUT_DIR="qiime2/taxonomy"
DB_DIR="qiime2/databases"  # Updated to match your structure

echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Database directory: $DB_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Verify input files exist
if [[ ! -f "$INPUT_DIR/CFM_16S_dada2_repseqs.qza" ]]; then
    echo "ERROR: 16S representative sequences not found: $INPUT_DIR/CFM_16S_dada2_repseqs.qza"
    exit 1
fi

if [[ ! -f "$INPUT_DIR/CFM_ITS1_dada2_repseqs.qza" ]]; then
    echo "ERROR: ITS1 representative sequences not found: $INPUT_DIR/CFM_ITS1_dada2_repseqs.qza"
    exit 1
fi

echo "Input files verified successfully"

# Taxonomic classification for 16S using SILVA
echo ""
echo "=== TAXONOMIC CLASSIFICATION: 16S rRNA V4 ==="
echo "Database: SILVA"
echo "Starting 16S classification at: $(date)"

# Define SILVA classifier path (diverse-weighted)
SILVA_CLASSIFIER="$DB_DIR/SILVA/silva-138-99-nb-diverse-weighted-classifier.qza"

# Check if SILVA classifier exists
if [[ ! -f "$SILVA_CLASSIFIER" ]]; then
    echo "ERROR: SILVA diverse-weighted classifier not found at: $SILVA_CLASSIFIER"
    echo "Please verify file location and try again"
    exit 1
else
    echo "SILVA classifier found: $SILVA_CLASSIFIER"
fi

echo "Classifying 16S sequences..."
qiime feature-classifier classify-sklearn \
    --i-classifier "$SILVA_CLASSIFIER" \
    --i-reads "$INPUT_DIR/CFM_16S_dada2_repseqs.qza" \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --o-classification "$OUTPUT_DIR/CFM_16S_taxonomy.qza" \
    --verbose

# Check 16S classification success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: 16S taxonomic classification completed at: $(date)"
    
    # Create 16S taxonomy visualization
    echo "Creating 16S taxonomy visualization..."
    qiime metadata tabulate \
        --m-input-file "$OUTPUT_DIR/CFM_16S_taxonomy.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_taxonomy_viz.qzv"
    
    # Create taxonomy barplot
    echo "Creating 16S taxonomy barplot..."
    qiime taxa barplot \
        --i-table "$INPUT_DIR/CFM_16S_dada2_table.qza" \
        --i-taxonomy "$OUTPUT_DIR/CFM_16S_taxonomy.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_taxonomy_barplot.qzv"
    
    echo "16S taxonomy visualizations created successfully"
else
    echo "ERROR: 16S taxonomic classification failed"
    echo "Check classifier path and database availability"
    exit 1
fi

echo ""

# Taxonomic classification for ITS1 using UNITE
echo "=== TAXONOMIC CLASSIFICATION: ITS1 Fungal ==="
echo "Database: UNITE"
echo "Starting ITS1 classification at: $(date)"

# Define UNITE classifier path (custom trained)
UNITE_CLASSIFIER="$DB_DIR/UNITE/unite_eukaryotes_dynamic_classifier.qza"

# Check if UNITE classifier exists
if [[ ! -f "$UNITE_CLASSIFIER" ]]; then
    echo "ERROR: UNITE classifier not found at: $UNITE_CLASSIFIER"
    echo "Please verify the custom-trained classifier location"
    exit 1
else
    echo "UNITE classifier found: $UNITE_CLASSIFIER"
fi

echo "Classifying ITS1 sequences..."
qiime feature-classifier classify-sklearn \
    --i-classifier "$UNITE_CLASSIFIER" \
    --i-reads "$INPUT_DIR/CFM_ITS1_dada2_repseqs.qza" \
    --p-n-jobs $SLURM_CPUS_PER_TASK \
    --o-classification "$OUTPUT_DIR/CFM_ITS1_taxonomy.qza" \
    --verbose

# Check ITS1 classification success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: ITS1 taxonomic classification completed at: $(date)"
    
    # Create ITS1 taxonomy visualization
    echo "Creating ITS1 taxonomy visualization..."
    qiime metadata tabulate \
        --m-input-file "$OUTPUT_DIR/CFM_ITS1_taxonomy.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_ITS1_taxonomy_viz.qzv"
    
    # Create taxonomy barplot
    echo "Creating ITS1 taxonomy barplot..."
    qiime taxa barplot \
        --i-table "$INPUT_DIR/CFM_ITS1_dada2_table.qza" \
        --i-taxonomy "$OUTPUT_DIR/CFM_ITS1_taxonomy.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_ITS1_taxonomy_barplot.qzv"
    
    echo "ITS1 taxonomy visualizations created successfully"
else
    echo "ERROR: ITS1 taxonomic classification failed"
    echo "Check classifier path and database availability"
    exit 1
fi

echo ""
echo "=== TAXONOMIC CLASSIFICATION COMPLETE ==="

# Final summary
echo "Checking output files..."

# List all created files
echo "Files created in $OUTPUT_DIR:"
ls -lah "$OUTPUT_DIR"/CFM_*_taxonomy*

echo ""
echo "=== CLASSIFICATION SUMMARY ==="

# 16S taxonomy summary
if [[ -f "$OUTPUT_DIR/CFM_16S_taxonomy.qza" ]]; then
    echo "16S Taxonomy Results:"
    qiime tools peek "$OUTPUT_DIR/CFM_16S_taxonomy.qza" | grep -E "features"
    
    size_16s_tax=$(du -h "$OUTPUT_DIR/CFM_16S_taxonomy.qza" | cut -f1)
    echo "16S taxonomy file size: $size_16s_tax"
fi

echo ""

# ITS1 taxonomy summary
if [[ -f "$OUTPUT_DIR/CFM_ITS1_taxonomy.qza" ]]; then
    echo "ITS1 Taxonomy Results:"
    qiime tools peek "$OUTPUT_DIR/CFM_ITS1_taxonomy.qza" | grep -E "features"
    
    size_its1_tax=$(du -h "$OUTPUT_DIR/CFM_ITS1_taxonomy.qza" | cut -f1)
    echo "ITS1 taxonomy file size: $size_its1_tax"
fi

echo ""
echo "Output files generated:"
echo "- CFM_16S_taxonomy.qza (bacterial/archaeal taxonomic assignments)"
echo "- CFM_16S_taxonomy_viz.qzv (16S taxonomy table visualization)"  
echo "- CFM_16S_taxonomy_barplot.qzv (16S interactive barplot)"
echo "- CFM_ITS1_taxonomy.qza (fungal/eukaryotic taxonomic assignments)"
echo "- CFM_ITS1_taxonomy_viz.qzv (ITS1 taxonomy table visualization)"
echo "- CFM_ITS1_taxonomy_barplot.qzv (ITS1 interactive barplot)"

echo ""
echo "Classifiers used:"
echo "- SILVA: Diverse-weighted classifier optimized for plant microbiomes"
echo "- UNITE: Custom-trained eukaryotes classifier with pathogen detection"

echo ""
echo "Next steps:"
echo "1. Download and examine taxonomy visualizations (.qzv files)"
echo "2. Check taxonomic assignments and coverage"
echo "3. Look for Phytophthora and other cacao pathogens in ITS1 data!"
echo "4. Filter taxa if needed (remove chloroplasts, mitochondria, etc.)"
echo "5. Proceed with alpha rarefaction and diversity analysis"

echo ""
echo "Taxonomic classification completed at: $(date)"
echo "Ready for alpha rarefaction and diversity analysis!"
