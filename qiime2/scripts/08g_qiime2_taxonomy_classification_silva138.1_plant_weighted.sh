#!/bin/bash
#SBATCH --job-name=taxonomy_plant_weighted
#SBATCH -p bigmem
#SBATCH --mem=500G
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --gres=tmp=100G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/08g_taxonomy_plant_weighted.out
#SBATCH --error=logs/08g_taxonomy_plant_weighted.err

# Script to classify 16S taxonomies using plant-weighted SILVA 138.1 classifiers
# Step 8g: Taxonomic classification using plant-corpus and plant-surface weighted classifiers
# Input: 16S DADA2 representative sequences (corrected data)
# Output: Taxonomic assignments and visualizations for plant-optimized classifiers

echo "Starting SILVA 138.1 plant-weighted classifier comparison for 16S taxonomy..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 500G"
echo "Purpose: Compare plant-corpus and plant-surface weighted SILVA classifiers"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
INPUT_DIR="qiime2/denoise"
OUTPUT_DIR="qiime2/taxonomy"
DB_DIR="qiime2/databases/SILVA"
METADATA_FILE="data/cfm_qiime2_metadata.txt"

echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Database directory: $DB_DIR"
echo "Metadata file: $METADATA_FILE"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Verify input files exist
if [[ ! -f "$INPUT_DIR/CFM_16S_dada2_repseqs.qza" ]]; then
    echo "ERROR: 16S representative sequences not found: $INPUT_DIR/CFM_16S_dada2_repseqs.qza"
    exit 1
fi

if [[ ! -f "$INPUT_DIR/CFM_16S_dada2_table.qza" ]]; then
    echo "ERROR: 16S feature table not found: $INPUT_DIR/CFM_16S_dada2_table.qza"
    exit 1
fi

# Check if metadata exists
if [[ -f "$METADATA_FILE" ]]; then
    echo "Metadata file found, will include in visualizations"
    METADATA_OPTION="--m-metadata-file $METADATA_FILE"
else
    echo "WARNING: Metadata file not found, proceeding without metadata"
    METADATA_OPTION=""
fi

echo "Input files verified successfully"

# Define plant-weighted classifiers and their short names
declare -A CLASSIFIERS
CLASSIFIERS["plant-corpus"]="$DB_DIR/silva138.1_plant-corpus_weighted_classifier.qza"
CLASSIFIERS["plant-surface"]="$DB_DIR/silva138.1_plant-surface_weighted_classifier.qza"

# Verify all classifiers exist
echo ""
echo "=== VERIFYING PLANT-WEIGHTED CLASSIFIERS ==="
for name in "${!CLASSIFIERS[@]}"; do
    classifier_path="${CLASSIFIERS[$name]}"
    if [[ -f "$classifier_path" ]]; then
        size=$(du -h "$classifier_path" | cut -f1)
        echo "✓ $name classifier found: $(basename $classifier_path) ($size)"
    else
        echo "✗ $name classifier NOT found: $classifier_path"
        echo "Please run 08e and 08f scripts to train plant-weighted classifiers"
        exit 1
    fi
done

echo ""
echo "=== RUNNING PLANT-WEIGHTED CLASSIFIER COMPARISON ==="

# Process each classifier
for classifier_name in "${!CLASSIFIERS[@]}"; do
    classifier_path="${CLASSIFIERS[$classifier_name]}"
    
    echo ""
    echo "--- PROCESSING: $classifier_name CLASSIFIER ---"
    echo "Classifier: $(basename $classifier_path)"
    echo "Starting classification at: $(date)"
    
    # Run taxonomic classification
    echo "Classifying 16S sequences with $classifier_name classifier..."
    qiime feature-classifier classify-sklearn \
        --i-classifier "$classifier_path" \
        --i-reads "$INPUT_DIR/CFM_16S_dada2_repseqs.qza" \
        --p-n-jobs $SLURM_CPUS_PER_TASK \
        --o-classification "$OUTPUT_DIR/CFM_16S_taxonomy_${classifier_name}.qza"

    # Check classification success
    if [[ $? -eq 0 ]]; then
        echo "SUCCESS: $classifier_name classification completed at: $(date)"

        # Create taxonomy visualization
        echo "Creating $classifier_name taxonomy visualization..."
        qiime metadata tabulate \
            --m-input-file "$OUTPUT_DIR/CFM_16S_taxonomy_${classifier_name}.qza" \
            --o-visualization "$OUTPUT_DIR/CFM_16S_taxonomy_${classifier_name}_viz.qzv"

        # Create taxonomy barplot with metadata
        if [[ -f "$METADATA_FILE" ]]; then
            echo "Creating $classifier_name taxonomy barplot with metadata..."
            qiime taxa barplot \
                --i-table "$INPUT_DIR/CFM_16S_dada2_table.qza" \
                --i-taxonomy "$OUTPUT_DIR/CFM_16S_taxonomy_${classifier_name}.qza" \
                $METADATA_OPTION \
                --o-visualization "$OUTPUT_DIR/CFM_16S_taxonomy_${classifier_name}_barplot_meta.qzv"
        else
            echo "WARNING: No metadata file found, skipping barplot creation for $classifier_name"
        fi

        echo "$classifier_name visualizations created successfully"
        
        # Quick stats
        size_taxonomy=$(du -h "$OUTPUT_DIR/CFM_16S_taxonomy_${classifier_name}.qza" | cut -f1)
        echo "$classifier_name taxonomy file size: $size_taxonomy"
        
    else
        echo "ERROR: $classifier_name classification failed"
        echo "Continuing with remaining classifiers..."
    fi
done

echo ""
echo "=== PLANT-WEIGHTED CLASSIFIER COMPARISON COMPLETE ==="

# Final summary
echo ""
echo "=== CLASSIFICATION SUMMARY ==="
echo "Files created in $OUTPUT_DIR:"
ls -lah "$OUTPUT_DIR"/CFM_16S_taxonomy_plant-*

echo ""
echo "=== CLASSIFIER COMPARISON RESULTS ==="

# Summary for each classifier
for classifier_name in "${!CLASSIFIERS[@]}"; do
    taxonomy_file="$OUTPUT_DIR/CFM_16S_taxonomy_${classifier_name}.qza"
    if [[ -f "$taxonomy_file" ]]; then
        echo ""
        echo "$classifier_name classifier results:"
        qiime tools peek "$taxonomy_file" | grep -E "features"
        
        # Get file size
        size=$(du -h "$taxonomy_file" | cut -f1)
        echo "$classifier_name taxonomy size: $size"
    else
        echo "$classifier_name: FAILED"
    fi
done

echo ""
echo "=== OUTPUT FILES GENERATED ==="
echo ""
echo "Taxonomy assignments (.qza files):"
echo "- CFM_16S_taxonomy_plant-corpus.qza (SILVA 138.1 plant-corpus weighted)"
echo "- CFM_16S_taxonomy_plant-surface.qza (SILVA 138.1 plant-surface weighted)"
echo ""
echo "Taxonomy visualizations (.qzv files):"
echo "- CFM_16S_taxonomy_[classifier]_viz.qzv (taxonomy tables)"
echo "- CFM_16S_taxonomy_[classifier]_barplot_meta.qzv (barplots with metadata)"

echo ""
echo "=== COMPLETE SILVA TAXONOMY COMPARISON ==="
echo ""
echo "You now have 4 SILVA taxonomies for comprehensive comparison:"
echo "1. CFM_16S_taxonomy_uniform.qza (uniform SILVA classifier)"
echo "2. CFM_16S_taxonomy_diverse-weighted.qza (diverse-weighted SILVA classifier)"
echo "3. CFM_16S_taxonomy_plant-corpus.qza (plant-corpus weighted classifier)"
echo "4. CFM_16S_taxonomy_plant-surface.qza (plant-surface weighted classifier)"

echo ""
echo "Plant-weighted classifier comparison completed at: $(date)"
echo "Ready for comprehensive SILVA taxonomy comparison!"
