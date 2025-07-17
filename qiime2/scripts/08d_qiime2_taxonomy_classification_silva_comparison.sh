#!/bin/bash
#SBATCH --job-name=taxonomy_silva_comparison
#SBATCH -p bigmem
#SBATCH --mem=500G
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --gres=tmp=100G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/08d_taxonomy_silva_comparison.out
#SBATCH --error=logs/08d_taxonomy_silva_comparison.err

# Script to compare 16S taxonomic classification using different SILVA classifiers
# Step 8d: Compare standard vs diverse-weighted vs plant-optimized SILVA classifiers
# Input: 16S DADA2 representative sequences (corrected data)
# Output: Taxonomic assignments and visualizations for each classifier

echo "Starting SILVA classifier comparison for 16S taxonomy..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 500G"
echo "Purpose: Compare different SILVA classifiers for optimal 16S classification"

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

# Define classifiers and their short names
declare -A CLASSIFIERS
CLASSIFIERS["uniform"]="$DB_DIR/silva-138-99-nb-classifier.qza"
CLASSIFIERS["diverse-weighted"]="$DB_DIR/silva-138-99-nb-diverse-weighted-classifier.qza"
CLASSIFIERS["plant-corpus"]="$DB_DIR/plant-corpus.qza"
CLASSIFIERS["plant-surface"]="$DB_DIR/plant-surface.qza"

# Verify all classifiers exist
echo ""
echo "=== VERIFYING SILVA CLASSIFIERS ==="
for name in "${!CLASSIFIERS[@]}"; do
    classifier_path="${CLASSIFIERS[$name]}"
    if [[ -f "$classifier_path" ]]; then
        echo "✓ $name classifier found: $(basename $classifier_path)"
    else
        echo "✗ $name classifier NOT found: $classifier_path"
        exit 1
    fi
done

echo ""
echo "=== RUNNING SILVA CLASSIFIER COMPARISON ==="

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
        --o-classification "$OUTPUT_DIR/CFM_16S_taxonomy_${classifier_name}.qza" \
        --verbose

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
echo "=== SILVA CLASSIFIER COMPARISON COMPLETE ==="

# Final summary
echo ""
echo "=== CLASSIFICATION SUMMARY ==="
echo "Files created in $OUTPUT_DIR:"
ls -lah "$OUTPUT_DIR"/CFM_16S_taxonomy_*

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
echo "- CFM_16S_taxonomy_uniform.qza (SILVA uniform classifier)"
echo "- CFM_16S_taxonomy_diverse-weighted.qza (SILVA diverse-weighted classifier)"
echo "- CFM_16S_taxonomy_plant-corpus.qza (Plant corpus classifier)"
echo "- CFM_16S_taxonomy_plant-surface.qza (Plant surface classifier)"
echo ""
echo "Taxonomy visualizations (.qzv files):"
echo "- CFM_16S_taxonomy_[classifier]_viz.qzv (taxonomy tables)"
echo "- CFM_16S_taxonomy_[classifier]_barplot_meta.qzv (barplots with metadata)"

echo ""
echo "=== COMPARISON GUIDELINES ==="
echo ""
echo "Classifier characteristics:"
echo "- uniform: General-purpose SILVA classifier"
echo "- diverse-weighted: Weighted for diverse environments"
echo "- plant-corpus: Optimized for plant-associated microbiomes"
echo "- plant-surface: Specialized for plant surface communities"
echo ""
echo "For cacao flower microbiome analysis, expect:"
echo "- plant-surface: Likely best for flower epiphytic communities"
echo "- plant-corpus: Good for plant-associated bacteria"
echo "- diverse: Good balance for mixed communities"
echo "- standard: Baseline reference"

echo ""
echo "Next steps:"
echo "1. Download and compare all barplot visualizations"
echo "2. Check classification rates and taxonomic resolution"
echo "3. Compare taxa detected by each classifier"
echo "4. Choose optimal classifier based on:"
echo "   - Classification success rate"
echo "   - Taxonomic resolution (species vs genus level)"
echo "   - Biological relevance for flower microbiomes"
echo "5. Use chosen classifier for downstream analysis"

echo ""
echo "SILVA classifier comparison completed at: $(date)"
echo "Ready to choose optimal classifier for cacao flower microbiome analysis!"
