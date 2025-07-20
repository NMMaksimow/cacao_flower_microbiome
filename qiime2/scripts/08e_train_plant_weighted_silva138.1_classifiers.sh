#!/bin/bash
#SBATCH --job-name=train_plant_classifiers
#SBATCH -p bigmem
#SBATCH --mem=240G
#SBATCH -c 48
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --gres=tmp:100G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/08e_train_plant_classifiers.out
#SBATCH --error=logs/08e_train_plant_classifiers.err

# Script to train plant-weighted SILVA 138.1 taxonomic classifiers
# Step 8e: Train custom plant-corpus and plant-surface weighted classifiers
# Input: SILVA 138.1 reference sequences, taxonomy, and plant-specific weights
# Output: Two custom-trained plant-optimized classifiers

echo "Starting plant-weighted SILVA 138.1 classifier training..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 240G"
echo "Purpose: Train plant-corpus and plant-surface weighted classifiers"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
OUTPUT_DIR="qiime2/databases/SILVA"

echo "Output directory: $OUTPUT_DIR"

# Create directory
mkdir -p "$OUTPUT_DIR"

echo "Downloading SILVA 138.1 reference data and plant weights to SILVA directory..."

# Download reference sequences (large file ~4GB)
echo ""
echo "=== DOWNLOADING SILVA 138.1 REFERENCE SEQUENCES ==="
echo "Downloading reference sequences from Zenodo..."
wget https://zenodo.org/record/6395539/files/ref-seqs.qza \
    -O "$OUTPUT_DIR/silva138.1_ref-seqs.qza"

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Reference sequences downloaded"
    ls -lah "$OUTPUT_DIR/silva138.1_ref-seqs.qza"
else
    echo "ERROR: Failed to download reference sequences"
    exit 1
fi

# Download reference taxonomy
echo ""
echo "=== DOWNLOADING SILVA 138.1 REFERENCE TAXONOMY ==="
echo "Downloading reference taxonomy..."
wget https://github.com/BenKaehler/readytowear/raw/master/data/silva_138_1/full_length/ref-tax.qza \
    -O "$OUTPUT_DIR/silva138.1_ref-tax.qza"

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Reference taxonomy downloaded"
    ls -lah "$OUTPUT_DIR/silva138.1_ref-tax.qza"
else
    echo "ERROR: Failed to download reference taxonomy"
    exit 1
fi

# Download plant-corpus weights
echo ""
echo "=== DOWNLOADING PLANT-CORPUS WEIGHTS ==="
echo "Downloading plant-corpus weights..."
wget https://github.com/BenKaehler/readytowear/raw/master/data/silva_138_1/full_length/plant-corpus.qza \
    -O "$OUTPUT_DIR/silva138.1_plant-corpus_weights.qza"

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Plant-corpus weights downloaded"
    ls -lah "$OUTPUT_DIR/silva138.1_plant-corpus_weights.qza"
else
    echo "ERROR: Failed to download plant-corpus weights"
    exit 1
fi

# Download plant-surface weights
echo ""
echo "=== DOWNLOADING PLANT-SURFACE WEIGHTS ==="
echo "Downloading plant-surface weights..."
wget https://github.com/BenKaehler/readytowear/raw/master/data/silva_138_1/full_length/plant-surface.qza \
    -O "$OUTPUT_DIR/silva138.1_plant-surface_weights.qza"

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Plant-surface weights downloaded"
    ls -lah "$OUTPUT_DIR/silva138.1_plant-surface_weights.qza"
else
    echo "ERROR: Failed to download plant-surface weights"
    exit 1
fi

# Verify downloaded files
echo ""
echo "=== VERIFYING DOWNLOADED FILES ==="
echo "Checking file types and integrity..."

echo "Reference sequences:"
qiime tools peek "$OUTPUT_DIR/silva138.1_ref-seqs.qza"

echo "Reference taxonomy:"
qiime tools peek "$OUTPUT_DIR/silva138.1_ref-tax.qza"

echo "Plant-corpus weights:"
qiime tools peek "$OUTPUT_DIR/silva138.1_plant-corpus_weights.qza"

echo "Plant-surface weights:"
qiime tools peek "$OUTPUT_DIR/silva138.1_plant-surface_weights.qza"

# Train plant-corpus classifier
echo ""
echo "=== TRAINING PLANT-CORPUS CLASSIFIER ==="
echo "Starting plant-corpus classifier training at: $(date)"
echo "This may take several hours with large reference database..."

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$OUTPUT_DIR/silva138.1_ref-seqs.qza" \
    --i-reference-taxonomy "$OUTPUT_DIR/silva138.1_ref-tax.qza" \
    --i-class-weight "$OUTPUT_DIR/silva138.1_plant-corpus_weights.qza" \
    --o-classifier "$OUTPUT_DIR/silva138.1_plant-corpus_weighted_classifier.qza" \
    --verbose

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Plant-corpus classifier training completed at: $(date)"
    
    # Get classifier info
    size_corpus=$(du -h "$OUTPUT_DIR/silva138.1_plant-corpus_weighted_classifier.qza" | cut -f1)
    echo "Plant-corpus classifier size: $size_corpus"
    
    # Verify classifier type
    echo "Plant-corpus classifier verification:"
    qiime tools peek "$OUTPUT_DIR/silva138.1_plant-corpus_weighted_classifier.qza"
else
    echo "ERROR: Plant-corpus classifier training failed"
    exit 1
fi

# Train plant-surface classifier
echo ""
echo "=== TRAINING PLANT-SURFACE CLASSIFIER ==="
echo "Starting plant-surface classifier training at: $(date)"
echo "This may take several hours with large reference database..."

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$OUTPUT_DIR/silva138.1_ref-seqs.qza" \
    --i-reference-taxonomy "$OUTPUT_DIR/silva138.1_ref-tax.qza" \
    --i-class-weight "$OUTPUT_DIR/silva138.1_plant-surface_weights.qza" \
    --o-classifier "$OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier.qza" \
    --verbose

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Plant-surface classifier training completed at: $(date)"
    
    # Get classifier info
    size_surface=$(du -h "$OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier.qza" | cut -f1)
    echo "Plant-surface classifier size: $size_surface"
    
    # Verify classifier type
    echo "Plant-surface classifier verification:"
    qiime tools peek "$OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier.qza"
else
    echo "ERROR: Plant-surface classifier training failed"
    exit 1
fi

echo ""
echo "=== PLANT CLASSIFIER TRAINING COMPLETE ==="

# Final summary
echo ""
echo "=== TRAINING SUMMARY ==="
echo "Files created in $OUTPUT_DIR:"
ls -lah "$OUTPUT_DIR"/silva138.1_plant-*

echo ""
echo "=== CLASSIFIER DETAILS ==="

if [[ -f "$OUTPUT_DIR/silva138.1_plant-corpus_weighted_classifier.qza" ]]; then
    echo "Plant-corpus classifier:"
    echo "- File: silva138.1_plant-corpus_weighted_classifier.qza"
    echo "- Size: $(du -h "$OUTPUT_DIR/silva138.1_plant-corpus_weighted_classifier.qza" | cut -f1)"
    echo "- Optimized for: Plant-associated microbiomes (rhizosphere, phyllosphere)"
fi

if [[ -f "$OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier.qza" ]]; then
    echo "Plant-surface classifier:"
    echo "- File: silva138.1_plant-surface_weighted_classifier.qza"
    echo "- Size: $(du -h "$OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier.qza" | cut -f1)"
    echo "- Optimized for: Plant surface microbiomes (leaf, flower, stem surfaces)"
fi

# Clean up temporary files
echo ""
echo "=== ALL FILES PRESERVED ==="
echo "All reference data and weights kept in: $OUTPUT_DIR"
echo "- silva138.1_ref-seqs.qza (reference sequences)"
echo "- silva138.1_ref-tax.qza (reference taxonomy)"  
echo "- silva138.1_plant-corpus_weights.qza (plant-corpus weights)"
echo "- silva138.1_plant-surface_weights.qza (plant-surface weights)"
echo "- silva138.1_plant-corpus_weighted_classifier.qza (trained classifier)"
echo "- silva138.1_plant-surface_weighted_classifier.qza (trained classifier)"

echo ""
echo "=== USAGE INSTRUCTIONS ==="
echo ""
echo "Your new plant-weighted classifiers are ready to use!"
echo ""
echo "For cacao flower microbiome analysis, use:"
echo "silva138.1_plant-surface_weighted_classifier.qza"
echo ""
echo "To test these classifiers, run:"
echo "qiime feature-classifier classify-sklearn \\"
echo "  --i-classifier $OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier.qza \\"
echo "  --i-reads your_sequences.qza \\"
echo "  --o-classification taxonomy_results.qza"

echo ""
echo "=== COMPARISON OPPORTUNITY ==="
echo ""
echo "You can now compare 4 different SILVA approaches:"
echo "1. silva-138-99-nb-classifier.qza (uniform weights)"
echo "2. silva-138-99-nb-diverse-weighted-classifier.qza (diverse weights)"
echo "3. silva138.1_plant-corpus_weighted_classifier.qza (plant-associated)"
echo "4. silva138.1_plant-surface_weighted_classifier.qza (plant surfaces)"
echo ""
echo "For flower microbiomes, expect plant-surface to perform best!"

echo ""
echo "Plant classifier training completed at: $(date)"
echo "Ready to classify cacao flower microbiomes with plant-optimized taxonomy!"
