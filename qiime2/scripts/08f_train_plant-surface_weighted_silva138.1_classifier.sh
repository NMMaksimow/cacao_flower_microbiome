#!/bin/bash
#SBATCH --job-name=train_plant_surface
#SBATCH -p bigmem
#SBATCH --mem=240G
#SBATCH -c 48
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --gres=tmp:50G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/08f_train_plant_surface_classifier.out
#SBATCH --error=logs/08f_train_plant_surface_classifier.err

# Script to train plant-surface weighted SILVA 138.1 taxonomic classifier
# Step 8f: Complete plant-surface classifier training (continuation of 08e)
# Input: SILVA 138.1 reference data and plant-surface weights (from 08e)
# Output: Plant-surface weighted classifier for flower microbiome analysis
# Note: This script repeats part of 08e, which hasn't been finished due to time limit (request time limit: 2-3 days next time when training multiple classifiers)

echo "Starting plant-surface weighted SILVA 138.1 classifier training..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 240G"
echo "Purpose: Complete plant-surface classifier training (continuation of 08e)"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
OUTPUT_DIR="qiime2/databases/SILVA"

echo "Output directory: $OUTPUT_DIR"

# Verify required input files exist (downloaded in 08e)
echo ""
echo "=== VERIFYING INPUT FILES ==="

if [[ ! -f "$OUTPUT_DIR/silva138.1_ref-seqs.qza" ]]; then
    echo "ERROR: Reference sequences not found: $OUTPUT_DIR/silva138.1_ref-seqs.qza"
    echo "Please run 08e script first to download required files"
    exit 1
fi

if [[ ! -f "$OUTPUT_DIR/silva138.1_ref-tax.qza" ]]; then
    echo "ERROR: Reference taxonomy not found: $OUTPUT_DIR/silva138.1_ref-tax.qza"
    echo "Please run 08e script first to download required files"
    exit 1
fi

if [[ ! -f "$OUTPUT_DIR/silva138.1_plant-surface_weights.qza" ]]; then
    echo "ERROR: Plant-surface weights not found: $OUTPUT_DIR/silva138.1_plant-surface_weights.qza"
    echo "Please run 08e script first to download required files"
    exit 1
fi

echo "✓ Reference sequences found: $(du -h "$OUTPUT_DIR/silva138.1_ref-seqs.qza" | cut -f1)"
echo "✓ Reference taxonomy found: $(du -h "$OUTPUT_DIR/silva138.1_ref-tax.qza" | cut -f1)"
echo "✓ Plant-surface weights found: $(du -h "$OUTPUT_DIR/silva138.1_plant-surface_weights.qza" | cut -f1)"

# Check if output already exists
if [[ -f "$OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier.qza" ]]; then
    echo ""
    echo "WARNING: Plant-surface classifier already exists!"
    echo "Existing file: $(du -h "$OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier.qza" | cut -f1)"
    echo "This script will overwrite the existing classifier"
    echo "Creating backup..."
    
    backup_timestamp=$(date +%Y%m%d_%H%M)
    cp "$OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier.qza" \
       "$OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier_backup_${backup_timestamp}.qza"
    echo "Backup created: silva138.1_plant-surface_weighted_classifier_backup_${backup_timestamp}.qza"
fi

# Verify file types
echo ""
echo "=== VERIFYING FILE TYPES ==="
echo "Reference sequences:"
qiime tools peek "$OUTPUT_DIR/silva138.1_ref-seqs.qza"

echo "Reference taxonomy:"
qiime tools peek "$OUTPUT_DIR/silva138.1_ref-tax.qza"

echo "Plant-surface weights:"
qiime tools peek "$OUTPUT_DIR/silva138.1_plant-surface_weights.qza"

# Train plant-surface classifier
echo ""
echo "=== TRAINING PLANT-SURFACE CLASSIFIER ==="
echo "Starting plant-surface classifier training at: $(date)"
echo "Expected training time: 12-24 hours"
echo "This classifier will be optimized for plant surface microbiomes (ideal for flowers)"

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$OUTPUT_DIR/silva138.1_ref-seqs.qza" \
    --i-reference-taxonomy "$OUTPUT_DIR/silva138.1_ref-tax.qza" \
    --i-class-weight "$OUTPUT_DIR/silva138.1_plant-surface_weights.qza" \
    --o-classifier "$OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier.qza" \
    --verbose

# Check training success
if [[ $? -eq 0 ]]; then
    echo ""
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
echo "=== PLANT-SURFACE CLASSIFIER TRAINING COMPLETE ==="

# Final summary
echo ""
echo "=== TRAINING SUMMARY ==="
echo "Plant-surface classifier successfully created:"
echo "- File: silva138.1_plant-surface_weighted_classifier.qza"
echo "- Size: $(du -h "$OUTPUT_DIR/silva138.1_plant-surface_weighted_classifier.qza" | cut -f1)"
echo "- Optimized for: Plant surface microbiomes (leaf, flower, stem surfaces)"
echo "- Ideal for: Cacao flower microbiome analysis"

echo ""
echo "=== COMPLETE SILVA CLASSIFIER COLLECTION ==="
echo "Available SILVA classifiers in $OUTPUT_DIR:"
ls -lah "$OUTPUT_DIR"/*.qza | grep classifier

echo ""
echo "You now have 4 SILVA classifiers for comparison:"
echo "1. silva-138-99-nb-classifier.qza (uniform weights)"
echo "2. silva-138-99-nb-diverse-weighted-classifier.qza (diverse weights)"
echo "3. silva138.1_plant-corpus_weighted_classifier.qza (plant-associated)"
echo "4. silva138.1_plant-surface_weighted_classifier.qza (plant surfaces) ← NEW!"

echo ""
echo "Plant-surface classifier training completed at: $(date)"
echo "Training completed successfully!"
