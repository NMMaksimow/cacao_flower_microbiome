#!/bin/bash
#SBATCH --job-name=train_unite_classifier
#SBATCH -p shared
#SBATCH --mem=240G
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -t 2-00:00:00
#SBATCH --gres=tmp:80G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/08a_train_unite_classifier.out
#SBATCH --error=logs/08a_train_unite_classifier.err

# Script to train custom UNITE eukaryotes classifier for ITS1
# Includes fungi AND other eukaryotes (especially Oomycetes incl. Phytophthora megakarya and P. palmivora)
# Input: UNITE raw sequences (.fasta) and taxonomy files (.txt)
# Output: Custom trained classifier (.qza)

echo "Starting UNITE classifier training..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 240G"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
DB_DIR="qiime2/databases/UNITE"
OUTPUT_DIR="qiime2/databases/UNITE"

echo "Database directory: $DB_DIR"
echo "Output directory: $OUTPUT_DIR"

# Define clustering approach (dynamic = expert-curated, most accurate)
CLUSTERING="dynamic"
FASTA_FILE="$DB_DIR/sh_refs_qiime_ver10_${CLUSTERING}_s_all_19.02.2025.fasta"
TAXONOMY_FILE="$DB_DIR/sh_taxonomy_qiime_ver10_${CLUSTERING}_s_all_19.02.2025.txt"

echo "Using clustering approach: $CLUSTERING"
echo "FASTA file: $FASTA_FILE"
echo "Taxonomy file: $TAXONOMY_FILE"

# Verify input files exist
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "ERROR: FASTA file not found: $FASTA_FILE"
    echo "Make sure UNITE archive is extracted"
    exit 1
fi

if [[ ! -f "$TAXONOMY_FILE" ]]; then
    echo "ERROR: Taxonomy file not found: $TAXONOMY_FILE"
    echo "Make sure UNITE archive is extracted"
    exit 1
fi

echo "Input files verified successfully"

# Check input file statistics
echo ""
echo "=== INPUT FILE STATISTICS ==="
sequences=$(grep -c "^>" "$FASTA_FILE")
taxonomy_lines=$(wc -l < "$TAXONOMY_FILE")
echo "Number of sequences: $sequences"
echo "Number of taxonomy entries: $taxonomy_lines"

echo ""
echo "=== STEP 1: IMPORT REFERENCE SEQUENCES ==="
echo "Starting sequence import at: $(date)"

# Import FASTA sequences
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "$FASTA_FILE" \
    --output-path "$OUTPUT_DIR/unite_eukaryotes_${CLUSTERING}_sequences.qza"

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Reference sequences imported"
else
    echo "ERROR: Failed to import reference sequences"
    exit 1
fi

echo ""
echo "=== STEP 2: IMPORT TAXONOMY ==="
echo "Starting taxonomy import at: $(date)"

# Import taxonomy file
qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path "$TAXONOMY_FILE" \
    --output-path "$OUTPUT_DIR/unite_eukaryotes_${CLUSTERING}_taxonomy.qza"

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Taxonomy imported"
else
    echo "ERROR: Failed to import taxonomy"
    exit 1
fi

echo ""
echo "=== STEP 3: TRAIN NAIVE BAYES CLASSIFIER ==="
echo "Starting classifier training at: $(date)"
echo "This may take several hours..."

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$OUTPUT_DIR/unite_eukaryotes_${CLUSTERING}_sequences.qza" \
    --i-reference-taxonomy "$OUTPUT_DIR/unite_eukaryotes_${CLUSTERING}_taxonomy.qza" \
    --o-classifier "$OUTPUT_DIR/unite_eukaryotes_${CLUSTERING}_classifier.qza" \
    --verbose

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Classifier training completed at: $(date)"
else
    echo "ERROR: Classifier training failed"
    exit 1
fi

echo ""
echo "=== TRAINING COMPLETE ==="

# Check output files
echo "Generated files:"
ls -lah "$OUTPUT_DIR"/unite_eukaryotes_${CLUSTERING}_*

# Get file sizes
classifier_size=$(du -h "$OUTPUT_DIR/unite_eukaryotes_${CLUSTERING}_classifier.qza" | cut -f1)
sequences_size=$(du -h "$OUTPUT_DIR/unite_eukaryotes_${CLUSTERING}_sequences.qza" | cut -f1)
taxonomy_size=$(du -h "$OUTPUT_DIR/unite_eukaryotes_${CLUSTERING}_taxonomy.qza" | cut -f1)

echo ""
echo "=== OUTPUT SUMMARY ==="
echo "Classifier: $classifier_size"
echo "Sequences: $sequences_size"  
echo "Taxonomy: $taxonomy_size"

echo ""
echo "Classifier file: unite_eukaryotes_${CLUSTERING}_classifier.qza"
echo "This classifier includes:"
echo "- All fungi from UNITE"
echo "- Phytophthora and other oomycetes"
echo "- Other eukaryotic microorganisms"
echo "- Expert-curated species boundaries (dynamic clustering)"

echo ""
echo "Next steps:"
echo "1. Update taxonomy classification script to use this classifier"
echo "2. Run taxonomic classification on your ITS1 data"
echo "3. Look for Phytophthora and other cacao pathogens!"

echo ""
echo "UNITE classifier training completed at: $(date)"
echo "Ready for pathogen detection in cacao flowers!"
