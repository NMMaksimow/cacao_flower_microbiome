#!/bin/bash
#SBATCH --job-name=qiime2_import_corrected
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --gres=tmp:20G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/06b_qiime2_import_16S_corrected_manifest.out
#SBATCH --error=logs/06b_qiime2_import_16S_corrected_manifest.err

# Script to re-import 16S sequences with corrected sample identities
# Step 6b: Import corrected 16S manifest (fixes cmt53/cmtPC lab mix-up)
# Input: Manually corrected manifest file for 16S amplicon
# Output: Overwrites existing 16S QIIME2 artifacts with corrected sample identities

echo "Starting QIIME2 import with corrected 16S manifest..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Purpose: Fix cmt53/cmtPC sample identity mix-up"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
OUTPUT_DIR="qiime2/denoise"
MANIFEST_DIR="qiime2/import"
MANIFEST_16S="$MANIFEST_DIR/manifest_16S.tsv"

echo "Output directory: $OUTPUT_DIR"
echo "Using corrected manifest: $MANIFEST_16S"

# Verify manifest file exists
if [[ ! -f "$MANIFEST_16S" ]]; then
    echo "ERROR: Corrected manifest file not found: $MANIFEST_16S"
    exit 1
fi

echo "Manifest file found and verified"

# Check manifest content
echo "=== MANIFEST VERIFICATION ==="
echo "Manifest header:"
head -1 "$MANIFEST_16S"
echo ""
echo "Sample count:"
wc -l < "$MANIFEST_16S"
echo ""
echo "First 5 samples:"
head -6 "$MANIFEST_16S"
echo ""

# Check for corrected samples specifically
echo "Verifying corrected sample mappings:"
echo "cmt53 mapping:"
grep "^cmt53" "$MANIFEST_16S" || echo "cmt53 not found in manifest"
echo "cmtPC mapping:"
grep "^cmtPC" "$MANIFEST_16S" || echo "cmtPC not found in manifest"
echo ""

# Backup existing files
echo "=== BACKING UP EXISTING FILES ==="
if [[ -f "$OUTPUT_DIR/CFM_16S_PE_import.qza" ]]; then
    echo "Backing up existing 16S import artifact..."
    cp "$OUTPUT_DIR/CFM_16S_PE_import.qza" "$OUTPUT_DIR/CFM_16S_PE_import_backup_$(date +%Y%m%d_%H%M).qza"
fi

if [[ -f "$OUTPUT_DIR/CFM_16S_PE_import_QC.qzv" ]]; then
    echo "Backing up existing 16S QC visualization..."
    cp "$OUTPUT_DIR/CFM_16S_PE_import_QC.qzv" "$OUTPUT_DIR/CFM_16S_PE_import_QC_backup_$(date +%Y%m%d_%H%M).qzv"
fi

# Import 16S sequences with corrected manifest
echo "=== IMPORTING CORRECTED 16S SEQUENCES ==="
echo "Starting corrected 16S import at: $(date)"

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST_16S" \
    --output-path "$OUTPUT_DIR/CFM_16S_PE_import.qza" \
    --input-format PairedEndFastqManifestPhred33V2

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: Corrected 16S sequences imported successfully"

    # Get basic statistics
    echo "Artifact information:"
    qiime tools peek "$OUTPUT_DIR/CFM_16S_PE_import.qza"

    # Create visualization
    echo "Creating corrected 16S visualization..."
    qiime demux summarize \
        --i-data "$OUTPUT_DIR/CFM_16S_PE_import.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_PE_import_QC.qzv"

    echo "Corrected 16S import completed successfully at: $(date)"
else
    echo "ERROR: Corrected 16S import failed"
    exit 1
fi

echo ""
echo "=== IMPORT SUMMARY ==="

# Check output files
if [[ -f "$OUTPUT_DIR/CFM_16S_PE_import.qza" ]]; then
    size_16s=$(du -h "$OUTPUT_DIR/CFM_16S_PE_import.qza" | cut -f1)
    echo "Corrected 16S artifact created: $size_16s"
else
    echo "ERROR: Corrected 16S artifact not found"
fi

if [[ -f "$OUTPUT_DIR/CFM_16S_PE_import_QC.qzv" ]]; then
    size_qc=$(du -h "$OUTPUT_DIR/CFM_16S_PE_import_QC.qzv" | cut -f1)
    echo "Corrected 16S QC visualization created: $size_qc"
else
    echo "ERROR: Corrected 16S QC visualization not found"
fi

# Verify sample correction was applied
echo ""
echo "=== VERIFICATION ==="
echo "Showing corrected sample mappings from manifest:"
echo "These samples had their file paths swapped to fix lab mix-up:"
echo ""
echo "cmt53 mapping (now pointing to original cmtPC files):"
grep "^cmt53" "$MANIFEST_16S"
echo ""
echo "cmtPC mapping (now pointing to original cmt53 files):"
grep "^cmtPC" "$MANIFEST_16S"
echo ""

echo "Extracting sample IDs from imported artifact..."
qiime tools export \
    --input-path "$OUTPUT_DIR/CFM_16S_PE_import.qza" \
    --output-path temp_export_check

if [[ -f "temp_export_check/MANIFEST" ]]; then
    echo "Verifying samples are present in imported data:"
    if grep -q "cmt53" temp_export_check/MANIFEST; then
        echo "✓ cmt53 found in imported data"
    else
        echo "✗ cmt53 NOT found in imported data"
    fi
    
    if grep -q "cmtPC" temp_export_check/MANIFEST; then
        echo "✓ cmtPC found in imported data"
    else
        echo "✗ cmtPC NOT found in imported data"
    fi
    
    # Clean up
    rm -rf temp_export_check
else
    echo "Could not verify sample correction (export failed)"
fi

echo ""
echo "Files in $OUTPUT_DIR:"
ls -lah "$OUTPUT_DIR"/CFM_16S*

echo ""
echo "=== CORRECTION COMPLETE ==="
echo "Lab mix-up correction applied:"
echo "- Sample cmt53 now correctly mapped to its actual data"
echo "- Sample cmtPC now correctly mapped to its actual data"
echo ""
echo "Next steps:"
echo "1. Download and examine CFM_16S_PE_import_QC.qzv at https://view.qiime2.org"
echo "2. Verify sample identities are correct in the visualization"
echo "3. Proceed with DADA2 denoising using corrected data"
echo "4. Update metadata file to match corrected sample identities"

echo ""
echo "Corrected QIIME2 import completed at: $(date)"
echo "Ready for downstream analysis with correct sample identities!"
