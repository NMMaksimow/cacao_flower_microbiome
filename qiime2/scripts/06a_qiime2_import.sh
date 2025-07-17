#!/bin/bash
#SBATCH --job-name=qiime2_import
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --gres=tmp:30G
#SBATCH --mail-type=ALL
#SBATCH --output=logs/06_qiime2_import.out
#SBATCH --error=logs/06_qiime2_import.err

# Script to import primer-trimmed sequences into QIIME2
# Step 6: Create QIIME2 artifacts for 16S and ITS1 amplicons
# Input: Adapter and primer-trimmed demultiplexed paired-end .fq.gz files
# Output: QIIME2 artifacts (.qza) ready for DADA2 denoising and visualisations for quality control (.qzv)

echo "Starting QIIME2 import..."
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Set directories
INPUT_16S="qiime2/import/primer_trimmed/16S"
INPUT_ITS1="qiime2/import/primer_trimmed/ITS1"
OUTPUT_DIR="qiime2/denoise"
MANIFEST_DIR="qiime2/import"

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$MANIFEST_DIR"

echo "Input 16S: $INPUT_16S"
echo "Input ITS1: $INPUT_ITS1"
echo "Output directory: $OUTPUT_DIR"
echo "Manifest directory: $MANIFEST_DIR"

# Function to create manifest file
# One TSV manifest for 16S and another for ITS1
# Manifest files map sample-id (without locus specifing suffixes) to absolute paths to forward and reverse read files .fq.gz
create_manifest() {
    local input_dir=$1
    local manifest_file=$2
    local amplicon=$3
    
    echo "Creating manifest for $amplicon amplicon..."
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$manifest_file"
    
    local sample_count=0
    local base_path=$(realpath "$input_dir")
    
    # Process each sample (look for .1.trimmed.fq.gz files)
    for r1_file in "$input_dir"/*.1.trimmed.fq.gz; do
        if [[ -f "$r1_file" ]]; then
            # Extract sample name and remove amplicon suffix
            local sample_name=$(basename "$r1_file" .1.trimmed.fq.gz)
            sample_name=$(echo "$sample_name" | sed 's/_16S$//' | sed 's/_ITS1$//')
            
            local r2_file="$input_dir"/$(basename "$r1_file" .1.trimmed.fq.gz).2.trimmed.fq.gz
            
            # Check if R2 file exists
            if [[ -f "$r2_file" ]]; then
                # Get absolute paths
                local r1_abs="$base_path/$(basename "$r1_file")"
                local r2_abs="$base_path/$(basename "$r2_file")"
                
                # Add to manifest (tab-separated)
                echo -e "$sample_name\t$r1_abs\t$r2_abs" >> "$manifest_file"
                sample_count=$((sample_count + 1))
            else
                echo "WARNING: R2 file not found for $sample_name"
            fi
        fi
    done
    
    echo "Manifest created with $sample_count samples for $amplicon"
    echo "First 5 entries:"
    head -6 "$manifest_file"
    echo ""
}

# Create manifest files
echo "=== CREATING MANIFEST FILES ==="
MANIFEST_16S="$MANIFEST_DIR/manifest_16S.tsv"
MANIFEST_ITS1="$MANIFEST_DIR/manifest_ITS1.tsv"

create_manifest "$INPUT_16S" "$MANIFEST_16S" "16S"
create_manifest "$INPUT_ITS1" "$MANIFEST_ITS1" "ITS1"

echo "Verifying manifest format:"
echo "16S manifest header:"
head -1 "$MANIFEST_16S"
echo "ITS1 manifest header:"
head -1 "$MANIFEST_ITS1"

# Import 16S sequences
echo "=== IMPORTING 16S SEQUENCES ==="
echo "Starting 16S import at: $(date)"

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST_16S" \
    --output-path "$OUTPUT_DIR/CFM_16S_PE_import.qza" \
    --input-format PairedEndFastqManifestPhred33V2

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: 16S sequences imported successfully"
    
    # Get basic statistics
    qiime tools peek "$OUTPUT_DIR/CFM_16S_PE_import.qza"
    
    # Create visualization
    echo "Creating 16S visualization..."
    qiime demux summarize \
        --i-data "$OUTPUT_DIR/CFM_16S_PE_import.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_16S_PE_import_QC.qzv"
    
    echo "16S import completed successfully at: $(date)"
else
    echo "ERROR: 16S import failed"
    exit 1
fi

echo ""

# Import ITS1 sequences
echo "=== IMPORTING ITS1 SEQUENCES ==="
echo "Starting ITS1 import at: $(date)"

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST_ITS1" \
    --output-path "$OUTPUT_DIR/CFM_ITS1_PE_import.qza" \
    --input-format PairedEndFastqManifestPhred33V2

if [[ $? -eq 0 ]]; then
    echo "SUCCESS: ITS1 sequences imported successfully"
    
    # Get basic statistics
    qiime tools peek "$OUTPUT_DIR/CFM_ITS1_PE_import.qza"
    
    # Create visualization
    echo "Creating ITS1 visualization..."
    qiime demux summarize \
        --i-data "$OUTPUT_DIR/CFM_ITS1_PE_import.qza" \
        --o-visualization "$OUTPUT_DIR/CFM_ITS1_PE_import_QC.qzv"
    
    echo "ITS1 import completed successfully at: $(date)"
else
    echo "ERROR: ITS1 import failed"
    exit 1
fi

echo ""
echo "=== IMPORT SUMMARY ==="

# Check output files
if [[ -f "$OUTPUT_DIR/CFM_16S_PE_import.qza" ]]; then
    size_16s=$(du -h "$OUTPUT_DIR/CFM_16S_PE_import.qza" | cut -f1)
    echo "16S artifact created: $size_16s"
else
    echo "ERROR: 16S artifact not found"
fi

if [[ -f "$OUTPUT_DIR/CFM_ITS1_PE_import.qza" ]]; then
    size_its1=$(du -h "$OUTPUT_DIR/CFM_ITS1_PE_import.qza" | cut -f1)
    echo "ITS1 artifact created: $size_its1"
else
    echo "ERROR: ITS1 artifact not found"
fi

# List all created files
echo ""
echo "Files created in $OUTPUT_DIR:"
ls -lah "$OUTPUT_DIR"

echo ""
echo "Next steps:"
echo "1. Download and examine .qzv files at https://view.qiime2.org"
echo "2. Determine DADA2 trimming parameters from quality plots"
echo "3. Run DADA2 denoising"

echo ""
echo "QIIME2 import completed at: $(date)"
echo "Ready for DADA2 denoising!"
