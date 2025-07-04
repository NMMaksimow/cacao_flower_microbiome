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
# Step 10: Alpha rarefaction curves for both 16S and ITS1 datasets
# Input: Filtered feature tables
# Output: Alpha rarefaction curves and depth recommendations

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

echo "Input denoise directory: $INPUT_DIR_DENOISE"
echo "Input filtered directory: $INPUT_DIR_FILTERED"
echo "Output directory: $OUTPUT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Verify input files exist
if [[ ! -f "$INPUT_DIR_DENOISE/CFM_16S_dada2_table.qza" ]]; then
    echo "ERROR: 16S unfiltered table not found: $INPUT_DIR_DENOISE/CFM_16S_dada2_table.qza"
    exit 1
fi

if [[ ! -f "$INPUT_DIR_FILTERED/CFM_16S_final_filtered.qza" ]]; then
    echo "ERROR: 16S filtered table (full) not found: $INPUT_DIR_FILTERED/CFM_16S_final_filtered.qza"
    exit 1
fi

if [[ ! -f "$INPUT_DIR_FILTERED/CFM_16S_prokaryotes_only.qza" ]]; then
    echo "WARNING: 16S prokaryotes-only table not found: $INPUT_DIR_FILTERED/CFM_16S_prokaryotes_only.qza"
    echo "Will skip prokaryotes-only rarefaction analysis"
    SKIP_PROKARYOTES=true
else
    SKIP_PROKARYOTES=false
fi

if [[ ! -f "$INPUT_DIR_DENOISE/CFM_ITS1_dada2_table.qza" ]]; then
    echo "ERROR: ITS1 unfiltered table not found: $INPUT_DIR_DENOISE/CFM_ITS1_dada2_table.qza"
    exit 1
fi

if [[ ! -f "$INPUT_DIR_FILTERED/CFM_ITS1_final_filtered.qza" ]]; then
    echo "ERROR: ITS1 filtered table not found: $INPUT_DIR_FILTERED/CFM_ITS1_final_filtered.qza"
    exit 1
fi

echo "Input files verified successfully"

# ============================================================================
# UNFILTERED DATA RAREFACTION (for comparison and quality assessment)
# ============================================================================

echo ""
echo "=== UNFILTERED DATA RAREFACTION ANALYSIS ==="
echo "Purpose: Quality assessment and comparison with filtered data"
echo "Starting unfiltered rarefaction analysis at: $(date)"

# 16S Unfiltered Alpha Rarefaction
echo ""
echo "--- 16S UNFILTERED RAREFACTION ---"
echo "Parameters:"
echo "  Max depth: 60,000 reads"
echo "  Steps: 50"
echo "  Iterations: 10"

qiime diversity alpha-rarefaction \
    --i-table "$INPUT_DIR_DENOISE/CFM_16S_dada2_table.qza" \
    --p-max-depth 60000 \
    --p-steps 50 \
    --p-iterations 10 \
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
echo "Parameters:"
echo "  Max depth: 40,000 reads"
echo "  Steps: 50"
echo "  Iterations: 10"

qiime diversity alpha-rarefaction \
    --i-table "$INPUT_DIR_DENOISE/CFM_ITS1_dada2_table.qza" \
    --p-max-depth 40000 \
    --p-steps 50 \
    --p-iterations 10 \
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
# FILTERED DATA RAREFACTION (for actual diversity analysis)
# ============================================================================

echo ""
echo "=== FILTERED DATA RAREFACTION ANALYSIS ==="
echo "Purpose: Determine optimal sampling depth for diversity analysis"
echo "Starting filtered rarefaction analysis at: $(date)"

# 16S Filtered Alpha Rarefaction
echo ""
echo "--- 16S FILTERED RAREFACTION ---"
echo "Starting 16S rarefaction analysis at: $(date)"
echo "Parameters:"
echo "  Max depth: 50,000 reads"
echo "  Steps: 50"
echo "  Iterations: 10"

# 16S Alpha Rarefaction
qiime diversity alpha-rarefaction \
    --i-table "$INPUT_DIR_FILTERED/CFM_16S_final_filtered.qza" \
    --p-max-depth 50000 \
    --p-steps 50 \
    --p-iterations 10 \
    --o-visualization "$OUTPUT_DIR/CFM_16S_alpha_rarefaction.qzv" \
    --verbose

# Check 16S rarefaction success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: 16S alpha rarefaction completed at: $(date)"
else
    echo "ERROR: 16S alpha rarefaction failed"
    exit 1
fi

echo ""

echo "=== ITS1 ALPHA RAREFACTION ==="
echo "Starting ITS1 rarefaction analysis at: $(date)"
echo "Parameters:"
echo "  Max depth: 30,000 reads"
echo "  Steps: 50"
echo "  Iterations: 10"

# ITS1 Alpha Rarefaction
qiime diversity alpha-rarefaction \
    --i-table "$INPUT_DIR_FILTERED/CFM_ITS1_final_filtered.qza" \
    --p-max-depth 30000 \
    --p-steps 50 \
    --p-iterations 10 \
    --o-visualization "$OUTPUT_DIR/CFM_ITS1_alpha_rarefaction.qzv" \
    --verbose

# Check ITS1 rarefaction success
if [[ $? -eq 0 ]]; then
    echo "SUCCESS: ITS1 alpha rarefaction completed at: $(date)"
else
    echo "ERROR: ITS1 alpha rarefaction failed"
    exit 1
fi

echo ""
echo "=== ALPHA RAREFACTION COMPLETE ==="

# ============================================================================
# GENERATE RAREFACTION SUMMARY REPORT
# ============================================================================

echo "Generating rarefaction summary report..."

# Create a Python script to analyze rarefaction depths
cat > "$OUTPUT_DIR/analyze_rarefaction_depth.py" << 'EOF'
#!/usr/bin/env python3
"""
Script to extract key information from QIIME2 feature table summaries
to help determine optimal rarefaction depths.
"""

import pandas as pd
import zipfile
import json
import sys
import os

def extract_sample_stats(qzv_file):
    """Extract sample frequency statistics from QIIME2 visualization file."""
    with zipfile.ZipFile(qzv_file, 'r') as zip_ref:
        # Find the data file
        file_list = zip_ref.namelist()
        data_file = None
        for file in file_list:
            if file.endswith('sample-frequency-detail.csv'):
                data_file = file
                break
        
        if data_file:
            with zip_ref.open(data_file) as f:
                df = pd.read_csv(f)
                return df
        else:
            print(f"Could not find sample frequency data in {qzv_file}")
            return None

def analyze_rarefaction_depth(df, dataset_name):
    """Analyze sample depths and suggest rarefaction parameters."""
    if df is None:
        return
    
    print(f"\n=== {dataset_name} RAREFACTION ANALYSIS ===")
    print(f"Total samples: {len(df)}")
    print(f"Min reads per sample: {df['frequency'].min():,}")
    print(f"Max reads per sample: {df['frequency'].max():,}")
    print(f"Mean reads per sample: {df['frequency'].mean():.0f}")
    print(f"Median reads per sample: {df['frequency'].median():.0f}")
    
    # Calculate percentiles
    percentiles = [10, 25, 50, 75, 90, 95, 99]
    print("\nPercentiles:")
    for p in percentiles:
        value = df['frequency'].quantile(p/100)
        print(f"  {p}th percentile: {value:,.0f}")
    
    # Suggest rarefaction depths
    print("\nSuggested rarefaction depths:")
    conservative = df['frequency'].quantile(0.10)  # 10th percentile
    moderate = df['frequency'].quantile(0.25)      # 25th percentile
    aggressive = df['frequency'].quantile(0.50)    # 50th percentile
    
    print(f"  Conservative (retain 90% samples): {conservative:,.0f}")
    print(f"  Moderate (retain 75% samples): {moderate:,.0f}")
    print(f"  Aggressive (retain 50% samples): {aggressive:,.0f}")
    
    # Count samples at different thresholds
    thresholds = [1000, 5000, 10000, 15000, 20000, 25000, 30000, 40000, 50000]
    print("\nSamples retained at different thresholds:")
    for threshold in thresholds:
        retained = (df['frequency'] >= threshold).sum()
        percent = (retained / len(df)) * 100
        print(f"  {threshold:,} reads: {retained} samples ({percent:.1f}%)")

if __name__ == "__main__":
    # Analyze 16S data
    qzv_16s = "qiime2/filtered/CFM_16S_final_filtered_summary.qzv"
    if os.path.exists(qzv_16s):
        df_16s = extract_sample_stats(qzv_16s)
        analyze_rarefaction_depth(df_16s, "16S BACTERIAL")
    else:
        print(f"16S summary file not found: {qzv_16s}")
    
    # Analyze ITS1 data
    qzv_its1 = "qiime2/filtered/CFM_ITS1_final_filtered_summary.qzv"
    if os.path.exists(qzv_its1):
        df_its1 = extract_sample_stats(qzv_its1)
        analyze_rarefaction_depth(df_its1, "ITS1 FUNGAL")
    else:
        print(f"ITS1 summary file not found: {qzv_its1}")
EOF

# Make the script executable
chmod +x "$OUTPUT_DIR/analyze_rarefaction_depth.py"

# Run the rarefaction analysis
echo "Running rarefaction depth analysis..."
python "$OUTPUT_DIR/analyze_rarefaction_depth.py" > "$OUTPUT_DIR/rarefaction_analysis_report.txt"

echo ""
echo "Checking output files..."

# List all created files
echo "Files created in $OUTPUT_DIR:"
ls -lah "$OUTPUT_DIR"/CFM_*

echo ""
echo "=== RAREFACTION ANALYSIS SUMMARY ==="

# Check file sizes for both unfiltered and filtered
echo "Unfiltered rarefaction curves:"
if [[ -f "$OUTPUT_DIR/CFM_16S_unfiltered_alpha_rarefaction.qzv" ]]; then
    size_16s_unfiltered=$(du -h "$OUTPUT_DIR/CFM_16S_unfiltered_alpha_rarefaction.qzv" | cut -f1)
    echo "- 16S unfiltered rarefaction: $size_16s_unfiltered"
fi

if [[ -f "$OUTPUT_DIR/CFM_ITS1_unfiltered_alpha_rarefaction.qzv" ]]; then
    size_its1_unfiltered=$(du -h "$OUTPUT_DIR/CFM_ITS1_unfiltered_alpha_rarefaction.qzv" | cut -f1)
    echo "- ITS1 unfiltered rarefaction: $size_its1_unfiltered"
fi

echo ""
echo "Filtered rarefaction curves:"
if [[ -f "$OUTPUT_DIR/CFM_16S_filtered_full_alpha_rarefaction.qzv" ]]; then
    size_16s_filtered_full=$(du -h "$OUTPUT_DIR/CFM_16S_filtered_full_alpha_rarefaction.qzv" | cut -f1)
    echo "- 16S filtered full dataset: $size_16s_filtered_full"
fi

if [[ -f "$OUTPUT_DIR/CFM_16S_prokaryotes_only_alpha_rarefaction.qzv" ]]; then
    size_16s_prokaryotes=$(du -h "$OUTPUT_DIR/CFM_16S_prokaryotes_only_alpha_rarefaction.qzv" | cut -f1)
    echo "- 16S prokaryotes-only dataset: $size_16s_prokaryotes"
fi

if [[ -f "$OUTPUT_DIR/CFM_ITS1_filtered_alpha_rarefaction.qzv" ]]; then
    size_its1_filtered=$(du -h "$OUTPUT_DIR/CFM_ITS1_filtered_alpha_rarefaction.qzv" | cut -f1)
    echo "- ITS1 filtered dataset: $size_its1_filtered"
fi

echo ""
echo "Output files generated:"
echo "UNFILTERED DATA (for comparison):"
echo "- CFM_16S_unfiltered_alpha_rarefaction.qzv"
echo "- CFM_ITS1_unfiltered_alpha_rarefaction.qzv"
echo ""
echo "FILTERED DATA (for analysis):"
echo "- CFM_16S_filtered_full_alpha_rarefaction.qzv (includes eukaryotic DNA)"
if [[ "$SKIP_PROKARYOTES" == "false" ]]; then
    echo "- CFM_16S_prokaryotes_only_alpha_rarefaction.qzv (bacterial microbiome only)"
fi
echo "- CFM_ITS1_filtered_alpha_rarefaction.qzv (fungal communities)"
echo ""
echo "ANALYSIS TOOLS:"
echo "- rarefaction_analysis_report.txt (comprehensive depth recommendations)"
echo "- analyze_rarefaction_depth.py (enhanced analysis script)"

echo ""
echo "Next steps:"
echo "1. Download and view ALL .qzv files at https://view.qiime2.org"
echo "2. Compare unfiltered vs filtered rarefaction curves"
echo "3. Compare 16S full vs prokaryotes-only datasets"
echo "4. Review the comprehensive text report for recommended depths"
echo "5. Choose your analysis approach:"
echo "   - Use 16S FULL for pathogen-inclusive analysis"
echo "   - Use 16S PROKARYOTES for standard bacterial microbiome"
echo "   - Use both approaches for comprehensive publication"

echo ""
echo "Important comparisons to make:"
echo "- Unfiltered vs filtered curves show filtering impact"
echo "- 16S full vs prokaryotes-only shows eukaryotic DNA contribution"
echo "- Use appropriate dataset recommendations for downstream analysis"
echo "- Document both approaches in your methods section"

echo ""
echo "Alpha rarefaction analysis completed at: $(date)"
echo "Ready for phylogenetic tree construction with optimal depths!"
