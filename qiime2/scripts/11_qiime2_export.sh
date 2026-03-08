#!/bin/bash
#SBATCH --job-name=qiime_export
#SBATCH -p shared
#SBATCH --mem=8G
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --output=logs/11_qiime_export.out
#SBATCH --error=logs/11_qiime_export.err
#SBATCH --mail-type=FAIL

echo "Starting QIIME2 export"
echo "Date: $(date)"
echo "Host: $(hostname)"

module load bioinformatics
module load qiime2/2024.10amplicon

# ------------------------------------------------------------------
# Create export directory
# ------------------------------------------------------------------
EXPORT_DIR="qiime2/export"
mkdir -p ${EXPORT_DIR}

# ------------------------------------------------------------------
# 16S exports
# ------------------------------------------------------------------

echo "Exporting 16S feature table..."
qiime tools export \
  --input-path qiime2/denoise/CFM_16S_dada2_table.qza \
  --output-path ${EXPORT_DIR}/CFM_16S_dada2_table

echo "Exporting 16S representative sequences..."
qiime tools export \
  --input-path qiime2/denoise/CFM_16S_dada2_repseqs.qza \
  --output-path ${EXPORT_DIR}/CFM_16S_dada2_repseqs

echo "Exporting 16S taxonomy (diverse-weighted)..."
qiime tools export \
  --input-path qiime2/taxonomy/CFM_16S_taxonomy_diverse-weighted.qza \
  --output-path ${EXPORT_DIR}/CFM_16S_taxonomy_diverse_weighted

# ------------------------------------------------------------------
# ITS1 exports
# ------------------------------------------------------------------

echo "Exporting ITS1 feature table..."
qiime tools export \
  --input-path qiime2/denoise/CFM_ITS1_dada2_table.qza \
  --output-path ${EXPORT_DIR}/CFM_ITS1_dada2_table

echo "Exporting ITS1 representative sequences..."
qiime tools export \
  --input-path qiime2/denoise/CFM_ITS1_dada2_repseqs.qza \
  --output-path ${EXPORT_DIR}/CFM_ITS1_dada2_repseqs

echo "Exporting ITS1 taxonomy..."
qiime tools export \
  --input-path qiime2/taxonomy/CFM_ITS1_taxonomy.qza \
  --output-path ${EXPORT_DIR}/CFM_ITS1_taxonomy

echo "QIIME2 export finished"
echo "Date: $(date)"
