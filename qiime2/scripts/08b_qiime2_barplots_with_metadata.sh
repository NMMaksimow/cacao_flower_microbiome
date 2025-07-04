#!/bin/bash

# This script generates QIIME2 visualizations coupled with metadata.

# Load required modules
module load bioinformatics
module load qiime2/2024.10amplicon

# Create metadata visualization within QIIME2
qiime metadata tabulate \
	--m-input-file data/cfm_qiime2_metadata.txt \
	--o-visualization data/cfm_metadata_summary.qzv

# Create 16S barplot visualization with metadata
qiime taxa barplot \
	--i-table qiime2/denoise/CFM_16S_dada2_table.qza \
	--i-taxonomy qiime2/taxonomy/CFM_16S_taxonomy.qza \
	--m-metadata-file data/cfm_qiime2_metadata.txt \
	--o-visualization qiime2/taxonomy/CFM_16S_taxonomy_barplot_with_metadata.qzv

# Create ITS1 barplot visualization with metadata
qiime taxa barplot \
	--i-table qiime2/denoise/CFM_ITS1_dada2_table.qza \
	--i-taxonomy qiime2/taxonomy/CFM_ITS1_taxonomy.qza \
	--m-metadata-file data/cfm_qiime2_metadata.txt \
	--o-visualization qiime2/taxonomy/CFM_ITS1_taxonomy_barplot_with_metadata.qzv


