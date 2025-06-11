# Cacao Flower Microbiome Analysis Log
This document tracks the progress of NGS data analysis for the cacao flower microbiome project.

=====================================================
 Background information on sampling and library prep
=====================================================

=====================================================
 Downloading, extracting and verifying the data
=====================================================

## Raw data downloaded from Novogene CSS, uploaded to Hamilton, extraction and integrity verification ##
# Extract the tar file
tar -xvf X204SC25042477-Z01-F001.tar

## Verify file integrity using MD5 checksums
md5sum -c MD5.txt 
X204SC25042477-Z01-F001.tar: OK

==========================
 Lane Merging
==========================

# Successfully merged L1 and L2 files for all 48 sub-libraries using cacao_flower_microbiome/qiime2/scripts/merge_PE_lanes.sh
# Output: 96 merged compressed R1.fq.gz and  R2.fq.gz files for each sublibrary were saved at cacao_flower_microbiome/qiime2/import/merged_files/
# Processing time: 29 seconds on cluster node cn032, log file stored at cacao_flower_microbiome/logs/slurm-13004770.out

==========================
 Adapter Trimming
==========================
## Illumina adapter removal using cutadapt ##

# Remove Illumina overhangs/adapters from multiplexed sublibrary lane merged files
# Script: cacao_flower_microbiome/qiime2/scripts/adapter_trimming_cutadapt.sh
# Input: 48 sublibrary pairs (96 files) from merged lane data
# Tool: cutadapt v4.9

# Adapter sequences removed:
# Forward (R1): TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG  
# Reverse (R2): GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

# Results:
# - Total read pairs processed: 890,963 across all sublibraries
# - Adapter contamination: <1% (20 adapters in R1, 7,133 in R2)
# - Data retention: 100% of reads preserved
# - Output: 96 adapter-trimmed files saved at qiime2/import/trimmed_reads/
# - Processing time: ~20 minutes, SLURM job 13163307
# - Log: logs/trim_adapters_13163307.out

# Status: ✅ COMPLETE - Ready for demultiplexing with internal barcodes

==========================
 Next: Demultiplexing
==========================
# TODO: Use Stacks process_radtags to separate samples by 8bp internal tags


