# Cacao Flower Microbiome Analysis Log
This document tracks the progress of NGS data analysis for the cacao flower microbiome project, examining the effects of pollinators on flower-associated microbial communities in West Africa. 
The microbial footprint of pollinators has been studied using a pollinator exclusion experimental design. Half of the flowers have been bagged before full anthesis, and the second half have been left to be openly pollinated. This experiment has been conducted across four management types: full sun, agroforest, near forest, and inside tropical forest.

## Background Information

### Sampling Design
- **Study system**: Cacao flowers (*Theobroma cacao*) across different farm management types
- **Locations**: 7 farms representing 4 management types: - Full sun (2 farms) - Agroforest (2 farms) - Near forest (2 farms) - Inside forest (1 farm)
- **Sampling strategy**: 7 trees per farm, 6 flowers per tree (3 bagged at bud stage + 3 openly pollinated) 
- **Total biological samples**: 294 flowers
- **Controls**: 14 extraction controls, 14 negative PCR controls, 14 mock communities
- **Research questions**: 1) Describe bacterial and fungal T. cacao flower communities for the first time, 2) Investigate effect of management type on microbiome, 3) Identify microbial footprint of pollinators, 4) Using microscopy data investigate correlation between microbiome structure and pollination intensity (PI) and pollen germination rate.

### Library Preparation
- **DNA extraction**: Separate extractions for bagged vs. open-pollinated flowers (14 batches total)
- **Amplicon targets**:
- 16S rRNA V4 region (bacterial communities; 515f_modified/806r_modified)
- ITS1 region (fungal communities; ITS1f/ITS2)
- **Primer design**: Combinatorial dual-indexing approach
- **External barcodes**: Sublibrary-level identification (48 sublibraries total)
- **Internal tags**: 8bp combinatorial tags for sample-level identification within sublibraries
- **Structure**: 5'-[Illumina overhang]-[8bp tag]-[locus-specific primer]-3'
- **Sequencing**: Illumina NovaSeq paired-end sequencing, 2×250 bp reads
- **Multiplexing**: Two-level approach enabling ~672 libraries (336 samples × 2 amplicons)

## Data Processing Workflow

### 1. Data Download and Verification
Raw data downloaded from Novogene CSS, uploaded to Hamilton cluster. Extract the tar file:

`bash`tar -xvf X204SC25042477-Z01-F001.tar`bash` 

Verify file integrity using MD5 checksums:

md5sum -c MD5.txt

**Result:** X204SC25042477-Z01-F001.tar: OK
 **Status**: ✅ **COMPLETE** - Data integrity verified ---

### 2. Lane Merging
**Combine L1 and L2 sequencing lanes for each sublibrary**
	- **Script**: `qiime2/scripts/merge_PE_lanes.sh`
	- **Input**: Raw .fq.gz files from 2 sequencing lanes
	- **Output**: 96 merged files (48 sublibraries × 2 read directions)
	- **Location**: `qiime2/import/merged_files/`
	- **Processing time**: 29 seconds on cluster node cn032
	- **Log file**: `logs/slurm-13004770.out`
**Status**: ✅ **COMPLETE** - All 48 sublibrary pairs merged successfully ---

### 3. Adapter Trimming
**Remove Illumina sequencing adapters using cutadapt**
	- **Script**: `qiime2/scripts/adapter_trimming_cutadapt.sh`
	- **Tool**: cutadapt v4.9
	- **Input**: 48 sublibrary pairs (96 files) from merged lane data
	- **Adapter sequences removed**:
	 - Forward (R1): `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`
	 - Reverse (R2): `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`
**Results**:
	- **Total read pairs processed**: 890,963 across all sublibraries
	- **Adapter contamination**: <1% (20 adapters in R1, 7,133 in R2)
	- **Data retention**: 100% of reads preserved
	- **Output**: 96 adapter-trimmed files at `qiime2/import/trimmed_reads/`
	- **Processing time**: ~20 minutes (SLURM job 13163307)
	- **Log file**: `logs/trim_adapters_13163307.out`
 **Status**: ✅ **COMPLETE** - Minimal adapter contamination, ready for demultiplexing ---

## Next Steps
	### 4. Demultiplexing (In Progress)
**Separate individual samples using internal barcodes**
	- **Tool**: Stacks process_radtags
	- **Strategy**: Use 8bp combinatorial internal tags to separate samples within sublibraries
	- **Expected output**: Individual sample files for each of 336 samples × 2 amplicons
	- **Challenge**: Combinatorial tags are reused across sublibraries, requiring sublibrary-specific processing

### 5. Amplicon Separation (Planned)
**Separate 16S and ITS1 amplicons**
	- **Tool**: cutadapt with primer-specific trimming
	- **Strategy**: Use locus-specific primer sequences to separate and trim amplicons
	- **Output**: Clean amplicon-specific files ready for QIIME2 import

### 6. QIIME2 Analysis (Planned)
**Standard microbiome analysis pipeline**
	- Import demultiplexed, amplicon-separated sequences
	- Quality filtering and denoising (DADA2)
	- Taxonomic classification
	- Diversity analysis
	- Statistical testing for pollinator effects

## File Organization
``` cacao_flower_microbiome/ ├── data/ │ ├── qiime2_cfm_metadata.txt # Sample metadata for QIIME2 │ └── raw_data/ # Original Novogene files ├── qiime2/ │ ├── import/ │ │ ├── merged_files/ # Lane-merged sublibraries │ │ 
└── trimmed_reads/ # Adapter-trimmed sublibraries │ └── scripts/ # Processing scripts ├── logs/ # SLURM job logs └── docs/ # Documentation ``` ---

*Last updated: June 11, 2025*
