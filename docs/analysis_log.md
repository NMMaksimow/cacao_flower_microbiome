# Cacao Flower Microbiome Analysis
This document tracks the progress of NGS data analysis for the cacao flower microbiome project, examining the effects of pollinators on flower-associated microbial communities in 7 cocoa farms in Ghana. 
The microbial footprint of pollinators has been studied using a pollinator exclusion experimental design. Half of the flowers have been bagged before full anthesis, and the second half have been left to be openly pollinated. This experiment has been conducted across four management types: full sun, agroforest, near forest, and inside tropical forest.

## Background Information

### Sampling Design
- **Study system**: Cacao flowers (*Theobroma cacao*) across different farm management types
- **Locations**: 7 farms representing 4 management types across deforestation gradient: - Full sun (2 farms) - Agroforest (2 farms) - Near forest (2 farms) - Inside forest (1 farm)
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
- **Internal tags**: 8 bp combinatorial tags for sample-level identification within sublibraries; all bagged/openly pollinated flowers from one farm and corresponding controls have the same tag combination; different sets of internal tags were used for ITS1 and 16S
- **Structure**: 5'-[Illumina overhang]-[8bp tag]-[locus-specific primer]-3'
- **Sequencing**: Illumina NovaSeq paired-end sequencing, 2×250 bp reads; sequencing was performed in two lanes (L1 & L2)
- **Multiplexing**: Two-level approach enabling ~672 libraries (336 samples × 2 amplicons)

## Data Processing Workflow

### 1. Data Download and Verification
Raw data downloaded from Novogene CSS, uploaded to Hamilton cluster. Extract the tar file:
```bat
tar -xvf X204SC25042477-Z01-F001.tar 
```

Verify file integrity using MD5 checksums:
```bat
md5sum -c MD5.txt
```
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
**Status**: **COMPLETE** - All 48 sublibrary pairs merged successfully ---  

### 3. Adapter Trimming
**Remove Illumina sequencing adapters using cutadapt**  
- **Script**: `qiime2/scripts/02_adapter_trimming_cutadapt.sh`  
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
 **Status**: **COMPLETE** - Minimal adapter contamination, ready for demultiplexing ---

### 4. Demultiplexing with Stacks
**Generating 48 TSV mapping files for each sublibrary**  
- **Script**: `qiime2/scripts/03_creating_separate_stacks_mapping_files.py`  
- **Tool**: python/3.12.6  
- **Input**: `qiime2/import/demux/stacks_sample_mapping_all_sublibraries.txt` is a TSV file generated in Excel from `qiime2/import/demux/stacks_sample_all_sublibraries.xlsx`  
- **Output**: 48 TSV mapping files in directory `qiime2/import/demux/internal_tag_mappings/`  with names corresponding to sublibraries names with addition suffix: _tags.tsv  
- **Logs**: the job has been performed on the login node, stdout saved at `logs/03_creating_mapping_files_stacks.out`  
- **Resulting file** has three columns (first col = forward 8 bp tag, 2nd col = reverse 8 bp tag, 3rd col = sample-id with _16S or _ITS1 suffix at the end, resulting in 14 rows (7 biological samples in each sublibrary for 2 genetic markers each)

**Separate individual samples using internal combinatorial tags and removing tags**  
- **Script**: `qiime2/scripts/04_demultiplex_sublib_to_samples_stacks.sh`  
- **Tool**: Stacks process_radtags v2.64  
- **Input**: 48 adapter-trimmed sublibrary pairs + 48 mapping files generated in the previous step, using a Python script parsing metadata  
- **Strategy**: Combinatorial internal tags (8bp forward + reverse)  
- **NB!** I used a different set of 8 bp internal tags for 16S and ITS1 primers. It allows me to separate reads by amplicons at this stage, omitting sorting by primer step  
- **Key commands**:  
```bat
process_radtags \
    	-P \
        -1 "$r1_file" \
       	-2 "$r2_file" \
        -b "$mapping_file" \
       	-o "$temp_output" \
        --inline_inline \
       	--disable-rad-check \
        --retain-header \
	-c -q -r
```

   
**Processing Results**: \
- **Success rate**: 100% (48/48 sublibraries processed successfully)  
- **Individual samples extracted**: 672 samples (336 × 2 amplicons)  
- **Total files created**: 2,688 (1,344 per amplicon directory)  
- **Read depth per sample**: 46,000-125,000 reads (excellent coverage)  
- **File organization**: Separated by amplicon (16S and ITS1 directories)  
- **Processing time**: ~35 minutes on cluster node cn057  
- **SLURM job**: 13187722

### 5. Locus-Specific Primer Trimming
**Remove primer sequences using cutadapt**
- **Script**: `qiime2/scripts/05_primer_trimming_cutadapt.sh`
- **Tool**: cutadapt v4.9 with IUPAC wildcard support
- **Input**: 672 demultiplexed sample files (336 samples × 2 amplicons)
- **Primer sequences trimmed**:
  - 16S V4: 515f_modified (`GTGYCAGCMGCCGCGGTAA`) / 806r_modified (`GGACTACNVGGGTWTCTAAT`)
  - ITS1: ITS1f (`CTTGGTCATTTAGAGGAAGTAA`) / ITS2 (`GCTGCGTTCTTCATCGATGC`)
- **Parameters**: `--match-read-wildcards`, `--minimum-length 50`, `--maximum-length 300`

**Processing Results**:
- **Success rate**: 100% (336/336 samples processed successfully for both amplicons)
- **Total files created**: 1,344 primer-trimmed files
- **Read retention**: High retention rates with clean primer removal
- **Output location**: `qiime2/import/primer_trimmed/16S/` and `qiime2/import/primer_trimmed/ITS1/`
- **Processing time**: ~25 minutes
- **SLURM job**: Job completed successfully with no failures

**Status**: **COMPLETE** - Primer sequences successfully removed, sequences ready for QIIME2 import

### 6. QIIME2 Data Import
**Import primer-trimmed sequences into QIIME2 format**
- **Script**: `qiime2/scripts/06_qiime2_import.sh`
- **Input**: Primer-trimmed paired-end FASTQ files
- **Sample ID strategy**: Removed amplicon suffixes for compatibility with a single metadata file
- **Manifest creation**: Tab-separated format with absolute file paths

**Import Results**:
- **16S import**: Successfully created `CFM_16S_PE_import.qza` (2.6GB)
- **ITS1 import**: Successfully created `CFM_ITS1_PE_import.qza` (1.6GB)
- **Quality visualisations**: Generated QC reports for both amplicons
- **Sample count**: 336 samples per amplicon successfully imported
- **Processing time**: ~13 minutes total
- **SLURM job**: Completed without errors

**Status**: **COMPLETE** - Sequences successfully imported into QIIME2 format

### 7. DADA2 Denoising and ASV Calling
**Quality filtering, denoising, and feature table generation**
- **Script**: `qiime2/scripts/07_qiime2_dada2_denoising.sh`
- **Tool**: DADA2 within QIIME2 2024.10amplicon
- **Truncation parameters**: Based on quality plot analysis
  - 16S: Forward 242bp (no truncation), Reverse 240bp
  - ITS1: Forward 230bp, Reverse 225bp
- **Algorithm**: Paired-end denoising with chimaera removal

**DADA2 Results**:
- **16S processing**: Successfully generated ASV table and representative sequences
- **ITS1 processing**: Successfully generated ASV table and representative sequences  
- **Chimaera removal**: Consensus method applied to both datasets
- **Output files**: Feature tables, representative sequences, and denoising statistics for both amplicons
- **Memory usage**: 64GB, 16 cores
- **Processing time**: Several hours for both amplicons

**Status**: **COMPLETE** - ASV tables and representative sequences generated for downstream analysis

### 8. Taxonomic Classification
**Assign taxonomy using reference databases**
- **Script**: `qiime2/scripts/08_qiime2_taxonomy_classification.sh`
- **Databases used**:
  - **16S**: SILVA 138 diverse-weighted classifier (plant microbiome optimised)
  - **ITS1**: Custom-trained UNITE eukaryotes classifier (February 2025 release)

**UNITE Classifier Training**:
- **Training script**: `qiime2/scripts/08a_train_unite_classifier.sh`
- **Source data**: UNITE eukaryotes database v10.0 (February 19, 2025)
- **Clustering approach**: Dynamic clustering with expert-curated species boundaries
- **Database scope**: Comprehensive eukaryotes including fungi, Phytophthora, and other oomycetes
- **Sequence count**: 266,589 reference sequences with corresponding taxonomy
- **Training parameters**: 240GB memory, 48-hour time limit on shared queue
- **Output classifier**: 545MB trained classifier optimised for pathogen detection
- **Training time**: ~8 hours
- **Pathogen coverage**: Includes 266+ Phytophthora entries and related plant pathogens

- **Memory requirements**: 500GB on bigmem queue due to large classifier size

**Classification Results**:
- **16S taxonomy**: 1.6MB taxonomy file, 24MB interactive barplot
- **ITS1 taxonomy**: 559KB taxonomy file, 33MB interactive barplot  
- **Processing time**: 16S (51 minutes), ITS1 (2h 9m)
- **Output visualisations**: Taxonomy tables and interactive barplots for both amplicons
- **Pathogen detection**: UNITE classifier includes Phytophthora and other cacao pathogens
- **Total runtime**: ~3 hours on high-memory nodes

**Status**: **COMPLETE** - Taxonomic assignments completed for bacterial and fungal communities

## Current Dataset Summary
- **Total samples**: 336 biological samples (294 flowers + 42 controls)
- **Amplicons**: 16S rRNA V4 (bacteria/archaea) and ITS1 (fungi/eukaryotes)
- **ASV calling**: Completed using DADA2 with quality-based parameters
- **Taxonomy**: Assigned using state-of-the-art reference databases
- **Pathogen detection**: Custom eukaryotes classifier includes known cacao pathogens
- **Data quality**: High-quality reads with successful processing through the full pipeline

## Next Steps
- Alpha rarefaction analysis to determine sampling depth
- Diversity analysis (alpha and beta diversity)
- Statistical testing for pollinator effects
- Pathogen identification in ITS1 data
- Integration with microscopy data (pollination intensity, pollen germination)

## Project Directory Structure
cacao_flower_microbiome/  
├── cacao_flower_microbiome.Rproj  
├── renv.lock  
├── .gitignore  
├── .git/  
├── data/  
│   └── raw_data/  
├── docs/  
│   └── analysis_log.md  
├── logs/  
│   ├── 01_lane_merging.out  
│   ├── 02_adapter_trimming.err  
│   ├── 02_adapter_trimming.out  
│   ├── 03_creating_mapping_files_stacks.out  
│   ├── 04b_demux_comparison_20250618_1315.log  
│   ├── 04b_demux_comparison_20250618_1330.log  
│   ├── 04_demultiplex_stacks.err  
│   ├── 04_demultiplex_stacks.out  
│   ├── 04_demultiplex_stacks_no_cqr.err  
│   ├── 04_demultiplex_stacks_no_cqr.out  
│   ├── 05_primer_trimming.err  
│   ├── 05_primer_trimming.out  
│   ├── 06_qiime2_import.err  
│   ├── 06_qiime2_import.out  
│   ├── 07_dada2_denoising.err  
│   ├── 07_dada2_denoising.out  
│   ├── 08a_train_unite_classifier.err  
│   ├── 08a_train_unite_classifier.out  
│   ├── 08_taxonomy_classification.err  
│   └── 08_taxonomy_classification.out  
├── qiime2/  
│   ├── databases/  
│   │   ├── SILVA/  
│   │   │   └── silva-138-99-nb-diverse-weighted-classifier.qza  
│   │   └── UNITE/  
│   │       ├── QIIME_ITS_readme_19.02.2025.pdf  
│   │       ├── sh_qiime_release_s_all_19.02.2025.tgz  
│   │       ├── sh_refs_qiime_ver10_dynamic_s_all_19.02.2025.fasta  
│   │       ├── sh_taxonomy_qiime_ver10_dynamic_s_all_19.02.2025.txt  
│   │       ├── unite_eukaryotes_dynamic_classifier.qza  
│   │       ├── unite_eukaryotes_dynamic_sequences.qza  
│   │       └── unite_eukaryotes_dynamic_taxonomy.qza  
│   ├── denoise/  
│   │   ├── CFM_16S_PE_import.qza  
│   │   ├── CFM_16S_PE_import_QC.qzv  
│   │   ├── CFM_16S_dada2_table.qza  
│   │   ├── CFM_16S_dada2_table_summary.qzv  
│   │   ├── CFM_16S_dada2_repseqs.qza  
│   │   ├── CFM_16S_dada2_repseqs_summary.qzv  
│   │   ├── CFM_16S_dada2_stats.qza  
│   │   ├── CFM_16S_dada2_stats_summary.qzv  
│   │   ├── CFM_ITS1_PE_import.qza  
│   │   ├── CFM_ITS1_PE_import_QC.qzv  
│   │   ├── CFM_ITS1_dada2_table.qza  
│   │   ├── CFM_ITS1_dada2_table_summary.qzv  
│   │   ├── CFM_ITS1_dada2_repseqs.qza  
│   │   ├── CFM_ITS1_dada2_repseqs_summary.qzv  
│   │   ├── CFM_ITS1_dada2_stats.qza  
│   │   └── CFM_ITS1_dada2_stats_summary.qzv  
│   ├── import/  
│   │   ├── manifest_16S.tsv  
│   │   ├── manifest_ITS1.tsv  
│   │   ├── merged_files/ (48 sublibrary pairs)  
│   │   ├── trimmed_reads/ (48 adapter-trimmed sublibrary pairs)  
│   │   ├── primer_trimmed/  
│   │   │   ├── 16S/ (672 primer-trimmed sample files)  
│   │   │   └── ITS1/ (672 primer-trimmed sample files)  
│   │   └── demux/
│   │       ├── stacks_sample_mapping_all_sublibraries.txt  
│   │       ├── internal_tag_mappings/ (48 mapping files)  
│   │       ├── demultiplexed_sample_files/  
│   │       │   ├── 16S/ (672 demultiplexed files + .rem files)  
│   │       │   └── ITS1/ (672 demultiplexed files + .rem files)  
│   │       └── demultiplexed_sample_files_no_cqr/  
│   │           ├── 16S/ (672 demultiplexed files)  
│   │           └── ITS1/ (672 demultiplexed files)  
│   ├── scripts/  
│   │   ├── 01_lane_merging_PE.sh  
│   │   ├── 02_adapter_trimming_cutadapt.sh  
│   │   ├── 03_creating_separate_stacks_mapping_files.py  
│   │   ├── 04_demultiplex_sublib_to_samples_stacks.sh  
│   │   ├── 04a_demultiplex_sublib_to_samples_stacks_no_cqr.sh  
│   │   ├── 04b_compare_demux_results.sh  
│   │   ├── 05_primer_trimming_cutadapt.sh  
│   │   ├── 06_qiime2_import.sh  
│   │   ├── 07_qiime2_dada2_denoising.sh  
│   │   ├── 08a_train_unite_classifier.sh  
│   │   └── 08_qiime2_taxonomy_classification.sh  
│   └── taxonomy/  
│       ├── CFM_16S_taxonomy.qza  
│       ├── CFM_16S_taxonomy_viz.qzv  
│       ├── CFM_16S_taxonomy_barplot.qzv  
│       ├── CFM_ITS1_taxonomy.qza  
│       ├── CFM_ITS1_taxonomy_viz.qzv  
│       └── CFM_ITS1_taxonomy_barplot.qzv  
├── reports/  
├── results/  
└── scripts/  

