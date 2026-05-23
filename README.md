## 🌺 Cacao Flower Microbiome Project Overview

This repository contains the bioinformatics analysis pipeline for investigating the effects of pollinators on flower-associated microbial communities in cacao (*Theobroma cacao* L.) across deforestation gradients in Ghana.

## 🔬 Research Questions

1. **Microbial Community Characterization**: Describe bacterial and fungal *T. cacao* flower communities for the first time with a focus on pathogens and mutualistic symbionts;
2. **Management Effects**: Investigate how farm management type affects flower microbiomes along the deforestation gradient;
3. **Pollinator Footprint**: Identify the microbial signature of flower-visitng animals;
4. **Functional Correlations**: Examine relationships between microbiome structure and the success of the early step of fertilisation.

## 🏞️ Study Design

### **Experimental Framework & Sampling Strategy**
- **Visitint insect exclusion experiment**: 3 bagged + 3 openly pollinated flowers per tree
- **Deforestation gradient**: 4 management types across 7 farms in Ghana at different stages of extensive farming:
  - Full sun (2 farms)
  - Agroforest (2 farms) 
  - Near forest (2 farms)
  - Inside tropical forest (1 farm)
- **Biological samples**: 294 flowers (7 trees/farm × 6 flowers/tree)
- **Controls**: 42 total (14 extraction + 14 PCR negative + 14 mock communities)

### **Molecular Methods**
- **Marker genes**: 16S rRNA V4 (515f/806r, bacteria) + ITS1 (ITS1f/ITS2, fungi)
- **Sequencing platform**: Illumina NovaSeq 2×250bp
- **Multiplexing and indexing**: Two-level approach enabling high sample throughput with combinatorial dual-indexing (48 sublibraries, 336 samples per marker)

## 📁 Repository Structure

The project keeps a separate top-level folder per dataset and per pipeline stage. Original CFM data and three re-analysed published datasets coexist; QIIME2 bioinformatics, R statistics and Python notebooks each live in their own area.

| Path | Contents |
|------|---------|
| `data/` | CFM metadata, microscopy data, QIIME2 metadata file; raw FASTQ (gitignored, on HPC). |
| `data_lewis/`, `data_schmidt/`, `data_wemheuer/` | Metadata + SRA `runinfo.csv` for re-analysed public datasets ([Lewis et al. 2024](https://doi.org/10.1094/PHYTOFR-08-23-0104-R), [Schmidt et al. 2023](https://doi.org/10.1128/msphere.00013-23), [Wemheuer et al. 2020](https://doi.org/10.3390/microorganisms8030405)). |
| `qiime2/` | QIIME2 pipeline for CFM original data: `scripts/`, `import/`, `denoise/`, `taxonomy/`, `taxonomy_nb_classifier_comparison/`, `filtered/`, `phylogeny/`, `rarefaction/`, `export/`, `databases/` (SILVA, UNITE). |
| `qiime2_lewis/`, `qiime2_schmidt/`, `qiime2_wemheuer/` | QIIME2 artefacts for re-analyses (`import/`, `denoise/`, `taxonomy/`, `phylogeny/`, `export/`, `scripts/`) — produced by the Snakemake workflow. |
| `workflow/` | Snakemake workflow (`Snakefile`) and per-dataset configs (`config_CFM.yaml`, `config_lewis.yaml`, `config_schmidt.yaml`, `config_wemheuer.yaml`). Covers import → DADA2 → taxonomy → phylogeny → export. CFM pre-processing (lane merging, adapter/primer trimming, demultiplexing) is handled by standalone scripts in `qiime2/scripts/`. |
| `scripts/` | R and Python scripts for downstream analysis: decontamination, off-target filtering, diversity, ordinations, ANCOMBC, indicator species, PhILR, cross-dataset comparisons, zeta diversity, microscopy data analysis. |
| `notebooks/` | Jupyter notebooks for exploratory analysis: dimensionality reduction, random forest, ANCOMBC, gamma diversity, TDA. |
| `results/` | Analysis outputs: `biom/` (feature tables), `figures/`, `tables/`, `rds/` (phyloseq intermediates), `zeta/`. |
| `sra_upload/` | Demultiplexed FASTQ files prepared for NCBI SRA upload + metadata. |
| `logs/` | SLURM stdout/stderr logs from HPC runs (preprocessing, DADA2, taxonomy, phylogeny, Snakemake). |
| `docs/` | `analysis_log.md` — detailed step-by-step log of the bioinformatics pipeline and draft methods text. |

### Reproducibility

- **R**: open the RStudio project (`cacao_flower_microbiome.Rproj`) and run `renv::restore()` to install all R package versions pinned in `renv.lock`.
- **Python**: `conda env create -f environment.yml` to recreate the `cfm_analysis` environment (Jupyter, scikit-learn, scikit-bio, biom-format, umap-learn).
- **QIIME2**: pipeline tested on QIIME2 2024.10 amplicon distribution; see `workflow/Snakefile` for parameters.

## 🔄 Current Status

**⚠️ This is an ongoing analysis repository**

The analysis pipeline is under active development. Raw sequencing data will be deposited in public repositories upon publication.

## 🤝 Open Science & Reproducibility

All analysis scripts are documented for reproducibility. Development includes AI assistance for code optimisation. Methodology and results will be made fully available upon publication.

---

**Keywords**: microbiome, anthosphere, pollination ecology, *Theobroma cacao*, deforestation, 16S rRNA, ITS1, QIIME2, metabarcoding, Ghana, pollen-pistil interaction
