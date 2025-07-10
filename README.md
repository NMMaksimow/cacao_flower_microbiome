## 🌺 Cacao Flower Microbiome Project Overview

This repository contains the bioinformatics analysis pipeline for investigating the effects of pollinators on flower-associated microbial communities in cacao (*Theobroma cacao* L.) across deforestation gradients in Ghana.

## 🔬 Research Questions

1. **Microbial Community Characterization**: Describe bacterial and fungal *T. cacao* flower communities for the first time with a focus on pathogens and mutualistic symbionts;
2. **Management Effects**: Investigate how farm management type affects flower microbiomes along the deforestation gradient;
3. **Pollinator Footprint**: Identify the microbial signature of flower-visitng animals;
4. **Functional Correlations**: Examine relationships between microbiome structure and the success of the early step of fertilisation.

## 🏞️ Study Design

### **Experimental Framework**
- **Pollinator Exclusion Experiment**: Bagged vs. openly pollinated flowers
- **Deforestation Gradient**: 4 management types across 7 farms in Ghana
  - Full sun (2 farms)
  - Agroforest (2 farms) 
  - Near forest (2 farms)
  - Inside tropical forest (1 farm)

### **Sampling Strategy**
- **Biological samples**: 294 flowers (7 trees/farm × 6 flowers/tree)
- **Controls**: 42 total (14 extraction + 14 PCR negative + 14 mock communities)
- **Design**: 3 bagged + 3 openly pollinated flowers per tree

### **Molecular Methods**
- **Targets**: 16S rRNA V4 (bacteria) + ITS1 (fungi)
- **Platform**: Illumina NovaSeq 2×250bp
- **Indexing**: Combinatorial dual-indexing (48 sublibraries, 672 total libraries)
- **Multiplexing**: Two-level approach enabling high sample throughput

## 📁 Repository Structure

```
cacao_flower_microbiome/
├── docs/                          # Analysis documentation
│   └── analysis_log.md            # Detailed progress tracking
├── data/                          # Metadata and small data files
├── qiime2/                        # QIIME2 analysis pipeline
│   ├── scripts/                   # Analysis scripts (01-10+)
│   ├── import/                    # Raw data import
│   ├── denoise/                   # DADA2 ASV calling
│   ├── taxonomy/                  # Taxonomic classification
│   ├── filtered/                  # Quality-filtered datasets
│   └── rarefaction/               # Alpha rarefaction analysis
├── logs/                          # SLURM job outputs
└── README.md                      # This file
```

## 🔄 Current Status

**⚠️ This is an ongoing analysis repository**

The analysis pipeline is under active development. Raw sequencing data will be deposited in public repositories upon publication.

## 🤝 Open Science & Reproducibility

All analysis scripts are documented for reproducibility. Development includes AI assistance for code optimisation. Methodology and results will be made fully available upon publication.

---

**Keywords**: microbiome, anthosphere, pollination ecology, *Theobroma cacao*, deforestation, 16S rRNA, ITS1, QIIME2, metabarcoding, Ghana, pollen-pistil interaction
