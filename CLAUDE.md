I am a "wet" biologist transitioning to bioinformatics.

Global goal: complete data analysis of the cacao flower microbiome project and
prepare a public repository to accompany a research article.

## Project in one sentence

16S (bacteria) and ITS1 (fungi) amplicon survey of *Theobroma cacao* flowers
across 7 farms in Ghana; core question — does insect access (bagged vs
unbagged flowers) and pollination intensity (phenotyped using microscopy) 
shift the flower microbiome?

## Setup

- **R**: Windows R 4.4.2 via RStudio; renv locks packages (`renv::restore()`).
  From WSL terminal: `"/mnt/c/Program Files/R/R-4.4.2/bin/Rscript.exe" script.R`
- **Python**: WSL Ubuntu; `ml` conda env (scikit-learn stack); `cfm_analysis`
  conda env (Jupyter, biom-format, umap-learn)
- **HPC Hamilton**: SSH alias `hamilton`; `jcnx71@login2.ham8`;
  project root `~/nobackup/cacao_flower_microbiome/`; **access expires 2026-08-23**

## Canonical R entry points

| File | Contents |
|------|----------|
| `results/rds/ps_16S_bacteria_biosamples.rds` | 16S phyloseq — post-decontam, bacteria only |
| `results/rds/ps_ITS1_otu97_fungi_biosamples.rds` | ITS1 OTU97 phyloseq — post-decontam, fungi only |
| `results/rds/16S_dist_list.rds` + `16S_meta_list.rds` | 16S distance matrices (08a output) |
| `results/rds/ITS1_otu97_dist_list.rds` + `ITS1_otu97_meta_list.rds` | ITS1 distance matrices (08b output) |

## Key analytical decisions

- **Beta diversity primary metric**: rCLR Aitchison —
  `vegdist(..., method = "robust.aitchison")`; non-rarefied, no zero imputation
- **Prevalence filter**: 10% of samples, applied before distance computation in 08a/08b
- **Alpha diversity (main text)**: Observed, Shannon, Berger-Parker, Faith PD (16S
  only); Chao1 computed but NOT reported (downward biased after DADA2 removes global
  singletons)
- **16S classifier**: SILVA 138 diverse-weighted NB (conservative; avoids over-annotation)
- **ITS1 classifier**: custom UNITE eukaryotes v10 (Feb 2025); eukaryotes release for
  oomycete / Phytophthora coverage
- **ITS1**: OTUs at 97% (vsearch de-novo) to collapse intragenomic length polymorphism;
  no ITS phylogeny (too variable for reliable inference)
- **Decontam**: prevalence method, threshold 0.55

## Script conventions

- Pattern: `NN_<marker>_<step>.R` (e.g. `08a_16S_pcoa_global.R`)
- Each script: one RDS in → one analysis → one RDS out; no monolithic scripts
- Commit style: lowercase imperative subject line

## Current active work

1. Stabilize modular R pipeline; beta diversity (08–10a/b/c scripts) recently restructured
2. Zeta diversity on HPC — failed Slurm runs need retry before 2026-08-23
3. Manuscript methods draft lives in `docs/analysis_log.md` (local only, not in git)

## Datasets

- **CFM**: this project (primary)
- **_lewis**, **_schmidt**, **_wemheuer**: re-analysed published datasets; taxonomy
  comparison only, not full re-analysis
