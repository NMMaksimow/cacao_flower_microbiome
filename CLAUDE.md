I am a "wet" biologist transitioning to bioinformatics.

Global goal: complete data analysis of the cacao flower microbiome project and
prepare a public repository to accompany a research article.

## Project in one sentence

16S (bacteria) and ITS1 (fungi) amplicon survey of *Theobroma cacao* flowers
across 7 farms in Ghana; core question — does insect access (bagged vs
unbagged flowers) and pollination intensity (phenotyped using microscopy) 
shift the flower microbiome?

## Setup

- **R**: Windows R 4.6.1 via RStudio; renv locks packages (`renv::restore()`).
  From WSL terminal: `"/mnt/c/Program Files/R/R-4.6.1/bin/Rscript.exe" script.R`
- **Python**: WSL Ubuntu, miniforge (conda/mamba). JupyterLab runs from `base`;
  notebook kernel `ml` (scikit-learn, biom-format, umap-learn; notebooks 01-03);
  host-genotype (04) and TDA (07a/07b, kernel `gamma`) notebooks are out of the
  paper's scope and live in the git-ignored `notebooks/exploratory/`
- **Notebooks in git**: outputs are stripped on commit by the nbstripout git
  filter (own env `nbstripout`); commit from the WSL terminal
- **HPC Hamilton**: SSH alias `hamilton`; `jcnx71@login2.ham8`;
  project root `~/nobackup/cacao_flower_microbiome/`;

## Canonical R entry points

Results are kept per processing version: `results/v1/` (v1 data, current R outputs)
and `results/v2/` (same layout, filled during the R switch-over). Paths below are
relative to `results/<version>/`.

| File | Contents |
|------|----------|
| `rds/ps_16S_bacteria_biosamples.rds` | 16S phyloseq — post-decontam, bacteria only |
| `rds/ps_ITS1_otu97_fungi_biosamples.rds` | ITS1 OTU97 phyloseq — post-decontam, fungi only |
| `rds/16S_dist_list.rds` + `16S_meta_list.rds` | 16S distance matrices (08a output) |
| `rds/ITS1_otu97_dist_list.rds` + `ITS1_otu97_meta_list.rds` | ITS1 distance matrices (08b output) |

## Key analytical decisions

- **Beta diversity primary metric**: CLR Aitchison with multiplicative zero
  replacement — `zCompositions::cmultRepl(method = "CZM", output = "p-counts")` →
  CLR → Euclidean (`aitchison_CLR` in 08a/08b); non-rarefied
- **rCLR** (`vegdist(..., method = "robust.aitchison")`): sensitivity metric only;
  it weights the zero/presence structure and inflated the ITS1 treatment R² (with
  significant PERMDISP) relative to CLR, Bray-Curtis and Jaccard
- **Prevalence filter**: 10% of samples, applied before distance computation in 08a/08b
- **Alpha diversity (main text)**: Observed, Shannon, Berger-Parker, Faith PD (16S
  only); Chao1 computed but NOT reported (downward biased after DADA2 removes global
  singletons)
- **16S classifier**: v1 SILVA 138 diverse-weighted NB (full-length); **v2 SILVA 144
  V4 515F/806R "average"** (official, DOI; habitat-averaged prior)
- **ITS1 classifier**: UNITE eukaryotes v10 (Feb 2025), eukaryotes release for
  oomycete / Phytophthora coverage; v1 custom-trained with singletons (s_all);
  **v2 community classifier without singletons** (more fungal reads to genus)
- **ITS1**: OTUs at 97% (vsearch de-novo) to collapse intragenomic length polymorphism;
  no ITS phylogeny (too variable for reliable inference)
- **Decontam**: prevalence method, threshold 0.55

## Script conventions

- Pattern: `NN_<marker>_<step>.R` (e.g. `08a_16S_pcoa_global.R`)
- Each script: one RDS in → one analysis → one RDS out; no monolithic scripts
- Commit style: lowercase imperative subject line

## Current active work

1. Stabilize modular R pipeline and results;
2. Stabilize ML approaches in notebooks;
3. Write manuscript methods and results; draft lives in `docs/cacao_flower_microbiome_nm.docx` (local only, not in git)

## v2 rerun: processing complete, R switch-over next

v1 left the 3' primer read-through untrimmed (~35% of ITS1 reads), merged the two
sequencing runs before DADA2 error learning and used an error model unsuited to
binned quality scores. v2 re-ran the amplicon processing (Phases 0-2 done,
2026-10-26); Phase 3 switches the R scripts and notebooks to v2. "Key analytical
decisions" above are otherwise unchanged.

- Pipeline: `workflow/Snakefile_v2` + `config_CFM_v2.yaml` (the v1
  `workflow/Snakefile` and `config_CFM.yaml` stay unchanged); v1 state is tag
  `v1-results`.
- v2 main choices: per-run DADA2 in R with the binned error model; 16S SILVA 144
  "average"; ITS1 UNITE no-singletons on 97% OTUs; three 16S trees, `bact` primary
  (`all`, `prok` are sensitivity analyses, decision rules fixed); ITS1 R1-only
  sensitivity branch (`CFM_ITS1_r1_*`).
- Exports for R: Hamilton `qiime2_v2/export/` (v1 folder names, plus
  `CFM_16S_taxonomy_average`, `CFM_16S_rooted_tree_{bact,all,prok}`, `RUN_INFO.txt`);
  QC tables and figures in `qiime2_v2/qc/`; SLURM logs in `logs/v2/`.
- For the R switch-over: SILVA 144 adds a `k__` rank (parse ranks by prefix) and
  renames genera (e.g. Aureimonas -> Mesocryomonas, Candidatus Tremblaya ->
  Tremblayella); attach a tree only after the taxonomy filter.
- QIIME 2 2026.4: `module load qiime2/2026.04qiime2`. Prefix every QIIME 2 and
  DADA2 command with `R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null`
  (`~/.Renviron` sets `R_LIBS_USER` for my own HPC R and breaks QIIME 2's R
  plugins; do not edit it).
- Snakemake 7.17.1: `module load snakemake/7.17.1`, SLURM via `--cluster`,
  profile `workflow/profiles/hamilton/`; code clone on Hamilton
  `~/nobackup/cfm_v2_code`.
- Decisions, results and phases: `docs/v2_rerun_plan.md` (local only, not in git).

## Datasets

- **CFM**: this project (primary)
- **_lewis**, **_schmidt**, **_wemheuer**: re-analysed published datasets; taxonomy
  comparison only, not full re-analysis
