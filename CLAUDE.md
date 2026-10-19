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
  notebook kernels: `ml` (scikit-learn, biom-format, umap-learn; notebooks
  01-04) and `gamma` (gudhi; TDA notebooks 07a/07b)
- **Notebooks in git**: outputs are stripped on commit by the nbstripout git
  filter (own env `nbstripout`); commit from the WSL terminal
- **HPC Hamilton**: SSH alias `hamilton`; `jcnx71@login2.ham8`;
  project root `~/nobackup/cacao_flower_microbiome/`;

## Canonical R entry points

| File | Contents |
|------|----------|
| `results/rds/ps_16S_bacteria_biosamples.rds` | 16S phyloseq — post-decontam, bacteria only |
| `results/rds/ps_ITS1_otu97_fungi_biosamples.rds` | ITS1 OTU97 phyloseq — post-decontam, fungi only |
| `results/rds/16S_dist_list.rds` + `16S_meta_list.rds` | 16S distance matrices (08a output) |
| `results/rds/ITS1_otu97_dist_list.rds` + `ITS1_otu97_meta_list.rds` | ITS1 distance matrices (08b output) |

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

1. Stabilize modular R pipeline and results;
2. Stabilize ML approaches in notebooks;
3. Write manuscript methods and results; draft lives in `docs/cacao_flower_microbiome_nm.docx` (local only, not in git)

## v2 rerun in progress

v1 ITS1 data carry untrimmed 3' primer read-through (step 05), so amplicon
processing is being re-run. Until v2 results replace them, "Key analytical
decisions" above describe v1.

- Branch `v2-rerun`; v1 state is tag `v1-results`. R analyses continue on
  `master` against v1 data.
- Main changes from v1: anchored 5' primer trimming plus 3' read-through
  removal; DADA2 in R with error models learned per sequencing run (v1 merged
  both runs before error learning) and an error function suited to NovaSeq
  binned quality scores; 16S taxonomy with SILVA 144 classifiers (choice
  pending); UNITE v10 classifier retrained in QIIME 2 2026.4. An ablation
  ladder separates the effects of these changes (see the plan).
- Hamilton outputs: `qiime2_v2/` (exports in `qiime2_v2/export/`, same layout as
  v1), SLURM logs in `logs/v2/`; `results/v2/` comes with the R switch-over.
- Pipeline code in `workflow/` with its own `config_CFM_v2.yaml`; the v1
  `workflow/Snakefile` and `config_CFM.yaml` stay unchanged.
- QIIME 2 2026.4: `module load qiime2/2026.04qiime2`. Prefix every QIIME 2 and
  DADA2 command with `R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null`
  (`~/.Renviron` sets `R_LIBS_USER` for my own HPC R and breaks QIIME 2's R
  plugins; do not edit it).
- Snakemake 7.17.1: `module load snakemake/7.17.1`, SLURM via `--cluster`,
  profile `workflow/profiles/hamilton/`.
- Full plan and phases: `docs/v2_rerun_plan.md` (local only, not in git).

## Datasets

- **CFM**: this project (primary)
- **_lewis**, **_schmidt**, **_wemheuer**: re-analysed published datasets; taxonomy
  comparison only, not full re-analysis
