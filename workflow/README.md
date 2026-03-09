# Snakemake QIIME2 Pipeline (from historical cluster scripts)

This folder provides a minimal Snakemake wrapper around the working QIIME2/Stacks/Cutadapt scripts in `qiime2/scripts/`, using the *same* file naming and locations as the historical Hamilton cluster run (see `logs/`).

## What it builds

Default target (`rule all`) builds the exported files used by downstream R scripts:

- `qiime2/export/CFM_16S_dada2_table/feature-table.biom`
- `qiime2/export/CFM_16S_dada2_repseqs/dna-sequences.fasta`
- `qiime2/export/CFM_16S_taxonomy_diverse_weighted/taxonomy.tsv`
- `qiime2/export/CFM_ITS1_dada2_table/feature-table.biom`
- `qiime2/export/CFM_ITS1_dada2_repseqs/dna-sequences.fasta`
- `qiime2/export/CFM_ITS1_taxonomy/taxonomy.tsv`

## Requirements / assumptions

- Designed for HPC (SLURM) where `module load` works, matching the original scripts (QIIME2 `2024.10amplicon`, `cutadapt/4.9`, `stacks/2.64`).
- Large database artifacts are **not stored** in this repo (see `.gitignore`). Provide them under `qiime2/databases/` or adjust script paths:
  - SILVA classifier(s) for 16S.
  - UNITE eukaryotes classifier for ITS1.
- The raw Novogene archive is expected at `data/raw_data/X204SC25042477-Z01-F001.tar` (configurable).

## Run

Example:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config.yaml --cores 1
```

On SLURM you will typically use a profile (not included here yet), e.g.:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config.yaml --profile slurm
```

## Notes

- The workflow is intentionally coarse-grained: it runs the existing scripts as-is to minimize drift from the historical run recorded in `logs/`.
- Logs from Snakemake rules are written to `logs/snakemake_*.log`.
