# Stacks barcode files for v2 demultiplexing

One file per sublibrary (48), used by `process_radtags -b` in
`workflow/rules/02_demultiplex.smk`. Columns: forward inline tag, reverse inline
tag, sample name with amplicon suffix (`_16S` or `_ITS1`).

The files are copies of `qiime2/import/demux/internal_tag_mappings/` (v1, made by
`qiime2/scripts/03_creating_separate_stacks_mapping_files.py`) with one deliberate
change.

## cmt53 / cmtPC correction (16S only)

A lab mix-up swapped the 16S libraries of sample `cmt53` and the mock community
`cmtPC`. Both use the tag pair `GAGTCACT` / `CGTAGGAA`, in two different
sublibraries. The correction swaps the two 16S sample names:

| file | v1 name | v2 name |
|---|---|---|
| `A3_13A_13_03_tags.tsv` | `cmt53_16S` | `cmtPC_16S` |
| `H3_20H_13_03_tags.tsv` | `cmtPC_16S` | `cmt53_16S` |

The ITS1 lines are unchanged. In v1 the same correction lived only in
`qiime2/import/manifest_16S.tsv` (swapped file paths, step 06b), so a rerun from
raw reads would have lost it. The SRA deposit already carries the corrected names
(`sra_upload/additional_files_for_sra/README.md.txt`).
