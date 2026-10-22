# =============================================================================
# Merge the runs' sequence tables and remove chimeras (one marker, one rung)
# =============================================================================
# Called by rule dada2_merge_chimeras (workflow/rules/07_dada2_merge_runs.smk):
#   Rscript dada2_merge_chimeras.R <threads> <seqtab_nochim.rds> \
#       <seqtab_nochim.tsv> <track.tsv> <seqtab_run1.rds> [<seqtab_run2.rds>]
#
# Read counts of the same sample and sequence are summed across runs. Chimeras are
# removed with the settings v1 used through q2-dada2 (consensus, minimum parental
# fold 1, no one-off parents). The TSV has one row per sequence and one column per
# sample.

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
threads    <- as.integer(args[1])
rds_out    <- args[2]
tsv_out    <- args[3]
track_out  <- args[4]
seqtab_in  <- args[-(1:4)]

tables <- lapply(seqtab_in, readRDS)
# One table (pooled runs, rung A) needs no merging
merged <- if (length(tables) > 1) mergeSequenceTables(tables = tables, repeats = "sum") else tables[[1]]
cat("Before chimera removal:", ncol(merged), "sequences,", sum(merged), "read pairs\n")

nochim <- removeBimeraDenovo(merged, method = "consensus", minFoldParentOverAbundance = 1,
                             allowOneOff = FALSE, multithread = threads, verbose = TRUE)
cat("After chimera removal:", ncol(nochim), "ASVs,", sum(nochim), "read pairs",
    sprintf("(%.2f%% of reads kept)", 100 * sum(nochim) / sum(merged)), "\n")

saveRDS(nochim, rds_out)
data.frame(sequence = colnames(nochim), t(nochim), check.names = FALSE) |>
  write.table(tsv_out, sep = "\t", quote = FALSE, row.names = FALSE)
data.frame(sample = rownames(merged), merged = rowSums(merged),
           nonchim = rowSums(nochim)[rownames(merged)]) |>
  write.table(track_out, sep = "\t", quote = FALSE, row.names = FALSE)
