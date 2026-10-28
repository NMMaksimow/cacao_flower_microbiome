## ============================================================================
## DADA2 denoising of the forward reads (R1) only: one run (lane), ITS1
## ============================================================================
## Called by rule dada2_denoise_r1 (workflow/rules/15_its1_r1_only.smk):
##   Rscript dada2_denoise_r1.R <filtered_dir> <errors_r1.rds> <threads> \
##       <max_length> <seqtab.rds> <track.tsv>
##
## ITS1 sensitivity branch: pairs whose amplicon is longer than ~430 nt (e.g.
## Saccharomyces, ~440 nt) cannot merge with 2x250 reads, so the main data lose
## them. Here R1 is denoised alone, from the same filtered reads and with the same
## error model as the main pipeline (dada2_denoise.R), pool = FALSE.
## Long amplicons give R1 reads of 220 or 221 nt for the same organism, which would
## become two ASVs; every sequence is therefore cut to max_length and the columns
## that become identical are summed. Shorter reads are complete short amplicons
## (their 3' primer was trimmed) and keep their length.

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
filtered_dir <- args[1]
err_r1       <- readRDS(args[2])
threads      <- as.integer(args[3])
max_length   <- as.integer(args[4])
seqtab_out   <- args[5]
track_out    <- args[6]

r1 <- sort(list.files(filtered_dir, pattern = "_R1\\.fq\\.gz$", full.names = TRUE))
samples <- sub("_R1\\.fq\\.gz$", "", basename(r1))
cat("Samples:", length(samples), "\n")

denoised <- dada(r1, err = err_r1, pool = FALSE, multithread = threads) |> setNames(samples)

# Sequence table (samples x sequences), then sequences cut to max_length and summed
seqtab <- makeSequenceTable(denoised) |>
  (\(tab) t(rowsum(t(tab), group = substr(colnames(tab), 1, max_length), reorder = FALSE)))()
storage.mode(seqtab) <- "integer"
saveRDS(seqtab, seqtab_out)
cat("R1 sequences after cutting to", max_length, "nt:", ncol(seqtab), "| reads:", sum(seqtab), "\n")

data.frame(
  sample      = samples,
  denoised_r1 = vapply(denoised, \(d) sum(getUniques(d)), numeric(1)),
  in_table    = rowSums(seqtab)[samples]
) |>
  write.table(track_out, sep = "\t", quote = FALSE, row.names = FALSE)
