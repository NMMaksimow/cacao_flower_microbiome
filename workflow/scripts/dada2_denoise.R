# =============================================================================
# DADA2 denoising and read merging of one run (lane), marker and error model
# =============================================================================
# Called by rule dada2_denoise (workflow/rules/06_dada2_denoise.smk):
#   Rscript dada2_denoise.R <filtered_dir> <errors_r1.rds> <errors_r2.rds> \
#       <threads> <seqtab.rds> <track.tsv>
#
# dada(pool = FALSE) denoises every sample on its own with the run's error model;
# mergePairs joins R1 and R2 (overlap >= 12 nt, no mismatch, as v1) and, with
# trimOverhang = TRUE, removes any part of a read that extends past the start of
# its mate (3' read-through left after primer trimming). Output: the sequence
# table (samples x merged sequences) and per-sample read counts.

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
filtered_dir <- args[1]
err_r1       <- readRDS(args[2])
err_r2       <- readRDS(args[3])
threads      <- as.integer(args[4])
seqtab_out   <- args[5]
track_out    <- args[6]

r1 <- sort(list.files(filtered_dir, pattern = "_R1\\.fq\\.gz$", full.names = TRUE))
r2 <- sub("_R1\\.fq\\.gz$", "_R2.fq.gz", r1)
samples <- sub("_R1\\.fq\\.gz$", "", basename(r1))
cat("Samples:", length(samples), "\n")

denoised_r1 <- dada(r1, err = err_r1, pool = FALSE, multithread = threads) |> setNames(samples)
denoised_r2 <- dada(r2, err = err_r2, pool = FALSE, multithread = threads) |> setNames(samples)
merged <- mergePairs(denoised_r1, r1, denoised_r2, r2,
                     minOverlap = 12, maxMismatch = 0, trimOverhang = TRUE) |>
  setNames(samples)

seqtab <- makeSequenceTable(merged)
saveRDS(seqtab, seqtab_out)
cat("Merged sequences:", ncol(seqtab), "| merged read pairs:", sum(seqtab), "\n")

reads_in <- function(denoised) sum(getUniques(denoised))
data.frame(
  sample      = samples,
  denoised_r1 = vapply(denoised_r1, reads_in, numeric(1)),
  denoised_r2 = vapply(denoised_r2, reads_in, numeric(1)),
  merged      = rowSums(seqtab)[samples]
) |>
  write.table(track_out, sep = "\t", quote = FALSE, row.names = FALSE)
