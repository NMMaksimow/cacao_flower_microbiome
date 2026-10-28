# =============================================================================
# DADA2 error model of one sequencing run (lane), marker and error function
# =============================================================================
# Called by rule dada2_learn_errors (workflow/rules/05_dada2_learn_errors.smk):
#   Rscript dada2_learn_errors.R <filtered_dir> <errfun> <qual_bins> <nbases> \
#       <randomize> <seed> <threads> <rds_r1> <rds_r2> <png_r1> <png_r2> <title>
#
# <errfun>: "loess"  = DADA2 default loessErrfun (smooth fit over all Q values);
#           "binned" = makeBinnedQualErrfun(qual_bins), which fits error rates only
#                      at the quality values the sequencer reports (NovaSeq bins).
# learnErrors draws reads from randomly chosen samples (fixed seed) until nbases
# bases are collected, then alternates denoising and fitting until the error rates
# no longer change.

suppressPackageStartupMessages({
  library(dada2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
filtered_dir <- args[1]
errfun       <- args[2]
qual_bins    <- as.integer(strsplit(args[3], ",")[[1]])
nbases       <- as.numeric(args[4])
randomize    <- as.logical(args[5])
seed         <- as.integer(args[6])
threads      <- as.integer(args[7])
rds_out      <- args[8:9]
png_out      <- args[10:11]
title        <- args[12]

error_function <- switch(errfun,
  loess  = loessErrfun,
  binned = makeBinnedQualErrfun(qual_bins)
)

for (i in 1:2) {
  read <- c("R1", "R2")[i]
  files <- sort(list.files(filtered_dir, pattern = sprintf("_%s\\.fq\\.gz$", read), full.names = TRUE))
  cat("==", read, ":", length(files), "files\n")

  set.seed(seed)
  err <- learnErrors(files, nbases = nbases, errorEstimationFunction = error_function,
                     randomize = randomize, multithread = threads, verbose = TRUE)
  saveRDS(err, rds_out[i])

  # plotErrors: one panel per substitution (e.g. A2C = A read as C);
  #   points       = observed error rates per quality score,
  #   black line   = rates fitted by the error function,
  #   red line     = rates expected from the nominal quality score (nominalQ = TRUE).
  # A good fit follows the points and falls with increasing quality.
  p_errors <- plotErrors(err, nominalQ = TRUE) +
    ggtitle(sprintf("%s | %s | %s", title, read, errfun)) +   # figure title: run | marker | read | function
    theme(plot.title = element_text(size = 11))                # title font size
  ggsave(png_out[i], plot = p_errors, width = 10, height = 8, dpi = 300)   # figure size in inches
}
