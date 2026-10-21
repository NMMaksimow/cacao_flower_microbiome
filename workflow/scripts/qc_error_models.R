# =============================================================================
# QC: compare DADA2 error functions on mock communities and controls
# =============================================================================
# Called by rule qc_error_models (workflow/rules/qc_04_error_models.smk):
#   Rscript qc_error_models.R <marker> <dada2_dir> <lanes> <errfuns> \
#       <expected_fasta> <min_pairs> <threads> <out_tsv> <out_monotonicity>
#
# For each run (lane) and error function, the mock communities (*PC), extraction
# controls (*EC) and PCR negatives (*NC) are denoised with that run's error model
# (dada, pool = FALSE), merged (mergePairs, trimOverhang = TRUE) and chimera-checked
# (removeBimeraDenovo, consensus). Each ASV is classed by its edit distance to the
# nearest expected mock sequence (ZymoBIOMICS D6305, workflow/resources/mock/):
#   expected      = identical to an expected sequence
#   near_expected = 1-3 differences (likely an uncorrected sequencing error)
#   other         = anything else (contaminant, cross-talk from other samples)
# Output 1: one row per sample, run and function. Output 2: whether the fitted
# error rate of every substitution falls with quality across the binned Q values.
# Sample types come from the name suffix, which matches metadata sample_type.

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
marker         <- args[1]
dada2_dir      <- args[2]
lanes          <- strsplit(args[3], ",")[[1]]
errfuns        <- strsplit(args[4], ",")[[1]]
expected_fasta <- args[5]
min_pairs      <- as.integer(args[6])
threads        <- as.integer(args[7])
out_tsv        <- args[8]
out_mono       <- args[9]

expected      <- readDNAStringSet(expected_fasta) |> as.character()
expected_org  <- sub("\\|.*", "", names(expected))
expected_main <- grepl("\\|variant1\\|", names(expected))

# Fitted error rate per substitution at the observed quality values; monotone if
# it never increases from one quality value to the next.
monotonicity <- function(err, lane, errfun, read) {
  rates <- getErrors(err)
  rates <- rates[!rownames(rates) %in% c("A2A", "C2C", "G2G", "T2T"), , drop = FALSE]
  q_cols <- colnames(rates)[colSums(err$trans) > 0]
  data.frame(lane = lane, errfun = errfun, read = read, transition = rownames(rates),
             rates[, q_cols, drop = FALSE] |> `colnames<-`(paste0("rate_Q", q_cols)),
             monotone = apply(rates[, q_cols, drop = FALSE], 1, \(x) all(diff(x) <= 0)),
             row.names = NULL)
}

samples_tsv <- list()
mono_tsv    <- list()
for (lane in lanes) {
  run_dir  <- file.path(dada2_dir, lane, marker)
  filtered <- read.delim(file.path(run_dir, "filter_counts.tsv"))
  r1 <- list.files(file.path(run_dir, "filtered"), pattern = "(PC|EC|NC)_R1\\.fq\\.gz$", full.names = TRUE) |> sort()
  r2 <- sub("_R1\\.fq\\.gz$", "_R2.fq.gz", r1)
  samples <- sub("_R1\\.fq\\.gz$", "", basename(r1))

  for (errfun in errfuns) {
    cat("==", marker, lane, errfun, ":", length(samples), "samples\n")
    err <- lapply(c("R1", "R2"), \(read) readRDS(file.path(run_dir, sprintf("errors_%s_%s.rds", errfun, read))))
    mono_tsv[[paste(lane, errfun)]] <- rbind(monotonicity(err[[1]], lane, errfun, "R1"),
                                             monotonicity(err[[2]], lane, errfun, "R2"))

    denoised_r1 <- dada(r1, err = err[[1]], pool = FALSE, multithread = threads)
    denoised_r2 <- dada(r2, err = err[[2]], pool = FALSE, multithread = threads)
    merged <- mergePairs(denoised_r1, r1, denoised_r2, r2, trimOverhang = TRUE) |> setNames(samples)
    seqtab <- makeSequenceTable(merged)
    nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads)

    distance  <- adist(colnames(nochim), expected)          # edit distance, ASV x expected
    asv_class <- cut(apply(distance, 1, min), c(-Inf, 0, 3, Inf),
                     labels = c("expected", "near_expected", "other"))
    asv_match <- match(colnames(nochim), expected)           # index of an identical expected sequence

    samples_tsv[[paste(lane, errfun)]] <- lapply(rownames(nochim), \(s) {
      present <- nochim[s, ] > 0
      found   <- unique(asv_match[present & !is.na(asv_match)])
      data.frame(
        marker = marker, lane = lane, errfun = errfun, sample = s,
        type = sub(".*(PC|EC|NC)$", "\\1", s),
        pairs_filtered = filtered$pairs_out[filtered$sample == s],
        pairs_merged = sum(seqtab[s, ]),
        pairs_nonchim = sum(nochim[s, ]),
        asv_expected = sum(present & asv_class == "expected"),
        asv_near_expected = sum(present & asv_class == "near_expected"),
        asv_other = sum(present & asv_class == "other"),
        reads_expected = sum(nochim[s, asv_class == "expected"]),
        reads_near_expected = sum(nochim[s, asv_class == "near_expected"]),
        reads_other = sum(nochim[s, asv_class == "other"]),
        organisms_main_variant_found = length(unique(expected_org[found[expected_main[found]]])),
        expected_variants_found = length(found)
      )
    }) |> do.call(what = rbind)
  }
}

result <- do.call(rbind, samples_tsv)
result$in_mock_summary <- result$type == "PC" & result$pairs_filtered >= min_pairs
write.table(result, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, mono_tsv), out_mono, sep = "\t", quote = FALSE, row.names = FALSE)

# Short summary in the SLURM log: mocks in the summary, per run and function
mocks <- result[result$in_mock_summary, ]
aggregate(cbind(asv_expected, asv_near_expected, asv_other, reads_near_expected, reads_other,
                organisms_main_variant_found) ~ lane + errfun, data = mocks, FUN = mean) |>
  print()
