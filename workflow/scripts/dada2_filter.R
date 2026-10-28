# =============================================================================
# DADA2 quality filtering of one sequencing run (lane) and marker
# =============================================================================
# Called by rule dada2_filter (workflow/rules/04_dada2_filter.smk):
#   Rscript dada2_filter.R <in_dir> <out_dir> <counts_tsv> <threads> \
#       <max_ee_r1> <max_ee_r2> <max_n> <min_len> <trunc_q> <qual_bins>
#
# filterAndTrim without truncation (ITS1 is variable in length), keeping read
# pairs with expected errors <= max_ee, no N and both mates >= min_len. Samples
# with no surviving pair get no output file. Afterwards the quality values present
# in the filtered reads are compared with the bins of the binned error function;
# the script stops if any other value occurs.

suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
})

args <- commandArgs(trailingOnly = TRUE)
in_dir     <- args[1]
out_dir    <- args[2]
counts_tsv <- args[3]
threads    <- as.integer(args[4])
max_ee     <- as.numeric(args[5:6])
max_n      <- as.integer(args[7])
min_len    <- as.integer(args[8])
trunc_q    <- as.integer(args[9])
qual_bins  <- as.integer(strsplit(args[10], ",")[[1]])

r1_in   <- sort(list.files(in_dir, pattern = "_R1\\.fq\\.gz$", full.names = TRUE))
r2_in   <- sub("_R1\\.fq\\.gz$", "_R2.fq.gz", r1_in)
samples <- sub("_R1\\.fq\\.gz$", "", basename(r1_in))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Samples:", length(samples), "\n")

counts <- filterAndTrim(
  fwd = r1_in, filt = file.path(out_dir, basename(r1_in)),
  rev = r2_in, filt.rev = file.path(out_dir, basename(r2_in)),
  maxEE = max_ee, maxN = max_n, minLen = min_len, truncQ = trunc_q,
  rm.phix = TRUE, compress = TRUE, multithread = threads
)

data.frame(sample = samples, pairs_in = counts[, "reads.in"], pairs_out = counts[, "reads.out"]) |>
  write.table(counts_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Read pairs in:", sum(counts[, "reads.in"]), " out:", sum(counts[, "reads.out"]),
    sprintf("(%.2f%%)", 100 * sum(counts[, "reads.out"]) / sum(counts[, "reads.in"])),
    "| samples without surviving pairs:", sum(counts[, "reads.out"] == 0), "\n")

# Quality values in the filtered reads (first 10,000 reads of every file)
observed_q <- list.files(out_dir, pattern = "\\.fq\\.gz$", full.names = TRUE) |>
  lapply(\(path) {
    stream <- FastqStreamer(path, n = 10000)
    on.exit(close(stream))
    as(quality(yield(stream)), "matrix") |> as.vector() |> unique()
  }) |>
  unlist() |>
  setdiff(NA) |>
  sort()
cat("Quality values in filtered reads:", observed_q, "| bins:", qual_bins, "\n")
if (!all(observed_q %in% qual_bins)) {
  stop("Quality values outside the bins: ", paste(setdiff(observed_q, qual_bins), collapse = ", "))
}
