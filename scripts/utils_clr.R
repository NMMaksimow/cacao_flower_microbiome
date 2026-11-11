## ============================================================================
## utils_clr — prevalence filter and CLR with CZM zero replacement
##
## Sourced by 08a, 08b (Aitchison CLR distance) and 08c (CLR matrices exported for
## the Python random forest notebooks 03, 03b, 03c). One function, so the PCoA,
## the PERMANOVA and the random forest all see the same CLR values.
##
## clr_czm(): zCompositions::cmultRepl(method = "CZM") replaces a zero in sample i
## with 0.65 * 0.5 / N_i (N_i = reads in sample i, about a third of one read),
## capped at 0.65 x the smallest observed proportion of that feature; non-zero
## values are rescaled so the row keeps its total. vegan::decostand("clr") then
## takes log(x) - mean(log(x)) per sample.
## ============================================================================

# Features present (count > 0) in at least `thresh` of the reference samples.
#   otu_t     samples x features count matrix
#   ref_rows  rows used to compute prevalence (NULL = all samples)
# Returns a logical vector over the columns of otu_t.
prevalence_keep <- function(otu_t, thresh, ref_rows = NULL) {
  ref <- if (is.null(ref_rows)) otu_t else otu_t[ref_rows, , drop = FALSE]
  colMeans(ref > 0) >= thresh
}

# CLR after CZM zero replacement. otu_t: samples x features counts, no all-zero
# columns. z.delete = FALSE keeps rows and columns with > 80% zeros: the prevalence
# filter, not cmultRepl, decides which features stay.
clr_czm <- function(otu_t) {
  otu_t |>
    zCompositions::cmultRepl(label = 0, method = "CZM", output = "p-counts",
                             z.delete = FALSE) |>
    vegan::decostand(method = "clr")
}
