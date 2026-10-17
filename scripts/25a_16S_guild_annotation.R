## ============================================================================
## 25a — 16S functional annotation with FAPROTAX (via microeco)
##
## FAPROTAX assigns functions by taxon NAME only: each function is a list of
## taxonomy patterns (mostly genus-level) whose cultured members all perform it.
## Abundances are not used for the assignment; they enter only when the per-ASV
## functions are summed into per-sample read shares. One ASV can carry several
## functions, so the function shares of a sample do not sum to 100%.
##
## Two inputs, both post-decontam, biological samples, NO prevalence filter:
##   bacteria     -- the canonical dataset used by every other analysis;
##   prokaryotes  -- Bacteria + Archaea (FAPROTAX covers archaeal functions such
##                   as ammonia oxidation and methanogenesis).
## Comparing the two is left for later.
##
## The first question is coverage: what share of ASVs and reads gets any
## annotation at all. It is reported at the end of the script.
##
## The taxonomy of the input RDS objects is never modified. microeco recognises
## prokaryotes by a "Kingdom" column, while SILVA calls the top rank "Domain";
## instead of renaming, for_what is set to "prok" explicitly.
##
## Input:
##   results/rds/ps_16S_bacteria_biosamples.rds
##   results/rds/ps_16S_prokaryotes_biosamples.rds   (created in 01a)
## Output:
##   results/rds/25a_16S_faprotax.rds
##   results/tables/25a_16S_faprotax_taxon_functions_<dataset>.csv
##   results/tables/25a_16S_faprotax_function_summary.csv
##   results/tables/25a_16S_faprotax_coverage.csv
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

INPUTS <- c(bacteria    = "ps_16S_bacteria_biosamples.rds",
            prokaryotes = "ps_16S_prokaryotes_biosamples.rds")

# Broad groups that almost every heterotroph matches. They inflate coverage
# without saying anything about ecology, so coverage is also reported without
# them. anaerobic_chemoheterotrophy is added by microeco, not by FAPROTAX itself.
GENERIC_FUNCTIONS <- c("chemoheterotrophy",
                       "aerobic_chemoheterotrophy",
                       "anaerobic_chemoheterotrophy")

# Treatment groups in sample_type
GROUPS <- c("bagged_flower", "unbagged_flower")

# ── 1. Libraries & data ───────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(microeco)

# Bundled FAPROTAX database; loaded here only to record its version
data("prok_func_FAPROTAX", package = "microeco")

# TAX_PLACEHOLDERS and INCERTAE_RE: names that carry no taxonomic information
source(here("scripts", "utils_tax_label.R"))

dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)

# Drop taxa with zero reads in biosamples (they were seen only in controls)
ps_list <- INPUTS |>
        map(\(f) readRDS(here("results", "rds", f))) |>
        map(\(ps) prune_taxa(taxa_sums(ps) > 0, ps))

iwalk(ps_list, \(ps, label)
      cat(sprintf("Loaded %-12s %d ASVs x %d samples\n", label, ntaxa(ps), nsamples(ps))))

# ── 2. Helper functions ───────────────────────────────────────────────────────

# Counts as taxa x samples, whichever orientation the object happens to use.
count_matrix <- function(physeq) {
        m <- as(otu_table(physeq), "matrix")
        if (taxa_are_rows(physeq) == FALSE) m <- t(m)
        m
}

# phyloseq -> microtable. Copies only; the phyloseq object is left as it is.
ps_to_microtable <- function(physeq) {
        microtable$new(
                otu_table    = as.data.frame(count_matrix(physeq)),
                tax_table    = as.data.frame(as(tax_table(physeq), "matrix")),
                sample_table = as(sample_data(physeq), "data.frame")
        )
}

# Genus is informative when it is assigned and is not a placeholder such as
# g__uncultured or g__Rhizobiales_Incertae_Sedis. FAPROTAX patterns are mostly
# genus-level, so this is the practical ceiling of the annotation.
has_informative_genus <- function(physeq) {
        bare <- as(tax_table(physeq), "matrix")[, "Genus"] |>
                sub(pattern = "^[a-z]__", replacement = "")
        !is.na(bare) & nchar(bare) > 0 &
                !bare %in% TAX_PLACEHOLDERS & !grepl(INCERTAE_RE, bare)
}

# Runs FAPROTAX; returns the per-ASV 0/1 table and per-sample % of reads.
run_faprotax <- function(physeq) {
        tf <- trans_func$new(ps_to_microtable(physeq))
        tf$for_what <- "prok"                       # SILVA top rank is "Domain"
        tf$cal_func(prok_database = "FAPROTAX")
        # abundance_weighted = TRUE: share of READS from ASVs with the function,
        # not share of ASVs; perc = TRUE: in percent; keep all-zero functions
        tf$cal_func_FR(abundance_weighted = TRUE, perc = TRUE, remove_zero = FALSE)
        list(taxon_functions      = tf$res_func,
             sample_function_perc = tf$res_func_FR)
}

# Coverage for a set of feature subsets: share of features, of pooled reads,
# and the per-sample distribution of read share (overall and per group).
coverage_table <- function(counts, feature_sets, groups) {
        sample_totals <- colSums(counts)
        imap_dfr(feature_sets, \(keep, measure) {
                per_sample <- 100 * colSums(counts[keep, , drop = FALSE]) / sample_totals
                tibble(
                        measure          = measure,
                        n_features       = sum(keep),
                        pct_features     = 100 * mean(keep),
                        pct_reads_pooled = 100 * sum(counts[keep, ]) / sum(counts),
                        sample_pct_median = median(per_sample),
                        sample_pct_q1     = unname(quantile(per_sample, 0.25)),
                        sample_pct_q3     = unname(quantile(per_sample, 0.75)),
                        sample_pct_median_bagged   = median(per_sample[groups == GROUPS[1]]),
                        sample_pct_median_unbagged = median(per_sample[groups == GROUPS[2]])
                )
        })
}

# ── 3. FAPROTAX annotation ────────────────────────────────────────────────────

faprotax <- map(ps_list, run_faprotax)

# Every input ASV has a row; read shares lie within 0-100
iwalk(faprotax, \(res, label) {
        stopifnot(identical(rownames(res$taxon_functions), taxa_names(ps_list[[label]])),
                  all(res$sample_function_perc >= 0 & res$sample_function_perc <= 100))
})

# The same ASV gets the same functions in both runs (bacteria is a subset)
shared <- rownames(faprotax$bacteria$taxon_functions)
stopifnot(all(shared %in% rownames(faprotax$prokaryotes$taxon_functions)),
          isTRUE(all.equal(faprotax$bacteria$taxon_functions,
                           faprotax$prokaryotes$taxon_functions[shared, ],
                           check.attributes = FALSE)))

# ── 4. Save per-ASV functions ─────────────────────────────────────────────────

# Long format: one row per ASV x function, annotated ASVs only
taxon_functions_long <- imap(faprotax, \(res, label)
        res$taxon_functions |>
                rownames_to_column("feature_id") |>
                pivot_longer(-feature_id, names_to = "func", values_to = "present") |>
                filter(present == 1) |>
                select(-present) |>
                left_join(as(tax_table(ps_list[[label]]), "matrix") |>
                                  as.data.frame() |>
                                  rownames_to_column("feature_id"),
                          by = "feature_id"))

iwalk(taxon_functions_long, \(df, label)
      write_csv(df, here("results", "tables",
                         sprintf("25a_16S_faprotax_taxon_functions_%s.csv", label))))

# ── 5. Coverage ───────────────────────────────────────────────────────────────

coverage <- imap_dfr(ps_list, \(ps, label) {
        fun <- faprotax[[label]]$taxon_functions
        specific <- setdiff(colnames(fun), GENERIC_FUNCTIONS)
        coverage_table(
                counts = count_matrix(ps),
                feature_sets = list(
                        all_asvs           = rep(TRUE, ntaxa(ps)),
                        informative_genus  = has_informative_genus(ps),
                        annotated_any      = rowSums(fun) > 0,
                        annotated_specific = rowSums(fun[, specific]) > 0
                ),
                groups = sample_data(ps)$sample_type
        ) |>
                mutate(dataset = label, .before = 1)
})

# Per function: number of ASVs and share of pooled reads
function_summary <- imap_dfr(ps_list, \(ps, label) {
        fun    <- as.matrix(faprotax[[label]]$taxon_functions)
        counts <- count_matrix(ps)
        tibble(dataset          = label,
               func             = colnames(fun),
               n_asvs           = colSums(fun),
               pct_reads_pooled = 100 * colSums(fun * rowSums(counts)) / sum(counts))
}) |>
        filter(n_asvs > 0) |>
        arrange(dataset, desc(pct_reads_pooled))

write_csv(coverage,         here("results", "tables", "25a_16S_faprotax_coverage.csv"))
write_csv(function_summary, here("results", "tables", "25a_16S_faprotax_function_summary.csv"))

saveRDS(list(faprotax         = faprotax,
             coverage         = coverage,
             function_summary = function_summary,
             faprotax_version = prok_func_FAPROTAX$ver),
        here("results", "rds", "25a_16S_faprotax.rds"),
        compress = "xz")

# ── 6. Report ─────────────────────────────────────────────────────────────────

cat("\n==== FAPROTAX coverage (percent) ====\n")
coverage |>
        mutate(across(where(is.double), \(x) round(x, 1))) |>
        print(width = Inf)

cat("\n==== Top 15 functions by pooled read share ====\n")
function_summary |>
        group_by(dataset) |>
        slice_head(n = 15) |>
        mutate(pct_reads_pooled = round(pct_reads_pooled, 2)) |>
        print(n = Inf)

# Spot checks: an insect endosymbiont, and the archaeal functions
cat("\n==== Spot checks ====\n")
taxon_functions_long$bacteria |>
        filter(Genus == "g__Wolbachia") |>
        distinct(Genus, func) |>
        print()
taxon_functions_long$prokaryotes |>
        filter(Domain == "d__Archaea") |>
        count(func, name = "n_archaeal_asvs") |>
        print()
