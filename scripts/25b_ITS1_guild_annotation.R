## ============================================================================
## 25b — ITS1 functional annotation with FungalTraits (via microeco)
##
## FungalTraits (Polme et al. 2020) is a genus-level table: one row per genus
## with its primary and secondary lifestyle, plant-pathogenic capacity,
## endophytic capability, animal-biotrophic capacity, etc. The assignment is a
## lookup of the Genus name only; abundances enter only when the per-OTU
## lifestyles are summed into per-sample read shares. One primary lifestyle per
## genus, so primary-lifestyle shares of a sample sum to 100% (with "unassigned").
##
## Input is the canonical fungi-only dataset: post-decontam, biological samples,
## NO prevalence filter. The non-fungal kingdoms in the ITS1 data (Viridiplantae,
## Rhizaria, Metazoa, ...) fall outside FungalTraits, and there are no
## Stramenopila OTUs, so nothing that FungalTraits could annotate is dropped.
##
## The first question is coverage: what share of OTUs and reads gets an
## annotation. OTUs are split into three classes:
##   annotated            -- genus found in FungalTraits;
##   genus_not_in_db      -- informative genus, but absent from FungalTraits;
##   no_informative_genus -- genus unassigned or a placeholder (g__unidentified,
##                           g__Fungi_gen_Incertae_sedis).
##
## The taxonomy of the input RDS object is never modified.
##
## Input:
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output:
##   results/rds/25b_ITS1_fungaltraits.rds
##   results/tables/25b_ITS1_fungaltraits_taxon_traits.csv
##   results/tables/25b_ITS1_fungaltraits_lifestyle_summary.csv
##   results/tables/25b_ITS1_fungaltraits_coverage.csv
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

INPUT <- "ps_ITS1_otu97_fungi_biosamples.rds"

# Treatment groups in sample_type
GROUPS <- c("bagged_flower", "unbagged_flower")

# ── 1. Libraries & data ───────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(microeco)

# TAX_PLACEHOLDERS and INCERTAE_RE: names that carry no taxonomic information
source(here("scripts", "utils_tax_label.R"))

dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)

# Drop OTUs with zero reads in biosamples (they were seen only in controls)
ps <- readRDS(here("results", "rds", INPUT)) |>
        (\(x) prune_taxa(taxa_sums(x) > 0, x))()

cat(sprintf("Loaded: %d OTUs x %d samples\n", ntaxa(ps), nsamples(ps)))

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
# g__unidentified or g__Fungi_gen_Incertae_sedis.
has_informative_genus <- function(physeq) {
        bare <- as(tax_table(physeq), "matrix")[, "Genus"] |>
                sub(pattern = "^[a-z]__", replacement = "")
        !is.na(bare) & nchar(bare) > 0 &
                !bare %in% TAX_PLACEHOLDERS & !grepl(INCERTAE_RE, bare)
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

# ── 3. FungalTraits annotation ────────────────────────────────────────────────

tf <- trans_func$new(ps_to_microtable(ps))
tf$for_what <- "fungi"
tf$cal_func(fungi_database = "FungalTraits")

# Raw traits per OTU (taxonomy columns + all FungalTraits columns; NA = no match)
taxon_traits <- tf$res_spe_func_raw_FungalTraits

stopifnot(identical(rownames(taxon_traits), taxa_names(ps)))

# Primary lifestyle per OTU, "unassigned" where the genus is not in FungalTraits
lifestyle <- taxon_traits$primary_lifestyle |>
        (\(x) if_else(is.na(x) | x == "", "unassigned", x))()

# Per-sample % of reads per primary lifestyle (lifestyles x samples).
# rowsum() adds up the rows (OTUs) that share a lifestyle.
counts <- count_matrix(ps)
sample_lifestyle_perc <- rowsum(counts, group = lifestyle) |>
        sweep(MARGIN = 2, STATS = colSums(counts), FUN = "/") |>
        (\(x) 100 * x)()

stopifnot(all(abs(colSums(sample_lifestyle_perc) - 100) < 1e-6))

# ── 4. Save per-OTU traits ────────────────────────────────────────────────────

taxon_traits |>
        rownames_to_column("feature_id") |>
        write_csv(here("results", "tables", "25b_ITS1_fungaltraits_taxon_traits.csv"))

# ── 5. Coverage ───────────────────────────────────────────────────────────────

annotated         <- lifestyle != "unassigned"
informative_genus <- has_informative_genus(ps)

coverage <- coverage_table(
        counts = counts,
        feature_sets = list(
                all_otus             = rep(TRUE, ntaxa(ps)),
                informative_genus    = informative_genus,
                annotated            = annotated,
                genus_not_in_db      = informative_genus & !annotated,
                no_informative_genus = !informative_genus & !annotated
        ),
        groups = sample_data(ps)$sample_type
) |>
        mutate(dataset = "fungi", .before = 1)

# Every OTU falls in exactly one of the three classes
stopifnot(sum(coverage$n_features[coverage$measure %in%
        c("annotated", "genus_not_in_db", "no_informative_genus")]) == ntaxa(ps))

# Per primary lifestyle: number of OTUs, pooled read share, per-group medians
groups <- sample_data(ps)$sample_type
lifestyle_summary <- tibble(
        primary_lifestyle = rownames(sample_lifestyle_perc),
        n_otus            = as.vector(table(lifestyle)[rownames(sample_lifestyle_perc)]),
        pct_reads_pooled  = 100 * rowSums(rowsum(counts, group = lifestyle)) / sum(counts),
        sample_pct_median_bagged   = apply(sample_lifestyle_perc[, groups == GROUPS[1]], 1, median),
        sample_pct_median_unbagged = apply(sample_lifestyle_perc[, groups == GROUPS[2]], 1, median)
) |>
        arrange(desc(pct_reads_pooled))

write_csv(coverage,          here("results", "tables", "25b_ITS1_fungaltraits_coverage.csv"))
write_csv(lifestyle_summary, here("results", "tables", "25b_ITS1_fungaltraits_lifestyle_summary.csv"))

saveRDS(list(taxon_traits          = taxon_traits,
             taxon_functions       = tf$res_func,     # microeco 0/1 trait table
             sample_lifestyle_perc = sample_lifestyle_perc,
             coverage              = coverage,
             lifestyle_summary     = lifestyle_summary),
        here("results", "rds", "25b_ITS1_fungaltraits.rds"),
        compress = "xz")

# ── 6. Report ─────────────────────────────────────────────────────────────────

cat("\n==== FungalTraits coverage (percent) ====\n")
coverage |>
        mutate(across(where(is.double), \(x) round(x, 1))) |>
        print(width = Inf)

cat("\n==== Primary lifestyles by pooled read share ====\n")
lifestyle_summary |>
        mutate(across(where(is.double), \(x) round(x, 2))) |>
        print(n = Inf, width = Inf)

# Most abundant informative genera missing from FungalTraits: candidates for
# the curated database
cat("\n==== Top 15 genera not in FungalTraits (by pooled reads) ====\n")
tibble(Genus = as(tax_table(ps), "matrix")[, "Genus"],
       reads = rowSums(counts)) |>
        filter(informative_genus & !annotated) |>
        group_by(Genus) |>
        summarise(n_otus = n(), pct_reads_pooled = 100 * sum(reads) / sum(counts)) |>
        arrange(desc(pct_reads_pooled)) |>
        slice_head(n = 15) |>
        print()

# Spot check: a common phylloplane genus
cat("\n==== Spot check ====\n")
taxon_traits |>
        filter(Genus %in% c("g__Cladosporium", "Cladosporium")) |>
        distinct(Genus, primary_lifestyle, Secondary_lifestyle) |>
        print()
