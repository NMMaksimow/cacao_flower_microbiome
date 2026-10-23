## ============================================================================
## Naive Bayes classifier comparison on the CFM ASVs or OTUs (16S or ITS1)
## Generalised from scripts/exploratory/01a_16S_compare_nb_classifiers_CFM.R
##
## Usage (from the project root, with the project's renv R):
##   Rscript workflow/scripts/plot_classifier_comparison.R <tax_prefix> \
##       <classifiers> <expected_mock.fasta> <title> <out.png>
##   <tax_prefix>   e.g. qiime2_v2/taxonomy/CFM_16S_taxonomy (files <prefix>_<name>.tsv,
##                  exported QIIME 2 taxonomy: Feature ID, Taxon, Confidence)
##   <classifiers>  comma-separated names, the first is the reference (e.g. uniform)
##
## Panels: 1. deepest rank per feature; 2. confidence distribution; 3. confidence
## bands; 4. pairwise genus agreement (confidence >= 0.7). Printed: KS tests,
## summary statistics, taxonomy of the expected mock features (IDs = md5 of the
## expected sequences in workflow/resources/mock/).
## ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

args        <- commandArgs(trailingOnly = TRUE)
tax_prefix  <- args[1]
classifiers <- strsplit(args[2], ",")[[1]]
expected_fa <- args[3]
title       <- args[4]
out_png     <- args[5]

# ── 0. Load and parse taxonomies ──────────────────────────────────────────────
# SILVA: d__;p__;c__;o__;f__;g__;s__   UNITE: k__;p__;c__;o__;f__;g__;s__
# A rank counts unless empty or uncultured/Unknown/metagenome/unidentified/Incertae_Sedis.
# SILVA 144 splits some taxa into "X" and "X--other"; both count as X.
uninformative <- "^(uncultured|Unknown|metagenome|unidentified|Incertae_Sedis)"
rank_levels <- c("unclassified", "domain", "phylum", "class", "order", "family", "genus", "species")

tax_parsed <- classifiers |>
  map(\(cl) read_tsv(sprintf("%s_%s.tsv", tax_prefix, cl), col_names = c("feature_id", "taxon", "confidence"),
                     skip = 1, show_col_types = FALSE) |> mutate(classifier = cl)) |>
  bind_rows() |>
  mutate(
    domain  = str_match(taxon, "(?:d|k)__([^;]+)")[, 2],
    phylum  = str_match(taxon, "p__([^;]+)")[, 2],
    class   = str_match(taxon, "c__([^;]+)")[, 2],
    order   = str_match(taxon, "o__([^;]+)")[, 2],
    family  = str_match(taxon, "f__([^;]+)")[, 2],
    genus   = str_match(taxon, "g__([^;]+)")[, 2],
    species = str_match(taxon, "s__([^;]+)")[, 2],
    across(domain:species, \(x) if_else(str_detect(x, uninformative), NA_character_, str_remove(x, "--other$"))),
    depth = case_when(
      !is.na(species) ~ "species", !is.na(genus) ~ "genus", !is.na(family) ~ "family",
      !is.na(order) ~ "order", !is.na(class) ~ "class", !is.na(phylum) ~ "phylum",
      !is.na(domain) ~ "domain", TRUE ~ "unclassified"
    ) |> factor(levels = rank_levels),
    classifier = factor(classifier, levels = classifiers)
  )

# ── 1. Resolution: deepest informative rank per feature ──────────────────────
p_resolution <- tax_parsed |>
  count(classifier, depth) |>
  group_by(classifier) |>
  mutate(pct = n / sum(n) * 100) |>
  ggplot(aes(x = pct, y = classifier, fill = depth)) +
  geom_col(position = "stack") +
  scale_fill_brewer(palette = "RdYlGn", direction = 1, name = "Deepest rank") +  # colour per rank
  labs(title = "Taxonomic resolution: deepest rank per feature", x = "% of features", y = NULL) +
  theme_bw(base_size = 12)

# ── 2. Confidence distribution + pairwise KS tests ───────────────────────────
combn(classifiers, 2, simplify = FALSE) |>
  map_dfr(\(pair) {
    ks <- ks.test(tax_parsed$confidence[tax_parsed$classifier == pair[1]],
                  tax_parsed$confidence[tax_parsed$classifier == pair[2]])
    tibble(pair = paste(pair, collapse = " vs "), D = unname(ks$statistic), p_value = ks$p.value)
  }) |>
  print()

p_confidence <- ggplot(tax_parsed, aes(x = classifier, y = confidence, fill = classifier)) +
  geom_violin(alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, show.legend = FALSE) +
  geom_hline(yintercept = 0.7, linetype = "dashed", colour = "grey40", linewidth = 0.4) +  # classify-sklearn cut-off
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Confidence score distribution", x = NULL, y = "Confidence") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

# ── 3. Confidence bands ───────────────────────────────────────────────────────
p_bands <- tax_parsed |>
  mutate(band = cut(confidence, c(-Inf, 0.7, 0.85, Inf), right = FALSE,
                    labels = c("low (<0.70)", "medium (0.70-0.85)", "high (>=0.85)"))) |>
  count(classifier, band) |>
  group_by(classifier) |>
  mutate(pct = n / sum(n) * 100) |>
  ggplot(aes(x = pct, y = classifier, fill = band)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("low (<0.70)" = "#d73027", "medium (0.70-0.85)" = "#fee090",
                               "high (>=0.85)" = "#1a9850"), name = "Confidence band") +
  labs(title = "Confidence bands per classifier", x = "% of features", y = NULL) +
  theme_bw(base_size = 12)

# ── 4. Pairwise genus agreement (features with a genus in both, conf >= 0.7) ──
genus_wide <- tax_parsed |>
  filter(confidence >= 0.7) |>
  select(feature_id, classifier, genus) |>
  pivot_wider(names_from = classifier, values_from = genus)

# % identical genus between two classifiers, over features with a genus in both
pair_agreement <- function(a, b) {
  both <- genus_wide[!is.na(genus_wide[[a]]) & !is.na(genus_wide[[b]]), ]
  tibble(agreement = mean(both[[a]] == both[[b]]) * 100, n = nrow(both))
}

p_agreement <- expand_grid(cl1 = classifiers, cl2 = classifiers) |>
  filter(match(cl1, classifiers) <= match(cl2, classifiers)) |>
  mutate(stats = map2(cl1, cl2, pair_agreement)) |>
  unnest(stats) |>
  mutate(cl1 = factor(cl1, levels = classifiers), cl2 = factor(cl2, levels = classifiers)) |>
  ggplot(aes(x = cl1, y = cl2, fill = agreement)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", agreement, scales::comma(n))), size = 3) +
  scale_fill_gradient(low = "white", high = "#2166ac", name = "% agreement") +  # tile colour scale
  coord_fixed() +
  labs(title = "Pairwise genus agreement (confidence >= 0.7)", x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())

# ── Printed summaries ─────────────────────────────────────────────────────────
tax_parsed |>
  group_by(classifier) |>
  summarise(n = n(), pct_genus = mean(depth %in% c("genus", "species")) * 100,
            pct_species = mean(depth == "species") * 100,
            pct_unclassified = mean(depth == "unclassified") * 100,
            median_conf = median(confidence), .groups = "drop") |>
  print()

# Expected mock features: their IDs are the md5 of the expected sequences
expected <- read_lines(expected_fa) |>
  (\(x) tibble(name = sub("^>", "", x[c(TRUE, FALSE)]), seq = x[c(FALSE, TRUE)]))() |>
  mutate(feature_id = map_chr(seq, \(s) digest::digest(s, algo = "md5", serialize = FALSE)))
tax_parsed |>
  inner_join(expected |> select(feature_id, name), by = "feature_id") |>
  transmute(name = sub("\\|copies.*", "", name), classifier,
            assigned = coalesce(species, genus, family, order, class, phylum, domain, "unclassified"),
            confidence = round(confidence, 3)) |>
  pivot_wider(names_from = classifier, values_from = c(assigned, confidence)) |>
  print(n = 50, width = Inf)

# ── Combine and save ──────────────────────────────────────────────────────────
p_combined <- ((p_resolution / p_bands / p_confidence) | p_agreement) +
  plot_annotation(title = title, theme = theme(plot.title = element_text(size = 13, face = "bold")))

ggsave(out_png, plot = p_combined, width = 18, height = 13, dpi = 300)   # size in inches
cat("Saved:", out_png, "\n")
