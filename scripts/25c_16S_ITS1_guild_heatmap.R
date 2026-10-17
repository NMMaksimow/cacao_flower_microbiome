## ============================================================================
## 25c — Preliminary guild figures: what is annotated, and where it sits
##
## A first look at the UNCURATED guild labels from 25a (FAPROTAX, 16S bacteria)
## and 25b (FungalTraits, ITS1 fungi), before any manual checking. Descriptive
## only; no statistics.
##
## Heatmaps (one figure per marker, two panels sharing the rows):
##   A "how much"            -- log10 of the median % of reads in each farm x
##                              treatment cell (21 flowers);
##   B "where more or less"  -- log2(cell median / median over all flowers of
##                              that guild). Guilds differ by 2-3 orders of
##                              magnitude, so on panel A alone colour mostly
##                              shows the row rank; panel B removes the row level
##                              and leaves the farm and treatment pattern.
## The cell value is the median over ALL flowers, not over the flowers where the
## guild was detected (as in 00b): a guild sums many taxa and is rarely zero.
##
## Each row label names the genus that contributes most of the guild's reads.
## Uncurated labels are easy to misread without it -- animal_parasite in ITS1 is
## mostly Hannaella, a phylloplane yeast -- and it doubles as a list of what to
## curate first.
##
## 16S: the three generic chemoheterotrophy groups are left out of the heatmap
## (they would fill the top rows and say nothing); their share is shown in the
## coverage figure. FAPROTAX functions are nested (methanol_oxidation inside
## methylotrophy, human_gut = mammal_gut), so some of the 20 rows repeat each other.
## ITS1: only 19 primary lifestyles occur, so all are shown, with "unassigned"
## as the last row.
##
## Input:
##   results/rds/25a_16S_faprotax.rds
##   results/rds/25b_ITS1_fungaltraits.rds
##   results/rds/ps_16S_bacteria_biosamples.rds
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output:
##   results/figures/25c_ITS1_guild_heatmap.png
##   results/figures/25c_16S_guild_heatmap.png
##   results/figures/25c_16S_ITS1_guild_coverage.png
##   results/figures/25c_ITS1_lifestyle_barplots.png
##   results/tables/25c_guild_dominant_genus.csv
##
## The lifestyle barplots show composition, which the medians of the heatmap
## cannot: medians of the parts do not add up to 100%. ITS1 only -- each OTU has
## one primary lifestyle, while FAPROTAX functions overlap and cannot be stacked.
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

N_TOP <- 20   # 16S: number of FAPROTAX functions shown (by pooled read share)

# Same generic groups as in 25a
GENERIC_FUNCTIONS <- c("chemoheterotrophy",
                       "aerobic_chemoheterotrophy",
                       "anaerobic_chemoheterotrophy")

# Canonical order, as in 00_curated_taxa_check.R
MGMT_LEVELS <- c("inside_forest", "near_forest", "agroforest", "full_sun")
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

# Panel B colour range in log2 units: +-3 = 8-fold above / below the guild median.
# Values further out are squished to the end colours.
LOG2_CAP <- 3

# ── 1. Libraries & data ───────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(ggh4x)       # facet_nested(): management strip above the farm strip
library(patchwork)   # panels A and B side by side

# TAX_PLACEHOLDERS and INCERTAE_RE: names that carry no taxonomic information
source(here("scripts", "utils_tax_label.R"))

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

# Zero-read taxa are dropped exactly as in 25a / 25b, so taxon order matches
ps_16S  <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds")) |>
        (\(x) prune_taxa(taxa_sums(x) > 0, x))()
ps_ITS1 <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds")) |>
        (\(x) prune_taxa(taxa_sums(x) > 0, x))()

fap <- readRDS(here("results", "rds", "25a_16S_faprotax.rds"))$faprotax$bacteria
ft  <- readRDS(here("results", "rds", "25b_ITS1_fungaltraits.rds"))

stopifnot(identical(rownames(fap$taxon_functions), taxa_names(ps_16S)),
          identical(rownames(ft$taxon_traits),    taxa_names(ps_ITS1)))

# ── 2. Helper functions ───────────────────────────────────────────────────────

# Counts as taxa x samples, whichever orientation the object happens to use.
count_matrix <- function(physeq) {
        m <- as(otu_table(physeq), "matrix")
        if (taxa_are_rows(physeq) == FALSE) m <- t(m)
        m
}

# Genus without its "g__" prefix
genus_bare <- function(physeq) {
        as(tax_table(physeq), "matrix")[, "Genus"] |>
                sub(pattern = "^[a-z]__", replacement = "")
}

# Same rule as in 25a / 25b
has_informative_genus <- function(physeq) {
        bare <- genus_bare(physeq)
        !is.na(bare) & nchar(bare) > 0 &
                !bare %in% TAX_PLACEHOLDERS & !grepl(INCERTAE_RE, bare)
}

# Sample metadata needed for the facets, with canonical factor levels
sample_meta <- function(physeq) {
        as(sample_data(physeq), "data.frame") |>
                as_tibble(rownames = "sample") |>
                select(sample, farm_id, management_type, sample_type, tree_id) |>
                mutate(tree_no         = readr::parse_number(as.character(tree_id)),  # "ib3" -> 3
                       farm_id         = factor(farm_id, levels = FARM_LEVELS),
                       management_type = factor(management_type, levels = MGMT_LEVELS),
                       treatment       = factor(sample_type,
                                                levels = c("bagged_flower", "unbagged_flower"),
                                                labels = c("bgd", "unbgd")))
}

# For each guild (column of a taxa x guild 0/1 matrix): number of taxa, pooled
# reads, and the genus contributing most of the guild's reads.
dominant_genus <- function(assign, counts, genus, unit) {
        reads <- rowSums(counts)
        genus <- replace_na(genus, "no genus")
        map_dfr(colnames(assign), \(g) {
                in_g     <- assign[, g] == 1
                by_genus <- tapply(reads[in_g], genus[in_g], sum)
                tibble(guild            = g,
                       n_features       = sum(in_g),
                       pct_reads_pooled = 100 * sum(reads[in_g]) / sum(counts),
                       top_genus        = names(which.max(by_genus)),
                       top_genus_pct    = 100 * max(by_genus) / sum(by_genus),
                       label = sprintf("%s (%d %s) · %s %.0f%%", g, sum(in_g), unit,
                                       names(which.max(by_genus)),
                                       100 * max(by_genus) / sum(by_genus)))
        })
}

# Median % reads per guild x farm x treatment, and the two fill values
cell_summary <- function(long) {
        row_median <- long |>
                group_by(guild) |>
                summarise(row_med = median(pct), .groups = "drop")
        long |>
                group_by(guild, label, management_type, farm_id, treatment) |>
                summarise(n_flowers = n(), med = median(pct), .groups = "drop") |>
                left_join(row_median, by = "guild") |>
                mutate(
                        # grey (NA) where the median is 0: log undefined
                        log10_med = if_else(med > 0, log10(med), NA_real_),
                        log2_rel  = if_else(med > 0 & row_med > 0,
                                            log2(med / row_med), NA_real_)
                )
}

# One heatmap panel; `fill_scale` decides what the colour means
heat_panel <- function(cells, fill_var, fill_scale, panel_title) {
        ggplot(cells, aes(x = treatment, y = label, fill = .data[[fill_var]])) +
                geom_tile(colour = "white", linewidth = 0.4) +
                # management strip above the farm strip; two cells per farm side by side
                ggh4x::facet_nested(~ management_type + farm_id) +
                fill_scale +
                labs(title = panel_title, x = NULL, y = NULL) +
                theme_bw(base_size = 10) +
                theme(
                        strip.background = element_rect(fill = "grey92"),
                        strip.text.x     = element_text(size = 7.5),
                        panel.grid       = element_blank(),
                        panel.spacing.x  = unit(0.25, "lines"),   # gap between farms
                        axis.text.x      = element_text(size = 7),
                        axis.text.y      = element_text(size = 8),
                        legend.position  = "bottom",              # under each panel
                        legend.key.width = unit(1.5, "cm"),
                        plot.title       = element_text(size = 10, face = "bold")
                )
}

# Fill scales shared by both markers
scale_level <- scale_fill_viridis_c(
        option    = "magma",
        direction = -1,                 # high share = dark, as in 00b
        na.value  = "grey90",           # median 0
        name      = expression(log[10]~"median % reads")
)
scale_relative <- scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",   # blue = less, red = more
        midpoint = 0,
        limits   = c(-LOG2_CAP, LOG2_CAP),
        oob      = scales::squish,      # beyond the cap -> end colours
        na.value = "grey90",
        name     = expression(log[2]~"(cell / guild median)")
)

# Two panels sharing the rows; panel B drops the repeated row labels
two_panel_heatmap <- function(cells, title, subtitle) {
        p_a <- heat_panel(cells, "log10_med", scale_level, "A  how much")
        p_b <- heat_panel(cells, "log2_rel", scale_relative, "B  where more or less") +
                theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        (p_a | p_b) +
                plot_annotation(title = title, subtitle = subtitle,
                                theme = theme(plot.title = element_text(size = 13)))
}

# ── 3. ITS1: per-sample shares and row labels ─────────────────────────────────

counts_ITS1 <- count_matrix(ps_ITS1)
meta_ITS1   <- sample_meta(ps_ITS1)

lifestyle <- ft$taxon_traits$primary_lifestyle |>
        (\(x) if_else(is.na(x) | x == "", "unassigned", x))()

# OTU x lifestyle 0/1 matrix (one lifestyle per OTU)
assign_ITS1 <- sapply(sort(unique(lifestyle)), \(l) as.numeric(lifestyle == l))

dom_ITS1 <- dominant_genus(assign_ITS1, counts_ITS1, genus_bare(ps_ITS1), "OTUs")

# Row order: highest pooled share at the top, "unassigned" always last.
# ggplot draws the first factor level at the bottom.
lab_order_ITS1 <- dom_ITS1 |>
        arrange(guild != "unassigned", pct_reads_pooled) |>
        pull(label)

long_ITS1 <- ft$sample_lifestyle_perc |>
        as.data.frame() |>
        rownames_to_column("guild") |>
        pivot_longer(-guild, names_to = "sample", values_to = "pct") |>
        left_join(meta_ITS1, by = "sample") |>
        left_join(dom_ITS1 |> select(guild, label), by = "guild") |>
        mutate(label = factor(label, levels = lab_order_ITS1))

# Cross-check against 25b: pooled bagged medians per lifestyle must agree
stopifnot(isTRUE(all.equal(
        long_ITS1 |> filter(treatment == "bgd") |> group_by(guild) |>
                summarise(m = median(pct)) |> arrange(guild) |> pull(m),
        ft$lifestyle_summary |> arrange(primary_lifestyle) |>
                pull(sample_pct_median_bagged),
        check.attributes = FALSE)))   # 25b stored the medians as a named vector

# ── 4. 16S: per-sample shares and row labels ──────────────────────────────────

counts_16S <- count_matrix(ps_16S)
meta_16S   <- sample_meta(ps_16S)
fun_16S    <- as.matrix(fap$taxon_functions)

dom_16S_all <- dominant_genus(fun_16S[, colSums(fun_16S) > 0],
                              counts_16S, genus_bare(ps_16S), "ASVs")

top_16S <- dom_16S_all |>
        filter(!guild %in% GENERIC_FUNCTIONS) |>
        slice_max(pct_reads_pooled, n = N_TOP, with_ties = FALSE)

long_16S <- fap$sample_function_perc |>
        rownames_to_column("sample") |>
        pivot_longer(-sample, names_to = "guild", values_to = "pct") |>
        filter(guild %in% top_16S$guild) |>
        left_join(meta_16S, by = "sample") |>
        left_join(top_16S |> select(guild, label), by = "guild") |>
        mutate(label = factor(label, levels = rev(top_16S$label)))  # top row = largest

# ── 5. Heatmaps ───────────────────────────────────────────────────────────────

cells_ITS1 <- cell_summary(long_ITS1)
cells_16S  <- cell_summary(long_16S)

cat("Flowers per farm x treatment cell:\n")
cells_ITS1 |> distinct(farm_id, treatment, n_flowers) |>
        pivot_wider(names_from = treatment, values_from = n_flowers) |> print()

p_heat_ITS1 <- two_panel_heatmap(
        cells_ITS1,
        title    = "ITS1 fungi — FungalTraits primary lifestyle (UNCURATED)",
        subtitle = paste0("median % of reads per farm x treatment (21 flowers); ",
                          "label: lifestyle (OTUs) · genus with most of its reads; ",
                          "grey = median 0"))

p_heat_16S <- two_panel_heatmap(
        cells_16S,
        title    = sprintf("16S bacteria — top %d FAPROTAX functions (UNCURATED)", N_TOP),
        subtitle = paste0("median % of reads per farm x treatment (21 flowers); ",
                          "generic chemoheterotrophy groups omitted; functions are nested, ",
                          "so some rows overlap; grey = median 0"))

# Height grows with the number of rows; width fits the long row labels
ggsave(here("results", "figures", "25c_ITS1_guild_heatmap.png"),
       plot = p_heat_ITS1, width = 18, height = 2.5 + 0.3 * nlevels(long_ITS1$label),
       dpi = 300)
ggsave(here("results", "figures", "25c_16S_guild_heatmap.png"),
       plot = p_heat_16S, width = 18, height = 2.5 + 0.3 * nlevels(long_16S$label),
       dpi = 300)

# ── 6. Coverage figure ────────────────────────────────────────────────────────

# Annotation classes shared by both markers, so one legend serves both rows.
# "generic only" exists for 16S alone (FungalTraits has no such catch-all).
COVERAGE_COLS <- c(
        "annotated"                        = "#1B7837",
        "generic chemoheterotrophy only"   = "#A6DBA0",
        "informative genus, not annotated" = "#9E9E9E",
        "no informative genus"             = "#E0E0E0"
)

class_ITS1 <- case_when(
        lifestyle != "unassigned"        ~ "annotated",
        has_informative_genus(ps_ITS1)   ~ "informative genus, not annotated",
        TRUE                             ~ "no informative genus")

specific_cols <- setdiff(colnames(fun_16S), GENERIC_FUNCTIONS)
class_16S <- case_when(
        rowSums(fun_16S[, specific_cols]) > 0 ~ "annotated",
        rowSums(fun_16S) > 0                  ~ "generic chemoheterotrophy only",
        has_informative_genus(ps_16S)         ~ "informative genus, not annotated",
        TRUE                                  ~ "no informative genus")

# Per-sample % of reads in each class: rowsum() adds up the taxa of one class
coverage_long <- function(counts, classes, meta, marker) {
        rowsum(counts, group = classes) |>
                sweep(MARGIN = 2, STATS = colSums(counts), FUN = "/") |>
                (\(x) 100 * x)() |>
                as.data.frame() |>
                rownames_to_column("class") |>
                pivot_longer(-class, names_to = "sample", values_to = "pct") |>
                left_join(meta, by = "sample") |>
                mutate(marker = marker)
}

cov_long <- bind_rows(
        coverage_long(counts_ITS1, class_ITS1, meta_ITS1, "ITS1"),
        coverage_long(counts_16S,  class_16S,  meta_16S,  "16S")
) |>
        mutate(class  = factor(class, levels = rev(names(COVERAGE_COLS))),
               marker = factor(marker, levels = c("ITS1", "16S")))

stopifnot(cov_long |> group_by(marker, sample) |>
                  summarise(s = sum(pct), .groups = "drop") |>
                  pull(s) |> (\(s) all(abs(s - 100) < 1e-6))())

# One plot per marker, stacked with patchwork. In a single faceted plot the two
# markers would share the x axis of each farm column, so every panel would carry
# the empty slots of the other marker's bars, and the bars could not be sorted
# separately per marker.
coverage_panel <- function(df, marker_title) {
        # bars sorted by annotated share within each panel, so the profile reads as a curve
        sample_order <- df |>
                filter(class == "annotated") |>
                arrange(desc(pct)) |>
                pull(sample)

        df |>
                mutate(sample = factor(sample, levels = sample_order)) |>
                ggplot(aes(x = sample, y = pct, fill = class)) +
                geom_col(width = 1) +                      # one bar per flower, no gaps
                ggh4x::facet_nested(~ management_type + farm_id + treatment,
                                    scales = "free_x", space = "free_x") +
                scale_fill_manual(values = COVERAGE_COLS, breaks = names(COVERAGE_COLS),
                                  name = NULL, drop = FALSE) +
                scale_y_continuous(expand = c(0, 0), breaks = c(0, 50, 100)) +
                labs(title = marker_title, x = NULL, y = "% of reads") +
                theme_bw(base_size = 10) +
                theme(
                        strip.background = element_rect(fill = "grey92"),
                        strip.text.x     = element_text(size = 7),
                        panel.grid       = element_blank(),
                        panel.spacing.x  = unit(0.15, "lines"),
                        axis.text.x      = element_blank(),  # 294 sample ids are unreadable
                        axis.ticks.x     = element_blank(),
                        plot.title       = element_text(size = 10, face = "bold")
                )
}

# The legend comes from the 16S panel only: it holds all four classes, while the
# ITS1 panel has no "generic" class and its legend would not merge with the other.
p_cov <- ((coverage_panel(filter(cov_long, marker == "ITS1"),
                          "ITS1 fungi — FungalTraits primary lifestyle") +
                   guides(fill = "none")) /
          coverage_panel(filter(cov_long, marker == "16S"),
                         "16S bacteria — FAPROTAX (annotated = at least one non-generic function)")) +
        plot_layout(guides = "collect") +             # one shared legend
        plot_annotation(
                title    = "Guild annotation coverage per flower (UNCURATED)",
                subtitle = "one bar per flower, sorted by annotated share within each panel",
                theme    = theme(plot.title = element_text(size = 13),
                                 legend.position = "bottom"))

ggsave(here("results", "figures", "25c_16S_ITS1_guild_coverage.png"),
       plot = p_cov, width = 18, height = 7, dpi = 300)

# ── 7. ITS1 primary lifestyle stacked bars: flower, tree, farm ───────────────

# Lifestyles with at least this pooled % of reads get their own colour; the
# smaller annotated ones are merged into "other annotated".
MIN_PCT_OWN_COLOUR <- 0.5

own_colour <- dom_ITS1 |>
        filter(guild != "unassigned", pct_reads_pooled >= MIN_PCT_OWN_COLOUR) |>
        arrange(desc(pct_reads_pooled)) |>
        pull(guild)

# Stack order = legend order: largest lifestyle first (bottom of the bar),
# "other annotated" and "unassigned" last (top of the bar)
LIFE_LEVELS <- c(own_colour, "other annotated", "unassigned")

# Okabe-Ito: colour-blind safe, 9 colours; the 9th (grey) is dropped so that grey
# is reserved for "other annotated" and "unassigned" (same light grey as
# "no informative genus" in the coverage figure).
LIFE_COLS <- c(
        setNames(unname(palette.colors(palette = "Okabe-Ito"))[c(2:8, 1)][seq_along(own_colour)],
                 own_colour),
        "other annotated" = "#9E9E9E",
        "unassigned"      = "#E0E0E0"
)

# Per-flower % in each group (sum of the lifestyles that share a group)
life_flower <- long_ITS1 |>
        mutate(group = case_when(guild %in% own_colour ~ guild,
                                 guild == "unassigned" ~ "unassigned",
                                 TRUE                  ~ "other annotated")) |>
        group_by(sample, farm_id, management_type, treatment, tree_no, group) |>
        summarise(pct = sum(pct), .groups = "drop")

# Flower position: trees 1..7 left to right, 3 flowers each, one empty slot
# between trees -> x = (tree_no - 1) * 4 + flower_no
flower_x <- life_flower |>
        distinct(sample, treatment, tree_no) |>
        arrange(treatment, tree_no, sample) |>
        group_by(treatment, tree_no) |>
        mutate(flower_no = row_number()) |>
        ungroup() |>
        mutate(x = (tree_no - 1) * 4 + flower_no)

life_flower <- life_flower |> left_join(select(flower_x, sample, x), by = "sample")

# Tree x treatment and farm x treatment: MEAN of the flowers' %. Means of the
# parts add up to 100% (medians do not), and every flower weighs the same
# (pooled reads would let deep libraries dominate).
life_tree <- life_flower |>
        group_by(farm_id, management_type, treatment, tree_no, group) |>
        summarise(n_flowers = n(), pct = mean(pct), .groups = "drop") |>
        mutate(x = tree_no)

life_farm <- life_flower |>
        group_by(farm_id, management_type, treatment, group) |>
        summarise(n_flowers = n(), pct = mean(pct), .groups = "drop") |>
        mutate(x = 1)

# Checks: 3 flowers per tree x treatment, 21 per farm x treatment; bars sum to 100
stopifnot(all(life_tree$n_flowers == 3), all(life_farm$n_flowers == 21))
walk(list(flower = life_flower |> group_by(sample),
          tree   = life_tree   |> group_by(farm_id, treatment, tree_no),
          farm   = life_farm   |> group_by(farm_id, treatment)),
     \(df) stopifnot(df |> summarise(s = sum(pct), .groups = "drop") |>
                             pull(s) |> (\(s) all(abs(s - 100) < 1e-6))()))

# One row of the figure; same strips as the coverage figure
lifestyle_panel <- function(df, row_title) {
        df |>
                mutate(group = factor(group, levels = rev(LIFE_LEVELS))) |>
                ggplot(aes(x = x, y = pct, fill = group)) +
                geom_col(width = 1) +                    # bars touch inside a tree
                ggh4x::facet_nested(~ management_type + farm_id + treatment,
                                    scales = "free_x", space = "free_x") +
                scale_fill_manual(values = LIFE_COLS, breaks = LIFE_LEVELS,
                                  name = NULL, drop = FALSE,
                                  guide = guide_legend(nrow = 2)) +  # wide, short legend

                scale_y_continuous(expand = c(0, 0), breaks = c(0, 50, 100)) +
                labs(title = row_title, x = NULL, y = "% of reads") +
                theme_bw(base_size = 10) +
                theme(
                        strip.background = element_rect(fill = "grey92"),
                        strip.text.x     = element_text(size = 7),
                        panel.grid       = element_blank(),
                        panel.spacing.x  = unit(0.15, "lines"),
                        axis.text.x      = element_blank(),
                        axis.ticks.x     = element_blank(),
                        plot.title       = element_text(size = 10, face = "bold")
                )
}

p_life <- (lifestyle_panel(life_flower, "flower (3 per tree, grouped by tree 1-7)") /
           lifestyle_panel(life_tree,   "tree (mean of 3 flowers)") /
           lifestyle_panel(life_farm,   "farm (mean of 21 flowers)")) +
        plot_layout(guides = "collect") +            # one legend for the three rows
        plot_annotation(
                title    = "ITS1 fungi — FungalTraits primary lifestyle (UNCURATED)",
                subtitle = sprintf(paste0(
                        "%% of reads; tree and farm bars are means of the flowers' %% ",
                        "(each flower weighs the same); lifestyles below %.1f%% of pooled ",
                        "reads merged into 'other annotated'"), MIN_PCT_OWN_COLOUR),
                theme    = theme(plot.title = element_text(size = 13))) &
        # `&` applies the theme to every row, so each legend is built horizontal
        # (2 rows, see guide_legend above) before patchwork collects them
        theme(legend.position = "bottom")

ggsave(here("results", "figures", "25c_ITS1_lifestyle_barplots.png"),
       plot = p_life, width = 18, height = 10, dpi = 300)

# ── 8. Dominant genus table ───────────────────────────────────────────────────

bind_rows(dom_ITS1    |> mutate(marker = "ITS1", .before = 1),
          dom_16S_all |> mutate(marker = "16S",  .before = 1)) |>
        select(-label) |>
        arrange(marker, desc(pct_reads_pooled)) |>
        write_csv(here("results", "tables", "25c_guild_dominant_genus.csv"))

cat("\n25c done: 3 figures, 1 table\n")
