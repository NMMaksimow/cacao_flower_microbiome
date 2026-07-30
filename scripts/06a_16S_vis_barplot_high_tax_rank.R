## ============================================================================
## 06a — 16S Phylum-level stacked barplots (CFM)
##
## Three aggregation levels (farm → tree → flower), all with the same layout:
##   Rows:    unbagged_flower (top) | bagged_flower (bottom)
##   Columns: farms ordered by management gradient; management strip at bottom
##   Fill:    Phylum, rare taxa (< RA_THRESHOLD_PHYLUM) collapsed to "Other"
##
## Key fix vs 06a_ (prev. 07a_): "Other" is now summed within each sample before any
## cross-sample averaging, so stacks correctly reach 100%.
##
## Input:  results/rds/ps_16S_bacteria_biosamples.rds
## Output: results/figures/06a_16S_barplot_phylum_*.png
## ============================================================================

library(here)
library(tidyverse)
library(phyloseq)
library(scales)
library(ggh4x)       # facet_nested() 

# ── 0. Parameters ─────────────────────────────────────────────────────────────

RA_THRESHOLD_PHYLUM <- 0.005   # Phyla with mean RA < 0.5% → "Other"

MGMT_LEVELS <- c("inside_forest", "near_forest", "agroforest", "full_sun")
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

# ── 1. Load data ──────────────────────────────────────────────────────────────

ps_16S_bacteria_biosamples <- readRDS(
        here("results", "rds", "ps_16S_bacteria_biosamples.rds")
)

cat(sprintf("Loaded: %d taxa × %d samples\n",
            ntaxa(ps_16S_bacteria_biosamples),
            nsamples(ps_16S_bacteria_biosamples)))

# ── 2. Base long table: RA at Phylum level, flower samples only ───────────────
# tax_glom() collapses ASVs to Phylum; transform_sample_counts() computes RA;
# psmelt() converts to long data.frame for ggplot.

fig_long <- ps_16S_bacteria_biosamples |>
        subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower")) |>
        tax_glom(taxrank = "Phylum", NArm = FALSE) |>
        transform_sample_counts(\(x) x / sum(x)) |>
        psmelt() |>
        mutate(
                # Replace NA taxonomy with explicit label
                Phylum          = if_else(is.na(Phylum), "p__Unassigned", Phylum),
                # unbagged first so it appears as the top facet row
                sample_type     = factor(sample_type,
                                         levels = c("unbagged_flower", "bagged_flower")),
                management_type = factor(management_type, levels = MGMT_LEVELS),
                farm_id         = factor(farm_id,         levels = FARM_LEVELS)
        )

# ── 3. Identify top Phyla ordered by descending mean RA ───────────────────────

top_phyla <- fig_long |>
        group_by(Phylum) |>
        summarise(mean_RA = mean(Abundance), .groups = "drop") |>
        filter(mean_RA >= RA_THRESHOLD_PHYLUM) |>
        arrange(desc(mean_RA)) |>   # most abundant first → drives legend order
        pull(Phylum)

cat(sprintf("Phyla shown individually (RA >= %.1f%%): %d\n",
            100 * RA_THRESHOLD_PHYLUM, length(top_phyla)))
print(top_phyla)

# Factor levels: desc RA order, p__Unassigned second-to-last, Other always last
phylum_factor_levels <- c(
        setdiff(top_phyla, "p__Unassigned"),
        if ("p__Unassigned" %in% top_phyla) "p__Unassigned",
        "Other"
)

# ── 4. Colour palette (kelly, drop white and near-black) ──────────────────────

named_phyla <- setdiff(phylum_factor_levels, "Other")
n_named     <- length(named_phyla)
base_pal    <- pals::kelly()[-c(1, 2)]   # 20 usable colours; starts with gold
phylum_pal  <- if (n_named <= length(base_pal)) {
        base_pal[seq_len(n_named)]
} else {
        colorRampPalette(base_pal)(n_named)
}
phylum_colors <- c(setNames(phylum_pal, named_phyla), "Other" = "grey80")

# ── 5. Shared theme helper ────────────────────────────────────────────────────

bar_theme <- function(x_angle = 45, x_size = 9) {
        theme_bw() +
        theme(
                strip.placement  = "outside",
                strip.background = element_blank(),
                strip.text.x     = element_text(size = 8),
                strip.text.y     = element_text(size = 9),
                panel.spacing.x  = unit(0.15, "lines"),
                legend.position  = "right",
                axis.text.x      = element_text(angle  = x_angle,
                                                hjust  = 1,
                                                vjust  = if (x_angle == 90) 0.5 else 1,
                                                size   = x_size)
        )
}

# ── Shared collapse helper ────────────────────────────────────────────────────
# Correctly sums all rare phyla into "Other" within each sample before any
# cross-sample aggregation.  Without the within-sample sum step, mean() over
# multiple rare phyla gives a tiny fraction instead of their correct total.

collapse_to_phylum_plot <- function(tbl) {
        tbl |>
                mutate(Phylum_plot = if_else(Phylum %in% top_phyla, Phylum, "Other")) |>
                # Step 1 — sum within sample: combines rare phyla into one "Other" row
                group_by(Sample, sample_type, management_type,
                         farm_id, tree_id, Phylum_plot) |>
                summarise(Abundance = sum(Abundance), .groups = "drop")
}


################################################################################
## Figure A — farm-level averages (7 bars × management type)                  ##
################################################################################
# Each bar = mean RA of all flowers in that farm.
# Two-step aggregation: sum within sample, then mean across samples in farm.

figA_tbl <- collapse_to_phylum_plot(fig_long) |>
        # Step 2 — average across all samples within farm
        group_by(sample_type, management_type, farm_id, Phylum_plot) |>
        summarise(Abundance = mean(Abundance), .groups = "drop") |>
        mutate(
                farm_id     = factor(farm_id,     levels = FARM_LEVELS),
                Phylum_plot = factor(Phylum_plot, levels = phylum_factor_levels)
        )

figA <- ggplot(figA_tbl, aes(x = farm_id, y = Abundance, fill = Phylum_plot)) +
        geom_col(width = 0.75, colour = NA) +
        ggh4x::facet_nested(
                sample_type ~ management_type,
                scales = "free_x", space = "free_x", switch = "x"
        ) +
        scale_y_continuous(
                labels = percent_format(),
                limits = c(0, 1.01),
                expand = expansion(mult = c(0, 0.01))
        ) +
        scale_fill_manual(
                values = phylum_colors,
                name   = sprintf("Phylum\n(mean RA \u2265 %.1f%%)", 100 * RA_THRESHOLD_PHYLUM)
        ) +
        bar_theme(x_angle = 45, x_size = 9) +
        labs(
                y        = "Relative abundance",
                x        = NULL,
                title    = "16S bacterial community — farm-level averages",
                subtitle = sprintf(
                        "Mean across all flowers per farm  |  Phyla < %.1f%% collapsed to Other",
                        100 * RA_THRESHOLD_PHYLUM
                )
        )

print(figA)

ggsave(
        here("results", "figures", "06a_16S_barplot_phylum_by_farm_mean.png"),
        plot = figA, width = 12, height = 8, dpi = 300
)
cat("Saved: 06a_16S_barplot_phylum_by_farm_mean.png\n")


################################################################################
## Figure B — tree-level averages (~21 bars per farm)                         ##
################################################################################
# Each bar = mean RA of all flowers in that tree.
# Two-step aggregation: sum within sample, then mean across samples in tree.

figB_tbl <- collapse_to_phylum_plot(fig_long) |>
        # Step 2 — average across flowers within tree
        group_by(sample_type, management_type, farm_id, tree_id, Phylum_plot) |>
        summarise(Abundance = mean(Abundance), .groups = "drop") |>
        mutate(
                tree_id     = factor(tree_id, levels = {
                        fig_long |>
                                distinct(farm_id, tree_id) |>
                                arrange(farm_id, tree_id) |>
                                pull(tree_id)
                }),
                Phylum_plot = factor(Phylum_plot, levels = phylum_factor_levels)
        )

figB <- ggplot(figB_tbl, aes(x = tree_id, y = Abundance, fill = Phylum_plot)) +
        geom_col(width = 0.9, colour = NA) +
        ggh4x::facet_nested(
                sample_type ~ management_type + farm_id,
                scales = "free_x", space = "free_x", switch = "x"
        ) +
        scale_x_discrete(drop = TRUE) +
        scale_y_continuous(
                labels = percent_format(),
                limits = c(0, 1.01),
                expand = expansion(mult = c(0, 0.01))
        ) +
        scale_fill_manual(
                values = phylum_colors,
                name   = sprintf("Phylum\n(mean RA \u2265 %.1f%%)", 100 * RA_THRESHOLD_PHYLUM)
        ) +
        bar_theme(x_angle = 90, x_size = 5) +
        labs(
                y        = "Relative abundance",
                x        = NULL,
                title    = "16S bacterial community — tree-level averages",
                subtitle = sprintf(
                        "Flowers averaged to tree level  |  Phyla < %.1f%% collapsed to Other",
                        100 * RA_THRESHOLD_PHYLUM
                )
        )

print(figB)

ggsave(
        here("results", "figures", "06a_16S_barplot_phylum_by_farm_tree_level.png"),
        plot = figB, width = 28, height = 8, dpi = 300
)
cat("Saved: 06a_16S_barplot_phylum_by_farm_tree_level.png\n")


################################################################################
## Figure C — flower level (individual samples)                               ##
################################################################################
# Each bar = one flower.
# X-axis: flower_id (= sample_id without leading letter — a display-only dummy
# variable that allows vertical alignment of bagged/unbagged rows across trees).
# Bars ordered: farm → tree → flower_id alphabetically within tree.

figC_tbl <- collapse_to_phylum_plot(fig_long) |>
        # flower_id strips the leading "c"/"p" prefix that encodes sample_type,
        # so the same physical flower has the same x label in both facet rows
        mutate(flower_id = str_sub(as.character(Sample), 2)) |>
        mutate(Phylum_plot = factor(Phylum_plot, levels = phylum_factor_levels))

# Determine x-axis order: farm → tree → flower_id
flower_order <- figC_tbl |>
        distinct(flower_id, farm_id, tree_id) |>
        arrange(farm_id, tree_id, flower_id) |>
        pull(flower_id) |>
        unique()

figC_tbl <- figC_tbl |>
        mutate(flower_id = factor(flower_id, levels = flower_order))

figC <- ggplot(figC_tbl, aes(x = flower_id, y = Abundance, fill = Phylum_plot)) +
        geom_col(width = 0.9, colour = NA) +
        # tree_id added to column facets so each tree cell has ~3 bars
        ggh4x::facet_nested(
                sample_type ~ management_type + farm_id + tree_id,
                scales = "free_x", space = "free_x", switch = "x"
        ) +
        scale_x_discrete(drop = TRUE) +
        scale_y_continuous(
                labels = percent_format(),
                limits = c(0, 1.01),
                expand = expansion(mult = c(0, 0.01))
        ) +
        scale_fill_manual(
                values = phylum_colors,
                name   = sprintf("Phylum\n(mean RA \u2265 %.1f%%)", 100 * RA_THRESHOLD_PHYLUM)
        ) +
        bar_theme(x_angle = 90, x_size = 4) +
        labs(
                y        = "Relative abundance",
                x        = NULL,
                title    = "16S bacterial community — flower level",
                subtitle = sprintf(
                        "Each bar = one flower; ordered by tree  |  Phyla < %.1f%% collapsed to Other",
                        100 * RA_THRESHOLD_PHYLUM
                )
        )

print(figC)

ggsave(
        here("results", "figures", "06a_16S_barplot_phylum_by_farm_flower_level.png"),
        plot = figC, width = 36, height = 8, dpi = 300
)
cat("Saved: 06a_16S_barplot_phylum_by_farm_flower_level.png\n")
