## ============================================================================
## 06c — ITS1 Class-level stacked barplots (trial, for comparison with 06b)
##
## Same three figures and layout as 06b (phylum), one rank deeper:
##   Rows:    unbagged_flower (top) | bagged_flower (bottom)
##   Columns: farms ordered by management gradient; management strip at bottom
##   Fill:    Class; classes with mean RA < RA_THRESHOLD go to "Other"
## At phylum rank the ITS1 bars are almost only Ascomycota and Basidiomycota;
## the class colours are shades of those two 06b colours (teal = Ascomycota,
## tan = Basidiomycota), so each 06b block reads as split into its classes.
## Class names "<taxon>_cls_Incertae_sedis" count as c__Unassigned.
##
## Input:  results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output: results/figures/06c_ITS1_barplot_class_by_farm_{mean,tree_level,flower_level}.png
## ============================================================================

library(here)
library(tidyverse)
library(phyloseq)
library(scales)
library(ggh4x)       # facet_nested()

# Data version and sample set: RESULTS_DIR, BIOSAMPLE_SUFFIX, SET_TAG
source(here("scripts", "utils_paths.R"))
source(here("scripts", "utils_barplot.R"))   # taxa_long(), top_taxa_levels(), colours, bar_theme()

# Under Rscript a plot device would write Rplots.pdf; figures go through ggsave
if (!interactive()) pdf(NULL)

# ── 0. Parameters ─────────────────────────────────────────────────────────────

RANK         <- "Class"
RA_THRESHOLD <- 0.005    # classes with mean RA < 0.5% -> "Other" (same rule as 06b)

MGMT_LEVELS <- c("inside_forest", "near_forest", "agroforest", "full_sun")
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

# ── 1. Data ───────────────────────────────────────────────────────────────────

ps <- readRDS(file.path(RESULTS_DIR, "rds", paste0("ps_ITS1_otu97_fungi_biosamples", BIOSAMPLE_SUFFIX, ".rds")))

meta <- data.frame(sample_data(ps)) |>
  transmute(Sample = as.character(sample_id),
            sample_type     = factor(sample_type, levels = c("unbagged_flower", "bagged_flower")),
            management_type = factor(management_type, levels = MGMT_LEVELS),
            farm_id         = factor(farm_id, levels = FARM_LEVELS),
            tree_id,
            sampling_deviation = as.character(sampling_deviation))

class_long <- taxa_long(ps, RANK)

# Shown classes and their colours (ordered by phylum, then by RA)
ranked  <- top_taxa_levels(class_long, RA_THRESHOLD)   # by decreasing mean RA
shown   <- its1_class_levels_colours(ranked, class_long, ps)
# Stacking order: most abundant class at the top, second most abundant at the
# bottom, each phylum kept as one block, Unassigned and Other in the middle
shown$levels <- anchor_levels(shown$levels, ranking = setdiff(ranked, c("c__Unassigned", "Other")),
                              group_of = class_to_phylum(ps))
flowers <- collapse_other(class_long, shown$levels) |>
  left_join(meta, by = "Sample")

cat("Classes shown:\n"); print(shown$levels)

# Every flower must sum to 100%
stopifnot(flowers |> summarise(s = sum(Abundance), .by = Sample) |> pull(s) |> (\(s) all(abs(s - 1) < 1e-9))())

# Shared scale and labels
fill_scale <- scale_fill_manual(
  values = shown$colours, drop = FALSE,
  name   = sprintf("Class (mean RA ≥ %.1f%%)\nteal = Ascomycota, tan = Basidiomycota", 100 * RA_THRESHOLD))
y_scale <- scale_y_continuous(labels = percent_format(), limits = c(0, 1.01),
                              expand = expansion(mult = c(0, 0.01)))
SUBTITLE_RULE <- sprintf("Classes < %.1f%% collapsed to Other", 100 * RA_THRESHOLD)

################################################################################
## Figure A — farm-level averages                                             ##
################################################################################

figA <- flowers |>
  summarise(Abundance = mean(Abundance), .by = c(sample_type, management_type, farm_id, taxon)) |>
  ggplot(aes(x = farm_id, y = Abundance, fill = taxon)) +
  geom_col(width = 0.75, colour = NA) +
  ggh4x::facet_nested(sample_type ~ management_type, scales = "free_x", space = "free_x", switch = "x") +
  y_scale + fill_scale +
  bar_theme(x_angle = 45, x_size = 9) +
  labs(y = "Relative abundance", x = NULL,
       title    = "ITS1 fungal community — farm-level averages (OTU97), Class",
       subtitle = paste("Mean across all flowers per farm  |", SUBTITLE_RULE))

ggsave(file.path(RESULTS_DIR, "figures", paste0("06c_ITS1_barplot_class_by_farm_mean", SET_TAG, ".png")),
       plot = figA, width = 12, height = 8, dpi = 300)

################################################################################
## Figure B — tree-level averages                                             ##
################################################################################

tree_order <- meta |> distinct(farm_id, tree_id) |> arrange(farm_id, tree_id) |> pull(tree_id)

figB <- flowers |>
  summarise(Abundance = mean(Abundance), .by = c(sample_type, management_type, farm_id, tree_id, taxon)) |>
  mutate(tree_id = factor(tree_id, levels = tree_order)) |>
  ggplot(aes(x = tree_id, y = Abundance, fill = taxon)) +
  geom_col(width = 0.9, colour = NA) +
  ggh4x::facet_nested(sample_type ~ management_type + farm_id, scales = "free_x", space = "free_x", switch = "x") +
  scale_x_discrete(drop = TRUE) +
  y_scale + fill_scale +
  bar_theme(x_angle = 90, x_size = 5) +
  labs(y = "Relative abundance", x = NULL,
       title    = "ITS1 fungal community — tree-level averages (OTU97), Class",
       subtitle = paste("Flowers averaged to tree level  |", SUBTITLE_RULE))

ggsave(file.path(RESULTS_DIR, "figures", paste0("06c_ITS1_barplot_class_by_farm_tree_level", SET_TAG, ".png")),
       plot = figB, width = 28, height = 8, dpi = 300)

################################################################################
## Figure C — flower level                                                    ##
################################################################################
# flower_id drops the c/p prefix that encodes sample_type, so the same position
# on a tree has the same x label in both rows (as 06b)

flowerC <- flowers |> mutate(flower_id = str_sub(Sample, 2))
flower_order <- flowerC |> distinct(flower_id, farm_id, tree_id) |>
  arrange(farm_id, tree_id, flower_id) |> pull(flower_id) |> unique()

figC <- flowerC |>
  mutate(flower_id = factor(flower_id, levels = flower_order)) |>
  ggplot(aes(x = flower_id, y = Abundance, fill = taxon)) +
  geom_col(width = 0.9, colour = NA) +
  ggh4x::facet_nested(sample_type ~ management_type + farm_id + tree_id,
                      scales = "free_x", space = "free_x", switch = "x") +
  scale_x_discrete(drop = TRUE) +
  y_scale + fill_scale +
  bar_theme(x_angle = 90, x_size = 4) +
  labs(y = "Relative abundance", x = NULL,
       title    = "ITS1 fungal community — flower level (OTU97), Class",
       subtitle = paste("Each bar = one flower; ordered by tree  |", SUBTITLE_RULE))

# 294 set: sampling deviations marked by label colour and a triangle above the bar
if (CFM_SAMPLE_SET == "with_deviations") {
  deviations <- meta |>
    filter(sampling_deviation != "none") |>
    mutate(flower_id = factor(str_sub(Sample, 2), levels = flower_order))
  figC <- mark_deviations(figC, deviations, label_size = 4)
}

ggsave(file.path(RESULTS_DIR, "figures", paste0("06c_ITS1_barplot_class_by_farm_flower_level", SET_TAG, ".png")),
       plot = figC, width = 36, height = 8, dpi = 300)

cat("06c done.\n")
