## ============================================================================
## 06d — Barplots ordered by microscopy phenotypes, 16S and ITS1
##
## One figure per phenotype (pollen grains pi_clsm, callose plug index, fungal
## colonisation score), only the flowers with a value. Three rows share the same
## flowers in the same order, so the bacteria and fungi of one flower line up:
##   1. the phenotype value of each flower (point, coloured by farm)
##   2. 16S community, Phylum (same taxa and colours as 06a)
##   3. ITS1 community, Class (same taxa and colours as 06c)
## Columns: unbagged flowers (left) and bagged flowers (right); within each,
## flowers from low to high value. Taxa shown by name are chosen on all flowers
## of the sample set (as in 06a/06c), not on the subset with a value.
## On the 294 set the labels of the flowers sampled outside the standard
## protocol are coloured by deviation type.
##
## Input:  results/rds/ps_16S_bacteria_biosamples.rds, ps_ITS1_otu97_fungi_biosamples.rds
## Output: results/figures/06d_16S_ITS1_barplot_by_{pi_clsm,callose_plug_index,fungal_colonization_score}.png
## ============================================================================

library(here)
library(tidyverse)
library(phyloseq)
library(scales)
library(patchwork)
library(ggtext)      # element_markdown(): coloured flower labels on the 294 set

# Data version and sample set: RESULTS_DIR, BIOSAMPLE_SUFFIX, SET_TAG, CFM_SAMPLE_SET
source(here("scripts", "utils_paths.R"))
source(here("scripts", "utils_barplot.R"))   # taxa_long(), top_taxa_levels(), colours, bar_theme()

# Under Rscript a plot device would write Rplots.pdf; figures go through ggsave
if (!interactive()) pdf(NULL)

# ── 0. Parameters ─────────────────────────────────────────────────────────────

RA_THRESHOLD <- 0.005     # taxa with mean RA < 0.5% -> "Other" (06a, 06b, 06c)
ITS1_RANK    <- "Class"   # "Phylum" to use the 06b rank instead of the 06c trial

VARIABLES <- c(
  pi_clsm                   = "Pollen grains (CLSM)",
  callose_plug_index        = "Callose plug index",
  fungal_colonization_score = "Fungal colonisation score"
)
LOG_VARIABLES <- c("pi_clsm", "callose_plug_index")   # skewed: log1p axis
LOG_BREAKS    <- c(0, 3, 10, 30, 100, 300)            # as in 05

FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
FARM_COLS   <- c(ib = "#1B9E77", vr = "#D95F02", sa = "#7570B3", kk = "#E7298A",
                 mt = "#66A61E", vi = "#E6AB02", yb = "#A6761D")   # Dark2, as 05 and 08a
MGMT <- c(ib = "inside forest", vr = "near forest", sa = "near forest",
          kk = "agroforest", mt = "agroforest", vi = "full sun", yb = "full sun")

# Label colours of the deviation flowers (294 set only)
DEVIATION_COLS <- c(no_buffer = "#E69F00", two_flowers = "#0072B2", insect_present = "#CC79A7")   # Okabe-Ito: orange, blue, reddish purple

TRACK_POINT_SIZE <- 1.6               # points of the phenotype row
LABEL_SIZE       <- 5                 # flower labels under the bottom row
PANEL_HEIGHTS    <- c(1.2, 3.5, 3.5)  # phenotype, 16S, ITS1 rows
FIG_WIDTH        <- 18                # inches; fixed, so the label size stays the same
FIG_HEIGHT       <- 10

LABEL_DEVIATIONS <- CFM_SAMPLE_SET == "with_deviations"

# ── 1. Data ───────────────────────────────────────────────────────────────────

ps_16S  <- readRDS(file.path(RESULTS_DIR, "rds", paste0("ps_16S_bacteria_biosamples", BIOSAMPLE_SUFFIX, ".rds")))
ps_ITS1 <- readRDS(file.path(RESULTS_DIR, "rds", paste0("ps_ITS1_otu97_fungi_biosamples", BIOSAMPLE_SUFFIX, ".rds")))

meta <- data.frame(sample_data(ps_16S)) |>
  transmute(Sample = as.character(sample_id),
            sample_type = factor(sample_type, levels = c("unbagged_flower", "bagged_flower"),
                                 labels = c("Unbagged", "Bagged")),
            farm_id = factor(farm_id, levels = FARM_LEVELS), tree_id,
            pi_clsm, callose_plug_index,
            fungal_colonization_score = as.numeric(as.character(fungal_colonization_score)),
            sampling_deviation = as.character(sampling_deviation))

# 16S: phylum, taxa and colours as 06a
long_16S   <- taxa_long(ps_16S, "Phylum")
levels_16S <- top_taxa_levels(long_16S, RA_THRESHOLD)
cols_16S   <- colours_16S_phylum(levels_16S)            # colours from the 06a order
# Stacking order: dominant phyla anchored to the top and bottom (anchor_levels())
levels_16S <- anchor_levels(levels_16S, ranking = setdiff(levels_16S, c("p__Unassigned", "Other")))
bars_16S   <- collapse_other(long_16S, levels_16S)

# ITS1: class (06c) or phylum (06b colours)
long_ITS1 <- taxa_long(ps_ITS1, ITS1_RANK)
if (ITS1_RANK == "Class") {
  shown_ITS1 <- its1_class_levels_colours(top_taxa_levels(long_ITS1, RA_THRESHOLD), long_ITS1, ps_ITS1)
} else {
  lv <- top_taxa_levels(long_ITS1, RA_THRESHOLD)
  named <- setdiff(lv, c("Other", "p__Unassigned"))
  shown_ITS1 <- list(levels = lv, colours = c(
    setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2")[-c(4, 5, 8)])(length(named)), named),
    p__Unassigned = "#E68FAC", Other = "grey80"))
}
# Stacking order: dominant taxa anchored to the top and bottom, phyla kept as blocks
ranked_ITS1 <- top_taxa_levels(long_ITS1, RA_THRESHOLD)
shown_ITS1$levels <- anchor_levels(shown_ITS1$levels,
                                   ranking  = setdiff(ranked_ITS1, c(str_subset(ranked_ITS1, "__Unassigned$"), "Other")),
                                   group_of = if (ITS1_RANK == "Class") class_to_phylum(ps_ITS1) else NULL)
bars_ITS1 <- collapse_other(long_ITS1, shown_ITS1$levels)

# Every flower must sum to 100% in both markers
walk(list(bars_16S, bars_ITS1), \(b)
  stopifnot(b |> summarise(s = sum(Abundance), .by = Sample) |> pull(s) |> (\(s) all(abs(s - 1) < 1e-9))()))

# ── 2. One figure per phenotype ───────────────────────────────────────────────

plot_variable <- function(variable) {
  # Flowers with a value, ordered within treatment from low to high value
  flowers <- meta |>
    filter(!is.na(.data[[variable]])) |>
    mutate(value = .data[[variable]]) |>
    arrange(sample_type, value, farm_id, tree_id, Sample)
  order <- flowers$Sample
  n_by  <- count(flowers, sample_type)
  strip_labels <- setNames(paste0(n_by$sample_type, " (n = ", n_by$n, ")"), n_by$sample_type)

  # Flower labels: black, or the deviation colour on the 294 set
  label_of <- flowers |>
    mutate(colour = if (LABEL_DEVIATIONS) coalesce(DEVIATION_COLS[sampling_deviation], "black") else "black",
           label  = paste0("<span style='color:", colour, "'>**", Sample, "**</span>")) |>
    select(Sample, label) |> deframe()

  facet <- facet_grid(. ~ sample_type, scales = "free_x", space = "free_x",
                      labeller = labeller(sample_type = strip_labels))
  x_scale <- scale_x_discrete(limits = \(x) order[order %in% x], labels = label_of)

  # Row 1: phenotype value of each flower
  y_value <- if (variable %in% LOG_VARIABLES) {
    scale_y_continuous(transform = "log1p", breaks = LOG_BREAKS)
  } else {
    scale_y_continuous(breaks = 0:5)
  }
  p_track <- ggplot(flowers, aes(x = Sample, y = value, colour = farm_id)) +
    geom_point(size = TRACK_POINT_SIZE) +
    facet + x_scale + y_value +
    scale_colour_manual(values = FARM_COLS, breaks = FARM_LEVELS, name = "Farm",
                        labels = paste0(FARM_LEVELS, " (", MGMT[FARM_LEVELS], ")")) +
    labs(x = NULL, y = VARIABLES[[variable]]) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = 10, face = "bold"),
          panel.grid.major.x = element_blank())

  # Rows 2 and 3: stacked bars of the same flowers
  bar_row <- function(bars, colours, legend_title, show_labels) {
    bars |>
      filter(Sample %in% order) |>
      left_join(select(flowers, Sample, sample_type), by = "Sample") |>
      ggplot(aes(x = Sample, y = Abundance, fill = taxon)) +
      geom_col(width = 0.9, colour = NA) +
      facet + x_scale +
      scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.01))) +
      scale_fill_manual(values = colours, drop = FALSE, name = legend_title) +
      labs(x = NULL, y = "Relative abundance") +
      bar_theme(x_angle = 90, x_size = LABEL_SIZE) +
      theme(strip.text.x = element_blank(),   # treatment strips only on the phenotype row
            axis.text.x = if (show_labels) element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = LABEL_SIZE)
                          else element_blank(),
            axis.ticks.x = if (show_labels) element_line() else element_blank())
  }
  p_16S  <- bar_row(bars_16S,  cols_16S,           "16S Phylum",            show_labels = FALSE)
  p_ITS1 <- bar_row(bars_ITS1, shown_ITS1$colours, paste("ITS1", ITS1_RANK), show_labels = TRUE)

  (p_track / p_16S / p_ITS1) +
    plot_layout(heights = PANEL_HEIGHTS, guides = "keep") +
    plot_annotation(
      title = sprintf("Community composition ordered by %s (%s set)", sub("^(.)", "\\L\\1", VARIABLES[[variable]], perl = TRUE), sub("_", "", SET_TAG)),
      caption = paste0("Only flowers with a value are shown. Taxa with mean RA < ", 100 * RA_THRESHOLD,
                       "% over all flowers are collapsed to Other.",
                       if (LABEL_DEVIATIONS) " Label colours: orange no buffer, blue two flowers, purple insect in the bag." else ""),
      theme = theme(plot.title = element_text(size = 12, face = "bold"),
                    plot.caption = element_text(size = 9, colour = "grey40", hjust = 0)))
}

# ── 3. Save ───────────────────────────────────────────────────────────────────

walk(names(VARIABLES), \(variable) {
  ggsave(file.path(RESULTS_DIR, "figures", paste0("06d_16S_ITS1_barplot_by_", variable, SET_TAG, ".png")),
         plot = plot_variable(variable), width = FIG_WIDTH, height = FIG_HEIGHT, dpi = 300)
  cat(sprintf("%s: %s\n", variable,
              meta |> filter(!is.na(.data[[variable]])) |> count(sample_type) |>
                mutate(s = paste(sample_type, n)) |> pull(s) |> paste(collapse = ", ")))
})
cat("06d done.\n")
