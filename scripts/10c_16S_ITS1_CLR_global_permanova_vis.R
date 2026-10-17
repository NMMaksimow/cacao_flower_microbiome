# ---- 0. Parameters -------------------------------------------------------
N_PERM <- 999     # permutations used in 10a/10b (caption only; no tests are run here)
ALPHA  <- 0.05    # significance threshold for "*"
DIST   <- "aitchison_CLR"   # primary beta diversity metric (CLR + CZM Aitchison)

# Marker key → RDS prefix (08a/08b) and table prefix (10a/10b)
MARKERS      <- c("16S" = "16S", "ITS1" = "ITS1_otu97")
TABLE_PREFIX <- c("16S" = "10a_16S", "ITS1" = "10b_ITS1")
MARKER_LABELS <- c("16S" = "16S bacteria", "ITS1" = "ITS1 fungi")

# Variance components: adonis2 term → label.
# Factor order controls the stack: the LAST level is drawn at the bottom of the bar,
# so tested terms (treatment first) sit on the axis and the thin treatment segment
# starts from zero. Legend order follows the stack top → bottom.
COMPONENTS <- c(
  "Residual"            = "Flowers (within tree x treatment)",
  "tree_id:sample_type" = "Tree x treatment",
  "tree_id"             = "Tree (within farm)",
  "farm_id"             = "Farm",
  "farm_id:sample_type" = "Farm x treatment",
  "sample_type"         = "Treatment"
)
# Tested terms in saturated colours, descriptive components in greys
COMPONENT_COLS <- c(
  "Treatment"                         = "#B2182B",
  "Farm x treatment"                  = "#7B3294",
  "Farm"                              = "#1B7837",
  "Tree (within farm)"                = "#9E9E9E",
  "Tree x treatment"                  = "#C8C8C8",
  "Flowers (within tree x treatment)" = "#EBEBEB"
)
# Segments below this share get their label outside the bar (ggrepel)
LABEL_MIN <- 0.03

# Bagged / unbagged colours match 10d_16S_ITS1_CLR_PCoA_by_farm.R
TREAT_COLS   <- c(bagged_flower = "#4C72B0", unbagged_flower = "#DD8452")
TREAT_LABELS <- c(bagged_flower = "Bagged", unbagged_flower = "Unbagged")

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(vegan)
library(ggrepel)    # non-overlapping labels for thin segments
library(patchwork)  # combine panels A and B

# ---- 2. Load CLR distances -----------------------------------------------
# Computed by 08a/08b. Only the primary metric is used.
clr <- names(MARKERS) |>
  set_names() |>
  lapply(\(mk) list(
    d    = readRDS(here("results", "rds", paste0(MARKERS[[mk]], "_dist_list.rds")))[[DIST]],
    meta = readRDS(here("results", "rds", paste0(MARKERS[[mk]], "_meta_list.rds")))[[DIST]]
  ))

# ---- 3. Load test results from 10a/10b -----------------------------------
# p-values come from the tree-level tests in 10a/10b, so the figure matches the
# manuscript table. SumOfSqs is kept for the consistency check in section 4.
tests <- names(MARKERS) |>
  lapply(\(mk) read_csv(here("results", "tables",
                             paste0(TABLE_PREFIX[[mk]], "_permanova_combined.csv")),
                        show_col_types = FALSE) |>
           filter(distance == DIST, level %in% c("tree", "tree_x_treatment")) |>
           transmute(marker = mk, term, level, SS_centroid = SumOfSqs, p = `Pr(>F)`)) |>
  bind_rows()

permdisp <- names(MARKERS) |>
  lapply(\(mk) read_csv(here("results", "tables",
                             paste0(TABLE_PREFIX[[mk]], "_permdisp.csv")),
                        show_col_types = FALSE) |>
           filter(distance == DIST) |>
           mutate(marker = mk)) |>
  bind_rows()

# ---- 4. Variance partition (descriptive) ---------------------------------
# One flower-level model per marker, without permutations (only sums of squares
# are needed): all six components share one denominator (total variation among
# the 294 flowers) and add up to 1. adonis2 p-values here are meaningless
# (permutations = 0) and are not used.

partition <- names(MARKERS) |>
  lapply(\(mk) {
    d <- clr[[mk]]$d   # adonis2 looks up the formula's left-hand side in this scope
    adonis2(d ~ farm_id + tree_id + sample_type +
              farm_id:sample_type + tree_id:sample_type,
            data = clr[[mk]]$meta, by = "terms", permutations = 0) |>
      as.data.frame() |>
      rownames_to_column("term") |>
      filter(term != "Total") |>
      transmute(marker = mk, term, df = Df, SS = SumOfSqs, proportion = R2)
  }) |>
  bind_rows() |>
  left_join(tests, by = c("marker", "term"))

# Consistency check: in the balanced design, flower-level SS equal the centroid
# SS from 10a/10b times the flowers per centroid (6 per tree, 3 per tree x treatment).
partition |>
  filter(!is.na(SS_centroid)) |>
  mutate(n_per_centroid = if_else(level == "tree", 6, 3),
         rel_diff       = abs(SS - n_per_centroid * SS_centroid) / SS) |>
  (\(x) if (any(x$rel_diff > 1e-6))
    stop("Variance partition does not match 10a/10b sums of squares; re-run 10a/10b"))()

cat("Variance partition (proportion of total variation among flowers):\n")
partition |>
  dplyr::select(marker, term, df, proportion, p) |>
  print()

# ---- 5. Dispersion -------------------------------------------------------
# Distance of each flower to the median of its farm x treatment group,
# same grouping as the global PERMDISP in 10a/10b.
disp <- names(MARKERS) |>
  lapply(\(mk) betadisper(clr[[mk]]$d,
                          droplevels(interaction(clr[[mk]]$meta$farm_id,
                                                 clr[[mk]]$meta$sample_type))) |>
           (\(bd) tibble(marker      = mk,
                         sample_type = clr[[mk]]$meta$sample_type,
                         distance    = bd$distances))()) |>
  bind_rows()

# ---- 6. Figure -----------------------------------------------------------
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

# Panel A data: segment midpoints computed by hand so that small segments can be
# labelled outside the bar. Stack starts from the bottom component (last level).
bar_data <- partition |>
  mutate(
    component = factor(COMPONENTS[term], levels = COMPONENTS),
    marker    = factor(MARKER_LABELS[marker], levels = MARKER_LABELS),
    label     = paste0(sprintf("%.1f%%", 100 * proportion),
                       if_else(!is.na(p) & p < ALPHA, " *", ""))
  ) |>
  group_by(marker) |>
  arrange(desc(component), .by_group = TRUE) |>
  mutate(ymid = cumsum(proportion) - proportion / 2) |>
  ungroup()

p_bar <- bar_data |>
  ggplot(aes(x = marker, y = proportion, fill = component)) +
  # width: bar width; colour: thin outline between segments
  geom_col(width = 0.55, colour = "white", linewidth = 0.3) +
  # labels inside segments that are large enough
  geom_text(data = filter(bar_data, proportion >= LABEL_MIN),
            aes(y = ymid, label = label), size = 3) +
  # labels for thin segments: placed to the right of the bar with a leader line;
  # nudge_x moves them right, direction = "y" lets them move only vertically
  geom_text_repel(data = filter(bar_data, proportion < LABEL_MIN),
                  aes(y = ymid, label = label),
                  size = 3, nudge_x = 0.45, direction = "y", hjust = 0,
                  segment.size = 0.3, min.segment.length = 0) +
  scale_fill_manual(values = COMPONENT_COLS, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Share of total variation among flowers") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right",
        panel.grid.major.x = element_blank())

# Panel B: distance to centroid, bagged vs unbagged; separate y scales because the
# two markers live in different Aitchison spaces
disp_labels <- permdisp |>
  transmute(marker = factor(MARKER_LABELS[marker], levels = MARKER_LABELS),
            label  = sprintf("F(%d, %d) = %.2f, p = %.3f",
                             df_group, df_resid, F_value, p_value))

p_disp <- disp |>
  mutate(marker = factor(MARKER_LABELS[marker], levels = MARKER_LABELS)) |>
  ggplot(aes(x = sample_type, y = distance, colour = sample_type)) +
  # outlier.shape = NA: outliers are already shown by the jittered points
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  # width: horizontal jitter; alpha: point transparency; size: point size
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  # PERMDISP result at the top of each panel; vjust pushes it below the edge
  geom_text(data = disp_labels, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 3) +
  facet_wrap(~ marker, scales = "free_y") +
  scale_colour_manual(values = TREAT_COLS, guide = "none") +
  scale_x_discrete(labels = TREAT_LABELS) +
  # 15% headroom so the PERMDISP label does not overlap the points
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "Distance to farm x treatment centroid") +
  theme_bw(base_size = 10)

caption_text <- paste(
  "A: shares of total variation among flowers (n = 294 per marker) from adonis2(d ~ farm + tree + treatment + farm:treatment + tree:treatment), Aitchison distance (CLR).",
  sprintf("   * p < %.2f from tree-level tests (%d permutations): farm on tree centroids, trees permuted freely between farms;", ALPHA, N_PERM),
  "   treatment and farm x treatment on tree x treatment centroids, permutations within tree (tree = unit of replication). Grey components are not tested.",
  "B: PERMDISP, distance of each flower to the median of its farm x treatment group; permutational test of treatment accounting for farm, permutations within tree.",
  sep = "\n")

# Panel widths: A narrower than B; tag_levels adds "A"/"B" labels
fig <- (p_bar + p_disp) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    tag_levels = "A",
    caption    = caption_text,
    # caption: left-aligned, small font
    theme      = theme(plot.caption = element_text(hjust = 0, size = 7))
  )

ggsave(here("results", "figures", "10c_16S_ITS1_CLR_global_permanova.png"),
       plot = fig, width = 12, height = 5.5, dpi = 300)

# ---- 7. Save outputs -----------------------------------------------------
partition |>
  mutate(component = COMPONENTS[term]) |>
  dplyr::select(marker, term, component, df, SS, proportion, p) |>
  write_csv(here("results", "tables", "10c_16S_ITS1_CLR_variance_partition.csv"))

message("Done. Figure and variance partition table saved.")
