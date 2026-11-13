## ============================================================================
## 05 — Rarefaction curves, 16S and ITS1, coloured by metadata
##
## Does each flower reach a plateau of observed taxa, and do the curves differ
## by treatment, farm or microscopy phenotype? Replaces 05a, 05b and 05c.
##
## Layout: 5 rows x 2 columns. Left column 16S (bacteria, ASVs), right column
## ITS1 (fungi, OTU97, the unit of every later ITS1 step). Rows: treatment, farm,
## pollen grains (pi_clsm), callose plug index, fungal colonisation score.
## One thin curve per flower; flowers without a value are left out of that row.
## Dashed lines: the rarefaction depths of the alpha-diversity steps 07a/07b.
## On the 294 set (CFM_SAMPLE_SET=with_deviations) the six flowers sampled
## outside the standard protocol are labelled at the end of their curve with the
## deviation type (insect_id for flowers with an insect inside the bag).
##
## Curves: vegan::rarecurve(), the analytic expectation of rarefy(), so the
## result is deterministic. About 15 s for both markers; nothing is cached.
## No Good's coverage table: in v2 no flower has a singleton taxon, so it is 1
## for every flower (same reason Chao1 is not reported).
##
## Input (01a, 01b):
##   results/rds/ps_16S_bacteria_biosamples.rds          (_with_deviations for 294)
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds      (_with_deviations for 294)
## Output:
##   results/figures/05_rarefaction_curves_288.png      (_294 for the 294 set)
##   results/figures/05_rarefaction_by_farm_288.png     (_294; section 6: the v1 05c
##     figure restored: 16S row / ITS1 row x 7 farms, curves bagged / unbagged)
## ============================================================================

# ---- 0. Parameters -------------------------------------------------------
STEP <- 500   # depth step of the curves, in reads (as in the old 05a-05c)

# Rarefaction depths of the alpha-diversity steps; must match
# 07a_hpc_16S_alpha_div_rarefaction.R and 07b_hpc_ITS1_alpha_div_rarefaction.R
ALPHA_DEPTHS <- list("16S" = c(2000, 10000), "ITS1" = c(5000, 15000))

FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")   # inside forest -> full sun
MGMT <- c(ib = "inside forest", vr = "near forest", sa = "near forest",
          kk = "agroforest", mt = "agroforest", vi = "full sun", yb = "full sun")

# Colours: treatment and farm as in 08a/08b (PALETTES)
TREATMENT_COLS <- c(bagged_flower = "#4C72B0", unbagged_flower = "#DD8452")
FARM_COLS      <- c(ib = "#1B9E77", vr = "#D95F02", sa = "#7570B3", kk = "#E7298A",
                    mt = "#66A61E", vi = "#E6AB02", yb = "#A6761D")   # Dark2

# Rows of the figure: metadata column -> panel subtitle
ROW_VARIABLES <- c(
  sample_type               = "Treatment",
  farm_id                   = "Farm",
  pi_clsm                   = "Pollen grains (CLSM), log scale",
  callose_plug_index        = "Callose plug index, log scale",
  fungal_colonization_score = "Fungal colonisation score"
)

MARKER_TITLES <- c("16S"  = "16S — Bacteria", "ITS1" = "ITS1 — Fungi")
Y_TITLES      <- c("16S"  = "Observed ASVs",  "ITS1" = "Observed OTUs (97%)")

LINE_WIDTH  <- 0.3    # curve thickness
LINE_ALPHA  <- 0.6    # curve transparency (lower = lighter where curves overlap)
X_TRANS     <- "identity"   # "sqrt" spreads out the low depths where the 07 depths sit
LABEL_SIZE  <- 2.2    # text size of the deviation labels (294 set only)
FIG_WIDTH   <- 12     # inches
FIG_HEIGHT  <- 18     # inches, about 3.4 per row

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(phyloseq)
library(vegan)
library(patchwork)
library(scales)
library(ggrepel)   # labels that do not overlap the curves or each other

# Data version and sample set: RESULTS_DIR, BIOSAMPLE_SUFFIX, SET_TAG, CFM_SAMPLE_SET
source(here("scripts", "utils_paths.R"))

LABEL_DEVIATIONS <- CFM_SAMPLE_SET == "with_deviations"

# Under Rscript, opening a plot device writes Rplots.pdf to the working directory;
# a null device prevents that stray file (the figure itself goes through ggsave)
if (!interactive()) pdf(NULL)

# ---- 2. Load -------------------------------------------------------------
ps_list <- list(
  "16S"  = readRDS(file.path(RESULTS_DIR, "rds", paste0("ps_16S_bacteria_biosamples", BIOSAMPLE_SUFFIX, ".rds"))),
  "ITS1" = readRDS(file.path(RESULTS_DIR, "rds", paste0("ps_ITS1_otu97_fungi_biosamples", BIOSAMPLE_SUFFIX, ".rds")))
)

# ---- 3. Curves -----------------------------------------------------------
# samples x taxa count matrix, as rarecurve() expects
count_matrix <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) t(m) else m
}

flower_metadata <- function(ps) {
  data.frame(sample_data(ps)) |>          # sample_data already holds sample_id (= row names)
    transmute(sample_id = as.character(sample_id), sample_type, farm_id, pi_clsm, callose_plug_index,
              fungal_colonization_score = factor(as.character(fungal_colonization_score), levels = 0:5),
              # Label of a deviation flower: its id and the deviation type, or the insect
              deviation_label = case_when(
                sampling_deviation == "none"           ~ NA_character_,
                sampling_deviation == "insect_present" ~ paste0(sample_id, " insect: ", insect_id),
                TRUE ~ paste0(sample_id, " ", str_replace_all(sampling_deviation, "_", " "))
              ))
}

# rarecurve() warns that the smallest count is 2: no flower has a singleton taxon
# in v2 (DADA2 drops global singletons). The warning is expected and harmless.
curves <- imap(ps_list, \(ps, marker) {
  count_matrix(ps) |>
    rarecurve(step = STEP, tidy = TRUE) |>
    as_tibble() |>
    transmute(sample_id = as.character(Site), depth = Sample, richness = Species, marker = marker) |>
    left_join(flower_metadata(ps), by = "sample_id")
}) |> list_rbind()

# Checks: one curve per flower, and each curve ends at the flower's depth and richness
iwalk(ps_list, \(ps, marker) {
  m    <- count_matrix(ps)
  ends <- curves |> filter(marker == !!marker) |> slice_max(depth, by = sample_id)
  stopifnot(nrow(ends) == nrow(m),
            all(ends$depth    == rowSums(m)[ends$sample_id]),
            all(ends$richness == rowSums(m > 0)[ends$sample_id]),
            all(rowSums(m) > 0))
})

# ---- 4. Plot -------------------------------------------------------------
# One colour scale per row. The two microscopy scales share their limits across
# markers, so a colour means the same value in both columns.
MAGMA_BREAKS <- c(0, 3, 10, 30, 100, 300)

ROW_SCALES <- list(
  sample_type = scale_colour_manual(values = TREATMENT_COLS, labels = c("Bagged", "Unbagged"),
                                    name = NULL),
  farm_id     = scale_colour_manual(values = FARM_COLS, breaks = FARM_LEVELS,
                                    labels = paste0(FARM_LEVELS, " (", MGMT[FARM_LEVELS], ")"),
                                    name = NULL),
  # magma, reversed so high values are dark; end = 0.85 keeps zero visible on white;
  # log1p because both variables are very skewed (median 0, max 349 and 443);
# MAGMA_BREAKS: legend ticks, spread evenly on the log scale
  pi_clsm            = scale_colour_viridis_c(option = "magma", direction = -1, end = 0.85,
                                              transform = "log1p", name = "Pollen\ngrains", breaks = MAGMA_BREAKS,
                                              limits = range(curves$pi_clsm, na.rm = TRUE)),
  callose_plug_index = scale_colour_viridis_c(option = "magma", direction = -1, end = 0.85,
                                              transform = "log1p", name = "Callose\nindex", breaks = MAGMA_BREAKS,
                                              limits = range(curves$callose_plug_index, na.rm = TRUE)),
  # Same colour family in steps: darker = more colonised
  fungal_colonization_score = scale_colour_viridis_d(option = "magma", direction = -1, begin = 0.1,
                                                     end = 0.85, drop = FALSE, name = "Score")
)

# Axis limits per marker, from all flowers, so the rows of a column compare
AXIS_LIMITS <- curves |>
  summarise(x_max = max(depth), y_max = max(richness), .by = marker)

plot_panel <- function(variable, marker) {
  df <- curves |>
    filter(marker == !!marker, !is.na(.data[[variable]])) |>
    mutate(colour_value = .data[[variable]])
  # Draw order: low values first, so flowers with high values stay on top
  if (variable %in% c("pi_clsm", "callose_plug_index", "fungal_colonization_score")) {
    df <- df |> mutate(sample_id = fct_reorder(sample_id, as.numeric(colour_value)))
  }
  lim      <- filter(AXIS_LIMITS, marker == !!marker)
  is_top   <- variable == names(ROW_VARIABLES)[1]
  is_last  <- variable == names(ROW_VARIABLES)[length(ROW_VARIABLES)]
  n_flower <- n_distinct(df$sample_id)

  p <- ggplot(df, aes(x = depth, y = richness, group = sample_id, colour = colour_value)) +
    geom_vline(xintercept = ALPHA_DEPTHS[[marker]], linetype = "dashed",
               colour = "grey55", linewidth = 0.4) +            # 07a/07b depths
    geom_line(linewidth = LINE_WIDTH, alpha = LINE_ALPHA) +
    ROW_SCALES[[variable]] +
    scale_x_continuous(transform = X_TRANS, labels = label_number(scale = 1e-3, suffix = "k")) +
    coord_cartesian(xlim = c(0, lim$x_max * 1.08), ylim = c(0, lim$y_max)) +   # 8% room for labels
    labs(title    = if (is_top) MARKER_TITLES[[marker]] else NULL,
         subtitle = sprintf("%s: %d flowers", ROW_VARIABLES[[variable]], n_flower),
         x = if (is_last) "Reads per flower" else NULL,
         y = Y_TITLES[[marker]]) +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          legend.position = if (marker == "ITS1") "right" else "none")   # one legend per row

  # Depth labels ("2k", "10k") above the dashed lines, top row only
  if (is_top) {
    p <- p + annotate("text", x = ALPHA_DEPTHS[[marker]], y = lim$y_max,
                      label = paste0(ALPHA_DEPTHS[[marker]] / 1000, "k"),
                      hjust = -0.15, vjust = 1, size = 2.8, colour = "grey40")
  }

  # 294 set: label the deviation flowers at the end of their curve
  if (LABEL_DEVIATIONS) {
    ends <- df |> filter(!is.na(deviation_label)) |> slice_max(depth, by = sample_id)
    p <- p +
      geom_point(data = ends, size = 0.9, colour = "black") +
      geom_text_repel(data = ends, aes(label = deviation_label), colour = "black",
                      size = LABEL_SIZE, direction = "y", hjust = 0, nudge_x = lim$x_max * 0.02,
                      segment.size = 0.2, segment.colour = "grey40", min.segment.length = 0,
                      box.padding = 0.15, max.overlaps = Inf, seed = 1)
  }
  p
}

p_rarefaction <- expand_grid(variable = names(ROW_VARIABLES), marker = names(ps_list)) |>
  pmap(plot_panel) |>
  wrap_plots(ncol = 2) +
  plot_annotation(
    tag_levels = "A",
    caption = sprintf(paste0("%s set. One curve per flower; flowers without a value are not shown. ",
                             "Dashed lines: rarefaction depths of the alpha-diversity steps (07a/07b).%s"),
                      sub("_", "", SET_TAG),
                      if (LABEL_DEVIATIONS) " Labelled: flowers sampled outside the standard protocol." else ""),
    theme = theme(plot.caption = element_text(size = 9, colour = "grey40", hjust = 0))
  )

# ---- 5. Save -------------------------------------------------------------
dir.create(file.path(RESULTS_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(RESULTS_DIR, "figures", paste0("05_rarefaction_curves", SET_TAG, ".png")),
       plot = p_rarefaction, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = 300)

curves |>
  distinct(marker, sample_id, .keep_all = TRUE) |>
  summarise(across(all_of(names(ROW_VARIABLES)), \(x) sum(!is.na(x))), .by = marker) |>
  print()

# ---- 6. Figure by farm (v1 05c restored) ---------------------------------
# Code of the v1 script 05c_rarefaction_by_farm.R (git show cded4c2^:scripts/...),
# base R graphics. v2 changes: ITS1 on OTU97 (v1 read the ASV object), objects and
# palette from this script, axis limits from AXIS_LIMITS (same maxima as v1's loops).
#
# Figure: 2 rows x 7 farms
#   Row 1 (top):    16S ASV bacteria
#   Row 2 (bottom): ITS1 OTU97 fungi
#
# Layout parameters:
#   oma[2] = 5 lines: left outer margin for vertical row labels ("16S", "ITS1")
#   oma[1] = 3 lines: bottom margin for shared legend strip
#   oma[3] = 2 lines: top margin for overall title
#   mar = per-panel margins; mar[3] tightened for row 2 (no panel title row)

# Count matrix of one farm (flowers x taxa)
farm_matrix <- function(ps, farm) {
  m    <- count_matrix(ps)
  meta <- data.frame(sample_data(ps))[rownames(m), , drop = FALSE]
  m[meta$farm_id == farm, , drop = FALSE]
}

# Shared limits per marker row: x = max reads, y = max richness + 10% padding
FARM_XLIM <- setNames(lapply(names(ps_list), \(mk) c(0, AXIS_LIMITS$x_max[AXIS_LIMITS$marker == mk])),
                      names(ps_list))
FARM_YLIM <- setNames(lapply(names(ps_list), \(mk) c(0, AXIS_LIMITS$y_max[AXIS_LIMITS$marker == mk] * 1.10)),
                      names(ps_list))
FARM_YLAB <- c("16S" = "Observed ASVs", "ITS1" = "Observed OTUs")

png(file.path(RESULTS_DIR, "figures", paste0("05_rarefaction_by_farm", SET_TAG, ".png")),
    width = 21, height = 10, units = "in", res = 300)

par(mfrow = c(2, length(FARM_LEVELS)),
    oma   = c(3, 5, 2, 0))

for (marker in names(ps_list)) {
  for (farm in FARM_LEVELS) {
    mat    <- farm_matrix(ps_list[[marker]], farm)
    meta_f <- data.frame(sample_data(ps_list[[marker]]))[rownames(mat), , drop = FALSE]
    cols   <- TREATMENT_COLS[as.character(meta_f$sample_type)]
    is_16S <- marker == "16S"

    # mar[3] (top): 2.5 lines for panel title only on top row; 0.4 on bottom row
    par(mar = c(2.5, 2.5, if (is_16S) 2.5 else 0.4, 0.4), mgp = c(1.5, 0.5, 0))

    vegan::rarecurve(
      mat,
      step  = STEP,
      col   = cols,
      label = FALSE,
      # x-label on the middle column of the bottom row only
      xlab  = if (!is_16S && farm == FARM_LEVELS[4]) "Sequencing depth" else "",
      ylab  = if (farm == FARM_LEVELS[1]) FARM_YLAB[[marker]] else "",   # leftmost panel only
      xlim  = FARM_XLIM[[marker]],   # shared x-limit across the panels of this marker
      ylim  = FARM_YLIM[[marker]],   # shared y-limit across the panels of this marker
      main  = if (is_16S) paste0(farm, " | ", MGMT[farm]) else ""
    )
    # vertical lines at the two rarefaction depths used in 07a / 07b alpha diversity
    abline(v = ALPHA_DEPTHS[[marker]], lty = 2, col = "grey40", lwd = 0.8)
  }
}

# ---- Row labels in outer left margin ----
# NDC positioning: centre of each row in the inner region
omi_h   <- par("omi")
din_h   <- par("din")[2]
inner_h <- din_h - omi_h[1] - omi_h[3]
row_labels <- c("16S — Bacteria", "ITS1 — Fungi")
for (i in seq_along(row_labels)) {
  # Centre of row i (1 = topmost)
  y_ndc <- (omi_h[1] + inner_h * (1 - (i - 0.5) / 2)) / din_h
  mtext(row_labels[i], side = 2, outer = TRUE,
        at = y_ndc, las = 0, cex = 0.80, line = 2)
}

mtext("Rarefaction curves by farm — 16S (top) and ITS1 (bottom)",
      outer = TRUE, cex = 1.0, line = 0.5, side = 3)

# ---- Shared legend in bottom outer margin ----
# fig = c(x1, x2, y1, y2) in normalised device coords.
# Increase 4th value (e.g. 0.06 -> 0.09) if legend appears cramped vertically.
par(fig = c(0, 1, 0, 0.06), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottom",
       legend = c("Bagged", "Unbagged"),
       col    = TREATMENT_COLS,
       lwd    = 2,
       horiz  = TRUE,
       cex    = 0.95,
       bty    = "n")

dev.off()

cat("05 done.\n")
