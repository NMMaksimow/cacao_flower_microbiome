# ---- 0. Parameters -------------------------------------------------------
# Distance metrics in canonical order (matches 20b_ITS1_pcoa_by_farms.R)
# ITS1 has no UniFrac (ITS too variable for reliable phylogenetic inference)
DIST_LEVELS <- c("aitchison_rCLR", "aitchison_CLR", "bray_curtis", "jaccard")

DIST_LABELS <- c(
  aitchison_rCLR = "Robust Aitchison (rCLR, vegan)",
  aitchison_CLR  = "Aitchison (cmultRepl, CLR)",
  bray_curtis    = "Bray-Curtis",
  jaccard        = "Jaccard"
)

# Coloring variables — one panel column each
COLOR_VARS <- c(
  "sample_type", "farm_id", "management_type",
  "is_pollination_clsm", "is_germination", "fungal_colonization_score"
)

VAR_LABELS <- c(
  sample_type              = "Sample type",
  farm_id                  = "Farm",
  management_type          = "Management",
  is_pollination_clsm      = "Pollination (CLSM)",
  is_germination           = "Germination",
  fungal_colonization_score = "Fungal colonization score"
)

# Variables for which ordiellipse is drawn (vegan silently skips NA rows)
ELLIPSE_VARS <- c(
  "sample_type", "farm_id", "management_type",
  "is_pollination_clsm", "is_germination",
  "fungal_colonization_score"
)

NA_COL <- "#BBBBBB"

# Columns stored as logical TRUE/FALSE in phyloseq sample_data (read from Excel
# with col_type = "logical"); recode to "1"/"0" to match palette keys.
LOGICAL_VARS <- c("is_pollination_clsm", "is_germination")

PALETTES <- list(
  sample_type = c(
    bagged_flower   = "#4C72B0",
    unbagged_flower = "#DD8452"
  ),
  farm_id = c(
    ib = "#1B9E77", vr = "#D95F02", sa = "#7570B3", kk = "#E7298A",
    mt = "#66A61E", vi = "#E6AB02", yb = "#A6761D"
  ),
  management_type = c(
    inside_forest = "#E41A1C",
    near_forest   = "#377EB8",
    agroforest    = "#4DAF4A",
    full_sun      = "#FF7F00"
  ),
  is_pollination_clsm = c("0" = "#4393C3", "1" = "#D6604D"),
  is_germination      = c("0" = "#4393C3", "1" = "#D6604D"),
  fungal_colonization_score = c(
    "0" = "#F7F7F7", "1" = "#FEE08B", "2" = "#FDAE61",
    "3" = "#F46D43", "4" = "#D73027", "5" = "#B30000"
  )
)

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(vegan)  # capscale(), scores(), ordiplot(), ordiellipse()

# ---- 2. Load data --------------------------------------------------------
dist_list <- readRDS(here("results", "rds", "19b_ITS1_dist_list.rds"))
meta_list <- readRDS(here("results", "rds", "19b_ITS1_meta_list.rds"))

missing_dist <- setdiff(DIST_LEVELS, names(dist_list))
if (length(missing_dist) > 0) {
  stop(sprintf(
    "dist_list missing: %s — re-run 19b_ITS1_beta_perm.R",
    paste(missing_dist, collapse = ", ")
  ))
}

# ---- 3. Global PCoA — all 294 samples across all distance metrics --------
# Layout: rows = distance metrics, columns = coloring variables.
# Algorithm:
#   - capscale(d ~ 1) computed ONCE per distance on the full matrix
#   - same ordination replotted 6 times with different point colours
#   - ordiellipse for categorical/binary variables (ELLIPSE_VARS); NA skipped by vegan
#   - points-only for fungal_colonization_score (ordinal, 7 levels incl. NA)
# Axis limits computed from the full ordination; consistent across colour panels in a row.

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

n_rows <- length(DIST_LEVELS)
n_cols <- length(COLOR_VARS)

# Figure: 4 distances × 6 colour variables
# Panel width ~4.2 in, panel height ~3.8 in + margins
# oma[2] = 5 lines: left outer margin for vertical distance-metric row labels
# oma[3] = 2.5 lines: top for overall title
png(here("results", "figures", "08b_ITS1_pcoa_global.png"),
    width = 26, height = 16, units = "in", res = 300)

par(mfrow = c(n_rows, n_cols),
    oma   = c(2, 5, 2.5, 0))

for (dname in DIST_LEVELS) {

  meta_d <- meta_list[[dname]]

  # Compute global ordination once for this distance metric
  ord <- capscale(dist_list[[dname]] ~ 1)
  sc  <- scores(ord, display = "sites")
  imp <- summary(ord)$cont$importance
  pct <- round(imp["Proportion Explained", 1:2] * 100, 1)

  # Shared axis limits (15% padding) — same across all 6 colour panels in this row
  xpad <- diff(range(sc[, 1])) * 0.15
  ypad <- diff(range(sc[, 2])) * 0.15
  xlim <- range(sc[, 1]) + c(-xpad, xpad)
  ylim <- range(sc[, 2]) + c(-ypad, ypad)

  for (col_var in COLOR_VARS) {

    # mar[1] (bottom): inter-row whitespace — increase (e.g. 2.5 → 4) for more spacing
    # mgp: axis title at 1.5 lines, tick labels at 0.5 lines
    mar_top <- if (dname == DIST_LEVELS[1]) 2.2 else 0.5
    par(mar = c(2.5, 2.5, mar_top, 0.5), mgp = c(1.5, 0.5, 0))

    # Map metadata values to colours; NA → NA_COL (grey)
    raw_vals <- as.character(meta_d[[col_var]])
    # Logical columns (TRUE/FALSE) must be recoded to "1"/"0" to match palette keys
    if (col_var %in% LOGICAL_VARS) {
      raw_vals[raw_vals == "TRUE"]  <- "1"
      raw_vals[raw_vals == "FALSE"] <- "0"
    }
    raw_vals[is.na(raw_vals) | raw_vals == "NA"] <- "NA_group"
    pal <- PALETTES[[col_var]]
    point_cols <- ifelse(raw_vals %in% names(pal), pal[raw_vals], NA_COL)

    # Column header shown only in the first distance row
    panel_main <- if (dname == DIST_LEVELS[1]) VAR_LABELS[col_var] else ""

    ordiplot(ord, display = "sites", type = "n",
             xlim = xlim, ylim = ylim,
             xlab = sprintf("PCoA1 (%.1f%%)", pct[1]),
             ylab = sprintf("PCoA2 (%.1f%%)", pct[2]),
             main = panel_main,
             cex.main = 0.85, cex.axis = 0.6, cex.lab = 0.70)

    # Confidence ellipses for categorical / binary variables
    # NA entries in groups_char are silently excluded by ordiellipse
    if (col_var %in% ELLIPSE_VARS) {
      groups_char <- as.character(meta_d[[col_var]])
      pal_for_ellipse <- PALETTES[[col_var]]
      ordiellipse(ord, groups = groups_char,
                  kind = "sd", conf = 0.865,
                  col  = pal_for_ellipse,
                  lwd  = 1.3, draw = "lines")
    }

    # Points drawn last (on top of ellipses)
    # cex: point size — increase (e.g. 0.55 → 0.8) if samples overlap too little
    points(sc, pch = 16, cex = 0.55, col = point_cols)

    # Group centroids — pch = 4 (x-cross), same approach as 20a/20b
    # tapply groups sample indices by raw_vals; colMeans gives the centroid coordinates.
    # valid_idx excludes "NA_group" so NA samples never produce a spurious centroid.
    if (col_var %in% ELLIPSE_VARS) {
      pal_c     <- PALETTES[[col_var]]
      valid_idx <- which(raw_vals %in% names(pal_c))
      ctrs      <- tapply(valid_idx, raw_vals[valid_idx],
                          \(i) colMeans(sc[i, , drop = FALSE]))
      for (grp in names(ctrs)) {
        points(t(ctrs[[grp]]), pch = 4, cex = 1.5, lwd = 2, col = pal_c[grp])
      }
    }
  }
}

# ---- Row labels (left outer margin) ----------------------------------------
# NDC positioning: same approach as 20b_ITS1_pcoa_by_farms.R
omi_h   <- par("omi")
din_h   <- par("din")[2]
inner_h <- din_h - omi_h[1] - omi_h[3]
for (i in seq_along(DIST_LEVELS)) {
  y_ndc <- (omi_h[1] + inner_h * (1 - (i - 0.5) / n_rows)) / din_h
  mtext(DIST_LABELS[DIST_LEVELS[i]], side = 2, outer = TRUE,
        at = y_ndc, las = 0, cex = 0.65, line = 2)
}

mtext("Global PCoA — ITS1 OTU97 fungi | prev >= 10% | BC/Jaccard rarefied 5k reads",
      outer = TRUE, cex = 1.0, line = 1, side = 3)

dev.off()

message("Done. Figure saved to results/figures/08b_ITS1_pcoa_global.png")
