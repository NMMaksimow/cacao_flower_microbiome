# ---- 0. Parameters -------------------------------------------------------
PREV_THRESHOLD <- 0.10   # must match 08b (prevalence filter to reproduce same samples)
FARM_LEVELS    <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
MGMT           <- c(ib = "inside_forest", vr = "near_forest",  sa = "near_forest",
                    kk = "agroforest",   mt = "agroforest",   vi = "full_sun",
                    yb = "full_sun")
# Colours match Python notebook (notebooks/01b_its1_pca_by_farm.ipynb)
COLS <- c(bagged_flower = "#4C72B0", unbagged_flower = "#DD8452")
# Canonical metric order (top → bottom rows): non-rarefied Aitchison first,
# then rarefied abundance-weighted, then rarefied binary.
# ITS1 has no UniFrac (ITS too variable for reliable phylogenetic inference).
DIST_LEVELS  <- c("aitchison_rCLR", "aitchison_CLR", "bray_curtis", "jaccard")
# Human-readable row labels for the outer left margin of the PCoA figure
DIST_LABELS <- c(
  aitchison_rCLR = "Robust Aitchison (rCLR, vegan)",
  aitchison_CLR  = "Aitchison (cmultRepl, CLR)",
  bray_curtis    = "Bray-Curtis",
  jaccard        = "Jaccard"
)

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(phyloseq)
library(vegan)   # capscale(), ordiplot(), ordiellipse(), scores()

# ---- 2. Load data --------------------------------------------------------
ps <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds")) |>
  subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))

ps_prev <- ps |>
  (\(x) prune_taxa(
    rowSums(otu_table(x) > 0) / nsamples(x) >= PREV_THRESHOLD, x
  ))()

dist_list <- readRDS(here("results", "rds", "ITS1_otu97_dist_list.rds"))
meta_list <- readRDS(here("results", "rds", "ITS1_otu97_meta_list.rds"))

missing_keys <- setdiff(DIST_LEVELS, names(dist_list))
if (length(missing_keys) > 0) {
  stop(sprintf(
    "dist_list missing: %s — re-run 08b_ITS1_pcoa_global.R to regenerate the RDS.",
    paste(missing_keys, collapse = ", ")
  ))
}
missing_meta <- setdiff(DIST_LEVELS, names(meta_list))
if (length(missing_meta) > 0) {
  stop(sprintf(
    "meta_list missing: %s — re-run 08b_ITS1_pcoa_global.R to regenerate the RDS.",
    paste(missing_meta, collapse = ", ")
  ))
}

# ---- 3. Within-farm PCoA (shared axis limits within each distance row) ---
# Layout: rows = 4 distances, columns = 7 farms (28 panels total).
# PCoA via vegan::capscale(d ~ 1): unconstrained ordination.
# Ellipses: ordiellipse(kind="sd", conf=0.865) → 2-sigma (matches Python n_std=2.0).
# Shared limits per row: all 7 farm panels share the same xlim/ylim so farm-to-farm
# community spread is directly comparable within each distance metric.
# Two-pass approach: pass 1 caches ordinations + finds row-level axis range;
# pass 2 draws all panels with the shared limits.
# Aitchison distances use all samples; rarefied distances (BC, Jaccard) use only
# samples that survived rarefaction — panel n may differ from the first row.

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

png(here("results", "figures", "09b_ITS1_pcoa_per_farm_shared_lim.png"),
    width = 20, height = 11, units = "in", res = 300)

# mfrow: 4 distance rows × 7 farm columns
# mar: bottom, left, top, right margins per panel (in lines)
# oma: outer margin: bottom (for legend) and top (for overall title)
par(mfrow = c(length(DIST_LEVELS), length(FARM_LEVELS)),
    mar   = c(2.8, 2.8, 2.2, 0.5),  # overridden per row inside loop
    oma   = c(3, 5, 2.5, 0))         # oma[2]=5: left margin for vertical row labels

for (dname in DIST_LEVELS) {
  meta_d <- meta_list[[dname]]

  # farm_n for panel titles: sample count per farm in this distance's metadata
  farm_n_d <- meta_d |> count(farm_id) |> deframe()

  # Pass 1: run capscale() for all farms, cache results, track axis range
  ord_cache  <- list()
  meta_cache <- list()
  row_xlim   <- c(Inf, -Inf)
  row_ylim   <- c(Inf, -Inf)

  for (farm in FARM_LEVELS) {
    farm_idx           <- which(!is.na(meta_d$farm_id) & meta_d$farm_id == farm)
    meta_cache[[farm]] <- meta_d[farm_idx, ]
    d_mat              <- as.matrix(dist_list[[dname]])
    d_farm             <- as.dist(d_mat[farm_idx, farm_idx])
    ord_cache[[farm]]  <- capscale(d_farm ~ 1)
    sc <- scores(ord_cache[[farm]], display = "sites")
    row_xlim <- c(min(row_xlim[1], min(sc[, 1])), max(row_xlim[2], max(sc[, 1])))
    row_ylim <- c(min(row_ylim[1], min(sc[, 2])), max(row_ylim[2], max(sc[, 2])))
  }
  # 15% padding so points/ellipses don't hug the axis limits
  xpad <- diff(row_xlim) * 0.15; row_xlim <- row_xlim + c(-xpad, xpad)
  ypad <- diff(row_ylim) * 0.15; row_ylim <- row_ylim + c(-ypad, ypad)

  # Pass 2: draw panels with the shared row-level limits
  for (farm in FARM_LEVELS) {
    ord       <- ord_cache[[farm]]
    meta_farm <- meta_cache[[farm]]
    # xpd=FALSE clips ALL drawing (incl. vegan's abline(h=0)/abline(v=0) origin lines)
    # to the current panel boundary; must be reset at start of each panel iteration.
    # mar_top: tight top margin for rows 2–4 since they have no panel title.
    mar_top <- if (dname == DIST_LEVELS[1]) 2.2 else 0.5
    # mar = c(bottom, left, top, right) in lines per panel.
    # mar[1] (bottom): controls the vertical gap between PCoA rows — increase to add
    #   more whitespace between rows (e.g. 3.5 → 4.5); decrease to tighten layout.
    # mar[2] (left): increase if Y-axis tick labels get clipped.
    # mgp = c(title_line, ticklabel_line, tickmark_line): axis title sits at 1.8 lines
    #   (must stay < mar[1] to avoid clipping); tick labels at 0.5 lines.
    par(xpd = FALSE, mar = c(3.5, 2.8, mar_top, 0.5), mgp = c(1.8, 0.5, 0))

    imp <- summary(ord)$cont$importance
    # "Proportion Explained" gives each axis's fraction of total inertia
    pct <- round(imp["Proportion Explained", 1:2] * 100, 1)

    # Farm + management title only for the top row; empty string for rows 2–4
    panel_main <- if (dname == DIST_LEVELS[1])
      paste0(farm, " | ", MGMT[farm], "\nn = ", farm_n_d[farm]) else ""

    ordiplot(ord, display = "sites", type = "n",
             xlab = sprintf("PCoA1 (%.1f%%)", pct[1]),
             ylab = sprintf("PCoA2 (%.1f%%)", pct[2]),
             main = panel_main,
             xlim = row_xlim, ylim = row_ylim,   # shared limits for this row
             cex.main = 0.75, cex.axis = 0.65, cex.lab = 0.70)

    # 2-sigma ellipses (conf=0.865 → sqrt(qchisq(0.865,2)) ≈ 2.0 σ, matching
    # Python notebooks confidence_ellipse(n_std=2.0); xpd=FALSE already active)
    ordiellipse(ord, groups = meta_farm$sample_type,
                kind = "sd", conf = 0.865,
                col = COLS, lwd = 1.5, draw = "lines")

    # Individual sample points; pch=16 = filled circle; cex = point size
    points(scores(ord, display = "sites"),
           pch = 16, cex = 0.7,
           col = COLS[meta_farm$sample_type])

    # Group centroids as X-cross (pch = 4); lwd controls cross arm thickness
    site_scores <- scores(ord, display = "sites")
    ctrs <- tapply(seq_len(nrow(site_scores)), meta_farm$sample_type,
                   \(i) colMeans(site_scores[i, , drop = FALSE]))
    for (grp in names(ctrs)) {
      points(t(ctrs[[grp]]), pch = 4, cex = 1.5, lwd = 2, col = COLS[grp])
    }
  }
}

# Row labels — drawn after all panels so par("omi") / par("din") are finalised.
# Explicit at= in normalised device coords (0=bottom, 1=top of full device) is
# required; without it mtext(outer=TRUE) ignores the outer margins and drifts.
omi_h   <- par("omi")    # outer margins in inches: [1]=bottom, [3]=top
din_h   <- par("din")[2] # device height in inches
inner_h <- din_h - omi_h[1] - omi_h[3]
n_rows  <- length(DIST_LEVELS)
for (i in seq_along(DIST_LEVELS)) {
  # Centre of row i (1=topmost) accounting for bottom + top outer margins
  y_ndc <- (omi_h[1] + inner_h * (1 - (i - 0.5) / n_rows)) / din_h
  # las=0: text parallel to axis = vertical (bottom-to-top) for side=2
  # line=2: centres label in the 5-line outer left margin
  mtext(DIST_LABELS[DIST_LEVELS[i]], side = 2, outer = TRUE,
        at = y_ndc, las = 0, cex = 0.70, line = 2)
}

mtext("Within-farm PCoA — ITS1 OTU97 fungi (shared axis limits within distance rows)",
      outer = TRUE, cex = 1.1, line = 1, side = 3)

# fig = c(x1, x2, y1, y2) in normalised device coords (0 = bottom/left, 1 = top/right).
# c(0, 1, 0, 0.07) reserves the bottom 7 % of the device height for the legend band.
# Increase the 4th value (e.g. 0.07 → 0.10) if the legend feels cramped vertically.
par(fig = c(0, 1, 0, 0.07), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
# Legend vertical position within the band: "bottom" anchors the box to the lower edge,
# "center" centres it, "top" raises it. Change the keyword to shift the legend up/down.
legend("bottom",
       legend = c("bagged", "unbagged"),
       col    = COLS,
       pch    = 16,
       lty    = 1,
       horiz  = TRUE,
       cex    = 0.95,
       bty    = "n",
       pt.cex = 1.3)

dev.off()

message("Done. PCoA saved to results/figures/09b_ITS1_pcoa_per_farm_shared_lim.png")
