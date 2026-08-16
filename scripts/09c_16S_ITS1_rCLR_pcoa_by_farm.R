# ---- 0. Parameters -------------------------------------------------------
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
MGMT        <- c(ib = "inside_forest", vr = "near_forest",  sa = "near_forest",
                 kk = "agroforest",   mt = "agroforest",   vi = "full_sun",
                 yb = "full_sun")
# Colours match Python notebooks (01a_16s_pca_by_farm.ipynb, 01b_its1_pca_by_farm.ipynb)
COLS <- c(bagged_flower = "#4C72B0", unbagged_flower = "#DD8452")
# Left-margin labels: one per row, top to bottom
ROW_LABELS <- c("16S bacteria (rCLR Aitchison)",
                "ITS1 OTU97 fungi (rCLR Aitchison)")

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(vegan)   # capscale(), ordiplot(), ordiellipse(), scores()

# ---- 2. Load data --------------------------------------------------------
# Only the rCLR slice is needed from each marker's RDS; no phyloseq loading required.
dist_list_16s <- readRDS(here("results", "rds", "16S_dist_list.rds"))
meta_list_16s <- readRDS(here("results", "rds", "16S_meta_list.rds"))

dist_list_its <- readRDS(here("results", "rds", "ITS1_otu97_dist_list.rds"))
meta_list_its <- readRDS(here("results", "rds", "ITS1_otu97_meta_list.rds"))

for (rds_name in c("16S_dist_list.rds", "ITS1_otu97_dist_list.rds")) {
  obj  <- if (grepl("16S", rds_name)) dist_list_16s else dist_list_its
  label <- if (grepl("16S", rds_name)) "08a" else "08b"
  if (!"aitchison_rCLR" %in% names(obj))
    stop(sprintf("%s missing 'aitchison_rCLR' — re-run %s to regenerate.", rds_name, label))
}

d_16s  <- dist_list_16s[["aitchison_rCLR"]]
m_16s  <- meta_list_16s[["aitchison_rCLR"]]
d_its  <- dist_list_its[["aitchison_rCLR"]]
m_its  <- meta_list_its[["aitchison_rCLR"]]

# ---- 3. Combined 2-row PCoA figure ---------------------------------------
# Layout: row 1 = 16S rCLR, row 2 = ITS1 OTU97 rCLR; columns = 7 farms.
# Each row has independent shared axis limits (16S and ITS1 live in different
# Aitchison spaces); limits computed in pass 1, applied in pass 2.
# Panel titles: row 1 shows farm code + management + 16S n; row 2 shows ITS1 n only
# (farm/management already visible in row 1 above each column).
# Ellipses: kind="sd", conf=0.865 → 2-sigma, matches Python confidence_ellipse(n_std=2).

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

png(here("results", "figures", "09c_16S_ITS1_rCLR_pcoa_by_farm.png"),
    width = 20, height = 9, units = "in", res = 300)

# mfrow: 2 marker rows × 7 farm columns
# oma[1]=3: bottom outer margin for legend band
# oma[2]=5: left outer margin for vertical row labels
# oma[3]=2.5: top outer margin for overall title
par(mfrow = c(2, length(FARM_LEVELS)),
    mar   = c(2.8, 2.8, 2.2, 0.5),
    oma   = c(3, 5, 2.5, 0))

markers <- list(
  list(d = d_16s, m = m_16s, row_idx = 1L),
  list(d = d_its, m = m_its, row_idx = 2L)
)

for (mk in markers) {
  dist_mat <- as.matrix(mk$d)
  meta_d   <- mk$m
  is_16s   <- mk$row_idx == 1L

  # sample counts per farm for panel titles
  farm_n <- meta_d |> count(farm_id) |> deframe()

  # Pass 1: cache ordinations and find shared xlim/ylim for this row
  ord_cache  <- list()
  meta_cache <- list()
  row_xlim   <- c(Inf, -Inf)
  row_ylim   <- c(Inf, -Inf)

  for (farm in FARM_LEVELS) {
    farm_idx           <- which(!is.na(meta_d$farm_id) & meta_d$farm_id == farm)
    meta_cache[[farm]] <- meta_d[farm_idx, ]
    d_farm             <- as.dist(dist_mat[farm_idx, farm_idx])
    ord_cache[[farm]]  <- capscale(d_farm ~ 1)
    sc <- scores(ord_cache[[farm]], display = "sites")
    row_xlim <- c(min(row_xlim[1], min(sc[, 1])), max(row_xlim[2], max(sc[, 1])))
    row_ylim <- c(min(row_ylim[1], min(sc[, 2])), max(row_ylim[2], max(sc[, 2])))
  }
  # 15% padding so points/ellipses don't hug the axis limits
  xpad <- diff(row_xlim) * 0.15; row_xlim <- row_xlim + c(-xpad, xpad)
  ypad <- diff(row_ylim) * 0.15; row_ylim <- row_ylim + c(-ypad, ypad)

  # Pass 2: draw panels with shared limits
  for (farm in FARM_LEVELS) {
    ord       <- ord_cache[[farm]]
    meta_farm <- meta_cache[[farm]]
    # xpd=FALSE clips drawing (including vegan's origin lines) to current panel boundary
    # mar_top: full top margin in row 1 (has title); tight in row 2 (no title)
    mar_top <- if (is_16s) 2.2 else 0.5
    # mar = c(bottom, left, top, right) in lines.
    # mar[1]: vertical gap between rows — increase (e.g. 3.5→4.5) to add whitespace.
    # mgp = c(axis_title, tick_labels, tick_marks) in lines from panel edge.
    par(xpd = FALSE, mar = c(3.5, 2.8, mar_top, 0.5), mgp = c(1.8, 0.5, 0))

    imp <- summary(ord)$cont$importance
    # "Proportion Explained": each axis's fraction of total inertia
    pct <- round(imp["Proportion Explained", 1:2] * 100, 1)

    # Row 1: full title with farm code, management type, and 16S n
    # Row 2: just ITS1 n (farm/mgmt already shown in row 1 above)
    panel_main <- if (is_16s)
      paste0(farm, " | ", MGMT[farm], "\nn = ", farm_n[farm])
    else
      ""

    ordiplot(ord, display = "sites", type = "n",
             xlab = sprintf("PCoA1 (%.1f%%)", pct[1]),
             ylab = sprintf("PCoA2 (%.1f%%)", pct[2]),
             main = panel_main,
             xlim = row_xlim, ylim = row_ylim,
             cex.main = 0.75, cex.axis = 0.65, cex.lab = 0.70)

    # 2-sigma ellipses (conf=0.865 → sqrt(qchisq(0.865,2)) ≈ 2.0 σ)
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

# Row labels in outer left margin.
# Explicit at= in normalised device coords required; mtext(outer=TRUE) without at=
# ignores outer margins and drifts. Formula mirrors 09a/09b.
omi_h   <- par("omi")    # outer margins in inches: [1]=bottom, [3]=top
din_h   <- par("din")[2] # device height in inches
inner_h <- din_h - omi_h[1] - omi_h[3]
for (i in seq_along(ROW_LABELS)) {
  # Centre of row i (1=topmost row) in normalised device coords
  y_ndc <- (omi_h[1] + inner_h * (1 - (i - 0.5) / length(ROW_LABELS))) / din_h
  # las=0: text parallel to axis = vertical (bottom-to-top) for side=2
  # line=2: centres label within the 5-line outer left margin
  mtext(ROW_LABELS[i], side = 2, outer = TRUE,
        at = y_ndc, las = 0, cex = 0.72, line = 2)
}

mtext("Within-farm PCoA — rCLR Aitchison (shared axis limits within marker rows)",
      outer = TRUE, cex = 1.1, line = 1, side = 3)

# Legend band: reserves the bottom 7% of device height.
# Increase the 4th value (0.07 → 0.10) if the legend feels cramped vertically.
par(fig = c(0, 1, 0, 0.07), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
# Legend vertical position: "bottom" anchors box to lower edge of the band.
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

message("Done. Figure saved to results/figures/09c_16S_ITS1_rCLR_pcoa_by_farm.png")
