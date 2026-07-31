# ---- Packages ----
library(here)
library(phyloseq)
library(vegan)

# ---- Parameters ----
# Farm order and management labels used throughout the project
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
MGMT <- c(
        ib = "inside_forest", vr = "near_forest",  sa = "near_forest",
        kk = "agroforest",   mt = "agroforest",   vi = "full_sun",
        yb = "full_sun"
)
# Canonical bagged/unbagged palette matching 20a_16S_pcoa_by_farms.R and 20b_ITS1_pcoa_by_farms.R
COLS <- c(bagged_flower = "#4C72B0", unbagged_flower = "#DD8452")

# Rarefaction depths used in 07a/07b alpha diversity analysis
DEPTH_16S_SHALLOW  <- 2000
DEPTH_16S_DEEP     <- 10000
DEPTH_ITS1_SHALLOW <- 5000
DEPTH_ITS1_DEEP    <- 15000

# ---- Load data ----
ps_16S  <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds"))
ps_ITS1 <- readRDS(here("results", "rds", "ps_ITS1_fungi_biosamples.rds"))

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

################################################################################
# Build per-farm matrices — pass 1: collect shared y-limits per marker row
################################################################################

build_farm_mat <- function(ps, farm) {
        ps_farm <- subset_samples(ps, farm_id == farm)
        mat     <- t(as(otu_table(ps_farm), "matrix"))
        mat[rowSums(mat) > 0, , drop = FALSE]
}

# Pass 1: find max observed richness across all farms per marker row
max_obs_16S  <- 0
max_obs_ITS1 <- 0
for (farm in FARM_LEVELS) {
        m16  <- build_farm_mat(ps_16S,  farm)
        mIT  <- build_farm_mat(ps_ITS1, farm)
        if (nrow(m16)  > 0) max_obs_16S  <- max(max_obs_16S,  max(rowSums(m16  > 0)))
        if (nrow(mIT)  > 0) max_obs_ITS1 <- max(max_obs_ITS1, max(rowSums(mIT  > 0)))
}
# 10% padding above max observed richness
ylim_16S  <- c(0, max_obs_16S  * 1.10)
ylim_ITS1 <- c(0, max_obs_ITS1 * 1.10)

# Pass 1b: shared x-limits per marker row (max total reads across all farms)
max_depth_16S  <- 0
max_depth_ITS1 <- 0
for (farm in FARM_LEVELS) {
        m16 <- build_farm_mat(ps_16S,  farm)
        mIT <- build_farm_mat(ps_ITS1, farm)
        if (nrow(m16) > 0) max_depth_16S  <- max(max_depth_16S,  max(rowSums(m16)))
        if (nrow(mIT) > 0) max_depth_ITS1 <- max(max_depth_ITS1, max(rowSums(mIT)))
}
xlim_16S  <- c(0, max_depth_16S)
xlim_ITS1 <- c(0, max_depth_ITS1)

################################################################################
# Figure: 2 rows × 7 farms
#   Row 1 (top):    16S ASV bacteria
#   Row 2 (bottom): ITS1 OTU97 fungi
#
# Layout parameters:
#   oma[2] = 5 lines: left outer margin for vertical row labels ("16S", "ITS1")
#   oma[1] = 3 lines: bottom margin for shared legend strip
#   oma[3] = 2 lines: top margin for overall title
#   mar = per-panel margins; mar[3] tightened for rows 2+ (no panel title row)
################################################################################

png(here("results", "figures", "05c_rarefaction_by_farm.png"),
    width = 21, height = 10, units = "in", res = 300)

par(mfrow = c(2, length(FARM_LEVELS)),
    oma   = c(3, 5, 2, 0))

# ---- Row 1: 16S ----
for (farm in FARM_LEVELS) {
        mat    <- build_farm_mat(ps_16S, farm)
        meta_f <- as.data.frame(sample_data(ps_16S))[rownames(mat), , drop = FALSE]
        cols   <- COLS[as.character(meta_f$sample_type)]

        # mar[3] (top): 2.5 lines for panel title only on top row; 0.4 on bottom row
        par(mar = c(2.5, 2.5, 2.5, 0.4), mgp = c(1.5, 0.5, 0))

        vegan::rarecurve(
                mat,
                step  = 500,
                col   = cols,
                label = FALSE,
                xlab  = "",                                                          # x-label only on bottom row
                ylab  = if (farm == FARM_LEVELS[1]) "Observed ASVs" else "",        # y-label on leftmost panel only
                xlim  = xlim_16S,   # shared x-limit across all 16S panels
                ylim  = ylim_16S,   # shared y-limit across all 16S panels
                main  = paste0(farm, " | ", MGMT[farm])
        )
        # vertical lines at the two rarefaction depths used in 07a alpha diversity
        abline(v = c(DEPTH_16S_SHALLOW, DEPTH_16S_DEEP), lty = 2, col = "grey40", lwd = 0.8)
}

# ---- Row 2: ITS1 ----
for (farm in FARM_LEVELS) {
        mat    <- build_farm_mat(ps_ITS1, farm)
        meta_f <- as.data.frame(sample_data(ps_ITS1))[rownames(mat), , drop = FALSE]
        cols   <- COLS[as.character(meta_f$sample_type)]

        # mar[3] (top): tight top margin — no panel title on bottom row
        par(mar = c(2.5, 2.5, 0.4, 0.4), mgp = c(1.5, 0.5, 0))

        vegan::rarecurve(
                mat,
                step  = 500,
                col   = cols,
                label = FALSE,
                xlab  = if (farm == FARM_LEVELS[4]) "Sequencing depth" else "",     # x-label on middle column only
                ylab  = if (farm == FARM_LEVELS[1]) "Observed OTUs" else "",        # y-label on leftmost panel only
                xlim  = xlim_ITS1,  # shared x-limit across all ITS1 panels
                ylim  = ylim_ITS1   # shared y-limit across all ITS1 panels
        )
        # vertical lines at the two rarefaction depths used in 07b alpha diversity
        abline(v = c(DEPTH_ITS1_SHALLOW, DEPTH_ITS1_DEEP), lty = 2, col = "grey40", lwd = 0.8)
}

# ---- Row labels in outer left margin ----
# NDC positioning: same approach as 20a_16S_pcoa_by_farms.R
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
# fig = c(x1, x2, y1, y2) in normalised device coords; same pattern as 20a/20b.
# Increase 4th value (e.g. 0.06 → 0.09) if legend appears cramped vertically.
par(fig = c(0, 1, 0, 0.06), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottom",
       legend = c("Bagged", "Unbagged"),
       col    = COLS,
       lwd    = 2,
       horiz  = TRUE,
       cex    = 0.95,
       bty    = "n")

dev.off()

message("Done. Figure saved to results/figures/05c_rarefaction_by_farm.png")
