# ---- 0. Parameters -------------------------------------------------------
PREV_THRESHOLD <- 0.10    # prevalence filter: features in >= 10% of all biosamples
RAREFY_DEPTH   <- 15000L  # rarefaction depth for depth-sensitive metrics (BC, Jaccard)

# Distance metrics in canonical order (matches 09b_ITS1_pcoa_by_farm.R)
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
  sample_type               = "Sample type",
  farm_id                   = "Farm",
  management_type           = "Management",
  is_pollination_clsm       = "Pollination (CLSM)",
  is_germination            = "Germination",
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
library(phyloseq)
library(vegan)          # capscale(), vegdist(), scores(), ordiplot(), ordiellipse()
library(zCompositions)  # cmultRepl() for multiplicative zero replacement (CLR pipeline)

# ---- 2. Load & filter ----------------------------------------------------
ps <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds")) |>
  subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))

# Global prevalence filter: keep features present in >= 10% of all biosamples.
# Mirrors filter_counts(counts_raw, 0.10) from notebooks/01b_its1_pca_by_farm.ipynb.
ps_prev <- ps |>
  (\(x) prune_taxa(
    rowSums(otu_table(x) > 0) / nsamples(x) >= PREV_THRESHOLD, x
  ))()

cat(sprintf("After prevalence filter: %d OTUs x %d samples\n",
            ntaxa(ps_prev), nsamples(ps_prev)))

# ---- 3. Distance matrices ------------------------------------------------
otu_mat <- as(otu_table(ps_prev), "matrix")          # features x samples
if (!taxa_are_rows(ps_prev)) otu_mat <- t(otu_mat)
otu_t   <- t(otu_mat)                                 # samples x features

# Create meta before distances so we can use its rownames as the canonical
# sample ID source — sample_data() preserves the original Excel IDs reliably.
meta    <- as(sample_data(ps_prev), "data.frame")
rownames(otu_t) <- rownames(meta)

# Robust Aitchison distance: rCLR + Euclidean, WITHOUT matrix completion.
# vegan >= 2.6-2 implements rCLR natively (geometric mean from non-zero values only;
# zeros receive 0). This differs from full RPCA (deicode, Python notebooks):
# deicode applies OptSpace matrix completion before rCLR. Use vegan rCLR for
# PERMANOVA — matrix completion is not designed for distance-matrix-based testing.
dist_ait <- vegdist(otu_t, method = "robust.aitchison")
# vegdist(robust.aitchison) sets Labels = "1","2",... instead of sample IDs
# (internal rCLR C routine strips rownames). Force correct Labels.
attr(dist_ait, "Labels") <- rownames(meta)

# CLR with multiplicative replacement (Martin-Fernandez et al. 2003).
# cmultRepl replaces zeros while preserving compositional totals (CZM method).
# decostand(method="clr") is safe here because cmultRepl has already removed all zeros.
otu_t_mr     <- zCompositions::cmultRepl(otu_t, label = 0, method = "CZM",
                                          output = "p-counts", z.delete = FALSE)
clr_mat      <- vegan::decostand(otu_t_mr, method = "clr")  # log(x) - mean(log(x)) per row
dist_ait_clr <- dist(clr_mat)
# dist() already sets Labels = rownames(clr_mat); no manual override needed

# Rarefaction for depth-sensitive metrics.
# Rationale: unbagged flowers have systematically more reads than bagged across all farms
# (confirmed in 04b depth boxplots). Without rarefaction, depth confounds the bagged vs
# unbagged comparison for Bray-Curtis and Jaccard.
# Aitchison metrics above are inherently depth-robust and stay non-rarefied.
# ITS1 uses 15000 reads, consistent with alpha diversity rarefaction depth; ITS1 libraries are deeper than 16S.
set.seed(42)
ps_rare <- suppressMessages(
  phyloseq::rarefy_even_depth(ps_prev, sample.size = RAREFY_DEPTH,
                               rngseed = 42L, replace = FALSE,
                               trimOTUs = TRUE, verbose = FALSE)
)
cat(sprintf("After rarefaction to %d reads: %d samples remaining (dropped %d)\n",
            RAREFY_DEPTH, nsamples(ps_rare),
            nsamples(ps_prev) - nsamples(ps_rare)))

meta_rare  <- as(sample_data(ps_rare), "data.frame")
otu_t_rare <- as(otu_table(ps_rare), "matrix") |>
  (\(m) if (taxa_are_rows(ps_rare)) t(m) else m)()
rownames(otu_t_rare) <- rownames(meta_rare)

# Bray-Curtis on rarefied counts (abundance-weighted)
dist_bc  <- vegdist(otu_t_rare, method = "bray")
attr(dist_bc,  "Labels") <- rownames(meta_rare)

# Jaccard on rarefied presence/absence (binary)
dist_jac <- vegdist(otu_t_rare, method = "jaccard", binary = TRUE)
attr(dist_jac, "Labels") <- rownames(meta_rare)

dist_list <- list(
  aitchison_rCLR = dist_ait,      # non-rarefied: robust rCLR (zeros as missing)
  aitchison_CLR  = dist_ait_clr,  # non-rarefied: CLR + multiplicative zero replacement
  bray_curtis    = dist_bc,        # rarefied 15k: abundance-weighted dissimilarity
  jaccard        = dist_jac        # rarefied 15k: binary presence/absence
)

# meta_list maps each distance key to the correct metadata frame for that distance.
# Rarefied distances use meta_rare (samples that survived rarefaction); non-rarefied
# distances use meta (all biological samples after prevalence filter).
meta_list <- list(
  aitchison_rCLR = meta,
  aitchison_CLR  = meta,
  bray_curtis    = meta_rare,
  jaccard        = meta_rare
)

saveRDS(dist_list, here("results", "rds", "ITS1_otu97_dist_list.rds"))
saveRDS(meta_list, here("results", "rds", "ITS1_otu97_meta_list.rds"))

# ---- 4. Global PCoA — all samples across all distance metrics ------------
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

# Figure: 4 distances x 6 colour variables
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

    # mar[1] (bottom): inter-row whitespace — increase (e.g. 2.5 -> 4) for more spacing
    # mgp: axis title at 1.5 lines, tick labels at 0.5 lines
    mar_top <- if (dname == DIST_LEVELS[1]) 2.2 else 0.5
    par(mar = c(2.5, 2.5, mar_top, 0.5), mgp = c(1.5, 0.5, 0))

    # Map metadata values to colours; NA -> NA_COL (grey)
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
    # cex: point size — increase (e.g. 0.55 -> 0.8) if samples overlap too little
    points(sc, pch = 16, cex = 0.55, col = point_cols)

    # Group centroids — pch = 4 (x-cross), same approach as 09b_ITS1_pcoa_by_farm.R
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
# NDC positioning: same approach as 09b_ITS1_pcoa_by_farm.R
omi_h   <- par("omi")
din_h   <- par("din")[2]
inner_h <- din_h - omi_h[1] - omi_h[3]
for (i in seq_along(DIST_LEVELS)) {
  y_ndc <- (omi_h[1] + inner_h * (1 - (i - 0.5) / n_rows)) / din_h
  mtext(DIST_LABELS[DIST_LEVELS[i]], side = 2, outer = TRUE,
        at = y_ndc, las = 0, cex = 0.65, line = 2)
}

mtext("Global PCoA — ITS1 OTU97 fungi | prev >= 10% | BC/Jaccard rarefied 15k reads",
      outer = TRUE, cex = 1.0, line = 1, side = 3)

dev.off()

message("Done. Figure saved to results/figures/08b_ITS1_pcoa_global.png")
