## ============================================================================
## 09d — Batch-effect check on processing controls (16S + ITS1), rCLR Aitchison
##
## Extraction ran in batches of 21 samples (centrifuge capacity), so one batch =
## half a farm of one treatment. Bagged flowers were extracted BEFORE unbagged on
## every one of the seven farms (gaps of 0-5 days). Any drift across that gap --
## reagent lot, kit age, ambient contamination -- would push all bagged one way and
## all unbagged the other, mimicking a treatment effect. Within-farm centering does
## not protect against this: it subtracts the mean of a farm's two batches and
## leaves the difference between them intact.
##
## Controls are biologically identical to each other, so any difference between
## those processed on bagged-days and those processed on unbagged-days is purely
## technical. Each control type is split exactly 7 vs 7.
##
## Input (BEFORE decontam -- decontam removes the very signal being measured, since
## the contaminant profile IS the batch signature):
##   results/rds/ps_16S_raw.rds
##   results/rds/ps_ITS1_otu97_raw.rds
## Output:
##   results/figures/09d_controls_rCLR_pcoa_by_day_type.png
##   results/figures/09d_controls_perm_heatmap.png
##   results/tables/09d_controls_permanova.csv
##   results/tables/09d_controls_depth_richness.csv
##
## This checks for process DRIFT. It does not remove the confounding: batches stay
## nested within farm x treatment, which no analysis can undo.
## ============================================================================

# ---- 0. Parameters ---------------------------------------------------------

N_PERM      <- 9999   # 7 vs 7 allows 3432 distinct splits, so the test is near-exact
PREV_THRESH <- 0.10   # within the 14-control subset: 10% of 14 means "seen in >= 2"
DEPTH_MIN   <- 1000   # sensitivity setting only; the primary panels keep all samples

# cibPC sequenced at 223 reads on 16S against 47,000 for the next lowest -- the
# library failed. On ITS1 the same sample is in the normal range (all mocks
# 34,707-56,319), so it is dropped from 16S only. Named explicitly rather than
# caught by DEPTH_MIN so the reason survives in the code.
FAILED_16S <- c("cibPC")

CONTROL_TYPES <- c(mock_community       = "Mock community",
                   extraction_control   = "Extraction blank",
                   negative_pcr_control = "PCR blank")

MARKERS <- c(`16S`  = "ps_16S_raw",
             ITS1   = "ps_ITS1_otu97_raw")

# Same hex values as 10d, where they mean bagged / unbagged flowers. Here they mean
# the day a control was processed alongside, which the legend must state explicitly.
DAY_COLS <- c(bagged_day = "#4C72B0", unbagged_day = "#DD8452")

# ---- 1. Libraries ----------------------------------------------------------

library(here)
library(tidyverse)
library(phyloseq)
library(vegan)      # vegdist(), capscale(), adonis2(), betadisper(), ordiplot()

# ---- 2. Load pre-decontam objects and label the processing day --------------

# Sample ID prefix encodes treatment: c = bagged, p = unbagged. Verified against
# the flower table -- all 147 bagged start with "c", all 147 unbagged with "p".
# extraction_batch is empty for mocks and PCR blanks, so the prefix is the only
# label that works for every control type.
ps_list <- MARKERS |>
  imap(\(obj, marker) {
    ps <- readRDS(here("results", "rds", paste0(obj, ".rds")))
    keep <- sample_names(ps)[
      as(sample_data(ps), "data.frame")$sample_type %in% names(CONTROL_TYPES)
    ]
    if (marker == "16S") keep <- setdiff(keep, FAILED_16S)
    ps <- prune_samples(keep, ps)
    sample_data(ps)$day_type <- ifelse(
      substr(sample_names(ps), 1, 1) == "c", "bagged_day", "unbagged_day")
    ps
  })

walk2(ps_list, names(ps_list), \(ps, marker) {
  md <- as(sample_data(ps), "data.frame")
  cat(sprintf("%s: %d controls\n", marker, nsamples(ps)))
  print(table(md$sample_type, md$day_type))
})

# ---- 3. One analysis: filter -> rCLR distance -> PERMANOVA + PERMDISP -------

analyse_controls <- function(ps, ctrl_type, use_filter, depth_min) {
  ps_sub <- prune_samples(
    sample_names(ps)[
      as(sample_data(ps), "data.frame")$sample_type == ctrl_type &
        sample_sums(ps) >= depth_min
    ], ps)

  # Prevalence is computed WITHIN this control subset. Inheriting the flower
  # analysis threshold would select taxa abundant in flowers, not in controls.
  if (use_filter) {
    ps_sub <- prune_taxa(
      rowSums(otu_table(ps_sub) > 0) / nsamples(ps_sub) >= PREV_THRESH, ps_sub)
  }
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)

  meta <- as(sample_data(ps_sub), "data.frame")
  # Carried so a degenerate panel can report why it is empty.
  meta$depth         <- as.integer(sample_sums(ps_sub)[rownames(meta)])
  meta$observed_taxa <- as.integer(colSums(otu_table(ps_sub) > 0)[rownames(meta)])
  out  <- tibble(n = nsamples(ps_sub), n_taxa = ntaxa(ps_sub),
                 n_bagged   = sum(meta$day_type == "bagged_day"),
                 n_unbagged = sum(meta$day_type == "unbagged_day"),
                 R2 = NA_real_, p_value = NA_real_, permdisp_p = NA_real_)

  # Degenerate cases: a low-biomass blank set can lose almost every taxon to the
  # prevalence filter. Return NA rather than a meaningless test.
  if (out$n_taxa < 3 || min(out$n_bagged, out$n_unbagged) < 3) {
    return(list(stats = out, dist = NULL, meta = meta))
  }

  # impute = FALSE is essential. vegan >= 2.7 defaults rclr to impute = TRUE, which
  # runs OptSpace matrix completion whenever zeros are present. That is full RPCA
  # (as in deicode), not the plain rCLR this project documents, and on small sparse
  # control sets it fails outright with a singular system. Passing impute = FALSE
  # keeps zeros as zeros and matches the method described in 08a.
  d <- tryCatch(
    vegdist(t(as(otu_table(ps_sub), "matrix")),
            method = "robust.aitchison", impute = FALSE),
    error = function(e) {
      message(sprintf("    distance failed (%s): %s", ctrl_type, conditionMessage(e)))
      NULL
    })
  if (is.null(d) || !all(is.finite(as.matrix(d)))) {
    return(list(stats = out, dist = NULL, meta = meta))
  }

  # Permutations are NOT blocked by farm, unlike 10a. Nothing from a farm enters a
  # blank -- it holds reagent and ambient contamination only -- so there is nothing
  # to block. Blocking would also cut the permutation space to 2^7 = 128 and raise
  # the p-value floor to 0.008.
  ad <- adonis2(d ~ day_type, data = meta, permutations = N_PERM)
  bd <- betadisper(d, factor(meta$day_type))

  out$R2         <- ad$R2[1]
  out$p_value    <- ad$`Pr(>F)`[1]
  out$permdisp_p <- permutest(bd, permutations = N_PERM)$tab$`Pr(>F)`[1]

  list(stats = out, dist = d, meta = meta)
}

# ---- 4. Sensitivity grid: filter on/off x depth threshold on/off ------------

grid_results <- expand_grid(
    marker    = names(ps_list),
    ctrl_type = names(CONTROL_TYPES),
    use_filter = c(TRUE, FALSE),
    depth_min  = c(0, DEPTH_MIN)
  ) |>
  pmap_dfr(\(marker, ctrl_type, use_filter, depth_min) {
    cat(sprintf("  %-5s %-22s filter=%-5s depth>=%-5d\n",
                marker, ctrl_type, use_filter, depth_min))
    analyse_controls(ps_list[[marker]], ctrl_type, use_filter, depth_min)$stats |>
      mutate(marker = marker, ctrl_type = ctrl_type,
             use_filter = use_filter, depth_min = depth_min, .before = 1)
  })

dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)
write_csv(grid_results, here("results", "tables", "09d_controls_permanova.csv"))

# Per-control depth and observed richness, so a reader can judge data quality.
# sample_data() already carries a column called sample_id, so the row names go into
# a differently named column to avoid a collision.
imap_dfr(ps_list, \(ps, marker) {
    md <- as(sample_data(ps), "data.frame")
    tibble(marker        = marker,
           control_id    = rownames(md),
           sample_type   = md$sample_type,
           day_type      = md$day_type,
           depth         = as.integer(sample_sums(ps)[rownames(md)]),
           observed_taxa = as.integer(colSums(otu_table(ps) > 0)[rownames(md)]))
  }) |>
  arrange(marker, sample_type, control_id) |>
  write_csv(here("results", "tables", "09d_controls_depth_richness.csv"))

# ---- 5. Figure 1: ordinations, 2 markers x 3 control types -----------------
# Primary setting: prevalence filter on, no depth threshold.

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

png(here("results", "figures", "09d_controls_rCLR_pcoa_by_day_type.png"),
    width = 14, height = 9, units = "in", res = 300)

# mfrow: 2 marker rows x 3 control-type columns
# oma[1]=4: bottom outer margin for the legend band
# oma[2]=4: left outer margin for the vertical row labels
# oma[3]=3: top outer margin for the overall title
par(mfrow = c(2, length(CONTROL_TYPES)),
    mar   = c(3.2, 3.2, 3.6, 0.6),
    oma   = c(4, 4, 3, 0))

walk(names(ps_list), \(marker) {
  walk(names(CONTROL_TYPES), \(ctrl) {
    res <- analyse_controls(ps_list[[marker]], ctrl, use_filter = TRUE, depth_min = 0)

    # A blank so clean that samples share almost no taxa cannot yield a distance
    # matrix. That is a result, not a failure: report the richness that caused it.
    if (is.null(res$dist)) {
      plot.new()
      box(col = "grey80")
      title(main = sprintf("%s | %s\nnot testable", CONTROL_TYPES[[ctrl]], marker),
            cex.main = 0.78)
      text(0.5, 0.55,
           sprintf("median %.0f taxa and %.0f reads per control\n(too clean to compare)",
                   median(res$meta$observed_taxa), median(res$meta$depth)),
           cex = 0.85, col = "grey30")
      return(invisible(NULL))
    }

    ord <- capscale(res$dist ~ 1)
    pct <- round(summary(ord)$cont$importance["Proportion Explained", 1:2] * 100, 1)
    st  <- res$stats

    # mgp = c(axis title, tick labels, tick marks) in lines from the panel edge.
    # The default c(3, 1, 0) puts the axis title outside the 3.2-line margin and it
    # gets clipped; 1.9 pulls it back in. Same fix as in 10d.
    # cex.main 0.78 keeps the four-line title inside the 3.6-line top margin.
    par(mgp = c(1.9, 0.5, 0))
    ordiplot(ord, display = "sites", type = "n",
             xlab = sprintf("PCoA1 (%.1f%%)", pct[1]),
             ylab = sprintf("PCoA2 (%.1f%%)", pct[2]),
             main = sprintf("%s | %s\nn = %d (%d vs %d), %d taxa\nPERMANOVA R2 = %.3f, p = %.4f\nPERMDISP p = %.3f",
                            CONTROL_TYPES[[ctrl]], marker, st$n, st$n_bagged,
                            st$n_unbagged, st$n_taxa, st$R2, st$p_value, st$permdisp_p),
             cex.main = 0.78, cex.axis = 0.7, cex.lab = 0.8)

    # 10d draws 2-sigma ellipses, but those assume ~21 points per group. On seven
    # they balloon far outside the panel and hide the data. A spider plus centroid
    # shows the same thing better here, and shows exactly what PERMANOVA tests:
    # the distance between group centroids relative to the within-group spread.
    ordispider(ord, groups = res$meta$day_type, col = DAY_COLS, lwd = 0.8,
               label = FALSE)
    points(scores(ord, display = "sites"), pch = 16, cex = 1.1,
           col = DAY_COLS[res$meta$day_type])
    # Centroids as an X-cross, matching the convention in 10d.
    local({
      sc <- scores(ord, display = "sites")
      ctr <- tapply(seq_len(nrow(sc)), res$meta$day_type,
                    \(i) colMeans(sc[i, , drop = FALSE]))
      iwalk(ctr, \(xy, grp) points(t(xy), pch = 4, cex = 1.8, lwd = 2.5,
                                   col = DAY_COLS[grp]))
    })
  })
})

# Row labels in the outer left margin. Explicit at= in normalised device coords is
# required; mtext(outer = TRUE) without at= ignores outer margins and drifts.
# Formula mirrors 10d.
local({
  omi_h   <- par("omi")
  din_h   <- par("din")[2]
  inner_h <- din_h - omi_h[1] - omi_h[3]
  walk2(names(ps_list), seq_along(ps_list), \(marker, i) {
    mtext(paste0(marker, " (rCLR Aitchison)"), side = 2, outer = TRUE, las = 0,
          at = (omi_h[1] + inner_h * (1 - (i - 0.5) / length(ps_list))) / din_h,
          cex = 0.8, line = 1.5)
  })
})

mtext("Processing controls by day type — is there batch drift between bagged-day and unbagged-day batches?",
      outer = TRUE, cex = 1.0, line = 0.8, side = 3)

# Legend band: reserves the bottom 8% of the device. Wording is deliberate -- these
# points are CONTROLS, and the colour is the day they were processed alongside, not
# the sample type. A short "bagged / unbagged" label would be read as flowers.
par(fig = c(0, 1, 0, 0.08), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottom",
       legend = c("processed alongside bagged flowers",
                  "processed alongside unbagged flowers"),
       col = DAY_COLS, pch = 16, lty = 1, horiz = TRUE,
       cex = 0.95, bty = "n", pt.cex = 1.3)

dev.off()

# ---- 6. Figure 2: PERMANOVA / PERMDISP heatmap, following 10a --------------
# Same visual grammar as 10a_16S_perm_farm_heatmap.png:
#   fill = PERMANOVA R2, label = PERMANOVA p, bold border = PERMDISP p < 0.05

p_heat <- grid_results |>
  mutate(
    row = factor(paste(marker, CONTROL_TYPES[ctrl_type], sep = " | "),
                 levels = rev(as.vector(outer(names(ps_list), CONTROL_TYPES,
                                              paste, sep = " | ")))),
    col = factor(sprintf("filter %s\ndepth >= %d",
                         ifelse(use_filter, "on", "off"), depth_min)),
    permdisp_sig = !is.na(permdisp_p) & permdisp_p < 0.05,
    p_label = case_when(
      is.na(p_value)   ~ "n/a",
      p_value < 0.001  ~ "<0.001 ***",
      p_value < 0.01   ~ paste0(sprintf("%.3f", p_value), " **"),
      p_value < 0.05   ~ paste0(sprintf("%.3f", p_value), " *"),
      TRUE             ~ sprintf("%.3f", p_value)
    )
  ) |>
  ggplot(aes(x = col, y = row)) +
  # Layer 1: tile fill by PERMANOVA R2 (thin grey grid between tiles)
  geom_tile(aes(fill = R2), colour = "grey80", linewidth = 0.3) +
  # Layer 2: bold border where dispersions differ, so a significant PERMANOVA
  # driven by spread rather than by a centroid shift is visible at a glance
  geom_tile(data = \(d) filter(d, permdisp_sig),
            fill = NA, colour = "black", linewidth = 1.2) +
  geom_text(aes(label = p_label), size = 2.6) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = expression(R^2), limits = c(0, NA),
                       na.value = "grey92") +
  labs(
    title    = "Batch drift between bagged-day and unbagged-day processing batches",
    subtitle = sprintf("adonis2(d ~ day_type, nperm = %d) on rCLR Aitchison distances\nFill = R²; label = PERMANOVA p; bold border = PERMDISP p < 0.05; grey = too few taxa to test",
                       N_PERM),
    x = "Sensitivity setting", y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(here("results", "figures", "09d_controls_perm_heatmap.png"),
       plot = p_heat, width = 9, height = 6, dpi = 300)

cat("\nDone. Wrote 2 figures and 2 tables.\n")
print(grid_results |> filter(use_filter, depth_min == 0) |>
        select(marker, ctrl_type, n, n_taxa, R2, p_value, permdisp_p))
