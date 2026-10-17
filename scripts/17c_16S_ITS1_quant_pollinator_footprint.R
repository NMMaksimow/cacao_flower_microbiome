## ============================================================================
## 17c — ANCOM-BC2: quantitative pollinator footprint (16S + ITS1 combined)
##
## Predictor: pi_clsm (continuous pollen intensity index, CLSM-based)
## Subset: unbagged flowers only; pi_clsm = 0 retained (true zero, no pollen
## detected on examination); NA excluded (pi_clsm not measured).
##
## fix_formula = "pi_clsm + farm_id": farm_id absorbs between-farm baseline;
## pi_clsm coefficient is the within-farm dose-response slope.
## rand_formula omitted: ~120 samples across 49 unique tree_id gives ~2.4
## obs/level; random intercept cannot be estimated.
##
## Complements 17a/17b (binary contrast at extremes) by using the full
## pollen gradient without discarding the grey zone (pi_clsm 1–39).
##
## Input:
##   results/rds/ps_16S_bacteria_biosamples.rds
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output:
##   results/figures/17c_quant_pollinator_volcano.png
##   results/figures/17c_quant_pollinator_scatter_16S.png
##   results/figures/17c_quant_pollinator_scatter_ITS1.png
##   results/figures/17c_quant_pollinator_sensitivity.png
##   results/tables/17c_quant_pollinator_footprint.csv
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

FDR_ALPHA   <- 0.05
PREV_THRESH <- 0.10
MIN_LIB     <- 500
N_SCATTER   <- 6      # top N significant taxa shown in scatter plots per marker

# Leverage sensitivity. pi_clsm is heavily right-skewed: of 102 unbagged flowers 38
# score zero, the median is 2, and the maximum is 349. In the design matrix
# (~ pi_clsm + farm_id, p = 8) the single most pollinated flower carries 34.5% of
# sum((x - xbar)^2), the top three carry 62.4% and the top five 76.3%; its leverage is
# 0.38 against a mean of p/n = 0.078. The slope is therefore estimated from a handful
# of flowers, so the whole model is refitted without them and the two fits compared.
PI_TRIM     <- 100    # sensitivity refit drops flowers at or above this pollen intensity

# ── 1. Libraries & data ───────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(ANCOMBC)
library(ggrepel)
library(patchwork)

# Taxonomy labels: deepest assigned rank, prefix kept, shared with the
# other ANCOM-BC2 scripts so every figure names a taxon the same way.
source(here("scripts", "utils_tax_label.R"))

ps_16s <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds"))
ps_its <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds"))

cat(sprintf("16S loaded:  %d taxa x %d samples\n", ntaxa(ps_16s), nsamples(ps_16s)))
cat(sprintf("ITS1 loaded: %d taxa x %d samples\n", ntaxa(ps_its), nsamples(ps_its)))

# ── 2. Helper functions ───────────────────────────────────────────────────────

filter_ps <- function(physeq, prev = PREV_THRESH, min_lib = MIN_LIB) {
        physeq <- prune_samples(sample_sums(physeq) >= min_lib, physeq)
        keep   <- rowSums(otu_table(physeq) > 0) / nsamples(physeq) >= prev
        prune_taxa(keep, physeq)
}

find_col <- function(res, pattern) grep(pattern, names(res), value = TRUE)[1]

run_ancombc2 <- function(physeq, fix_formula, rand_formula = NULL) {
        args <- list(
                data         = physeq,
                fix_formula  = fix_formula,
                p_adj_method = "BH",
                alpha        = FDR_ALPHA,
                verbose      = FALSE
        )
        if (!is.null(rand_formula)) args$rand_formula <- rand_formula
        tryCatch(do.call(ancombc2, args), error = function(e) {
                message("  ancombc2() failed: ", conditionMessage(e)); NULL
        })
}

# CLR with pseudocount 0.5 for visualization only.
# True rCLR (used in PCoA) cannot be computed per-taxon when zero counts are
# present (zero → -Inf after log). Returns taxa × samples matrix.
clr_vis <- function(ps_sub, taxon_ids) {
        otu_mat <- as(otu_table(ps_sub), "matrix")
        if (!taxa_are_rows(ps_sub)) otu_mat <- t(otu_mat)
        all_log    <- log(otu_mat + 0.5)
        gmeans     <- colMeans(all_log)             # geometric mean per sample
        focal_log  <- log(otu_mat[taxon_ids, , drop = FALSE] + 0.5)
        sweep(focal_log, 2, gmeans, "-")            # taxa × samples
}

# Slope volcano for one marker.
# Uses static legend labels so guides can be collected across the 2-panel patchwork.
plot_slope_volcano <- function(res, lfc_col, q_col, diff_col, label_fn, ps_raw,
                                title = "", fdr_col = NULL) {
        if (is.null(res)) {
                return(
                        ggplot() + theme_void() +
                                annotate("text", x = 0.5, y = 0.5, hjust = 0.5,
                                         label = "Analysis failed", size = 5) +
                                labs(title = title)
                )
        }

        res$tax_label <- label_fn(ps_raw, res$taxon)

        df <- res |>
                dplyr::mutate(
                        lfc        = .data[[lfc_col]],
                        q          = pmax(.data[[q_col]], 1e-300),
                        signif     = as.logical(.data[[diff_col]]),
                        neg_log10q = -log10(q),
                        direction  = dplyr::case_when(
                                signif & lfc > 0 ~ "pos_slope",
                                signif & lfc < 0 ~ "neg_slope",
                                TRUE             ~ "ns"
                        )
                )

        n_sig <- sum(df$signif, na.rm = TRUE)
        # Taxa passing FDR alone, shown alongside so the drop caused by the
        # sensitivity filter is visible rather than looking like missing points.
        n_fdr <- if (is.null(fdr_col)) NA_integer_
                 else sum(as.logical(res[[fdr_col]]), na.rm = TRUE)
        label_df <- df |> dplyr::filter(signif) |> dplyr::slice_max(abs(lfc), n = 15)

        ggplot(df, aes(x = lfc, y = neg_log10q)) +
                geom_point(aes(colour = direction), size = 1.8, alpha = 0.75) +
                geom_hline(yintercept = -log10(FDR_ALPHA),
                           linetype = "dashed", colour = "grey50", linewidth = 0.4) +
                geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
                ggrepel::geom_text_repel(
                        data = label_df, aes(label = tax_label),
                        size = 2.8, max.overlaps = 20, segment.colour = "grey60"
                ) +
                scale_colour_manual(
                        values = c(pos_slope = "#D6604D",
                                   neg_slope = "#4393C3",
                                   ns        = "#CCCCCC"),
                        # static labels so guides can be collected in the patchwork
                        labels = c(pos_slope = "higher with pollen",
                                   neg_slope = "lower with pollen",
                                   ns        = "not significant"),
                        breaks = c("pos_slope", "neg_slope", "ns"),
                        drop   = FALSE,
                        name   = NULL
                ) +
                labs(
                        title    = title,
                        subtitle = sprintf(
                                "ANCOM-BC2  |  FDR < %.2f  |  prevalence ≥ %.0f%%  |  sig = %d of %d passing FDR",
                                FDR_ALPHA, 100 * PREV_THRESH, n_sig, n_fdr
                        ),
                        x = "regression slope (log2 per pi_clsm unit)",
                        y = expression(-log[10](q[BH]))
                ) +
                theme_bw(base_size = 11) +
                theme(
                        legend.position = "bottom",
                        plot.title      = element_text(size = 12, face = "bold"),
                        plot.subtitle   = element_text(size = 9,  colour = "grey40")
                )
}

# Build long data frame for scatter plots: top N significant taxa by |slope|.
build_scatter_df <- function(ps_sub, res, lfc_col, q_col, diff_col,
                              label_fn, n = N_SCATTER) {
        sig_res <- res |>
                dplyr::filter(as.logical(.data[[diff_col]])) |>
                dplyr::slice_max(abs(.data[[lfc_col]]), n = n)

        if (nrow(sig_res) == 0) return(NULL)

        taxon_ids <- sig_res$taxon
        clr_mat   <- clr_vis(ps_sub, taxon_ids)    # taxa × samples

        sdata <- as(sample_data(ps_sub), "data.frame")[colnames(clr_mat),
                                                        c("pi_clsm", "farm_id")]
        sdata$farm_id <- factor(sdata$farm_id)

        # Partial (added-variable) coordinates. The model estimates the pi_clsm effect
        # WITHIN farms, so plotting raw pi_clsm against raw CLR shows the marginal
        # trend instead -- which for several taxa runs opposite to the fitted slope.
        # Regressing both axes on farm_id first removes the between-farm component;
        # the OLS slope through the residuals then equals, exactly, the pi_clsm
        # coefficient of lm(clr ~ pi_clsm + farm_id). It still will not equal the
        # ANCOM-BC2 slope numerically (that model works on bias-corrected log
        # abundances, not CLR with a pseudocount), but the sign now agrees.
        pi_resid <- residuals(lm(pi_clsm ~ farm_id, data = sdata))

        dplyr::bind_rows(lapply(seq_along(taxon_ids), function(i) {
                clr_i <- as.numeric(clr_mat[i, ])
                dplyr::tibble(
                        taxon     = taxon_ids[i],
                        tax_label = label_fn(ps_sub, taxon_ids[i]),
                        slope     = sig_res[[lfc_col]][i],
                        q_val     = sig_res[[q_col]][i],
                        clr       = clr_i,
                        clr_resid = residuals(lm(clr_i ~ sdata$farm_id)),
                        pi_clsm   = sdata$pi_clsm,
                        pi_resid  = pi_resid,
                        farm_id   = as.character(sdata$farm_id)
                )
        }))
}

# Assemble 3-column patchwork of partial-regression (added-variable) plots,
# coloured by farm_id. Both axes carry residuals after regressing on farm_id, so the
# fitted line reproduces the within-farm slope the model actually estimates. Zero on
# either axis means "average for that farm", not "no pollen" / "no abundance".
plot_scatter_panels <- function(scatter_df) {
        if (is.null(scatter_df) || nrow(scatter_df) == 0) return(NULL)

        # Order panels: highest |slope| first
        tax_order <- scatter_df |>
                dplyr::distinct(tax_label, slope) |>
                dplyr::arrange(dplyr::desc(abs(slope))) |>
                dplyr::pull(tax_label)

        panels <- lapply(tax_order, function(lbl) {
                d <- dplyr::filter(scatter_df, tax_label == lbl)
                ggplot(d, aes(x = pi_resid, y = clr_resid, colour = farm_id)) +
                        geom_point(size = 1.5, alpha = 0.7) +
                        geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
                        geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
                        geom_smooth(method = "lm", se = TRUE,
                                    colour = "grey30", fill = "grey70",
                                    linewidth = 0.7, alpha = 0.3,
                                    aes(group = 1)) +  # single line across all farms
                        labs(
                                title    = lbl,
                                subtitle = sprintf("slope = %.4f  |  q = %.3g",
                                                   unique(d$slope), unique(d$q_val)),
                                x        = "Pollen intensity | farm (residual)",
                                y        = "CLR abundance | farm (residual)",
                                colour   = "Farm"
                        ) +
                        theme_bw(base_size = 9) +
                        theme(
                                plot.title    = element_text(size = 9, face = "bold.italic"),
                                plot.subtitle = element_text(size = 7, colour = "grey40"),
                                legend.position = "none"  # collected at patchwork level
                        )
        })

        n_panels <- length(panels)
        ncols    <- min(3L, n_panels)

        wrap_plots(panels, ncol = ncols) +
                plot_layout(guides = "collect") &
                theme(legend.position = "bottom")
}

# ── 3. Build subset phyloseq ──────────────────────────────────────────────────

cat("\n── 3. Subsetting to unbagged flowers with measured pi_clsm ─────────────\n")

build_subset <- function(ps_raw, label) {
        sdata <- as(sample_data(ps_raw), "data.frame")
        keep  <- sdata$sample_type == "unbagged_flower" & !is.na(sdata$pi_clsm)
        ps_sub <- prune_samples(keep, ps_raw) |> filter_ps()
        pi_rng <- range(sample_data(ps_sub)$pi_clsm, na.rm = TRUE)
        cat(sprintf("  %s: %d samples  pi_clsm [%.0f, %.0f]  taxa_after_filter=%d\n",
                    label, nsamples(ps_sub), pi_rng[1], pi_rng[2], ntaxa(ps_sub)))
        ps_sub
}

ps_16s_sub <- build_subset(ps_16s, "16S")
ps_its_sub <- build_subset(ps_its, "ITS1")

# ── 4. Run ANCOMBC2 ───────────────────────────────────────────────────────────

cat("\n── 4. Running ANCOMBC2 ──────────────────────────────────────────────────\n")
cat("  16S ...\n")
out_16s <- run_ancombc2(ps_16s_sub, fix_formula = "pi_clsm + farm_id")

cat("  ITS1 ...\n")
out_its <- run_ancombc2(ps_its_sub, fix_formula = "pi_clsm + farm_id")

# ── 5. Extract results ────────────────────────────────────────────────────────

res_16s <- if (!is.null(out_16s)) out_16s$res else NULL
res_its <- if (!is.null(out_its)) out_its$res else NULL

# Column names for continuous pi_clsm: lfc_pi_clsm, se_pi_clsm, q_pi_clsm, diff_pi_clsm
#
# Significance requires BOTH the FDR call and ANCOM-BC2's own sensitivity check.
# res carries three parallel columns per model term: diff_*, passed_ss_* and
# diff_robust_*, the last being exactly the conjunction of the first two (verified on
# this project's data). passed_ss marks taxa whose call survives adding pseudo-counts
# (0.1, 0.5, 1) to zero counts; the package vignette treats taxa that fail it as not
# significant. This matters most here: the unfiltered run called half the bacterial
# community differentially abundant along a covariate that the binary contrast (17a)
# and the random forest cannot detect at all.
#
# find_col() returns the first match and ANCOM-BC2 emits the diff_* block before
# diff_robust_*, so "^diff" still selects the plain FDR column; the stopifnot() below
# turns that ordering assumption into a loud failure rather than a silent one.
lfc_col_16s  <- find_col(res_16s, "^lfc.*pi_clsm")
se_col_16s   <- find_col(res_16s, "^se.*pi_clsm")
q_col_16s    <- find_col(res_16s, "^q.*pi_clsm")
fdr_col_16s  <- find_col(res_16s, "^diff.*pi_clsm")
ss_col_16s   <- find_col(res_16s, "^passed_ss.*pi_clsm")
diff_col_16s <- find_col(res_16s, "^diff_robust.*pi_clsm")

lfc_col_its  <- find_col(res_its, "^lfc.*pi_clsm")
se_col_its   <- find_col(res_its, "^se.*pi_clsm")
q_col_its    <- find_col(res_its, "^q.*pi_clsm")
fdr_col_its  <- find_col(res_its, "^diff.*pi_clsm")
ss_col_its   <- find_col(res_its, "^passed_ss.*pi_clsm")
diff_col_its <- find_col(res_its, "^diff_robust.*pi_clsm")

if (!is.null(res_16s)) {
        stopifnot(fdr_col_16s != diff_col_16s, is.na(ss_col_16s) == FALSE)
        cat(sprintf("  16S:  significant taxa = %d (FDR) -> %d after sensitivity filter\n",
                    sum(as.logical(res_16s[[fdr_col_16s]]), na.rm = TRUE),
                    sum(as.logical(res_16s[[diff_col_16s]]), na.rm = TRUE)))
}
if (!is.null(res_its)) {
        stopifnot(fdr_col_its != diff_col_its, is.na(ss_col_its) == FALSE)
        cat(sprintf("  ITS1: significant taxa = %d (FDR) -> %d after sensitivity filter\n",
                    sum(as.logical(res_its[[fdr_col_its]]), na.rm = TRUE),
                    sum(as.logical(res_its[[diff_col_its]]), na.rm = TRUE)))
}

# ── 5b. Leverage sensitivity: refit without the most pollinated flowers ───────

cat(sprintf("\n── 5b. Sensitivity refit, dropping flowers with pi_clsm >= %d ──────────\n",
            PI_TRIM))

# Flowers are pruned from the ALREADY filtered subset rather than the analysis being
# rebuilt from scratch, so the same taxa enter both fits and the comparison is about
# the removed flowers, not about a changed feature set.
trim_subset <- function(ps_sub, label) {
        ps_out <- prune_samples(sample_data(ps_sub)$pi_clsm < PI_TRIM, ps_sub)
        cat(sprintf("  %s: %d of %d flowers kept\n",
                    label, nsamples(ps_out), nsamples(ps_sub)))
        ps_out
}

out_16s_trim <- run_ancombc2(trim_subset(ps_16s_sub, "16S"),
                             fix_formula = "pi_clsm + farm_id")
out_its_trim <- run_ancombc2(trim_subset(ps_its_sub, "ITS1"),
                             fix_formula = "pi_clsm + farm_id")

# One row per taxon: full fit beside trimmed fit. Taxa that ANCOM-BC2's internal
# prevalence cut drops from the trimmed run come back as NA, which counts as "lost".
sens_table <- function(res_full, res_trim, marker) {
        if (is.null(res_full)) return(NULL)
        full <- dplyr::tibble(
                marker = marker,
                taxon  = res_full$taxon,
                slope  = res_full[[find_col(res_full, "^lfc.*pi_clsm")]],
                diff   = as.logical(res_full[[find_col(res_full, "^diff_robust.*pi_clsm")]])
        )
        if (is.null(res_trim)) {
                return(dplyr::mutate(full, slope_trim = NA_real_,
                                     q_trim = NA_real_, diff_trim = NA))
        }
        full |> dplyr::left_join(
                dplyr::tibble(
                        taxon      = res_trim$taxon,
                        slope_trim = res_trim[[find_col(res_trim, "^lfc.*pi_clsm")]],
                        q_trim     = res_trim[[find_col(res_trim, "^q.*pi_clsm")]],
                        diff_trim  = as.logical(
                                res_trim[[find_col(res_trim, "^diff_robust.*pi_clsm")]])
                ), by = "taxon")
}

# `%in% TRUE` rather than isTRUE(): vectorised, and it maps NA to FALSE, which is the
# reading we want for taxa missing from the trimmed fit.
sens <- dplyr::bind_rows(
        sens_table(res_16s, if (is.null(out_16s_trim)) NULL else out_16s_trim$res, "16S"),
        sens_table(res_its, if (is.null(out_its_trim)) NULL else out_its_trim$res, "ITS1")
) |>
        dplyr::mutate(
                still = diff_trim %in% TRUE,
                outcome = dplyr::case_when(
                        diff & still & sign(slope) != sign(slope_trim) ~ "sign flip",
                        diff & still                                   ~ "robust to trimming",
                        diff                                           ~ "lost after trimming",
                        still                                          ~ "gained after trimming",
                        TRUE                                           ~ "not significant"
                )
        )

for (mk in unique(sens$marker)) {
        s <- dplyr::filter(sens, marker == mk)
        cat(sprintf("  %s: %d significant in the full fit -> %s\n", mk, sum(s$diff),
                    paste(sprintf("%s %d", names(table(s$outcome[s$diff])),
                                  as.integer(table(s$outcome[s$diff]))), collapse = ", ")))
}

# Facet labels carry the counts, so the figure reads without the log.
sens <- sens |>
        dplyr::mutate(marker_lab = sprintf(
                "%s — %d significant, %d survive trimming",
                marker,
                ave(diff, marker, FUN = sum),
                ave(diff & still, marker, FUN = sum)))

OUTCOME_COLS <- c("robust to trimming"    = "#1B7837",
                  "lost after trimming"   = "#D6604D",
                  "gained after trimming" = "#4393C3",
                  "sign flip"             = "#762A83",
                  "not significant"       = "#CCCCCC")

fig_sens <- ggplot(sens, aes(x = slope, y = slope_trim, colour = outcome)) +
        # y = x: a point on the diagonal is unchanged by dropping the flowers
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    colour = "grey50", linewidth = 0.4) +
        geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
        geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
        geom_point(size = 1.6, alpha = 0.8) +
        scale_colour_manual(values = OUTCOME_COLS, breaks = names(OUTCOME_COLS),
                            drop = FALSE, name = NULL) +
        facet_wrap(~ marker_lab, scales = "free") +   # slopes differ in scale by marker
        labs(
                title    = "Leverage check: slope with and without the most pollinated flowers",
                subtitle = sprintf(
                        "refit after dropping flowers with pi_clsm >= %d  |  the top three flowers carry 62%% of the predictor variance",
                        PI_TRIM),
                x = "regression slope, all flowers",
                y = sprintf("regression slope, pi_clsm < %d", PI_TRIM)
        ) +
        theme_bw(base_size = 11) +
        theme(
                legend.position = "bottom",
                plot.title      = element_text(size = 12, face = "bold"),
                plot.subtitle   = element_text(size = 9, colour = "grey40"),
                strip.text      = element_text(size = 9, face = "bold")
        )

# ── 6. Slope volcano (2-panel patchwork) ──────────────────────────────────────

cat("\n── 6. Building slope volcano ────────────────────────────────────────────\n")

p_vol_16s <- plot_slope_volcano(
        res_16s, lfc_col_16s, q_col_16s, diff_col_16s,
        make_tax_label, ps_16s,
        title = "16S Bacteria", fdr_col = fdr_col_16s
)

p_vol_its <- plot_slope_volcano(
        res_its, lfc_col_its, q_col_its, diff_col_its,
        make_tax_label, ps_its,
        title = "ITS1 OTU97 Fungi", fdr_col = fdr_col_its
)

fig_volcano <- (p_vol_16s / p_vol_its) +
        plot_layout(guides = "collect") +
        plot_annotation(
                title   = "ANCOMBC2: Quantitative pollinator footprint (pi_clsm continuous)",
                caption = sprintf(
                        "FDR < %.2f  |  prevalence ≥ %.0f%%  |  fix_formula: pi_clsm + farm_id  |  unbagged flowers only",
                        FDR_ALPHA, 100 * PREV_THRESH
                ),
                theme = theme(
                        plot.title   = element_text(size = 13, face = "bold"),
                        plot.caption = element_text(size = 7,  colour = "grey40")
                )
        ) &
        theme(legend.position = "bottom")

# ── 7. Scatter plots for top N taxa per marker ────────────────────────────────

cat(sprintf("\n── 7. Building scatter plots (top %d taxa per marker) ───────────────\n",
            N_SCATTER))

scatter_16s <- build_scatter_df(
        ps_16s_sub, res_16s, lfc_col_16s, q_col_16s, diff_col_16s,
        make_tax_label
)
fig_scatter_16s <- plot_scatter_panels(scatter_16s)

scatter_its <- build_scatter_df(
        ps_its_sub, res_its, lfc_col_its, q_col_its, diff_col_its,
        make_tax_label
)
fig_scatter_its <- plot_scatter_panels(scatter_its)

# ── 8. Save outputs ───────────────────────────────────────────────────────────

cat("\n── 8. Saving outputs ────────────────────────────────────────────────────\n")

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

ggsave(
        here("results", "figures", "17c_quant_pollinator_volcano.png"),
        plot   = fig_volcano,
        width  = 10,
        height = 12,   # 2 panels × 6 in each
        dpi    = 300
)
cat("  Saved: 17c_quant_pollinator_volcano.png\n")

if (!is.null(fig_scatter_16s)) {
        ggsave(
                here("results", "figures", "17c_quant_pollinator_scatter_16S.png"),
                plot   = fig_scatter_16s,
                width  = 12,
                height = 8,
                dpi    = 300
        )
        cat("  Saved: 17c_quant_pollinator_scatter_16S.png\n")
} else {
        cat("  No significant 16S taxa for scatter plots.\n")
}

if (!is.null(fig_scatter_its)) {
        ggsave(
                here("results", "figures", "17c_quant_pollinator_scatter_ITS1.png"),
                plot   = fig_scatter_its,
                width  = 12,
                height = 8,
                dpi    = 300
        )
        cat("  Saved: 17c_quant_pollinator_scatter_ITS1.png\n")
} else {
        cat("  No significant ITS1 taxa for scatter plots.\n")
}

ggsave(
        here("results", "figures", "17c_quant_pollinator_sensitivity.png"),
        plot   = fig_sens,
        width  = 11,
        height = 6,
        dpi    = 300
)
cat("  Saved: 17c_quant_pollinator_sensitivity.png\n")

# ── 9. Save CSV ───────────────────────────────────────────────────────────────

csv_rows <- list()

if (!is.null(res_16s) && !is.na(lfc_col_16s)) {
        csv_rows[[1]] <- dplyr::tibble(
                marker    = "16S",
                taxon     = res_16s$taxon,
                tax_label = make_tax_label(ps_16s, res_16s$taxon),
                slope     = res_16s[[lfc_col_16s]],
                se        = res_16s[[se_col_16s]],
                q         = res_16s[[q_col_16s]],
                diff_fdr  = as.logical(res_16s[[fdr_col_16s]]),   # FDR alone
                passed_ss = as.logical(res_16s[[ss_col_16s]]),    # pseudo-count check
                diff      = as.logical(res_16s[[diff_col_16s]])   # both
        )
}

if (!is.null(res_its) && !is.na(lfc_col_its)) {
        csv_rows[[2]] <- dplyr::tibble(
                marker    = "ITS1",
                taxon     = res_its$taxon,
                tax_label = make_tax_label(ps_its, res_its$taxon),
                slope     = res_its[[lfc_col_its]],
                se        = res_its[[se_col_its]],
                q         = res_its[[q_col_its]],
                diff_fdr  = as.logical(res_its[[fdr_col_its]]),
                passed_ss = as.logical(res_its[[ss_col_its]]),
                diff      = as.logical(res_its[[diff_col_its]])
        )
}

if (length(csv_rows) > 0) {
        dplyr::bind_rows(csv_rows) |>
                # the trimmed refit travels in the same table, so a reader can compare
                # the two fits taxon by taxon without joining files
                dplyr::left_join(dplyr::select(sens, marker, taxon,
                                               slope_trim, q_trim, diff_trim, outcome),
                                 by = c("marker", "taxon")) |>
                dplyr::arrange(marker, q) |>
                readr::write_csv(here("results", "tables",
                                      "17c_quant_pollinator_footprint.csv"))
        cat("  Saved: 17c_quant_pollinator_footprint.csv\n")
} else {
        cat("  No results to save.\n")
}

cat("\n17c done.\n")
