## ============================================================================
## 16b — ANCOM-BC2: ITS1 OTU97 fungi — insect footprint
##
## Contrast: bagged_flower (reference) vs unbagged_flower
## All flowers; farm_id as fixed covariate; tree_id as random blocking factor.
##
## ANCOM-BC2 complements ISA (13b): ISA finds taxon-level indicators;
## ANCOM-BC2 estimates log2 fold-change with uncertainty on the same contrast.
##
## Input:
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output:
##   results/figures/16b_ancombc2_insect_footprint.png
##   results/tables/16b_ancombc2_insect_footprint.csv
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

FDR_ALPHA   <- 0.05
PREV_THRESH <- 0.10
MIN_LIB     <- 500

RANK_PLOT <- "Genus"

# ── 1. Libraries & data ───────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(ANCOMBC)
library(ggrepel)

# Taxonomy labels: deepest assigned rank, prefix kept, shared with the
# other ANCOM-BC2 scripts so every figure names a taxon the same way.
source(here("scripts", "utils_tax_label.R"))

ps_raw <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds"))

cat(sprintf("Loaded: %d taxa x %d samples\n", ntaxa(ps_raw), nsamples(ps_raw)))

# ── 2. Helper functions ───────────────────────────────────────────────────────

filter_ps <- function(physeq,
                      prev = PREV_THRESH,
                      min_lib = MIN_LIB) {
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

# Volcano: log2FC (x) vs -log10(q_BH) (y).
# res must have a tax_label column added before calling.
# x_label format "log2FC (numerator / denominator)" drives legend text.
plot_volcano <- function(res, lfc_col, q_col, diff_col, fdr_col = NULL,
                         x_label = "log2FC",
                         title   = "") {

        if (is.null(res)) {
                return(ggplot() + theme_void() +
                               annotate("text", x = 0.5, y = 0.5, hjust = 0.5,
                                        label = "Analysis failed", size = 5) +
                               labs(title = title))
        }

        df <- res |>
                dplyr::mutate(
                        lfc        = .data[[lfc_col]],
                        q          = pmax(.data[[q_col]], 1e-300),
                        signif     = as.logical(.data[[diff_col]]),
                        neg_log10q = -log10(q),
                        direction  = dplyr::case_when(
                                signif & lfc > 0 ~ "higher_num",
                                signif & lfc < 0 ~ "higher_den",
                                TRUE             ~ "ns"
                        )
                )

        n_num <- sum(df$direction == "higher_num", na.rm = TRUE)
        n_den <- sum(df$direction == "higher_den", na.rm = TRUE)

        # Taxa passing FDR alone, before the sensitivity filter -- shown so that the
        # drop is visible on the figure rather than looking like missing points.
        n_fdr <- if (is.null(fdr_col)) NA_integer_
                 else sum(as.logical(res[[fdr_col]]), na.rm = TRUE)

        m         <- regmatches(x_label, regexec("\\((.+?)\\s*/\\s*(.+?)\\)", x_label))[[1]]
        num_label <- if (length(m) >= 3) trimws(m[2]) else "numerator"
        den_label <- if (length(m) >= 3) trimws(m[3]) else "denominator"

        label_df <- df |> dplyr::filter(signif) |>
                dplyr::slice_max(abs(lfc), n = 15)

        ggplot(df, aes(x = lfc, y = neg_log10q)) +
                geom_point(aes(colour = direction), size = 1.8, alpha = 0.75) +
                geom_hline(yintercept = -log10(FDR_ALPHA),
                           linetype = "dashed", colour = "grey50", linewidth = 0.4) +
                geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
                ggrepel::geom_text_repel(
                        data           = label_df,
                        aes(label = tax_label),
                        size           = 2.8, max.overlaps = 20,
                        segment.colour = "grey60"
                ) +
                scale_colour_manual(
                        values = c(higher_num = "#D6604D",
                                   higher_den = "#4393C3",
                                   ns         = "#CCCCCC"),
                        labels = c(
                                higher_num = sprintf("↑ in %s (n=%d)", num_label, n_num),
                                higher_den = sprintf("↑ in %s (n=%d)", den_label, n_den),
                                ns         = sprintf("Not significant (n=%d)",
                                                     sum(df$direction == "ns"))
                        ),
                        name = NULL
                ) +
                labs(
                        title    = title,
                        subtitle = sprintf(
                                "ANCOM-BC2  |  FDR < %.2f  |  prevalence ≥ %.0f%%  |  sig = %d of %d passing FDR",
                                FDR_ALPHA, 100 * PREV_THRESH, n_num + n_den, n_fdr
                        ),
                        x = x_label,
                        y = expression(-log[10](q[BH]))
                ) +
                theme_bw(base_size = 11) +
                theme(
                        legend.position = "bottom",
                        plot.title      = element_text(size = 12, face = "bold"),
                        plot.subtitle   = element_text(size = 9,  colour = "grey40")
                )
}

# ── 3. Build contrast phyloseq ────────────────────────────────────────────────

cat("\n── 3. Insect footprint contrast (bagged vs unbagged) ───────────────────\n")

ps_f <- subset_samples(
        ps_raw,
        sample_data(ps_raw)$sample_type %in% c("bagged_flower", "unbagged_flower")
)

sample_data(ps_f)$sample_type <- factor(
        sample_data(ps_f)$sample_type,
        levels = c("bagged_flower", "unbagged_flower")
)

n_bagged   <- sum(sample_data(ps_f)$sample_type == "bagged_flower")
n_unbagged <- sum(sample_data(ps_f)$sample_type == "unbagged_flower")
cat(sprintf("  bagged: %d  |  unbagged: %d\n", n_bagged, n_unbagged))

ps_f <- filter_ps(ps_f)
cat(sprintf("  After prevalence/library filter: %d taxa x %d samples\n",
            ntaxa(ps_f), nsamples(ps_f)))

# ── 4. Run ANCOM-BC2 ──────────────────────────────────────────────────────────

cat("  Running ancombc2()...\n")
out <- run_ancombc2(
        ps_f,
        fix_formula  = "sample_type + farm_id",
        rand_formula = "(1|tree_id)"
)

# ── 5. Build volcano plot ─────────────────────────────────────────────────────

res <- if (!is.null(out)) out$res else NULL

lfc_col  <- find_col(res, "^lfc.*unbagged")
q_col    <- find_col(res, "^q.*unbagged")
# Significance requires BOTH the FDR call and ANCOM-BC2's own sensitivity check.
# res carries three parallel columns per model term: diff_*, passed_ss_* and
# diff_robust_*, the last being exactly the conjunction of the first two (verified on
# this project's data). passed_ss marks taxa whose call survives adding pseudo-counts
# (0.1, 0.5, 1) to zero counts; the package vignette treats taxa that fail it as not
# significant, and doing so lowers the false discovery rate.
#
# find_col() returns the first match and ANCOM-BC2 emits the diff_* block before
# diff_robust_*, so "^diff" still selects the plain FDR column; the stopifnot() below
# turns that ordering assumption into a loud failure rather than a silent one.
fdr_col  <- find_col(res, "^diff.*unbagged")
ss_col   <- find_col(res, "^passed_ss.*unbagged")
diff_col <- find_col(res, "^diff_robust.*unbagged")

if (!is.null(res)) {
        cat(sprintf("  LFC column: %s\n", lfc_col))
        stopifnot(fdr_col != diff_col, is.na(ss_col) == FALSE)
        cat(sprintf("  Significant taxa: %d (FDR) -> %d after sensitivity filter\n",
                    sum(as.logical(res[[fdr_col]]), na.rm = TRUE),
                    sum(as.logical(res[[diff_col]]), na.rm = TRUE)))
        res$tax_label <- make_tax_label(ps_raw, res$taxon)
}

p <- plot_volcano(
        res, lfc_col, q_col, diff_col, fdr_col,
        x_label = "log2FC (unbagged / bagged)",
        title   = "ITS1 OTU97 Fungi — Insect footprint (bagged vs unbagged)"
)

# ── 6. Save outputs ───────────────────────────────────────────────────────────

cat("\n── 6. Saving outputs ────────────────────────────────────────────────────\n")

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

ggsave(
        here("results", "figures", "16b_ancombc2_insect_footprint.png"),
        plot   = p,
        width  = 10,
        height = 7,
        dpi    = 300
)
cat("  Saved: 16b_ancombc2_insect_footprint.png\n")

res |>
        dplyr::select(taxon, tax_label,
                      lfc       = dplyr::all_of(lfc_col),
                      q         = dplyr::all_of(q_col),
                      diff_fdr  = dplyr::all_of(fdr_col),   # FDR alone
                      passed_ss = dplyr::all_of(ss_col),    # survived pseudo-count check
                      diff      = dplyr::all_of(diff_col)) |>  # both -- used for figures
        dplyr::arrange(q) |>
        readr::write_csv(here("results", "tables",
                              "16b_ancombc2_insect_footprint.csv"))
cat("  Saved: 16b_ancombc2_insect_footprint.csv\n")

cat(sprintf("\n16b done. %d significant taxa (FDR < %.2f + sensitivity filter).\n",
            sum(as.logical(res[[diff_col]]), na.rm = TRUE), FDR_ALPHA))
