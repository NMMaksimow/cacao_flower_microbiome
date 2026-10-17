## ============================================================================
## 16c — ANCOM-BC2: insect footprint by farm (16S + ITS1 combined)
##
## Contrast: bagged_flower (reference) vs unbagged_flower
## Analysis run within each farm independently; no farm_id covariate.
##
## rand_formula omitted: singular fits across farms indicate between-tree variance
## is effectively 0 (flowers from the same tree do not cluster together in beta
## diversity). Including tree_id reduces power without correcting real pseudoreplication.
##
## Result: 2-row x 7-column patchwork (row 1 = 16S, row 2 = ITS1; columns =
## farms in canonical order), mirroring the PCoA layout used in 09a, 09b and 10d.
##
## Input:
##   results/rds/ps_16S_bacteria_biosamples.rds
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output:
##   results/figures/16c_ancombc2_insect_footprint_by_farm.png
##   results/tables/16c_ancombc2_insect_footprint_by_farm.csv
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

FDR_ALPHA   <- 0.05
PREV_THRESH <- 0.10    # applied to each farm's subset independently
MIN_LIB     <- 500
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

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

# Compact volcano for one farm panel.
# label_fn:     make_tax_label (kept as an argument so the panel builder stays
#               independent of where the label map comes from)
# show_x_label: draw x-axis title (bottom row only)
# show_y_label: draw y-axis title (leftmost column only)
plot_farm_volcano <- function(out, ps_raw, label_fn, farm,
                               show_x_label = TRUE, show_y_label = TRUE) {

        title <- toupper(farm)

        if (is.null(out)) {
                return(
                        ggplot() + theme_void(base_size = 8) +
                                annotate("text", x = 0.5, y = 0.5, hjust = 0.5,
                                         label = "Failed", size = 2.5, colour = "grey50") +
                                labs(title = title) +
                                theme(plot.title = element_text(size = 8, face = "bold",
                                                                hjust = 0.5))
                )
        }

        res      <- out$res
        lfc_col  <- find_col(res, "^lfc.*unbagged")
        q_col    <- find_col(res, "^q.*unbagged")
        # diff_robust_* is ANCOM-BC2's conjunction of diff_* (FDR) and passed_ss_*
        # (survives adding pseudo-counts to zero counts). Colouring by it keeps the
        # per-farm panels consistent with the pooled analyses in 16a/16b.
        diff_col <- find_col(res, "^diff_robust.*unbagged")

        res$tax_label <- label_fn(ps_raw, res$taxon)

        df <- res |>
                dplyr::mutate(
                        lfc        = .data[[lfc_col]],
                        q          = pmax(.data[[q_col]], 1e-300),
                        signif     = as.logical(.data[[diff_col]]),
                        neg_log10q = -log10(q),
                        direction  = dplyr::case_when(
                                signif & lfc > 0 ~ "higher_unbagged",
                                signif & lfc < 0 ~ "higher_bagged",
                                TRUE             ~ "ns"
                        )
                )

        # label top 5 significant taxa by |lfc|; fewer labels keep compact panels readable
        label_df <- df |> dplyr::filter(signif) |> dplyr::slice_max(abs(lfc), n = 5)

        ggplot(df, aes(x = lfc, y = neg_log10q)) +
                geom_point(aes(colour = direction), size = 1.2, alpha = 0.75) +
                geom_hline(yintercept = -log10(FDR_ALPHA),
                           linetype = "dashed", colour = "grey50", linewidth = 0.3) +
                geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
                ggrepel::geom_text_repel(
                        data = label_df, aes(label = tax_label),
                        size = 2, max.overlaps = 10, segment.colour = "grey60"
                ) +
                scale_colour_manual(
                        values = c(higher_unbagged = "#D6604D",
                                   higher_bagged   = "#4393C3",
                                   ns              = "#CCCCCC"),
                        # static labels (no per-panel counts) so guides can be collected
                        labels = c(higher_unbagged = "higher in unbagged",
                                   higher_bagged   = "higher in bagged",
                                   ns              = "not significant"),
                        # explicit breaks ensure consistent legend order and all levels shown
                        breaks = c("higher_unbagged", "higher_bagged", "ns"),
                        drop   = FALSE,
                        name   = NULL
                ) +
                labs(
                        title = title,
                        x     = if (show_x_label) "log2FC (unbagged / bagged)" else NULL,
                        y     = if (show_y_label) expression(-log[10](q[BH]))   else NULL
                ) +
                theme_bw(base_size = 8) +
                theme(
                        legend.position = "bottom",
                        plot.title      = element_text(size = 8, face = "bold", hjust = 0.5),
                        axis.title      = element_text(size = 7)
                )
}

# ── 3. Per-farm analysis ──────────────────────────────────────────────────────

cat("\n── 3. Per-farm ANCOMBC2 ─────────────────────────────────────────────────\n")

# Subsets ps_raw to one farm's bagged + unbagged flowers, applies prevalence
# filter, and runs ancombc2 with fix_formula = "sample_type".
run_farm <- function(ps_raw, farm) {
        sdata <- as(sample_data(ps_raw), "data.frame")
        keep  <- sdata$farm_id == farm &
                 sdata$sample_type %in% c("bagged_flower", "unbagged_flower")
        ps_f  <- prune_samples(keep, ps_raw)

        n_bag <- sum(sample_data(ps_f)$sample_type == "bagged_flower")
        n_unb <- sum(sample_data(ps_f)$sample_type == "unbagged_flower")

        sample_data(ps_f)$sample_type <- factor(
                sample_data(ps_f)$sample_type,
                levels = c("bagged_flower", "unbagged_flower")
        )
        ps_f <- filter_ps(ps_f)
        cat(sprintf("  %s: bagged=%d  unbagged=%d  taxa_after_filter=%d\n",
                    farm, n_bag, n_unb, ntaxa(ps_f)))

        if (ntaxa(ps_f) == 0 || nsamples(ps_f) < 10) {
                message("  Skipping ", farm, ": too few taxa/samples after filter")
                return(NULL)
        }

        tryCatch(
                ancombc2(data = ps_f, fix_formula = "sample_type",
                         p_adj_method = "BH", alpha = FDR_ALPHA, verbose = FALSE),
                error = function(e) {
                        message("  ancombc2() failed for ", farm, ": ", conditionMessage(e))
                        NULL
                }
        )
}

plots_16s <- vector("list", length(FARM_LEVELS))
plots_its <- vector("list", length(FARM_LEVELS))
csv_rows  <- list()

cat("  --- 16S ---\n")
for (i in seq_along(FARM_LEVELS)) {
        farm <- FARM_LEVELS[i]
        out  <- run_farm(ps_16s, farm)
        plots_16s[[i]] <- plot_farm_volcano(
                out, ps_16s, make_tax_label, farm,
                show_x_label = FALSE,   # x label on ITS1 (bottom) row only
                show_y_label = (i == 1) # y label on leftmost panel only
        )
        if (!is.null(out)) {
                res      <- out$res
                lfc_col  <- find_col(res, "^lfc.*unbagged")
                q_col    <- find_col(res, "^q.*unbagged")
                fdr_col  <- find_col(res, "^diff.*unbagged")          # FDR alone
                ss_col   <- find_col(res, "^passed_ss.*unbagged")     # sensitivity check
                diff_col <- find_col(res, "^diff_robust.*unbagged")   # both
                if (!is.na(lfc_col) && !is.na(q_col) && !is.na(diff_col)) {
                        csv_rows[[length(csv_rows) + 1]] <- dplyr::tibble(
                                marker    = "16S",
                                farm      = farm,
                                taxon     = res$taxon,
                                tax_label = make_tax_label(ps_16s, res$taxon),
                                lfc       = res[[lfc_col]],
                                q         = res[[q_col]],
                                diff_fdr  = as.logical(res[[fdr_col]]),
                                passed_ss = as.logical(res[[ss_col]]),
                                diff      = as.logical(res[[diff_col]])
                        )
                }
        }
}

cat("  --- ITS1 ---\n")
for (i in seq_along(FARM_LEVELS)) {
        farm <- FARM_LEVELS[i]
        out  <- run_farm(ps_its, farm)
        plots_its[[i]] <- plot_farm_volcano(
                out, ps_its, make_tax_label, farm,
                show_x_label = TRUE,    # x label on bottom row
                show_y_label = (i == 1)
        )
        if (!is.null(out)) {
                res      <- out$res
                lfc_col  <- find_col(res, "^lfc.*unbagged")
                q_col    <- find_col(res, "^q.*unbagged")
                fdr_col  <- find_col(res, "^diff.*unbagged")          # FDR alone
                ss_col   <- find_col(res, "^passed_ss.*unbagged")     # sensitivity check
                diff_col <- find_col(res, "^diff_robust.*unbagged")   # both
                if (!is.na(lfc_col) && !is.na(q_col) && !is.na(diff_col)) {
                        csv_rows[[length(csv_rows) + 1]] <- dplyr::tibble(
                                marker    = "ITS1",
                                farm      = farm,
                                taxon     = res$taxon,
                                tax_label = make_tax_label(ps_its, res$taxon),
                                lfc       = res[[lfc_col]],
                                q         = res[[q_col]],
                                diff_fdr  = as.logical(res[[fdr_col]]),
                                passed_ss = as.logical(res[[ss_col]]),
                                diff      = as.logical(res[[diff_col]])
                        )
                }
        }
}

# ── 4. Assemble 2 x 7 patchwork ──────────────────────────────────────────────

cat("\n── 4. Assembling 2 x 7 patchwork ──────────────────────────────────────\n")

# Row-level patchworks with row label as title
row_16s <- wrap_plots(plots_16s, nrow = 1) +
        plot_annotation(
                title = "16S Bacteria",
                theme = theme(plot.title = element_text(size = 9, face = "bold"))
        )

row_its <- wrap_plots(plots_its, nrow = 1) +
        plot_annotation(
                title = "ITS1 OTU97 Fungi",
                theme = theme(plot.title = element_text(size = 9, face = "bold"))
        )

fig <- (row_16s / row_its) +
        plot_layout(guides = "collect") +
        plot_annotation(
                title   = "ANCOMBC2: Insect footprint by farm (bagged vs unbagged)",
                caption = sprintf(
                        "FDR < %.2f + sensitivity filter  |  prevalence >= %.0f%%  |  fix_formula: sample_type  |  no random effects",
                        FDR_ALPHA, 100 * PREV_THRESH
                ),
                theme = theme(
                        plot.title   = element_text(size = 12, face = "bold"),
                        plot.caption = element_text(size = 7,  colour = "grey40")
                )
        ) &
        theme(legend.position = "bottom")

# ── 5. Save outputs ───────────────────────────────────────────────────────────

cat("\n── 5. Saving outputs ────────────────────────────────────────────────────\n")

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

ggsave(
        here("results", "figures", "16c_ancombc2_insect_footprint_by_farm.png"),
        plot   = fig,
        width  = 20,   # 7 columns x ~2.9 inches each
        height = 10,   # 2 rows x 5 inches each
        dpi    = 300
)
cat("  Saved: 16c_ancombc2_insect_footprint_by_farm.png\n")

if (length(csv_rows) > 0) {
        dplyr::bind_rows(csv_rows) |>
                dplyr::arrange(marker, farm, q) |>
                readr::write_csv(here("results", "tables",
                                      "16c_ancombc2_insect_footprint_by_farm.csv"))
        cat("  Saved: 16c_ancombc2_insect_footprint_by_farm.csv\n")
} else {
        cat("  No results to save (all farms failed).\n")
}

cat(sprintf("\n16c done. %d farms x 2 markers analysed.\n", length(FARM_LEVELS)))
