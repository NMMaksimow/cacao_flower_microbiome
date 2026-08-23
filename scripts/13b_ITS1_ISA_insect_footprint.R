## ============================================================================
## 13b — Indicator Species Analysis (ISA) — ITS1 Fungi (CFM)
##
## Presence/absence ISA (IndVal.g) on ITS1 OTU97 data.
## Contrast: bagged vs unbagged flowers — insect-vectored fungi
##
## Farm-confounding test flags any indicator taxon where >= FARM_CONC of
## all detections come from a single farm (stored in CSV; visible in figure
## as columns with colour in only one farm).
##
## Input:
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output:
##   results/figures/13b_ISA_insect_footprint.png
##   results/tables/13b_ISA_insect_footprint.csv
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

P_PERM      <- 0.05   # permutation p-value threshold (per multipatt call)
FDR_ALPHA   <- 0.05   # BH-corrected p-value threshold
STAT_THRESH <- 0.30   # minimum IndVal.g statistic
N_PERM      <- 999    # permutations per multipatt() call
FARM_CONC   <- 0.70   # flag if >= this fraction of detections are in one farm

RANK_PLOT   <- "Genus"  # taxonomic rank for heatmap aggregation

FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

# ── 1. Libraries & data ───────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(indicspecies)
library(patchwork)
library(vegan)

ps_ITS1 <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds"))

cat(sprintf("Loaded phyloseq: %d taxa x %d samples\n",
            ntaxa(ps_ITS1), nsamples(ps_ITS1)))

# ── Helper functions ──────────────────────────────────────────────────────────

prune_for_isa <- function(physeq) {
        physeq |>
                (\(x) prune_taxa(taxa_sums(x) > 1, x))() |>
                (\(x) prune_taxa(apply(otu_table(x) > 0, 1, sum) > 1, x))()
}

run_isa <- function(physeq, group_var, n_perm = N_PERM) {
        grp <- sample_data(physeq)[[group_var]]
        otu <- as(otu_table(physeq), "matrix")
        pa  <- if (taxa_are_rows(physeq)) t(otu) > 0 else otu > 0
        multipatt(pa, grp, func = "IndVal.g", control = how(nperm = n_perm))
}

extract_isa_hits <- function(mp_result, physeq, group_var,
                              p_perm      = P_PERM,
                              fdr_alpha   = FDR_ALPHA,
                              stat_thresh = STAT_THRESH) {

        group_levels <- levels(factor(sample_data(physeq)[[group_var]]))

        sign_df <- as.data.frame(mp_result$sign) |>
                tibble::rownames_to_column("ASV") |>
                dplyr::mutate(p_adj = p.adjust(p.value, method = "BH"))

        hits <- sign_df |>
                dplyr::filter(
                        !is.na(p.value),
                        p.value <= p_perm,
                        p_adj   <= fdr_alpha,
                        stat    >= stat_thresh
                ) |>
                dplyr::mutate(
                        indicator_group = dplyr::case_when(
                                index <= length(group_levels) ~ group_levels[index],
                                TRUE ~ "both"
                        )
                )

        if (nrow(hits) == 0) {
                cat("  No significant indicators found.\n")
                return(tibble())
        }

        otu <- as(otu_table(physeq), "matrix")
        pa  <- if (taxa_are_rows(physeq)) otu > 0 else t(otu) > 0
        grp <- sample_data(physeq)[[group_var]]

        prev_tbl <- purrr::map_dfr(group_levels, function(lvl) {
                idx    <- which(grp == lvl)
                pa_sub <- pa[hits$ASV, idx, drop = FALSE]
                tibble::tibble(ASV = hits$ASV, group = lvl,
                               prevalence = rowMeans(pa_sub))
        }) |>
                tidyr::pivot_wider(
                        names_from   = group,
                        values_from  = prevalence,
                        names_prefix = "prev_"
                )

        tax_df <- as.data.frame(tax_table(physeq)) |>
                tibble::rownames_to_column("ASV")

        hits |>
                dplyr::left_join(prev_tbl, by = "ASV") |>
                dplyr::left_join(tax_df,   by = "ASV")
}

test_farm_confounding <- function(hits, physeq, farm_conc = FARM_CONC) {
        if (nrow(hits) == 0) return(tibble())

        otu <- as(otu_table(physeq), "matrix")
        pa  <- if (taxa_are_rows(physeq)) otu > 0 else t(otu) > 0
        pa  <- pa[hits$ASV, , drop = FALSE]

        meta <- as(sample_data(physeq), "data.frame") |>
                tibble::rownames_to_column("Sample") |>
                dplyr::select(Sample, farm_id)

        purrr::map_dfr(hits$ASV, function(asv) {
                detected <- colnames(pa)[pa[asv, ]]
                if (length(detected) < 2) {
                        return(tibble::tibble(
                                ASV               = asv,
                                chi_p             = NA_real_,
                                pct_dominant_farm = 1.0,
                                farm_confounded   = TRUE
                        ))
                }
                fc <- meta |>
                        dplyr::filter(Sample %in% detected) |>
                        dplyr::count(farm_id)

                chi_p   <- tryCatch(chisq.test(fc$n)$p.value,
                                    error = function(e) NA_real_)
                pct_dom <- max(fc$n) / sum(fc$n)

                tibble::tibble(
                        ASV               = asv,
                        chi_p             = chi_p,
                        pct_dominant_farm = pct_dom,
                        farm_confounded   = pct_dom >= farm_conc
                )
        })
}

# Combined heatmap: 2 treatment columns + 7 farm columns in one figure.
# Layout: facet_grid(. ~ section, scales = "free_x", space = "free_x") gives
# panel widths proportional to number of columns (2:7 ratio).
# Farm columns show prevalence within the focal treatment (indicator_grp) per farm,
# making it visible whether an indicator taxon is present across all farms or
# concentrated in one (farm-confounded taxa show colour in one column only).
plot_indicator_heatmap_with_farms <- function(
        hits, physeq, indicator_grp, group_cols,
        group_var     = "sample_type",
        rank          = RANK_PLOT,
        min_asv       = 1,
        feature_label = "OTUs",   # "OTUs" for ITS1, "ASVs" for 16S
        title         = NULL,
        show_x_labels = TRUE,
        farm_levels   = FARM_LEVELS
) {
        if (nrow(hits) == 0 || !"indicator_group" %in% colnames(hits)) {
                return(ggplot() + theme_void() +
                               labs(title = if (!is.null(title)) title else
                                       paste("No indicators for", indicator_grp)))
        }

        rank_sym <- rlang::sym(rank)

        # ── Genus-level aggregation (preserves existing ordering logic) ───────────

        hits_genus <- hits |>
                dplyr::filter(indicator_group == indicator_grp,
                              !is.na(!!rank_sym)) |>
                dplyr::group_by(!!rank_sym, Order, Class, Phylum) |>
                dplyr::summarise(
                        n_feature = dplyr::n(),
                        dplyr::across(dplyr::starts_with("prev_"),
                                      ~ mean(.x, na.rm = TRUE)),
                        ordering  = mean(
                                .data[[paste0("prev_", group_cols[2])]] -
                                        .data[[paste0("prev_", group_cols[1])]],
                                na.rm = TRUE
                        ),
                        .groups = "drop"
                ) |>
                dplyr::filter(n_feature >= min_asv) |>
                dplyr::arrange(ordering) |>
                dplyr::rowwise() |>
                dplyr::mutate(
                        Taxon_label = {
                                val <- !!rank_sym
                                # UNITE Incertae_sedis entries are uninformative at genus level;
                                # fall back to higher ranks in that case
                                if (!is.na(val) && !grepl("Incertae_sedis", val))
                                        as.character(val)
                                else if (!is.na(Order) && !grepl("Incertae_sedis", Order))
                                        as.character(Order)
                                else if (!is.na(Class) && !grepl("Incertae_sedis", Class))
                                        as.character(Class)
                                else if (!is.na(Phylum) && !grepl("Incertae_sedis", Phylum))
                                        as.character(Phylum)
                                else as.character(val)
                        },
                        Taxon_label = paste0(Taxon_label, " (", n_feature, " ", feature_label, ")")
                ) |>
                dplyr::ungroup() |>
                dplyr::mutate(Taxon_label = factor(Taxon_label, unique(Taxon_label)))

        if (nrow(hits_genus) == 0) {
                return(ggplot() + theme_void() +
                               labs(title = if (!is.null(title)) title else
                                       paste("No indicators for", indicator_grp)))
        }

        taxon_levels   <- levels(hits_genus$Taxon_label)
        genus_to_label <- hits_genus |>
                dplyr::select(!!rank_sym, Taxon_label) |>
                dplyr::mutate(Taxon_label = as.character(Taxon_label))

        # ── Treatment prevalence columns (section = "Treatment") ─────────────────

        treatment_data <- hits_genus |>
                tidyr::pivot_longer(
                        cols         = dplyr::starts_with("prev_"),
                        names_to     = "col_label",
                        values_to    = "prevalence",
                        names_prefix = "prev_"
                ) |>
                dplyr::select(Taxon_label, col_label, prevalence) |>
                dplyr::mutate(
                        section     = "Treatment",
                        Taxon_label = as.character(Taxon_label)
                )

        # ── Farm prevalence columns (section = "Farm", focal treatment only) ─────

        ps_focal <- prune_samples(
                sample_data(physeq)[[group_var]] == indicator_grp, physeq
        )
        otu <- as(otu_table(ps_focal), "matrix")
        pa  <- if (taxa_are_rows(ps_focal)) otu > 0 else t(otu) > 0

        meta   <- as(sample_data(ps_focal), "data.frame") |>
                tibble::rownames_to_column("Sample") |>
                dplyr::select(Sample, farm_id)
        farm_n <- meta |> dplyr::count(farm_id, name = "n_farm")

        indicator_asvs <- hits |>
                dplyr::filter(indicator_group == indicator_grp) |>
                dplyr::pull(ASV)
        present_asvs <- intersect(indicator_asvs, rownames(pa))

        if (length(present_asvs) > 0) {
                farm_raw <- as.data.frame(t(pa[present_asvs, , drop = FALSE])) |>
                        tibble::rownames_to_column("Sample") |>
                        tidyr::pivot_longer(-Sample, names_to = "ASV",
                                            values_to = "detected") |>
                        dplyr::filter(detected) |>
                        dplyr::left_join(meta, by = "Sample") |>
                        dplyr::left_join(
                                hits |> dplyr::select(ASV, !!rank_sym),
                                by = "ASV"
                        ) |>
                        # collapse: a sample is positive for a genus if any OTU detected
                        dplyr::group_by(Sample, farm_id, !!rank_sym) |>
                        dplyr::summarise(.groups = "drop") |>
                        dplyr::group_by(farm_id, !!rank_sym) |>
                        dplyr::summarise(n_detected = dplyr::n(), .groups = "drop") |>
                        dplyr::left_join(farm_n, by = "farm_id") |>
                        dplyr::mutate(prevalence = n_detected / n_farm) |>
                        dplyr::left_join(genus_to_label,
                                         by = as.character(rank_sym)) |>
                        dplyr::filter(!is.na(Taxon_label))
        } else {
                farm_raw <- tibble::tibble(
                        farm_id = character(), Taxon_label = character(),
                        prevalence = numeric()
                )
        }

        # Complete the grid: zero prevalence for undetected genus × farm combos
        farm_data <- tidyr::expand_grid(
                Taxon_label = taxon_levels,
                col_label   = farm_levels
        ) |>
                dplyr::left_join(
                        farm_raw |>
                                dplyr::mutate(col_label = as.character(farm_id)) |>
                                dplyr::select(Taxon_label, col_label, prevalence),
                        by = c("Taxon_label", "col_label")
                ) |>
                dplyr::mutate(
                        prevalence = dplyr::if_else(is.na(prevalence), 0, prevalence),
                        section    = "Farm"
                )

        # ── Combine and plot ──────────────────────────────────────────────────────

        plot_df <- dplyr::bind_rows(
                treatment_data |> dplyr::select(Taxon_label, section, col_label, prevalence),
                farm_data      |> dplyr::select(Taxon_label, section, col_label, prevalence)
        ) |>
                dplyr::mutate(
                        section     = factor(section, c("Treatment", "Farm")),
                        col_label   = factor(col_label, c(group_cols, farm_levels)),
                        Taxon_label = factor(Taxon_label, taxon_levels)
                )

        plot_title <- if (!is.null(title)) title else
                paste0("Indicators for ", indicator_grp,
                       "  (", rank, " level,  IndVal.g >= ", STAT_THRESH,
                       ",  BH p <= ", FDR_ALPHA, ")")

        ggplot(plot_df, aes(x = col_label, y = Taxon_label, fill = prevalence)) +
                geom_tile(colour = "grey90") +
                # facet_grid separates treatment (2 cols) from farm (7 cols);
                # space = "free_x" makes panel widths proportional to column count (2:7)
                facet_grid(. ~ section, scales = "free_x", space = "free_x") +
                scale_fill_viridis_c(
                        option = "magma", limits = c(0, 1),
                        name   = "Prevalence"  # fraction of group/farm samples with taxon detected
                ) +
                labs(x = NULL, y = NULL, title = plot_title) +
                theme_minimal(base_size = 10) +
                theme(
                        panel.grid       = element_blank(),
                        strip.text       = element_text(size = 8, face = "bold"),
                        strip.background = element_blank(),
                        axis.text.x      = if (show_x_labels)
                                element_text(angle = 30, hjust = 1)
                        else
                                element_blank(),
                        axis.ticks.x     = if (show_x_labels) element_line()
                                           else element_blank()
                )
}

# ── 2. ISA: bagged vs unbagged ────────────────────────────────────────────────

cat("\n── 2. Bagged vs unbagged ISA (insect footprint) ────────────────────────\n")

ps_bvu <- ps_ITS1 |>
        subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower")) |>
        prune_for_isa()

cat(sprintf("  Dataset: %d taxa x %d samples  (bagged: %d, unbagged: %d)\n",
            ntaxa(ps_bvu), nsamples(ps_bvu),
            sum(sample_data(ps_bvu)$sample_type == "bagged_flower"),
            sum(sample_data(ps_bvu)$sample_type == "unbagged_flower")))

set.seed(42)
cat(sprintf("  Running multipatt (N_PERM = %d)...\n", N_PERM))
mp_bvu <- run_isa(ps_bvu, "sample_type")

hits_bvu <- extract_isa_hits(mp_bvu, ps_bvu, "sample_type")
if (nrow(hits_bvu) > 0) {
        conf_bvu <- test_farm_confounding(hits_bvu, ps_bvu)
        hits_bvu <- dplyr::left_join(hits_bvu, conf_bvu, by = "ASV")
}

cat(sprintf("  Indicators: %d total  (%d unbagged, %d bagged)  |  %d farm-confounded\n",
            nrow(hits_bvu),
            if (nrow(hits_bvu) > 0) sum(hits_bvu$indicator_group == "unbagged_flower") else 0L,
            if (nrow(hits_bvu) > 0) sum(hits_bvu$indicator_group == "bagged_flower")   else 0L,
            if (nrow(hits_bvu) > 0) sum(hits_bvu$farm_confounded, na.rm = TRUE)        else 0L))

# ── 3. Figure ─────────────────────────────────────────────────────────────────

grp_cols_bvu <- c("bagged_flower", "unbagged_flower")

# Top panel: x-labels suppressed (shared with bagged panel below)
p_unbagged <- plot_indicator_heatmap_with_farms(
        hits_bvu, ps_bvu, "unbagged_flower", grp_cols_bvu,
        title         = "Unbagged-flower indicators (insect-vectored fungi)",
        show_x_labels = FALSE,
        feature_label = "OTUs"
)

# Bottom panel: x-labels visible (treatment names + farm IDs)
p_bagged <- plot_indicator_heatmap_with_farms(
        hits_bvu, ps_bvu, "bagged_flower", grp_cols_bvu,
        title         = "Bagged-flower indicators",
        show_x_labels = TRUE,
        feature_label = "OTUs"
)

n_rows_ub <- if (nrow(hits_bvu) > 0) {
        hits_bvu |>
                dplyr::filter(indicator_group == "unbagged_flower",
                              !is.na(.data[[RANK_PLOT]])) |>
                dplyr::pull(RANK_PLOT) |> unique() |> length()
} else 0L

n_rows_b <- if (nrow(hits_bvu) > 0) {
        hits_bvu |>
                dplyr::filter(indicator_group == "bagged_flower",
                              !is.na(.data[[RANK_PLOT]])) |>
                dplyr::pull(RANK_PLOT) |> unique() |> length()
} else 0L

fig <- (p_unbagged / p_bagged) +
        plot_layout(heights = c(max(2, n_rows_ub), max(2, n_rows_b)),
                    guides  = "collect") +
        plot_annotation(
                title   = "ISA: Bagged vs unbagged flowers — ITS1 OTU97 fungi",
                caption = sprintf(
                        "IndVal.g >= %.2f  |  perm p <= %.2f  |  BH p <= %.2f  |  N_PERM = %d\nFarm columns: prevalence within focal treatment per farm",
                        STAT_THRESH, P_PERM, FDR_ALPHA, N_PERM
                ),
                theme = theme(plot.title = element_text(face = "bold"))
        ) &
        theme(legend.position = "right")

# ── 4. Save outputs ───────────────────────────────────────────────────────────

cat("\n── 4. Saving outputs ────────────────────────────────────────────────────\n")

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

ggsave(
        here("results", "figures", "13b_ISA_insect_footprint.png"),
        plot   = fig,
        width  = 12,
        height = max(5, (n_rows_ub + n_rows_b) * 0.6 + 3),
        dpi    = 300
)
cat("  Saved: 13b_ISA_insect_footprint.png\n")

if (nrow(hits_bvu) > 0) {
        hits_bvu |>
                dplyr::select(ASV, indicator_group, stat, p.value, p_adj,
                              dplyr::starts_with("prev_"),
                              Kingdom:Species,
                              chi_p, pct_dominant_farm, farm_confounded) |>
                readr::write_csv(here("results", "tables",
                                      "13b_ISA_insect_footprint.csv"))
        cat("  Saved: 13b_ISA_insect_footprint.csv\n")
}

cat("\n13b done.\n")
