## ============================================================================
## 15b — Indicator Species Analysis (ISA) — ITS1 Fungi OTU97 (CFM)
##        Fungal colonization footprint: no colonization vs maximum colonization
##
## Contrast: unbagged flowers only; fungal_colonization_score 0 vs 5.
## Intermediate scores (1–4) excluded (small n per score, not suitable for ISA).
##   score_0 — fungal_colonization_score == 0 (n ≈ 43, no visible colonization)
##   score_5 — fungal_colonization_score == 5 (n ≈ 20, maximum colonization)
##
## Binary ISA using IndVal.g. IndVal.g normalises specificity by within-group
## means, so it is robust to the ~43:20 imbalance; no downsampling needed.
##
## Farm-confounding test flags any indicator taxon where >= FARM_CONC of all
## detections come from a single farm.
##
## Input:
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output:
##   results/figures/15b_ISA_fungal_colonization.png
##   results/tables/15b_ISA_fungal_colonization.csv
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

P_PERM      <- 0.05   # permutation p-value threshold
FDR_ALPHA   <- 0.05   # BH-corrected p-value threshold
STAT_THRESH <- 0.30   # minimum IndVal.g statistic
N_PERM      <- 999    # permutations per multipatt() call
FARM_CONC   <- 0.70   # flag if >= this fraction of detections are in one farm

RANK_PLOT   <- "Genus"

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
                                TRUE ~ "combination"
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

# Combined heatmap: treatment columns + 7 farm columns in one figure.
# Layout: facet_grid(. ~ section, scales = "free_x", space = "free_x") gives
# panel widths proportional to number of columns (2:7 ratio for 2 groups).
# Farm columns show prevalence within the focal group (indicator_grp) per farm.
#
# group_var: name of the sample_data column holding the group factor
#            ("colonization_class" for 15a/15b)
# n_per_group: optional named numeric vector — appends "(n=N)" to treatment labels
plot_indicator_heatmap_with_farms <- function(
        hits, physeq, indicator_grp, group_cols,
        group_var     = "colonization_class",
        rank          = RANK_PLOT,
        min_asv       = 1,
        feature_label = "OTUs",
        title         = NULL,
        show_x_labels = TRUE,
        farm_levels   = FARM_LEVELS,
        n_per_group   = NULL
) {
        if (nrow(hits) == 0 || !"indicator_group" %in% colnames(hits)) {
                return(ggplot() + theme_void() +
                               labs(title = if (!is.null(title)) title else
                                       paste("No indicators for", indicator_grp)))
        }

        rank_sym <- rlang::sym(rank)

        # Pre-compute column name vectors for 2-group ordering
        focal_prev_col  <- paste0("prev_", indicator_grp)
        other_prev_col  <- paste0("prev_", setdiff(group_cols, indicator_grp))

        # ── Genus-level aggregation ────────────────────────────────────────────────

        hits_genus <- hits |>
                dplyr::filter(indicator_group == indicator_grp,
                              !is.na(!!rank_sym)) |>
                dplyr::group_by(!!rank_sym, Order, Class, Phylum) |>
                dplyr::summarise(
                        n_feature = dplyr::n(),
                        dplyr::across(dplyr::starts_with("prev_"),
                                      ~ mean(.x, na.rm = TRUE)),
                        .groups = "drop"
                ) |>
                dplyr::filter(n_feature >= min_asv) |>
                # 2-group ordering: focal prevalence minus the other group's prevalence
                dplyr::mutate(
                        ordering = .data[[focal_prev_col]] - .data[[other_prev_col]]
                ) |>
                dplyr::arrange(ordering) |>
                dplyr::rowwise() |>
                dplyr::mutate(
                        Taxon_label = {
                                val <- !!rank_sym
                                # Incertae_sedis is an uninformative placeholder; fall back up the tree
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

        # ── Farm prevalence columns (section = "Farm", focal group only) ─────────

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

        # Build custom x-axis labels: add "(n=N)" to treatment columns if provided
        x_labels <- if (!is.null(n_per_group)) {
                c(
                        setNames(
                                paste0(names(n_per_group), "\n(n=", n_per_group, ")"),
                                names(n_per_group)
                        ),
                        setNames(farm_levels, farm_levels)
                )
        } else NULL

        p <- ggplot(plot_df,
                    aes(x = col_label, y = Taxon_label, fill = prevalence)) +
                geom_tile(colour = "grey90") +
                # facet_grid separates treatment from farm columns;
                # space = "free_x" makes panel widths proportional to column count
                facet_grid(. ~ section, scales = "free_x", space = "free_x") +
                scale_fill_viridis_c(
                        option = "magma", limits = c(0, 1),
                        name   = "Prevalence"
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

        if (!is.null(x_labels)) {
                p <- p + scale_x_discrete(labels = x_labels)
        }
        p
}

# ── 2. Build contrast phyloseq ────────────────────────────────────────────────

cat("\n── 2. Fungal colonization ISA (score 0 vs 5) ───────────────────────────\n")

sdata <- as(sample_data(ps_ITS1), "data.frame")
keep  <- sdata$sample_type == "unbagged_flower" &
         !is.na(sdata$fungal_colonization_score) &
         sdata$fungal_colonization_score %in% c(0L, 5L)
ps_col <- prune_samples(keep, ps_ITS1)

sample_data(ps_col)$colonization_class <- factor(
        ifelse(sample_data(ps_col)$fungal_colonization_score == 0L,
               "score_0", "score_5"),
        levels = c("score_0", "score_5")
)

n_per_group <- table(sample_data(ps_col)$colonization_class) |>
        as.integer() |>
        setNames(levels(sample_data(ps_col)$colonization_class))

cat(sprintf("  Groups: score_0 n=%d, score_5 n=%d\n",
            n_per_group["score_0"], n_per_group["score_5"]))

ps_bvu <- ps_col |> prune_for_isa()

cat(sprintf("  ISA dataset: %d taxa x %d samples\n",
            ntaxa(ps_bvu), nsamples(ps_bvu)))

# ── 3. Run ISA ────────────────────────────────────────────────────────────────

set.seed(42)
cat(sprintf("  Running multipatt (N_PERM = %d)...\n", N_PERM))
mp_bvu <- run_isa(ps_bvu, "colonization_class")

hits_bvu <- extract_isa_hits(mp_bvu, ps_bvu, "colonization_class")
if (nrow(hits_bvu) > 0) {
        conf_bvu <- test_farm_confounding(hits_bvu, ps_bvu)
        hits_bvu <- dplyr::left_join(hits_bvu, conf_bvu, by = "ASV")

        cat(sprintf(
                "  Indicators: %d total  (score_0: %d, score_5: %d, combination: %d)  |  %d farm-confounded\n",
                nrow(hits_bvu),
                if (nrow(hits_bvu) > 0) sum(hits_bvu$indicator_group == "score_0") else 0L,
                if (nrow(hits_bvu) > 0) sum(hits_bvu$indicator_group == "score_5") else 0L,
                if (nrow(hits_bvu) > 0) sum(hits_bvu$indicator_group == "combination") else 0L,
                sum(hits_bvu$farm_confounded, na.rm = TRUE)
        ))
}

# ── 4. Figure ─────────────────────────────────────────────────────────────────

grp_cols_bvu <- c("score_0", "score_5")

# Top panel: high-colonization indicators, x-labels suppressed
p_high <- plot_indicator_heatmap_with_farms(
        hits_bvu, ps_bvu, "score_5", grp_cols_bvu,
        title         = "High-colonization indicators (score = 5)",
        show_x_labels = FALSE,
        feature_label = "OTUs",
        n_per_group   = n_per_group
)

# Bottom panel: no-colonization indicators, x-labels visible
p_low <- plot_indicator_heatmap_with_farms(
        hits_bvu, ps_bvu, "score_0", grp_cols_bvu,
        title         = "No-colonization indicators (score = 0)",
        show_x_labels = TRUE,
        feature_label = "OTUs",
        n_per_group   = n_per_group
)

n_rows <- function(indicator_grp) {
        if (nrow(hits_bvu) == 0) return(1L)
        hits_bvu |>
                dplyr::filter(indicator_group == indicator_grp,
                              !is.na(.data[[RANK_PLOT]])) |>
                dplyr::pull(RANK_PLOT) |> unique() |> length() |>
                max(1L)
}

n_high <- n_rows("score_5")
n_low  <- n_rows("score_0")

fig <- (p_high / p_low) +
        plot_layout(heights = c(max(2, n_high), max(2, n_low)),
                    guides  = "collect") +
        plot_annotation(
                title   = "ISA: Fungal colonization footprint — ITS1 OTU97 fungi  (unbagged flowers)",
                caption = sprintf(
                        "IndVal.g >= %.2f  |  perm p <= %.2f  |  BH p <= %.2f  |  N_PERM = %d\nFarm columns: prevalence within focal group per farm",
                        STAT_THRESH, P_PERM, FDR_ALPHA, N_PERM
                ),
                theme = theme(plot.title = element_text(face = "bold"))
        ) &
        theme(legend.position = "right")

# ── 5. Save outputs ───────────────────────────────────────────────────────────

cat("\n── 5. Saving outputs ────────────────────────────────────────────────────\n")

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

ggsave(
        here("results", "figures", "15b_ISA_fungal_colonization.png"),
        plot   = fig,
        width  = 12,
        height = max(6, (n_high + n_low) * 0.6 + 3),
        dpi    = 300
)
cat("  Saved: 15b_ISA_fungal_colonization.png\n")

if (nrow(hits_bvu) > 0) {
        hits_bvu |>
                dplyr::select(ASV, indicator_group, stat, p.value, p_adj,
                              dplyr::starts_with("prev_"),
                              Kingdom:Species,
                              chi_p, pct_dominant_farm, farm_confounded) |>
                readr::write_csv(here("results", "tables",
                                      "15b_ISA_fungal_colonization.csv"))
        cat("  Saved: 15b_ISA_fungal_colonization.csv\n")
}

cat(sprintf("\n15b done. %d indicators found.\n", nrow(hits_bvu)))
