## ============================================================================
## 13c — Indicator Species Analysis (ISA) by farm — 16S + ITS1 combined
##
## Presence/absence ISA (IndVal.g) run within each of 7 farms independently.
## Contrast: bagged_flower (reference) vs unbagged_flower.
##
## Notes:
##   - farm_confounding test omitted: within a single farm 100% of detections
##     come from that farm by definition; the test is meaningless here
##   - BH correction is per-farm; results are exploratory
##   - prune_for_isa() is more stringent on ~42 within-farm samples than on
##     ~200 globally — expect fewer testable taxa per farm than in 13a/13b
##   - runtime ~10–20 min; intermediate RDS saved before figure assembly
##
## Input:
##   results/rds/ps_16S_bacteria_biosamples.rds
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output:
##   results/rds/13c_ISA_by_farm_results.rds
##   results/figures/13c_ISA_insect_footprint_by_farm.png
##   results/tables/13c_ISA_insect_footprint_by_farm.csv
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

P_PERM      <- 0.05
FDR_ALPHA   <- 0.05
STAT_THRESH <- 0.30
N_PERM      <- 999

RANK_PLOT   <- "Genus"
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

# ── 1. Libraries & data ───────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(indicspecies)
library(patchwork)
library(vegan)

ps_16s <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds"))
ps_its <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds"))

cat(sprintf("16S loaded:  %d taxa x %d samples\n", ntaxa(ps_16s), nsamples(ps_16s)))
cat(sprintf("ITS1 loaded: %d taxa x %d samples\n", ntaxa(ps_its), nsamples(ps_its)))

# ── 2. Helper functions ───────────────────────────────────────────────────────

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
                return(tibble::tibble())
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

# Compact 2-column prevalence heatmap (bagged | unbagged) for one farm panel.
# No farm columns: within-farm ISA; farm breakdown is not meaningful.
# label_type:    "16S" uses nchar fallback; "ITS1" uses Incertae_sedis fallback.
# show_x_labels: suppress x-axis labels on top row (16S).
plot_farm_isa_panel <- function(hits, farm,
                                 label_type    = c("16S", "ITS1"),
                                 feature_label = "features",
                                 rank          = RANK_PLOT,
                                 show_x_labels = TRUE) {

        label_type <- match.arg(label_type)
        title      <- toupper(farm)

        void_panel <- function() {
                ggplot() + theme_void(base_size = 8) +
                        annotate("text", x = 0.5, y = 0.5, hjust = 0.5,
                                 label = "No hits", size = 2.5, colour = "grey50") +
                        labs(title = title) +
                        theme(plot.title      = element_text(size = 8, face = "bold",
                                                             hjust = 0.5),
                              legend.position = "bottom")
        }

        if (nrow(hits) == 0) return(void_panel())

        grp_cols <- c("bagged_flower", "unbagged_flower")
        rank_sym <- rlang::sym(rank)

        # uninformative genus check differs by classifier
        if (label_type == "16S") {
                bad <- function(x) is.na(x) || nchar(as.character(x)) <= 3
        } else {
                bad <- function(x) is.na(x) || grepl("Incertae_sedis", as.character(x))
        }

        hits_genus <- hits |>
                dplyr::filter(!is.na(!!rank_sym)) |>
                dplyr::group_by(!!rank_sym, Order, Class, Phylum) |>
                dplyr::summarise(
                        n_feature = dplyr::n(),
                        dplyr::across(dplyr::starts_with("prev_"),
                                      ~ mean(.x, na.rm = TRUE)),
                        # unbagged indicators (positive ordering) sorted to top of y-axis
                        ordering  = mean(
                                .data[["prev_unbagged_flower"]] -
                                        .data[["prev_bagged_flower"]],
                                na.rm = TRUE
                        ),
                        .groups = "drop"
                ) |>
                dplyr::arrange(ordering) |>
                dplyr::rowwise() |>
                dplyr::mutate(
                        Taxon_label = {
                                val <- !!rank_sym
                                # fall back to higher ranks when genus is uninformative
                                if      (!bad(val))     as.character(val)
                                else if (!bad(Order))  as.character(Order)
                                else if (!bad(Class))  as.character(Class)
                                else if (!bad(Phylum)) as.character(Phylum)
                                else as.character(val)
                        },
                        Taxon_label = paste0(Taxon_label, " (", n_feature, " ",
                                             feature_label, ")")
                ) |>
                dplyr::ungroup() |>
                dplyr::mutate(Taxon_label = factor(Taxon_label, unique(Taxon_label)))

        if (nrow(hits_genus) == 0) return(void_panel())

        plot_df <- hits_genus |>
                tidyr::pivot_longer(
                        cols         = dplyr::starts_with("prev_"),
                        names_to     = "col_label",
                        values_to    = "prevalence",
                        names_prefix = "prev_"
                ) |>
                dplyr::filter(col_label %in% grp_cols) |>
                dplyr::mutate(
                        col_label   = factor(col_label, grp_cols),
                        Taxon_label = factor(Taxon_label, levels(hits_genus$Taxon_label))
                )

        ggplot(plot_df, aes(x = col_label, y = Taxon_label, fill = prevalence)) +
                geom_tile(colour = "grey90") +
                scale_fill_viridis_c(
                        option = "magma", limits = c(0, 1),
                        name   = "Prevalence"  # fraction of farm group's samples with taxon detected
                ) +
                scale_x_discrete(labels = c(bagged_flower   = "bagged",
                                            unbagged_flower = "unbagged")) +
                labs(x = NULL, y = NULL, title = title) +
                theme_minimal(base_size = 8) +
                theme(
                        panel.grid       = element_blank(),
                        plot.title       = element_text(size = 8, face = "bold",
                                                        hjust = 0.5),
                        axis.text.y      = element_text(size = 6),
                        axis.text.x      = if (show_x_labels)
                                element_text(angle = 30, hjust = 1, size = 7)
                        else
                                element_blank(),
                        axis.ticks.x     = if (show_x_labels) element_line()
                                           else element_blank(),
                        legend.position  = "bottom"
                )
}

# Subsets one farm's bagged + unbagged flowers, prunes for ISA, runs multipatt,
# and extracts hits. Returns list(hits, ps) or list(hits = empty tibble, ps = NULL).
run_farm_isa <- function(ps_raw, farm, group_var = "sample_type") {
        sdata <- as(sample_data(ps_raw), "data.frame")
        keep  <- sdata$farm_id == farm &
                 sdata[[group_var]] %in% c("bagged_flower", "unbagged_flower")
        ps_f  <- prune_samples(keep, ps_raw)
        ps_f  <- prune_for_isa(ps_f)

        n_bag <- sum(sample_data(ps_f)[[group_var]] == "bagged_flower")
        n_unb <- sum(sample_data(ps_f)[[group_var]] == "unbagged_flower")
        n_tax <- ntaxa(ps_f)
        cat(sprintf("  %s: bagged=%d  unbagged=%d  testable_taxa=%d\n",
                    farm, n_bag, n_unb, n_tax))

        # guard: multipatt() errors with < 5 taxa or any group having only 1 sample
        if (n_tax < 5 || n_bag < 2 || n_unb < 2) {
                message("  Skipping ", farm, ": insufficient taxa or samples after filter")
                return(list(hits = tibble::tibble(), ps = NULL))
        }

        mp <- tryCatch(
                run_isa(ps_f, group_var),
                error = function(e) {
                        message("  multipatt() failed for ", farm, ": ",
                                conditionMessage(e))
                        NULL
                }
        )
        if (is.null(mp)) return(list(hits = tibble::tibble(), ps = NULL))

        hits <- extract_isa_hits(mp, ps_f, group_var)
        list(hits = hits, ps = ps_f)
}

# ── 3. Per-farm ISA loop ──────────────────────────────────────────────────────

cat("\n── 3. Per-farm ISA (both markers) ──────────────────────────────────────\n")
set.seed(42)   # before the loop for reproducible permutation sequences

results <- list()

cat("  --- 16S ---\n")
for (farm in FARM_LEVELS) {
        r <- run_farm_isa(ps_16s, farm)
        results[[paste0("16S_", farm)]] <- r
        cat(sprintf("    -> %d significant indicators\n", nrow(r$hits)))
}

cat("  --- ITS1 ---\n")
for (farm in FARM_LEVELS) {
        r <- run_farm_isa(ps_its, farm)
        results[[paste0("ITS1_", farm)]] <- r
        cat(sprintf("    -> %d significant indicators\n", nrow(r$hits)))
}

# ── 4. Save intermediate RDS before figure assembly ───────────────────────────

cat("\n── 4. Saving intermediate RDS ───────────────────────────────────────────\n")
dir.create(here("results", "rds"), showWarnings = FALSE, recursive = TRUE)
saveRDS(results, here("results", "rds", "13c_ISA_by_farm_results.rds"))
cat("  Saved: 13c_ISA_by_farm_results.rds\n")

# ── 5. Build per-farm panels ──────────────────────────────────────────────────

cat("\n── 5. Building panels ───────────────────────────────────────────────────\n")

panels_16s <- lapply(FARM_LEVELS, function(farm) {
        plot_farm_isa_panel(
                results[[paste0("16S_", farm)]]$hits,
                farm          = farm,
                label_type    = "16S",
                feature_label = "ASVs",
                show_x_labels = FALSE   # x labels on ITS1 (bottom) row only
        )
})

panels_its <- lapply(FARM_LEVELS, function(farm) {
        plot_farm_isa_panel(
                results[[paste0("ITS1_", farm)]]$hits,
                farm          = farm,
                label_type    = "ITS1",
                feature_label = "OTUs",
                show_x_labels = TRUE
        )
})

# ── 6. Assemble 2 × 7 patchwork ──────────────────────────────────────────────

cat("\n── 6. Assembling 2 x 7 patchwork ──────────────────────────────────────\n")

row_16s <- wrap_plots(panels_16s, nrow = 1) +
        plot_annotation(
                title = "16S Bacteria",
                theme = theme(plot.title = element_text(size = 9, face = "bold"))
        )

row_its <- wrap_plots(panels_its, nrow = 1) +
        plot_annotation(
                title = "ITS1 OTU97 Fungi",
                theme = theme(plot.title = element_text(size = 9, face = "bold"))
        )

fig <- (row_16s / row_its) +
        plot_layout(guides = "collect") +
        plot_annotation(
                title   = "ISA by farm: Insect footprint (bagged vs unbagged)",
                caption = sprintf(
                        "IndVal.g >= %.2f  |  perm p <= %.2f  |  BH p <= %.2f  |  N_PERM = %d  |  per-farm BH correction",
                        STAT_THRESH, P_PERM, FDR_ALPHA, N_PERM
                ),
                theme = theme(
                        plot.title   = element_text(size = 12, face = "bold"),
                        plot.caption = element_text(size = 7,  colour = "grey40")
                )
        ) &
        theme(legend.position = "bottom")

# ── 7. Save outputs ───────────────────────────────────────────────────────────

cat("\n── 7. Saving outputs ────────────────────────────────────────────────────\n")

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

ggsave(
        here("results", "figures", "13c_ISA_insect_footprint_by_farm.png"),
        plot   = fig,
        width  = 20,   # 7 columns x ~2.9 inches each
        height = 10,   # 2 rows x 5 inches each
        dpi    = 300
)
cat("  Saved: 13c_ISA_insect_footprint_by_farm.png\n")

csv_rows <- dplyr::bind_rows(lapply(FARM_LEVELS, function(farm) {
        hits_16s <- results[[paste0("16S_",  farm)]]$hits
        hits_its <- results[[paste0("ITS1_", farm)]]$hits
        dplyr::bind_rows(
                if (nrow(hits_16s) > 0) dplyr::mutate(hits_16s, marker = "16S",  farm = farm) else NULL,
                if (nrow(hits_its) > 0) dplyr::mutate(hits_its, marker = "ITS1", farm = farm) else NULL
        )
}))

if (nrow(csv_rows) > 0) {
        csv_rows |>
                dplyr::rename(feature_id = ASV) |>  # OTU97 hashes for ITS1, ASV IDs for 16S
                dplyr::select(
                        marker, farm, feature_id, indicator_group,
                        stat, p.value, p_adj,
                        dplyr::any_of(c("prev_bagged_flower", "prev_unbagged_flower")),
                        dplyr::any_of(c("Kingdom", "Phylum", "Class", "Order",
                                        "Family", "Genus", "Species"))
                ) |>
                dplyr::arrange(marker, farm, p_adj) |>
                readr::write_csv(here("results", "tables",
                                      "13c_ISA_insect_footprint_by_farm.csv"))
        cat("  Saved: 13c_ISA_insect_footprint_by_farm.csv\n")
} else {
        cat("  No significant indicators across all farms and markers.\n")
}

cat(sprintf("\n13c done. %d farms x 2 markers analysed.\n", length(FARM_LEVELS)))
