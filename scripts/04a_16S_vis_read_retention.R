# ---- Packages ----
library(here)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(scales)

# ---- Load intermediate objects (saved by 01a) ----
ps_16S_raw              <- readRDS(here("results", "rds", "ps_16S_raw.rds"))
ps_16S_decontam         <- readRDS(here("results", "rds", "ps_16S_decontam.rds"))
ps_16S_bacteria         <- readRDS(here("results", "rds", "ps_16S_bacteria.rds"))
ps_16S_bacteria_biosamples <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds"))
decontam_16S_results          <- readRDS(here("results", "rds", "decontam_prev055_16S.rds"))

# Metadata as plain data.frame (needed for joins by sample_id)
metadata <- sample_data(ps_16S_raw) |>
        as.data.frame() |>
        (\(x) { class(x) <- "data.frame"; x })()

# ---- Output directories ----
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)


################################################################################
# 1. Decontam diagnostic: ASV prevalence in controls vs biological samples
################################################################################
# ---- Sequencing retention after decontamination and removal off-targe sequences ----
read_retention_loss_table_16S <-
        dplyr::bind_rows(
                ps_16S_raw |>
                        sample_sums() |>
                        tibble::enframe(name = "sample_id",
                                        value = "n_reads") |>
                        dplyr::mutate(stage = "Raw"),
                
                ps_16S_decontam |>
                        sample_sums() |>
                        tibble::enframe(name = "sample_id",
                                        value = "n_reads") |>
                        dplyr::mutate(stage = "After decontam"),
                
                ps_16S_bacteria |>
                        sample_sums() |>
                        tibble::enframe(name = "sample_id",
                                        value = "n_reads") |>
                        dplyr::mutate(stage = "Bacteria only")
        ) |>
        dplyr::mutate(sample_id = stringr::str_trim(as.character(sample_id))) |>
        tidyr::complete(
                sample_id,
                stage,
                fill = list(n_reads = 0)
        ) |> 
        dplyr::group_by(sample_id) |>
        dplyr::mutate(
                raw_reads = n_reads[stage == "Raw"],
                after_decontam_reads = n_reads[stage == "After decontam"],
                bacteria_only_reads = n_reads[stage == "Bacteria only"],
                
                fraction_removed_decontam =
                        (raw_reads - after_decontam_reads) / raw_reads,
                
                fraction_removed_off_target =
                        (after_decontam_reads - bacteria_only_reads) / raw_reads,
                
                fraction_retained =
                        bacteria_only_reads / raw_reads
        ) |> 
        dplyr::ungroup() |>
        dplyr::left_join(
                ps_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        dplyr::mutate(sample_id = stringr::str_trim(as.character(sample_id))),
                by = "sample_id"
        ) |>
        tidyr::pivot_longer(
                cols = c(
                        fraction_removed_decontam,
                        fraction_removed_off_target,
                        fraction_retained
                ),
                names_to = "component",
                values_to = "fraction_of_raw"
        ) |>
        dplyr::mutate(
                component = factor(
                        component,
                        levels = c(
                                "fraction_removed_decontam",
                                "fraction_removed_off_target",
                                "fraction_retained"
                        ),
                        labels = c(
                                "Removed by decontam",
                                "Removed as off-target",
                                "Retained (Bacteria only)"
                        )
                )
        )

# Helper function to plot decontamination and off-target filtering results
plot_read_retention_stacked_16S <-
        function(grouping_variable,
                 grouping_levels = NULL,
                 facet_variable  = NULL,
                 facet_levels    = NULL,
                 filter_function = NULL) {

                group_vars <- c(grouping_variable, "component")
                if (!is.null(facet_variable)) group_vars <- c(facet_variable, group_vars)

                p <- read_retention_loss_table_16S |>
                        (\(x) {
                                if (!is.null(filter_function)) {
                                        filter_function(x)
                                } else {
                                        x
                                }
                        })() |>
                        (\(x) {
                                if (!is.null(grouping_levels)) {
                                        dplyr::mutate(
                                                x,
                                                !!grouping_variable :=
                                                        factor(
                                                                .data[[grouping_variable]],
                                                                levels = grouping_levels
                                                        )
                                        )
                                } else {
                                        x
                                }
                        })() |>
                        (\(x) {
                                if (!is.null(facet_variable) && !is.null(facet_levels)) {
                                        dplyr::mutate(
                                                x,
                                                !!facet_variable :=
                                                        factor(
                                                                .data[[facet_variable]],
                                                                levels = facet_levels
                                                        )
                                        )
                                } else {
                                        x
                                }
                        })() |>
                        dplyr::group_by(across(all_of(group_vars))) |>
                        dplyr::summarise(
                                mean_fraction_of_raw =
                                        mean(fraction_of_raw, na.rm = TRUE),
                                .groups = "drop"
                        ) |>
                        ggplot(
                                aes(
                                        x = .data[[grouping_variable]],
                                        y = mean_fraction_of_raw,
                                        fill = component
                                )
                        ) +
                        geom_bar(stat = "identity") +
                        labs(
                                x = if (!is.null(facet_variable))
                                        paste0(grouping_variable, " and ", facet_variable)
                                else
                                        grouping_variable,
                                y = "Mean fraction of raw reads",
                                fill = "Processing outcome"
                        ) +
                        theme_bw() +
                        theme(
                                axis.text.x = element_text(
                                        angle = 45,
                                        hjust = 1
                                ),
                                legend.position = "right"
                        )

                if (!is.null(facet_variable)) {
                        p <- p +
                                facet_grid(
                                        . ~ .data[[facet_variable]],
                                        scales = "free_x",
                                        space  = "free_x",
                                        switch = "x"
                                ) +
                                theme(
                                        strip.placement  = "outside",
                                        strip.background = element_blank(),
                                        axis.text.x      = element_text(
                                                angle = 45, hjust = 1,
                                                margin = margin(t = 5)
                                        )
                                )
                }

                p
        }

# ---- Do decontamination and off-target sequences removal affect experimentla groups and phenotypes differently? ----
# Plot decontamination and off-target filtering results by sublibrary_id
p_ret_sublibrary <- plot_read_retention_stacked_16S(grouping_variable = "sublibrary_id")
ggsave(here("results", "figures", "04a_16S_read_retention_by_sublibrary.png"),
       plot = p_ret_sublibrary, width = 10, height = 5, dpi = 300)

# Plot decontamination and off-target filtering results by sample type
p_ret_sample_type <- plot_read_retention_stacked_16S(
        grouping_variable = "sample_type",
        grouping_levels = c("extraction_control",
                            "negative_pcr_control",
                            "mock_community",
                            "bagged_flower",
                            "unbagged_flower"
        )
)
ggsave(here("results", "figures", "04a_16S_read_retention_by_sample_type.png"),
       plot = p_ret_sample_type, width = 8, height = 5, dpi = 300)

# Plot decontamination and off-target filtering results by farm_id with management_type facets (only biological samples)
p_ret_farm_mgmt <- plot_read_retention_stacked_16S(
        grouping_variable = "farm_id",
        grouping_levels   = c("ib", "vr", "sa", "kk", "mt", "vi", "yb"),
        facet_variable    = "management_type",
        facet_levels      = c("inside_forest", "near_forest", "agroforest", "full_sun"),
        filter_function   =
                \(x) dplyr::filter(
                        x,
                        sample_type %in% c("bagged_flower", "unbagged_flower") &
                                !is.na(sample_type)
                )
)
ggsave(here("results", "figures", "04a_16S_read_retention_by_farm_management.png"),
       plot = p_ret_farm_mgmt, width = 10, height = 5, dpi = 300)

# Plot decontamination and off-target filtering results by is_pollination_clsm (only phenotyped flowers)
p_ret_pollination <- plot_read_retention_stacked_16S(
        grouping_variable = "is_pollination_clsm",
        grouping_levels = c(FALSE, TRUE),
        filter_function =
                \(x) dplyr::filter(
                        x,
                        !is.na(is_pollination_clsm)
                )
)
ggsave(here("results", "figures", "04a_16S_read_retention_by_pollination.png"),
       plot = p_ret_pollination, width = 6, height = 5, dpi = 300)

# Plot decontamination and off-target filtering results by is_germination (only phenotyped flowers)
p_ret_germination <- plot_read_retention_stacked_16S(
        grouping_variable = "is_germination",
        grouping_levels = c(FALSE, TRUE),
        filter_function =
                \(x) dplyr::filter(
                        x,
                        !is.na(is_germination)
                )
)
ggsave(here("results", "figures", "04a_16S_read_retention_by_germination.png"),
       plot = p_ret_germination, width = 6, height = 5, dpi = 300)

# Does sequencing depth confound pollination intensity? No
p_scatter_pi <- read_retention_loss_table_16S |>
        dplyr::filter(
                !is.na(pi_clsm),
                pi_clsm > 0,
                component == "Retained (Bacteria only)"
        ) |>
        ggplot(
                aes(
                        x = pi_clsm,
                        y = fraction_of_raw
                )
        ) +
        geom_point(alpha = 0.6) +
        # geom_smooth(method = "loess", se = TRUE) +
        labs(
                x = "Number of pollen grains (pi_clsm)",
                y = "Fraction of raw reads retained"
        ) +
        theme_bw()
ggsave(here("results", "figures", "04a_16S_depth_vs_pi_clsm.png"),
       plot = p_scatter_pi, width = 6, height = 5, dpi = 300)

# Does sequencing depth confound callose_plug_index? No
p_scatter_callose <- read_retention_loss_table_16S |>
        dplyr::filter(
                !is.na(callose_plug_index),
                callose_plug_index > 0,
                component == "Retained (Bacteria only)"
        ) |>
        ggplot(
                aes(
                        x = callose_plug_index,
                        y = fraction_of_raw
                )
        ) +
        geom_point(alpha = 0.6) +
        # geom_smooth(method = "loess", se = TRUE) +
        labs(
                x = "Callose plug index",
                y = "Fraction of raw reads retained"
        ) +
        theme_bw()
ggsave(here("results", "figures", "04a_16S_depth_vs_callose_plug_index.png"),
       plot = p_scatter_callose, width = 6, height = 5, dpi = 300)


# ---- Are my experimental and phenotype groups confounded by sequencing depth? ----
# sublibrary_id
p_box_sublibrary <- dplyr::bind_rows(
        ps_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |>
        dplyr::left_join(
                ps_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })(),
                by = "sample_id"
        ) |>
        ggplot(aes(x = sublibrary_id, y = n_reads, colour = sample_type)) +
        geom_boxplot() +
        theme_bw() +
        theme(
                axis.text.x = element_text(
                        angle = 45,
                        hjust = 1
                ),
                legend.position = "bottom"
        )
ggsave(here("results", "figures", "04a_16S_depth_boxplot_by_sublibrary.png"),
       plot = p_box_sublibrary, width = 10, height = 5, dpi = 300)

# sample_type
p_box_sample_type <- dplyr::bind_rows(
        ps_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |>
        dplyr::left_join(
                ps_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })(),
                by = "sample_id"
        ) |>
        dplyr::mutate(sample_type = factor(sample_type,
                                           levels = c("extraction_control",
                                                      "negative_pcr_control",
                                                      "mock_community",
                                                      "bagged_flower",
                                                      "unbagged_flower"
                                           )
        )) |>
        ggplot(aes(x = sample_type, y = n_reads)) +
        geom_boxplot() +
        theme_bw()
ggsave(here("results", "figures", "04a_16S_depth_boxplot_by_sample_type.png"),
       plot = p_box_sample_type, width = 7, height = 5, dpi = 300)

# bagged Vs unbagged
p_box_treatment <- dplyr::bind_rows(
        ps_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |>
        dplyr::left_join(
                ps_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })(),
                by = "sample_id"
        ) |>
        filter(sample_type %in% c("bagged_flower", "unbagged_flower")) |>
        ggplot(aes(x = sample_type, y = n_reads)) +
        geom_boxplot() +
        theme_bw()
ggsave(here("results", "figures", "04a_16S_depth_boxplot_by_flower_treatment.png"),
       plot = p_box_treatment, width = 5, height = 5, dpi = 300)

# farm_id by management_type (faceted, no legend)
p_box_farm <- dplyr::bind_rows(
        ps_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |>
        dplyr::left_join(
                ps_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })(),
                by = "sample_id"
        ) |>
        filter(sample_type %in% c("bagged_flower", "unbagged_flower")) |>
        mutate(
                farm_id = factor(farm_id,
                                 levels = c("ib", "vr", "sa", "kk", "mt", "vi", "yb")),
                management_type = factor(management_type,
                                         levels = c("inside_forest", "near_forest",
                                                    "agroforest", "full_sun"))
        ) |>
        ggplot(aes(x = farm_id, y = n_reads)) +
        geom_boxplot() +
        facet_grid(
                . ~ management_type,
                scales = "free_x",
                space  = "free_x",
                switch = "x"
        ) +
        labs(x = "Farm ID and management type", y = "n_reads") +
        theme_bw() +
        theme(
                strip.placement  = "outside",
                strip.background = element_blank(),
                axis.text.x      = element_text(margin = margin(t = 5))
        )
ggsave(here("results", "figures", "04a_16S_depth_boxplot_by_farm.png"),
       plot = p_box_farm, width = 9, height = 5, dpi = 300)

# sequencing depth by farm_id and sample_type (bagged vs unbagged within each farm)
p_depth_farm_type <- dplyr::bind_rows(
        ps_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id", value = "n_reads")
) |>
        dplyr::left_join(
                ps_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })(),
                by = "sample_id"
        ) |>
        dplyr::filter(sample_type %in% c("bagged_flower", "unbagged_flower")) |>
        dplyr::mutate(
                farm_id         = factor(farm_id,
                                         levels = c("ib", "vr", "sa", "kk", "mt", "vi", "yb")),
                management_type = factor(management_type,
                                         levels = c("inside_forest", "near_forest",
                                                    "agroforest", "full_sun")),
                sample_type     = factor(sample_type,
                                         levels = c("bagged_flower", "unbagged_flower"))
        ) |>
        ggplot(aes(x = farm_id, y = n_reads, fill = sample_type)) +
        geom_boxplot(position = position_dodge(width = 0.8)) +
        facet_grid(
                . ~ management_type,
                scales = "free_x",
                space  = "free_x",
                switch = "x"
        ) +
        labs(x = "Farm ID and management type", y = "Reads (bacteria only)",
             fill = "Sample type") +
        theme_bw() +
        theme(
                strip.placement  = "outside",
                strip.background = element_blank(),
                axis.text.x      = element_text(margin = margin(t = 5)),
                legend.position  = "bottom"
        )
ggsave(here("results", "figures", "04a_16S_depth_boxplot_by_farm_and_sample_type.png"),
       plot = p_depth_farm_type, width = 9, height = 5, dpi = 300)

# is_pollination_clsm
p_box_pollination <- dplyr::bind_rows(
        ps_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |>
        dplyr::left_join(
                ps_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })(),
                by = "sample_id"
        ) |>
        filter(!is.na(is_pollination_clsm)) |>
        ggplot(aes(x = is_pollination_clsm, y = n_reads)) +
        geom_boxplot() +
        theme_bw()
ggsave(here("results", "figures", "04a_16S_depth_boxplot_by_pollination.png"),
       plot = p_box_pollination, width = 5, height = 5, dpi = 300)

# is_germination
p_box_germination <- dplyr::bind_rows(
        ps_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |>
        dplyr::left_join(
                ps_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })(),
                by = "sample_id"
        ) |>
        filter(!is.na(is_germination)) |>
        ggplot(aes(x = is_germination, y = n_reads)) +
        geom_boxplot() +
        theme_bw()
ggsave(here("results", "figures", "04a_16S_depth_boxplot_by_germination.png"),
       plot = p_box_germination, width = 5, height = 5, dpi = 300)

# ---- Is there correlation between DNA concentration and sequencing depth? ----
# DNA was measured using Qubit after locus-specific PCR, bead cleaning & pooling triplicates
# Sequencing depth here is the number of reads per sample in data set after decontamination and off-targer sequences filtering
p_scatter_dna <- dplyr::bind_rows(
        ps_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |>
        dplyr::left_join(
                ps_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })(),
                by = "sample_id"
        ) |>
        ggplot(aes(x = dna_concentration_hs_16s , y = n_reads, color = sample_type)) +
        geom_point() +
        theme_bw()
ggsave(here("results", "figures", "04a_16S_depth_vs_dna_concentration.png"),
       plot = p_scatter_dna, width = 7, height = 5, dpi = 300)

message("04a done")
