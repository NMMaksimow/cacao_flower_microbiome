# ---- Packages ----
library(here)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(scales)

# Data version switch: CFM_VERSION, EXPORT_DIR, RESULTS_DIR
source(here("scripts", "utils_paths.R"))

# ---- Load objects ----
# The raw objects hold every sequenced sample (controls and all 294 flowers),
# so this QC covers the sampling deviations as well
ps_16S_raw      <- readRDS(file.path(RESULTS_DIR, "rds", "ps_16S_raw.rds"))
ps_16S_decontam <- readRDS(file.path(RESULTS_DIR, "rds", "ps_16S_decontam.rds"))
ps_16S_bacteria <- readRDS(file.path(RESULTS_DIR, "rds", "ps_16S_bacteria.rds"))

# ITS1 at OTU97, the level used by all later ITS1 steps
ps_ITS1_raw      <- readRDS(file.path(RESULTS_DIR, "rds", "ps_ITS1_otu97_raw.rds"))
ps_ITS1_decontam <- readRDS(file.path(RESULTS_DIR, "rds", "ps_ITS1_otu97_decontam.rds"))
ps_ITS1_fungi    <- readRDS(file.path(RESULTS_DIR, "rds", "ps_ITS1_otu97_fungi.rds"))

dir.create(file.path(RESULTS_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

################################################################################
# Shared constants
################################################################################

# ggplot2 default hue palette for 5 sample types (alphabetical order):
# bagged_flower, extraction_control, mock_community, negative_pcr_control, unbagged_flower
SAMPLE_TYPE_PALETTE <- c(
        "bagged_flower"        = "#F8766D",
        "unbagged_flower"      = "#E76BF3",
        "mock_community"       = "#00BF7D",
        "extraction_control"   = "#A3A500",
        "negative_pcr_control" = "#00B0F6")

# Applied to ITS1 (right-side) panels so y-axis tick labels are not duplicated.
# Every figure gives 16S and ITS1 the same y range (0-100% for fractions, a shared
# read-depth ceiling for depths), so the 16S (left) axis is valid for both.
HIDE_Y_AXIS <- theme(axis.text.y  = element_blank(),
                     axis.ticks.y = element_blank())

COMPONENT_LEVELS <- c("fraction_removed_decontam",
                       "fraction_removed_off_target",
                       "fraction_retained")
COMPONENT_LABELS <- c("Removed by decontam",
                       "Removed as off-target",
                       "Retained")

################################################################################
# Data-building helpers
################################################################################

make_retention_table <- function(ps_raw, ps_decontam, ps_final) {
        bind_rows(
                ps_raw      |> sample_sums() |> enframe("sample_id", "n_reads") |> mutate(stage = "Raw"),
                ps_decontam |> sample_sums() |> enframe("sample_id", "n_reads") |> mutate(stage = "After decontam"),
                ps_final    |> sample_sums() |> enframe("sample_id", "n_reads") |> mutate(stage = "Final")
        ) |>
                mutate(sample_id = str_trim(as.character(sample_id))) |>
                complete(sample_id, stage, fill = list(n_reads = 0)) |>
                group_by(sample_id) |>
                mutate(
                        raw_reads            = n_reads[stage == "Raw"],
                        after_decontam_reads = n_reads[stage == "After decontam"],
                        final_reads          = n_reads[stage == "Final"],
                        fraction_removed_decontam   = (raw_reads - after_decontam_reads) / raw_reads,
                        fraction_removed_off_target = (after_decontam_reads - final_reads) / raw_reads,
                        fraction_retained           = final_reads / raw_reads
                ) |>
                ungroup() |>
                left_join(
                        ps_raw |> sample_data() |> as.data.frame() |>
                                (\(x) { class(x) <- "data.frame"; x })() |>
                                mutate(sample_id = str_trim(as.character(sample_id))),
                        by = "sample_id"
                ) |>
                pivot_longer(
                        cols      = all_of(COMPONENT_LEVELS),
                        names_to  = "component",
                        values_to = "fraction_of_raw"
                ) |>
                mutate(component = factor(component,
                                          levels = COMPONENT_LEVELS,
                                          labels = COMPONENT_LABELS))
}

make_depth_table <- function(ps_final, ps_raw) {
        ps_final |> sample_sums() |>
                enframe("sample_id", "n_reads") |>
                left_join(
                        ps_raw |> sample_data() |> as.data.frame() |>
                                (\(x) { class(x) <- "data.frame"; x })(),
                        by = "sample_id"
                )
}

################################################################################
# Plot helpers
#
# make_stack parameters:
#   hide_x_axis  – removes tick labels, ticks, and axis title (use for top rows)
#   x_angle      – text rotation (0 = horizontal)
#   x_label      – x-axis title
#   x_scale_labels – passed to scale_x_discrete(labels=)
#
# Y scale is always 0–100% so panels are directly comparable across figures.
#
# make_box parameters:
#   depth_ylim   – numeric: upper y limit (coord_cartesian); pass the same value
#                  to the 16S and ITS1 panels so their depths share one scale
#   hline        – numeric: draw horizontal red dashed reference line
#   show_n       – TRUE: print "n=X" above each box
#
# make_scatter parameters:
#   y_ylim        – numeric c(min, max): y limits shared by the 16S and ITS1 panels
#   colour_levels – sample types listed in the legend (all panels of a figure
#                   must use the same set, so patchwork merges their legends)
#   jitter_width  – horizontal jitter of the points (for an ordinal x, e.g. a score)
################################################################################

make_stack <- function(tbl,
                       grouping_variable,
                       grouping_levels  = NULL,
                       facet_variable   = NULL,
                       facet_levels     = NULL,
                       filter_fn        = NULL,
                       x_angle          = 45,
                       x_label          = NULL,
                       x_scale_labels   = waiver(),
                       hide_x_axis      = FALSE) {

        group_vars <- c(grouping_variable, "component")
        if (!is.null(facet_variable)) group_vars <- c(facet_variable, group_vars)
        hjust_val <- if (x_angle > 0) 1 else 0.5

        p <- tbl |>
                (\(x) if (!is.null(filter_fn)) filter_fn(x) else x)() |>
                (\(x) if (!is.null(grouping_levels))
                        mutate(x, !!grouping_variable :=
                                       factor(.data[[grouping_variable]], levels = grouping_levels))
                else x)() |>
                (\(x) if (!is.null(facet_variable) && !is.null(facet_levels))
                        mutate(x, !!facet_variable :=
                                       factor(.data[[facet_variable]], levels = facet_levels))
                else x)() |>
                group_by(across(all_of(group_vars))) |>
                summarise(mean_fraction = mean(fraction_of_raw, na.rm = TRUE),
                          .groups = "drop") |>
                ggplot(aes(x    = .data[[grouping_variable]],
                           y    = mean_fraction,
                           fill = component)) +
                geom_bar(stat = "identity") +
                scale_x_discrete(labels = x_scale_labels) +
                scale_y_continuous(
                        labels = percent_format(),
                        limits = c(0, 1),
                        expand = expansion(mult = c(0, 0.04))
                ) +
                labs(x    = if (hide_x_axis) NULL else (if (is.null(x_label)) grouping_variable else x_label),
                     y    = "Mean fraction of raw reads",
                     fill = "Processing outcome") +
                theme_bw() +
                theme(axis.text.x     = element_text(angle = x_angle, hjust = hjust_val),
                      legend.position = "right")

        if (hide_x_axis) {
                p <- p + theme(axis.text.x  = element_blank(),
                               axis.ticks.x = element_blank())
        }

        if (!is.null(facet_variable)) {
                p <- p +
                        facet_grid(. ~ .data[[facet_variable]],
                                   scales = "free_x", space = "free_x", switch = "x") +
                        theme(strip.placement  = "outside",
                              strip.background = element_blank())
                if (!hide_x_axis) {
                        p <- p + theme(axis.text.x = element_text(angle  = x_angle,
                                                                   hjust  = hjust_val,
                                                                   margin = margin(t = 5)))
                }
        }
        p
}

make_box <- function(depth_tbl,
                     grouping_variable,
                     depth_ylim,
                     grouping_levels     = NULL,
                     facet_variable      = NULL,
                     facet_levels        = NULL,
                     filter_fn           = NULL,
                     x_angle             = 45,
                     x_label             = NULL,
                     x_scale_labels      = waiver(),
                     hline               = NULL,
                     hline_colour        = "#D55E00",
                     show_n              = FALSE,
                     fill_by_sample_type = FALSE) {
                     # fill_by_sample_type = TRUE: boxes colored by sample_type using
                     # SAMPLE_TYPE_PALETTE; each x group then shows one box per sample
                     # type present (e.g. bagged and unbagged side by side per farm).

        hjust_val <- if (x_angle > 0) 1 else 0.5

        p <- depth_tbl |>
                (\(x) if (!is.null(filter_fn)) filter_fn(x) else x)() |>
                (\(x) if (!is.null(grouping_levels))
                        mutate(x, !!grouping_variable :=
                                       factor(.data[[grouping_variable]], levels = grouping_levels))
                else x)() |>
                (\(x) if (!is.null(facet_variable) && !is.null(facet_levels))
                        mutate(x, !!facet_variable :=
                                       factor(.data[[facet_variable]], levels = facet_levels))
                else x)() |>
                ggplot(aes(x = .data[[grouping_variable]], y = n_reads)) +
                geom_boxplot() +
                scale_x_discrete(labels = x_scale_labels) +
                # coord_cartesian zooms without dropping points from the box statistics
                coord_cartesian(ylim = c(0, depth_ylim)) +
                labs(x = if (is.null(x_label)) grouping_variable else x_label,
                     y = "Read depth (final)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = x_angle, hjust = hjust_val))

        # Color boxes by sample_type and add to shared legend
        if (fill_by_sample_type) {
                p <- p +
                        aes(fill = sample_type) +
                        scale_fill_manual(values = SAMPLE_TYPE_PALETTE,
                                          name   = "Sample type")
        }

        if (!is.null(hline)) {
                p <- p + geom_hline(yintercept  = hline,
                                    colour      = hline_colour,
                                    linetype    = "dashed",
                                    linewidth   = 0.7)
        }

        if (show_n) {
                p <- p + stat_summary(
                        geom     = "text",
                        fun.data = function(y) data.frame(
                                y     = Inf,
                                label = paste0("n=", sum(!is.na(y)))
                        ),
                        vjust = 1.4, size = 2.4
                )
        }

        if (!is.null(facet_variable)) {
                p <- p +
                        facet_grid(. ~ .data[[facet_variable]],
                                   scales = "free_x", space = "free_x", switch = "x") +
                        theme(strip.placement  = "outside",
                              strip.background = element_blank(),
                              axis.text.x      = element_text(angle  = x_angle,
                                                              hjust  = hjust_val,
                                                              margin = margin(t = 5)))
        }
        p
}

# Points colored by sample_type. The legend lists every type in colour_levels
# (limits =, and show.legend = TRUE draws a key even for types absent from a
# panel), so patchwork merges the legends of all panels into one.
make_scatter <- function(tbl,
                         x_variable,
                         y_variable,
                         y_ylim,
                         filter_fn     = NULL,
                         x_label       = x_variable,
                         y_label       = y_variable,
                         y_labels      = waiver(),
                         colour_levels = names(SAMPLE_TYPE_PALETTE),
                         jitter_width  = 0) {
        tbl |>
                (\(x) if (!is.null(filter_fn)) filter_fn(x) else x)() |>
                ggplot(aes(x      = .data[[x_variable]],
                           y      = .data[[y_variable]],
                           colour = sample_type)) +
                geom_point(alpha = 0.6, na.rm = TRUE,      # alpha: point transparency
                           # seed: the same jitter on every run
                           position    = position_jitter(width = jitter_width, height = 0,
                                                         seed = 1),
                           show.legend = TRUE) +
                scale_colour_manual(values = SAMPLE_TYPE_PALETTE,
                                    limits = colour_levels,
                                    name   = "Sample type") +
                scale_y_continuous(labels = y_labels) +
                coord_cartesian(ylim = y_ylim) +
                labs(x = x_label, y = y_label) +
                theme_bw()
}

# Shared read-depth ceiling for the 16S and ITS1 panels of one figure:
# the deepest sample of either marker (after filter_fn) plus 10% headroom
# for the "n=" labels.
shared_depth_ylim <- function(filter_fn = identity) {
        c(filter_fn(depth_16S)$n_reads, filter_fn(depth_ITS1)$n_reads) |>
                max(na.rm = TRUE) * 1.10
}

# Standard 2×2 composite: stacked bars (top) / boxplots (bottom),
# 16S (left) / ITS1 (right).
# Stacked bar x-axis is hidden; y-label on left panels only.
# Panel tags A-D run left to right, top to bottom (plot_annotation tag_levels).
make_composite <- function(stack_16S, stack_ITS1, box_16S, box_ITS1,
                           filename, width = 14, height = 10) {
        fig <- (
                (stack_16S + ggtitle("16S — Bacteria")) |
                (stack_ITS1 + ggtitle("ITS1 — Fungi") + labs(y = NULL) + HIDE_Y_AXIS)
        ) /
        (
                box_16S |
                (box_ITS1 + labs(y = NULL) + HIDE_Y_AXIS)
        ) +
        plot_layout(guides = "collect") +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "right")    # collected legends to the right of all panels

        print(fig)
        ggsave(file.path(RESULTS_DIR, "figures", filename),
               plot = fig, width = width, height = height, dpi = 300)
        invisible(fig)
}

################################################################################
# Build data tables
################################################################################

tbl_16S  <- make_retention_table(ps_16S_raw,  ps_16S_decontam,  ps_16S_bacteria)
tbl_ITS1 <- make_retention_table(ps_ITS1_raw, ps_ITS1_decontam, ps_ITS1_fungi)

depth_16S  <- make_depth_table(ps_16S_bacteria, ps_16S_raw)
depth_ITS1 <- make_depth_table(ps_ITS1_fungi,   ps_ITS1_raw)

################################################################################
# Figure 1: By sequencing sublibrary
#   Stack  — mean fraction retained/removed per sublibrary; no x-axis labels
#            (hidden to save space; labels shown on boxplot panel below)
#   Boxplot — read depth per sample colored by sample_type; shared y-ceiling;
#             red dashed line at 10 000 reads (minimum acceptable depth)
################################################################################

subl_ylim <- shared_depth_ylim()

# Sublibrary labels are long and many: vertical and small so they do not overlap
SUBLIBRARY_X_TEXT <- theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                                      vjust = 0.5, size = 6))

make_composite(
        stack_16S  = make_stack(tbl_16S,  "sublibrary_id", hide_x_axis = TRUE),
        stack_ITS1 = make_stack(tbl_ITS1, "sublibrary_id", hide_x_axis = TRUE),
        # Bottom row: boxplots colored by sample_type; each sublibrary may
        # contain multiple sample types shown as separate colored boxes
        box_16S    = make_box(depth_16S,  "sublibrary_id", depth_ylim = subl_ylim,
                              x_label = "Sublibrary", hline = 10000,
                              fill_by_sample_type = TRUE) + SUBLIBRARY_X_TEXT,
        box_ITS1   = make_box(depth_ITS1, "sublibrary_id", depth_ylim = subl_ylim,
                              x_label = "Sublibrary", hline = 10000,
                              fill_by_sample_type = TRUE) + SUBLIBRARY_X_TEXT,
        filename   = "04_read_retention_sublibrary.png",
        width = 16, height = 10
)

################################################################################
# Figure 2: By sample type
#   Stack  — mean retention profile per sample type category
#   Boxplot — read depth distribution; shows whether controls vs. flowers differ
################################################################################

sample_type_levels <- c("extraction_control", "negative_pcr_control",
                        "mock_community", "bagged_flower", "unbagged_flower")
sample_type_ylim   <- shared_depth_ylim()

make_composite(
        stack_16S  = make_stack(tbl_16S,  "sample_type",
                                grouping_levels = sample_type_levels,
                                x_angle = 0, hide_x_axis = TRUE),
        stack_ITS1 = make_stack(tbl_ITS1, "sample_type",
                                grouping_levels = sample_type_levels,
                                x_angle = 0, hide_x_axis = TRUE),
        box_16S    = make_box(depth_16S,  "sample_type",
                              depth_ylim      = sample_type_ylim,
                              grouping_levels = sample_type_levels,
                              x_angle = 0, x_label = "Sample type",
                              fill_by_sample_type = TRUE),
        box_ITS1   = make_box(depth_ITS1, "sample_type",
                              depth_ylim      = sample_type_ylim,
                              grouping_levels = sample_type_levels,
                              x_angle = 0, x_label = "Sample type",
                              fill_by_sample_type = TRUE),
        filename   = "04_read_retention_sample_type.png",
        width = 14, height = 10
)

################################################################################
# Figure 3: Farm × management type
#   Stack  — retention profile per farm; faceted by management type (strip at bottom)
#   Boxplot — read depth per farm, bagged and unbagged side by side; checks for
#             sequencing depth confounds with treatment, management type or geography
################################################################################

farm_levels <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
mgmt_levels <- c("inside_forest", "near_forest", "agroforest", "full_sun")
bio_fn      <- \(x) filter(x, sample_type %in% c("bagged_flower", "unbagged_flower"))
farm_ylim   <- shared_depth_ylim(bio_fn)

make_composite(
        stack_16S  = make_stack(tbl_16S, "farm_id",
                                grouping_levels = farm_levels,
                                facet_variable  = "management_type",
                                facet_levels    = mgmt_levels,
                                filter_fn       = bio_fn,
                                x_angle = 0, hide_x_axis = TRUE),
        stack_ITS1 = make_stack(tbl_ITS1, "farm_id",
                                grouping_levels = farm_levels,
                                facet_variable  = "management_type",
                                facet_levels    = mgmt_levels,
                                filter_fn       = bio_fn,
                                x_angle = 0, hide_x_axis = TRUE),
        box_16S    = make_box(depth_16S, "farm_id",
                              depth_ylim      = farm_ylim,
                              grouping_levels = farm_levels,
                              facet_variable  = "management_type",
                              facet_levels    = mgmt_levels,
                              filter_fn       = bio_fn,
                              x_angle = 0, x_label = "Farm",
                              fill_by_sample_type = TRUE),
        box_ITS1   = make_box(depth_ITS1, "farm_id",
                              depth_ylim      = farm_ylim,
                              grouping_levels = farm_levels,
                              facet_variable  = "management_type",
                              facet_levels    = mgmt_levels,
                              filter_fn       = bio_fn,
                              x_angle = 0, x_label = "Farm",
                              fill_by_sample_type = TRUE),
        filename   = "04_read_retention_farm_management.png",
        width = 18, height = 10
)

################################################################################
# Figure 4: Phenotype variables
#   Columns (L→R): Flower treatment | Pollination (CLSM) | Germination | Fungal hyphae
#   All 16S panels left (cols 1–4), all ITS1 panels right (cols 5–8)
#   Top row  — stacked bars: retention profile per phenotype group; x hidden
#   Bottom row — boxplots: read depth by group on one shared scale; x-title
#              shown on all panels so "Yes / No" labels are interpretable;
#              n= per group shown above each box
################################################################################

yn_labels   <- c("FALSE" = "No", "TRUE" = "Yes")
bool_levels <- c("FALSE", "TRUE")

flower_fn_ph <- \(x) filter(x, sample_type %in% c("bagged_flower", "unbagged_flower"))
poll_fn  <- \(x) filter(x, !is.na(is_pollination_clsm)) |>
        mutate(is_pollination_clsm = as.character(is_pollination_clsm))
germ_fn  <- \(x) filter(x, !is.na(is_germination)) |>
        mutate(is_germination = as.character(is_germination))
hyph_fn  <- \(x) filter(x, !is.na(is_hyphae_present)) |>
        mutate(is_hyphae_present = as.character(is_hyphae_present))

flower_lvl <- c("bagged_flower", "unbagged_flower")
flower_lbl <- c("bagged_flower" = "Bagged", "unbagged_flower" = "Unbagged")

# All phenotyped samples are flowers, so the flower ceiling covers every panel
pheno_ylim <- shared_depth_ylim(flower_fn_ph)

# ── Stacked bars (hide x-axis on all panels) ──────────────────────────────────
s16_flow <- make_stack(tbl_16S,  "sample_type", grouping_levels = flower_lvl,
                       filter_fn = flower_fn_ph, x_angle = 0,
                       x_scale_labels = flower_lbl, hide_x_axis = TRUE)
sIT_flow <- make_stack(tbl_ITS1, "sample_type", grouping_levels = flower_lvl,
                       filter_fn = flower_fn_ph, x_angle = 0,
                       x_scale_labels = flower_lbl, hide_x_axis = TRUE)

s16_poll <- make_stack(tbl_16S,  "is_pollination_clsm", grouping_levels = bool_levels,
                       filter_fn = poll_fn, x_angle = 0,
                       x_scale_labels = yn_labels, hide_x_axis = TRUE)
sIT_poll <- make_stack(tbl_ITS1, "is_pollination_clsm", grouping_levels = bool_levels,
                       filter_fn = poll_fn, x_angle = 0,
                       x_scale_labels = yn_labels, hide_x_axis = TRUE)

s16_germ <- make_stack(tbl_16S,  "is_germination", grouping_levels = bool_levels,
                       filter_fn = germ_fn, x_angle = 0,
                       x_scale_labels = yn_labels, hide_x_axis = TRUE)
sIT_germ <- make_stack(tbl_ITS1, "is_germination", grouping_levels = bool_levels,
                       filter_fn = germ_fn, x_angle = 0,
                       x_scale_labels = yn_labels, hide_x_axis = TRUE)

s16_hyph <- make_stack(tbl_16S,  "is_hyphae_present", grouping_levels = bool_levels,
                       filter_fn = hyph_fn, x_angle = 0,
                       x_scale_labels = yn_labels, hide_x_axis = TRUE)
sIT_hyph <- make_stack(tbl_ITS1, "is_hyphae_present", grouping_levels = bool_levels,
                       filter_fn = hyph_fn, x_angle = 0,
                       x_scale_labels = yn_labels, hide_x_axis = TRUE)

# ── Boxplots (x-title on every panel; n= per group; shared depth scale) ──────
b16_flow <- make_box(depth_16S,  "sample_type", depth_ylim = pheno_ylim,
                     grouping_levels = flower_lvl,
                     filter_fn = flower_fn_ph, x_angle = 0, show_n = TRUE,
                     x_scale_labels = flower_lbl, x_label = "Flower treatment")
bIT_flow <- make_box(depth_ITS1, "sample_type", depth_ylim = pheno_ylim,
                     grouping_levels = flower_lvl,
                     filter_fn = flower_fn_ph, x_angle = 0, show_n = TRUE,
                     x_scale_labels = flower_lbl, x_label = "Flower treatment")

b16_poll <- make_box(depth_16S,  "is_pollination_clsm", depth_ylim = pheno_ylim,
                     grouping_levels = bool_levels,
                     filter_fn = poll_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Pollination (CLSM)")
bIT_poll <- make_box(depth_ITS1, "is_pollination_clsm", depth_ylim = pheno_ylim,
                     grouping_levels = bool_levels,
                     filter_fn = poll_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Pollination (CLSM)")

b16_germ <- make_box(depth_16S,  "is_germination", depth_ylim = pheno_ylim,
                     grouping_levels = bool_levels,
                     filter_fn = germ_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Pollen germination")
bIT_germ <- make_box(depth_ITS1, "is_germination", depth_ylim = pheno_ylim,
                     grouping_levels = bool_levels,
                     filter_fn = germ_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Pollen germination")

b16_hyph <- make_box(depth_16S,  "is_hyphae_present", depth_ylim = pheno_ylim,
                     grouping_levels = bool_levels,
                     filter_fn = hyph_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Fungal hyphae present")
bIT_hyph <- make_box(depth_ITS1, "is_hyphae_present", depth_ylim = pheno_ylim,
                     grouping_levels = bool_levels,
                     filter_fn = hyph_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Fungal hyphae present")

# ── Assemble: all 16S left (cols 1–4), all ITS1 right (cols 5–8) ─────────────
# One wrap_plots per row (not nested blocks) so patchwork allocates equal
# pixel width to every column regardless of y-axis presence on col 1.
# y-label on col 1 (s16_flow / b16_flow) only; all other panels hide y-axis.
# Panel tags are set by hand (labs(tag =)) on the first panel of each 4-panel
# block: A 16S retention, B ITS1 retention, C 16S depth, D ITS1 depth.

top_row <- wrap_plots(
        s16_flow + ggtitle("16S — Bacteria") + labs(tag = "A"),   # col 1: y-label shown
        s16_poll + labs(y = NULL) + HIDE_Y_AXIS,
        s16_germ + labs(y = NULL) + HIDE_Y_AXIS,
        s16_hyph + labs(y = NULL) + HIDE_Y_AXIS,
        sIT_flow + ggtitle("ITS1 — Fungi") + labs(y = NULL, tag = "B") + HIDE_Y_AXIS,
        sIT_poll + labs(y = NULL) + HIDE_Y_AXIS,
        sIT_germ + labs(y = NULL) + HIDE_Y_AXIS,
        sIT_hyph + labs(y = NULL) + HIDE_Y_AXIS,
        nrow = 1
)

bot_row <- wrap_plots(
        b16_flow + labs(tag = "C"),                                     # col 1: y-label shown
        b16_poll + labs(y = NULL) + HIDE_Y_AXIS,
        b16_germ + labs(y = NULL) + HIDE_Y_AXIS,
        b16_hyph + labs(y = NULL) + HIDE_Y_AXIS,
        bIT_flow + labs(y = NULL, tag = "D") + HIDE_Y_AXIS,
        bIT_poll + labs(y = NULL) + HIDE_Y_AXIS,
        bIT_germ + labs(y = NULL) + HIDE_Y_AXIS,
        bIT_hyph + labs(y = NULL) + HIDE_Y_AXIS,
        nrow = 1
)

fig4 <- (top_row / bot_row) +
        plot_layout(guides = "collect") &
        theme(legend.position = "right")

print(fig4)
ggsave(file.path(RESULTS_DIR, "figures", "04_read_retention_phenotypes.png"),
       plot = fig4, width = 28, height = 10, dpi = 300)

################################################################################
# Figure 5 + table: read depth and DNA yield vs microscopy (flowers only)
#   Rows: number of pollen grains (CLSM) | callose plug index | fungal
#         colonization score (ordinal 0-5, points jittered horizontally)
#   Columns: 16S | ITS1; y = final read depth on one scale for all panels
#   Subtitles: Spearman rho for all flowers, bagged only and unbagged only.
#   Each farm x treatment group was extracted in its own batch, so a correlation
#   within one treatment cannot come from the extraction batch.
#   The table holds the same correlations for DNA concentration as well.
################################################################################

dir.create(file.path(RESULTS_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

MICROSCOPY_VARIABLES <- c(pi_clsm                   = "Number of pollen grains (CLSM)",
                          callose_plug_index        = "Callose plug index",
                          fungal_colonization_score = "Fungal colonization score")

# Microscopy values as numbers (the fungal score is an ordered factor 0-5
# in the phyloseq sample data)
as_number <- \(v) as.numeric(as.character(v))

# One row per flower, marker, response (DNA, depth) and microscopy variable;
# zeros are kept (0 pollen grains and score 0 are real observations)
microscopy_flowers <- bind_rows(
        depth_16S  |> transmute(marker = "16S",  sample_type,
                                dna = dna_concentration_hs_16s, n_reads,
                                across(all_of(names(MICROSCOPY_VARIABLES)), as_number)),
        depth_ITS1 |> transmute(marker = "ITS1", sample_type,
                                dna = dna_concentration_hs_its1, n_reads,
                                across(all_of(names(MICROSCOPY_VARIABLES)), as_number))
) |>
        filter(sample_type %in% flower_lvl) |>
        pivot_longer(c(dna, n_reads), names_to = "response", values_to = "response_value") |>
        pivot_longer(all_of(names(MICROSCOPY_VARIABLES)),
                     names_to = "microscopy", values_to = "microscopy_value") |>
        filter(!is.na(response_value), !is.na(microscopy_value))

# exact = FALSE: normal approximation, as ties make the exact p-value unavailable
microscopy_spearman <- bind_rows(
        microscopy_flowers |> mutate(subset = "all"),
        microscopy_flowers |> mutate(subset = str_remove(sample_type, "_flower"))
) |>
        mutate(subset = factor(subset, levels = c("all", "bagged", "unbagged"))) |>
        group_by(marker, response, microscopy, subset) |>
        group_modify(\(x, key) {
                test <- cor.test(x$microscopy_value, x$response_value,
                                 method = "spearman", exact = FALSE)
                tibble(n = nrow(x), rho = unname(test$estimate), p = test$p.value)
        }) |>
        ungroup()

print(microscopy_spearman, n = Inf)
write_csv(microscopy_spearman,
          file.path(RESULTS_DIR, "tables", "04_microscopy_dna_depth_spearman.csv"))

# Subtitle text: one line per subset, read depth as the response
spearman_label <- \(marker_name, variable_name) microscopy_spearman |>
        filter(marker == marker_name, response == "n_reads", microscopy == variable_name) |>
        arrange(subset) |>
        (\(x) sprintf("%s: rho = %.2f, p = %.1e, n = %d", x$subset, x$rho, x$p, x$n))() |>
        paste(collapse = "\n") |>
        (\(s) paste0("Spearman (read depth)\n", s))()

make_microscopy_panel <- function(depth_tbl, marker_name, variable_name, y_ylim) {
        make_scatter(depth_tbl, variable_name, "n_reads",
                     y_ylim        = y_ylim,
                     # flowers with a value of this variable, as a number (continuous x)
                     filter_fn     = \(x) flower_fn_ph(x) |>
                             filter(!is.na(.data[[variable_name]])) |>
                             mutate(across(all_of(variable_name), as_number)),
                     x_label       = MICROSCOPY_VARIABLES[[variable_name]],
                     y_label       = "Read depth (final)",
                     y_labels      = label_number(big.mark = ""),   # 100000, not 1e+05
                     colour_levels = flower_lvl,
                     # jitter only the ordinal score; 0.15 keeps the score columns apart
                     jitter_width  = if (variable_name == "fungal_colonization_score") 0.15 else 0) +
                labs(subtitle = spearman_label(marker_name, variable_name)) +
                theme(plot.subtitle = element_text(size = 8))   # size of the Spearman summary
}

microscopy_ylim <- c(0, shared_depth_ylim(flower_fn_ph))

# Panel tags A-F run left to right, top to bottom
fig5 <- wrap_plots(
        make_microscopy_panel(depth_16S,  "16S",  "pi_clsm", microscopy_ylim) +
                ggtitle("16S — Bacteria"),
        make_microscopy_panel(depth_ITS1, "ITS1", "pi_clsm", microscopy_ylim) +
                ggtitle("ITS1 — Fungi") + labs(y = NULL) + HIDE_Y_AXIS,
        make_microscopy_panel(depth_16S,  "16S",  "callose_plug_index", microscopy_ylim),
        make_microscopy_panel(depth_ITS1, "ITS1", "callose_plug_index", microscopy_ylim) +
                labs(y = NULL) + HIDE_Y_AXIS,
        make_microscopy_panel(depth_16S,  "16S",  "fungal_colonization_score", microscopy_ylim),
        make_microscopy_panel(depth_ITS1, "ITS1", "fungal_colonization_score", microscopy_ylim) +
                labs(y = NULL) + HIDE_Y_AXIS,
        ncol = 2
) +
        plot_layout(guides = "collect") +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "right")

print(fig5)
ggsave(file.path(RESULTS_DIR, "figures", "04_microscopy_depth.png"),
       plot = fig5, width = 12, height = 14, dpi = 300)

################################################################################
# Figure 6 + table: DNA yield and read depth, and bagged vs unbagged flowers
#   Row 1 — final read depth vs DNA concentration after PCR (Qubit HS, measured
#           per marker after bead cleaning and pooling the triplicates); all samples
#   Rows 2-3 — is amplicon DNA yield (row 2) or sequencing depth (row 3)
#   confounded with the treatment? Three checks per marker and variable:
#     - flowers: Mann-Whitney test (unpaired Wilcoxon), bagged vs unbagged
#     - trees:   mean of the bagged vs mean of the unbagged flowers of each tree
#                (3 + 3 by design), paired Wilcoxon test
#     - farms:   number of farms where the unbagged median is higher
#   Panels: one grey line per tree joins its bagged and unbagged means; boxplots
#   summarise the tree means. The y range of each row is shared by 16S and ITS1.
################################################################################

dna_depth_ylim <- c(0, shared_depth_ylim())

c16_dna <- make_scatter(depth_16S,  "dna_concentration_hs_16s",  "n_reads",
                        y_ylim   = dna_depth_ylim,
                        x_label  = "DNA concentration after PCR (Qubit HS)",
                        y_label  = "Read depth (final)",
                        y_labels = label_number(big.mark = ""))
cIT_dna <- make_scatter(depth_ITS1, "dna_concentration_hs_its1", "n_reads",
                        y_ylim   = dna_depth_ylim,
                        x_label  = "DNA concentration after PCR (Qubit HS)",
                        y_label  = "Read depth (final)",
                        y_labels = label_number(big.mark = ""))

TREATMENT_VARIABLES <- c(dna     = "DNA concentration after PCR (Qubit HS)",
                         n_reads = "Read depth (final)")

# One row per flower, marker and variable; the DNA column is marker-specific
treatment_flowers <- bind_rows(
        depth_16S  |> transmute(marker = "16S",  sample_type, farm_id, tree_id,
                                dna = dna_concentration_hs_16s, n_reads),
        depth_ITS1 |> transmute(marker = "ITS1", sample_type, farm_id, tree_id,
                                dna = dna_concentration_hs_its1, n_reads)
) |>
        filter(sample_type %in% flower_lvl) |>
        pivot_longer(c(dna, n_reads), names_to = "variable", values_to = "value") |>
        filter(!is.na(value))

# One row per tree: mean of its bagged and mean of its unbagged flowers
tree_means <- treatment_flowers |>
        group_by(marker, variable, farm_id, tree_id, sample_type) |>
        summarise(value = mean(value), .groups = "drop") |>
        pivot_wider(names_from = sample_type, values_from = value) |>
        filter(!is.na(bagged_flower), !is.na(unbagged_flower))

# exact = FALSE: normal approximation, as ties make the exact p-value unavailable
treatment_tests <- treatment_flowers |>
        group_by(marker, variable) |>
        summarise(
                n_bagged        = sum(sample_type == "bagged_flower"),
                n_unbagged      = sum(sample_type == "unbagged_flower"),
                median_bagged   = median(value[sample_type == "bagged_flower"]),
                median_unbagged = median(value[sample_type == "unbagged_flower"]),
                p_mann_whitney  = wilcox.test(value ~ sample_type, exact = FALSE)$p.value,
                .groups = "drop"
        ) |>
        left_join(
                tree_means |>
                        group_by(marker, variable) |>
                        summarise(
                                n_trees               = n(),
                                trees_unbagged_higher = sum(unbagged_flower > bagged_flower),
                                p_paired_wilcoxon     = wilcox.test(unbagged_flower, bagged_flower,
                                                                    paired = TRUE,
                                                                    exact  = FALSE)$p.value,
                                .groups = "drop"
                        ),
                by = c("marker", "variable")
        ) |>
        left_join(
                treatment_flowers |>
                        group_by(marker, variable, farm_id, sample_type) |>
                        summarise(value = median(value), .groups = "drop") |>
                        pivot_wider(names_from = sample_type, values_from = value) |>
                        group_by(marker, variable) |>
                        summarise(n_farms               = n(),
                                  farms_unbagged_higher = sum(unbagged_flower > bagged_flower),
                                  .groups = "drop"),
                by = c("marker", "variable")
        )

print(treatment_tests, width = Inf)
write_csv(treatment_tests,
          file.path(RESULTS_DIR, "tables", "04_treatment_dna_depth_tests.csv"))

make_paired <- function(marker_name, variable_name, y_ylim) {
        stats <- treatment_tests |> filter(marker == marker_name, variable == variable_name)

        tree_means |>
                filter(marker == marker_name, variable == variable_name) |>
                mutate(tree_id = as.character(tree_id)) |>
                pivot_longer(c(bagged_flower, unbagged_flower),
                             names_to = "sample_type", values_to = "tree_mean") |>
                ggplot(aes(x = sample_type, y = tree_mean)) +
                # One line per tree; colour/alpha set how strongly the lines show
                geom_line(aes(group = tree_id), colour = "grey60", alpha = 0.5) +
                # width: box width; alpha: box fill transparency (lines show through);
                # no legend: the x-axis already names the treatment
                geom_boxplot(aes(fill = sample_type), width = 0.4,
                             outlier.shape = NA, alpha = 0.7, show.legend = FALSE) +
                scale_fill_manual(values = SAMPLE_TYPE_PALETTE, name = "Sample type") +
                scale_x_discrete(labels = flower_lbl) +
                scale_y_continuous(labels = label_number(big.mark = "")) +  # 100000, not 1e+05
                coord_cartesian(ylim = y_ylim) +
                labs(x        = NULL,
                     y        = paste0(TREATMENT_VARIABLES[[variable_name]], "\n(tree mean)"),
                     subtitle = sprintf(paste0("Mann-Whitney p = %.1e (flowers)\n",
                                               "paired Wilcoxon p = %.1e (trees)\n",
                                               "unbagged higher: %d/%d trees, %d/%d farms"),
                                        stats$p_mann_whitney, stats$p_paired_wilcoxon,
                                        stats$trees_unbagged_higher, stats$n_trees,
                                        stats$farms_unbagged_higher, stats$n_farms)) +
                theme_bw() +
                theme(plot.subtitle = element_text(size = 8))   # size of the test summary
}

# Row y ceiling: the highest tree mean of either marker plus 5% headroom
row_ylim <- \(variable_name) tree_means |>
        filter(variable == variable_name) |>
        select(bagged_flower, unbagged_flower) |>
        max() |>
        (\(m) c(0, m * 1.05))()

# Panel tags A-F run left to right, top to bottom
fig6 <- wrap_plots(
        c16_dna + ggtitle("16S — Bacteria"),
        cIT_dna + ggtitle("ITS1 — Fungi") + labs(y = NULL) + HIDE_Y_AXIS,
        make_paired("16S",  "dna",     row_ylim("dna")),
        make_paired("ITS1", "dna",     row_ylim("dna")) + labs(y = NULL) + HIDE_Y_AXIS,
        make_paired("16S",  "n_reads", row_ylim("n_reads")),
        make_paired("ITS1", "n_reads", row_ylim("n_reads")) + labs(y = NULL) + HIDE_Y_AXIS,
        ncol = 2
) +
        plot_layout(guides = "collect") +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "right")

print(fig6)
ggsave(file.path(RESULTS_DIR, "figures", "04_treatment_dna_depth.png"),
       plot = fig6, width = 12, height = 15, dpi = 300)

################################################################################
# Figure 7 + table: batch checks of DNA yield and read depth (flowers only)
#   Within a farm, bagged and unbagged flowers were extracted (and amplified for
#   ITS1) in separate batches, so treatment and batch are aligned. Two parts of
#   the design break this alignment; both are found from the metadata below:
#   (a) the one extraction batch that holds both treatments (rows 1-2).
#       Does unbagged > bagged hold inside this batch? Boxes: the farm's flowers
#       from its own single-treatment batches; points: the mixed-batch flowers.
#   (b) farms whose bagged and unbagged flowers share one 16S PCR batch (row 3).
#       Is the 16S DNA gap similar in these farms and in the farms with
#       separate 16S PCR batches? (No ITS1 PCR batch holds both treatments.)
#   Tests are descriptive: few flowers in (a), few farms in (b).
################################################################################

flowers_16S <- depth_16S |> filter(sample_type %in% flower_lvl)

# Batches that hold both treatments of the same farm
mixed_extraction_batch <- flowers_16S |>
        mutate(extraction_batch = as.character(extraction_batch)) |>
        group_by(extraction_batch) |>
        filter(n_distinct(sample_type) == 2) |>
        pull(extraction_batch) |>
        unique()
shared_16S_pcr_farms <- flowers_16S |>
        group_by(farm_id, pcr_batch_16s) |>
        filter(n_distinct(sample_type) == 2) |>
        pull(farm_id) |>
        unique() |>
        as.character()

message("Mixed extraction batch: ", paste(mixed_extraction_batch, collapse = ", "),
        " | farms with a shared 16S PCR batch: ", paste(shared_16S_pcr_farms, collapse = ", "))

# One row per flower, marker and variable
batch_flowers <- bind_rows(
        depth_16S  |> transmute(marker = "16S",  sample_type, farm_id,
                                extraction_batch = as.character(extraction_batch),
                                dna = dna_concentration_hs_16s, n_reads),
        depth_ITS1 |> transmute(marker = "ITS1", sample_type, farm_id,
                                extraction_batch = as.character(extraction_batch),
                                dna = dna_concentration_hs_its1, n_reads)
) |>
        filter(sample_type %in% flower_lvl) |>
        mutate(mixed_batch = extraction_batch %in% mixed_extraction_batch) |>
        pivot_longer(c(dna, n_reads), names_to = "variable", values_to = "value") |>
        filter(!is.na(value))

mixed_batch_farms <- batch_flowers |> filter(mixed_batch) |> pull(farm_id) |> unique()

# Bagged vs unbagged summary of grouped flowers; exact = FALSE as above
summarise_treatment <- \(x) x |>
        summarise(
                n_bagged              = sum(sample_type == "bagged_flower"),
                n_unbagged            = sum(sample_type == "unbagged_flower"),
                median_bagged         = median(value[sample_type == "bagged_flower"]),
                median_unbagged       = median(value[sample_type == "unbagged_flower"]),
                ratio_unbagged_bagged = median_unbagged / median_bagged,
                p_mann_whitney        = wilcox.test(value ~ sample_type, exact = FALSE)$p.value,
                .groups = "drop"
        )

batch_tests <- bind_rows(
        # (a) mixed extraction batch: per farm and both farms together
        batch_flowers |>
                filter(mixed_batch) |>
                (\(x) bind_rows(x, mutate(x, farm_id = "both")))() |>
                group_by(marker, variable, farm_id) |>
                summarise_treatment() |>
                mutate(check = "mixed extraction batch"),
        # (b) 16S DNA: per farm, then per PCR group (farms pooled)
        batch_flowers |>
                filter(marker == "16S", variable == "dna") |>
                mutate(pcr_group = if_else(farm_id %in% shared_16S_pcr_farms,
                                           "shared 16S PCR batch",
                                           "separate 16S PCR batches")) |>
                (\(x) bind_rows(x, mutate(x, farm_id = paste("all:", pcr_group))))() |>
                group_by(marker, variable, pcr_group, farm_id) |>
                summarise_treatment() |>
                mutate(check = "16S PCR batch")
) |>
        relocate(check)

print(batch_tests, n = Inf, width = Inf)
write_csv(batch_tests, file.path(RESULTS_DIR, "tables", "04_batch_checks_dna_depth.csv"))

# (a) panel: boxes of the farm's own-batch flowers, points of the mixed-batch flowers
make_mixed_batch_panel <- function(marker_name, variable_name, y_ylim) {
        stats <- batch_tests |>
                filter(check == "mixed extraction batch", marker == marker_name,
                       variable == variable_name, farm_id == "both")

        batch_flowers |>
                filter(marker == marker_name, variable == variable_name,
                       farm_id %in% mixed_batch_farms) |>
                ggplot(aes(x = farm_id, y = value, fill = sample_type)) +
                # alpha: box fill transparency, so the mixed-batch points stand out
                geom_boxplot(data = \(d) filter(d, !mixed_batch),
                             position = position_dodge(width = 0.8), width = 0.7,
                             outlier.shape = NA, alpha = 0.4) +
                # size: mixed-batch point size; jitter.width spreads the 3 points of a group
                geom_point(data = \(d) filter(d, mixed_batch),
                           position = position_jitterdodge(jitter.width = 0.15,
                                                           dodge.width  = 0.8, seed = 1),
                           shape = 21, size = 3, colour = "black") +
                scale_fill_manual(values = SAMPLE_TYPE_PALETTE, limits = flower_lvl,
                                  name = "Sample type") +
                scale_y_continuous(labels = label_number(big.mark = "")) +
                coord_cartesian(ylim = y_ylim) +
                labs(x        = "Farm",
                     y        = TREATMENT_VARIABLES[[variable_name]],
                     subtitle = sprintf(paste0("Mixed batch %s, both farms: unbagged / bagged ",
                                               "median = %.2f,\nMann-Whitney p = %.2f ",
                                               "(n = %d + %d; descriptive)"),
                                        paste(mixed_extraction_batch, collapse = ", "),
                                        stats$ratio_unbagged_bagged, stats$p_mann_whitney,
                                        stats$n_bagged, stats$n_unbagged)) +
                theme_bw() +
                theme(plot.subtitle = element_text(size = 8))   # size of the test summary
}

# Row y ceiling for (a): the highest value of either marker in the mixed-batch farms
mixed_row_ylim <- \(variable_name) batch_flowers |>
        filter(variable == variable_name, farm_id %in% mixed_batch_farms) |>
        pull(value) |>
        max() |>
        (\(m) c(0, m * 1.05))()

# (b) panel: 16S DNA by farm and treatment, farms grouped by 16S PCR batch layout
pcr_group_stats <- batch_tests |> filter(check == "16S PCR batch", str_starts(farm_id, "all:"))

p_pcr_16S <- batch_flowers |>
        filter(marker == "16S", variable == "dna") |>
        mutate(farm_id   = factor(farm_id, levels = farm_levels),
               pcr_group = if_else(farm_id %in% shared_16S_pcr_farms,
                                   "shared 16S PCR batch", "separate 16S PCR batches") |>
                       factor(levels = c("shared 16S PCR batch", "separate 16S PCR batches"))) |>
        ggplot(aes(x = farm_id, y = value, fill = sample_type)) +
        # same alpha as panels A-D; no legend: the colours match the legend of A-D
        geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
                     alpha = 0.4, show.legend = FALSE) +
        scale_fill_manual(values = SAMPLE_TYPE_PALETTE, limits = flower_lvl,
                          name = "Sample type") +
        # strip at the bottom, below the farm labels (as in figure 3)
        facet_grid(. ~ pcr_group, scales = "free_x", space = "free_x", switch = "x") +
        labs(x        = "Farm",
             y        = TREATMENT_VARIABLES[["dna"]],
             title    = "16S — DNA by 16S PCR batch layout",
             subtitle = paste(sprintf("%s: unbagged / bagged median = %.2f",
                                      str_remove(pcr_group_stats$farm_id, "all: "),
                                      pcr_group_stats$ratio_unbagged_bagged),
                              collapse = "\n")) +
        theme_bw() +
        theme(strip.placement  = "outside",
              strip.background = element_blank(),
              plot.subtitle    = element_text(size = 8))

# Panel tags A-E run left to right, top to bottom
fig7 <- (
        (make_mixed_batch_panel("16S",  "dna", mixed_row_ylim("dna")) +
                 ggtitle("16S — Bacteria")) |
        (make_mixed_batch_panel("ITS1", "dna", mixed_row_ylim("dna")) +
                 ggtitle("ITS1 — Fungi") + labs(y = NULL) + HIDE_Y_AXIS)
) /
        (
        make_mixed_batch_panel("16S",  "n_reads", mixed_row_ylim("n_reads")) |
        (make_mixed_batch_panel("ITS1", "n_reads", mixed_row_ylim("n_reads")) +
                 labs(y = NULL) + HIDE_Y_AXIS)
) /
        p_pcr_16S +
        plot_layout(guides = "collect") +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "right")

print(fig7)
ggsave(file.path(RESULTS_DIR, "figures", "04_batch_checks_dna_depth.png"),
       plot = fig7, width = 12, height = 15, dpi = 300)

message("04 done")
