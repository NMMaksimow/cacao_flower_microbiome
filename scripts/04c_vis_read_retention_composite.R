# ---- Packages ----
library(here)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(scales)
library(ggtext)

# ---- Load objects ----
ps_16S_raw      <- readRDS(here("results", "rds", "ps_16S_raw.rds"))
ps_16S_decontam <- readRDS(here("results", "rds", "ps_16S_decontam.rds"))
ps_16S_bacteria <- readRDS(here("results", "rds", "ps_16S_bacteria.rds"))

ps_ITS1_raw      <- readRDS(here("results", "rds", "ps_ITS1_raw.rds"))
ps_ITS1_decontam <- readRDS(here("results", "rds", "ps_ITS1_decontam.rds"))
ps_ITS1_fungi    <- readRDS(here("results", "rds", "ps_ITS1_fungi.rds"))

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

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
# The 16S (left) panel already carries the scale; ITS1 shares it via coord_cartesian
# or scale_y_continuous(limits = c(0,1)), so the ticks would be redundant.
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
#   hline        – numeric: draw horizontal red dashed reference line
#   show_n       – TRUE: print "n=X" above each box
#   markdown_x   – TRUE: render x-axis labels as markdown/HTML (for coloured labels)
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
                     markdown_x          = FALSE,
                     fill_by_sample_type = FALSE) {
                     # fill_by_sample_type = TRUE: boxes colored by sample_type using
                     # SAMPLE_TYPE_PALETTE; useful when grouping_variable is sublibrary_id
                     # so each sublibrary shows a distinct box per sample type present.

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
                labs(x = if (is.null(x_label)) grouping_variable else x_label,
                     y = "Read depth (final)") +
                theme_bw() +
                theme(axis.text.x = if (markdown_x)
                        element_markdown(angle = x_angle, hjust = hjust_val)
                      else
                        element_text(angle = x_angle, hjust = hjust_val))

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

# Build colored x-axis labels for sublibrary plots.
# Each sublibrary label is colored by the dominant sample type it contains.
make_colored_sublibrary_labels <- function(depth_tbl, palette) {
        dominant <- depth_tbl |>
                group_by(sublibrary_id, sample_type) |>
                summarise(n = n(), .groups = "drop") |>
                group_by(sublibrary_id) |>
                slice_max(n, n = 1, with_ties = FALSE) |>
                ungroup()
        col <- palette[dominant$sample_type]
        col[is.na(col)] <- "black"
        setNames(
                sprintf("<span style='color:%s'>%s</span>", col, dominant$sublibrary_id),
                dominant$sublibrary_id
        )
}

# Standard 2×2 composite: stacked bars (top) / boxplots (bottom),
# 16S (left) / ITS1 (right).
# Stacked bar x-axis is hidden; y-label on left panels only.
make_composite <- function(stack_16S, stack_ITS1, box_16S, box_ITS1,
                           filename, width = 14, height = 10) {
        fig <- (
                (stack_16S + ggtitle("16S \u2014 Bacteria")) |
                (stack_ITS1 + ggtitle("ITS1 \u2014 Fungi") + labs(y = NULL) + HIDE_Y_AXIS)
        ) /
        (
                box_16S |
                (box_ITS1 + labs(y = NULL) + HIDE_Y_AXIS)
        ) +
        plot_layout(guides = "collect") &
        theme(legend.position = "right")

        print(fig)
        ggsave(here("results", "figures", filename),
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

# Shared y ceiling: both markers on the same absolute read-count scale
subl_ylim <- max(depth_16S$n_reads, depth_ITS1$n_reads, na.rm = TRUE) * 1.10

fig1 <- (
        # Top row: stacked bars — fraction of reads at each pipeline stage
        (make_stack(tbl_16S,  "sublibrary_id", hide_x_axis = TRUE) +
                 ggtitle("16S \u2014 Bacteria")) |
        (make_stack(tbl_ITS1, "sublibrary_id", hide_x_axis = TRUE) +
                 ggtitle("ITS1 \u2014 Fungi") + labs(y = NULL) + HIDE_Y_AXIS)
) /
(
        # Bottom row: boxplots colored by sample_type; each sublibrary may
        # contain multiple sample types shown as separate colored boxes
        (make_box(depth_16S, "sublibrary_id", x_angle = 45, x_label = "Sublibrary",
                  hline = 10000, fill_by_sample_type = TRUE) +
                 coord_cartesian(ylim = c(0, subl_ylim))) |
        (make_box(depth_ITS1, "sublibrary_id", x_angle = 45, x_label = "Sublibrary",
                  hline = 10000, fill_by_sample_type = TRUE) +
                 coord_cartesian(ylim = c(0, subl_ylim)) + labs(y = NULL) + HIDE_Y_AXIS)
) +
plot_layout(guides = "collect") &
theme(legend.position = "right")

print(fig1)
ggsave(here("results", "figures", "04c_composite_sublibrary.png"),
       plot = fig1, width = 16, height = 10, dpi = 300)

################################################################################
# Figure 2: By sample type
#   Stack  — mean retention profile per sample type category
#   Boxplot — read depth distribution; shows whether controls vs. flowers differ
################################################################################

sample_type_levels <- c("extraction_control", "negative_pcr_control",
                        "mock_community", "bagged_flower", "unbagged_flower")

make_composite(
        stack_16S  = make_stack(tbl_16S,  "sample_type",
                                grouping_levels = sample_type_levels,
                                x_angle = 0, hide_x_axis = TRUE),
        stack_ITS1 = make_stack(tbl_ITS1, "sample_type",
                                grouping_levels = sample_type_levels,
                                x_angle = 0, hide_x_axis = TRUE),
        box_16S    = make_box(depth_16S,  "sample_type",
                              grouping_levels = sample_type_levels,
                              x_angle = 0, x_label = "Sample type",
                              fill_by_sample_type = TRUE),
        box_ITS1   = make_box(depth_ITS1, "sample_type",
                              grouping_levels = sample_type_levels,
                              x_angle = 0, x_label = "Sample type",
                              fill_by_sample_type = TRUE),
        filename   = "04c_composite_sample_type.png",
        width = 14, height = 10
)

################################################################################
# Figure 3: Farm × management type
#   Stack  — retention profile per farm; faceted by management type (strip at bottom)
#   Boxplot — read depth per farm; checks for sequencing depth confounds
#             with management type or geography
################################################################################

farm_levels <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
mgmt_levels <- c("inside_forest", "near_forest", "agroforest", "full_sun")
bio_fn      <- \(x) filter(x, sample_type %in% c("bagged_flower", "unbagged_flower"))

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
                              grouping_levels = farm_levels,
                              facet_variable  = "management_type",
                              facet_levels    = mgmt_levels,
                              filter_fn       = bio_fn,
                              x_angle = 0, x_label = "Farm"),
        box_ITS1   = make_box(depth_ITS1, "farm_id",
                              grouping_levels = farm_levels,
                              facet_variable  = "management_type",
                              facet_levels    = mgmt_levels,
                              filter_fn       = bio_fn,
                              x_angle = 0, x_label = "Farm"),
        filename   = "04c_composite_farm_management.png",
        width = 18, height = 10
)

################################################################################
# Figure 4: Phenotype variables
#   Columns (L→R): Flower treatment | Pollination (CLSM) | Germination | Fungal hyphae
#   Each column: 16S left panel, ITS1 right panel (y-axis ticks hidden)
#   Top row  — stacked bars: retention profile per phenotype group; x hidden
#   Bottom row — boxplots: read depth by group; x-title shown on both markers
#              so "Yes / No" labels are interpretable without knowing the variable;
#              n= per group shown above each box
#
# NOTE: replace "is_hyphae_present" with the actual column name in your metadata.
################################################################################

yn_labels   <- c("FALSE" = "No", "TRUE" = "Yes")
bool_levels <- c("FALSE", "TRUE")

flower_fn_ph <- \(x) filter(x, sample_type %in% c("bagged_flower", "unbagged_flower"))
poll_fn  <- \(x) filter(x, !is.na(is_pollination_clsm)) |>
        mutate(is_pollination_clsm = as.character(is_pollination_clsm))
germ_fn  <- \(x) filter(x, !is.na(is_germination)) |>
        mutate(is_germination = as.character(is_germination))
hyph_fn  <- \(x) filter(x, !is.na(is_hyphae_present)) |>   # FIXME column name
        mutate(is_hyphae_present = as.character(is_hyphae_present))

flower_lvl <- c("bagged_flower", "unbagged_flower")
flower_lbl <- c("bagged_flower" = "Bagged", "unbagged_flower" = "Unbagged")

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

s16_hyph <- make_stack(tbl_16S,  "is_hyphae_present", grouping_levels = bool_levels,   # FIXME
                       filter_fn = hyph_fn, x_angle = 0,
                       x_scale_labels = yn_labels, hide_x_axis = TRUE)
sIT_hyph <- make_stack(tbl_ITS1, "is_hyphae_present", grouping_levels = bool_levels,   # FIXME
                       filter_fn = hyph_fn, x_angle = 0,
                       x_scale_labels = yn_labels, hide_x_axis = TRUE)

# ── Boxplots (x-title on 16S panel only per variable; n= per group) ───────────
b16_flow <- make_box(depth_16S,  "sample_type", grouping_levels = flower_lvl,
                     filter_fn = flower_fn_ph, x_angle = 0, show_n = TRUE,
                     x_scale_labels = flower_lbl, x_label = "Flower treatment")
bIT_flow <- make_box(depth_ITS1, "sample_type", grouping_levels = flower_lvl,
                     filter_fn = flower_fn_ph, x_angle = 0, show_n = TRUE,
                     x_scale_labels = flower_lbl, x_label = "Flower treatment")

b16_poll <- make_box(depth_16S,  "is_pollination_clsm", grouping_levels = bool_levels,
                     filter_fn = poll_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Pollination (CLSM)")
bIT_poll <- make_box(depth_ITS1, "is_pollination_clsm", grouping_levels = bool_levels,
                     filter_fn = poll_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Pollination (CLSM)")

b16_germ <- make_box(depth_16S,  "is_germination", grouping_levels = bool_levels,
                     filter_fn = germ_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Pollen germination")
bIT_germ <- make_box(depth_ITS1, "is_germination", grouping_levels = bool_levels,
                     filter_fn = germ_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Pollen germination")

b16_hyph <- make_box(depth_16S,  "is_hyphae_present", grouping_levels = bool_levels,
                     filter_fn = hyph_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Fungal hyphae present")
bIT_hyph <- make_box(depth_ITS1, "is_hyphae_present", grouping_levels = bool_levels,
                     filter_fn = hyph_fn, x_angle = 0, show_n = TRUE,
                     x_scale_labels = yn_labels, x_label = "Fungal hyphae present")

# ── Assemble: all 16S left (cols 1–4), all ITS1 right (cols 5–8) ─────────────
# One wrap_plots per row (not nested blocks) so patchwork allocates equal
# pixel width to every column regardless of y-axis presence on col 1.
# y-label on col 1 (s16_flow / b16_flow) only; all other panels hide y-axis.
# x-title kept on all boxplot panels so "Yes / No" labels are interpretable.

top_row <- wrap_plots(
        s16_flow + ggtitle("16S \u2014 Bacteria"),          # col 1: y-label shown
        s16_poll + labs(y = NULL) + HIDE_Y_AXIS,
        s16_germ + labs(y = NULL) + HIDE_Y_AXIS,
        s16_hyph + labs(y = NULL) + HIDE_Y_AXIS,
        sIT_flow + ggtitle("ITS1 \u2014 Fungi") + labs(y = NULL) + HIDE_Y_AXIS,
        sIT_poll + labs(y = NULL) + HIDE_Y_AXIS,
        sIT_germ + labs(y = NULL) + HIDE_Y_AXIS,
        sIT_hyph + labs(y = NULL) + HIDE_Y_AXIS,
        nrow = 1
)

bot_row <- wrap_plots(
        b16_flow,                                            # col 1: y-label shown
        b16_poll + labs(y = NULL) + HIDE_Y_AXIS,
        b16_germ + labs(y = NULL) + HIDE_Y_AXIS,
        b16_hyph + labs(y = NULL) + HIDE_Y_AXIS,
        bIT_flow + labs(y = NULL) + HIDE_Y_AXIS,
        bIT_poll + labs(y = NULL) + HIDE_Y_AXIS,
        bIT_germ + labs(y = NULL) + HIDE_Y_AXIS,
        bIT_hyph + labs(y = NULL) + HIDE_Y_AXIS,
        nrow = 1
)

fig4 <- (top_row / bot_row) +
        plot_layout(guides = "collect") &
        theme(legend.position = "right")

print(fig4)
ggsave(here("results", "figures", "04c_composite_phenotypes.png"),
       plot = fig4, width = 28, height = 10, dpi = 300)

message("04c done")
