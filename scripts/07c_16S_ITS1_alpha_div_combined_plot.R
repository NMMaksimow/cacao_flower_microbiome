# Combined 16S | ITS1 alpha diversity figures.
# Inputs (from 07a / 07b): rarefaction summary RDS files, global treatment + farm table,
# per-farm treatment tables. Primary depth: 16S 2k, ITS1 5k; 10k / 15k = sensitivity.
# Figures:
#   07c_hill_profiles_16S_ITS1   Hill profiles, columns 16S | ITS1, rows low | high depth
#   07c_forest_treatment_16S_ITS1  treatment effect (unbagged - bagged) with 95% CI,
#                                  all farms (global model) and per farm, primary depth
#   07c_panel1_lmm_16S_ITS1      per-farm boxplots, bagged vs unbagged (primary depth)
#   07c_panel2_pooled_16S_ITS1   per-farm boxplots, treatments pooled (primary depth)

# ---- Libraries ----
library(here)
library(tidyverse)
library(ggh4x)
library(patchwork)

# Data version and sample set: RESULTS_DIR, BIOSAMPLE_SUFFIX, SET_TAG
source(here("scripts", "utils_paths.R"))

# ---- Parameters ----
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

# Treatment colours (same as 07a / 07b)
SAMPLE_TYPE_COLOURS <- c(bagged_flower = "#4C72B0", unbagged_flower = "#DD8452")
SAMPLE_TYPE_LABELS  <- c(bagged_flower = "bagged", unbagged_flower = "unbagged")

# Single colour for panel 2 (treatments pooled)
POOL_COL <- "#4A7AB5"

# Management type strip labeller — shared between plot_pooled_one calls
MGMT_LABELLER <- as_labeller(c(
    inside_forest = "inside\nforest",
    near_forest   = "near\nforest",
    agroforest    = "agroforest",
    full_sun      = "full\nsun"
))

# Main-text metrics (InvSimpson dropped; it appears only as Hill q = 2)
METRICS_16S  <- c("Observed_median", "Shannon_median", "BergerParker_median", "FaithPD_median")
METRICS_ITS1 <- c("Observed_median", "Shannon_median", "BergerParker_median")

# ---- Load pre-computed data (produced by 07a / 07b) ----
read_summary <- function(name) {
    readRDS(file.path(RESULTS_DIR, "rds", paste0(name, BIOSAMPLE_SUFFIX, ".rds")))
}
alpha_16S_2k   <- read_summary("alpha_16S_2k_999iters_summary")
alpha_16S_10k  <- read_summary("alpha_16S_10k_999iters_summary")
alpha_ITS1_5k  <- read_summary("alpha_ITS1_otu97_5k_999iters_summary")
alpha_ITS1_15k <- read_summary("alpha_ITS1_otu97_15k_999iters_summary")

read_table <- function(name, col_types = NULL) {
    read_csv(file.path(RESULTS_DIR, "tables", paste0(name, SET_TAG, ".csv")),
             col_types = col_types, show_col_types = FALSE)
}
# Global models: marker, depth, metric, estimate, conf_low, conf_high, p, q, farm_* (LRT)
# marker and depth read as text ("16S" would otherwise be guessed as a number)
GLOBAL_COL_TYPES <- cols(marker = col_character(), depth = col_character(), .default = col_guess())
effects_global <- bind_rows(read_table("07a_treatment_farm_effects", GLOBAL_COL_TYPES),
                            read_table("07b_treatment_farm_effects", GLOBAL_COL_TYPES))
# Per-farm models: farm, metric, estimate, conf_low, conf_high, p, q, stars
lmm_16S_2k   <- read_table("07a_2k_lmm_per_farm")
lmm_16S_10k  <- read_table("07a_10k_lmm_per_farm")
lmm_ITS1_5k  <- read_table("07b_5k_lmm_per_farm")
lmm_ITS1_15k <- read_table("07b_15k_lmm_per_farm")

dir.create(file.path(RESULTS_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

################################################################################
#### 1. Hill diversity profiles: 16S | ITS1, low depth (top) / high depth ######
#### Thin line = one flower; thick line + points = group median; ##############
#### ribbon = interquartile range. Log y axis: q = 0 (richness) is an order of #
#### magnitude above q >= 1, which a linear axis would flatten. ################
################################################################################

# x positions of the Hill orders; q = Inf drawn at 6 to keep a continuous axis
HILL_Q_POSITIONS <- c(`0` = 0, `0.5` = 0.5, `1` = 1, `1.5` = 1.5,
                      `2` = 2, `3` = 3, `5` = 5, `Inf` = 6)

# Long table of Hill numbers per flower (q = 2 is InvSimpson, q = Inf is 1 / Berger-Parker)
hill_long <- function(alpha_sum) {
    alpha_sum |>
        transmute(
            sample_id, sample_type,
            `0`   = Observed_median,
            `0.5` = Hill_q0_5_median,
            `1`   = Hill_q1_median,
            `1.5` = Hill_q1_5_median,
            `2`   = InvSimpson_median,
            `3`   = Hill_q3_median,
            `5`   = Hill_q5_median,
            `Inf` = 1 / BergerParker_median
        ) |>
        pivot_longer(-c(sample_id, sample_type), names_to = "q", values_to = "diversity") |>
        mutate(q_plot = HILL_Q_POSITIONS[q])
}

# One Hill panel. ylim: shared within a marker so the two depths compare directly.
# Only flowers that reached the depth are in alpha_sum, so n differs between rows;
# n per treatment is added to the panel title.
plot_hill <- function(alpha_sum, title, y_label, ylim) {
    long  <- hill_long(alpha_sum)
    n_trt <- table(alpha_sum$sample_type)
    title <- sprintf("%s (n = %d bagged, %d unbagged)", title,
                     n_trt[["bagged_flower"]], n_trt[["unbagged_flower"]])
    # median and quartiles: stay positive on the log axis (mean - SD does not)
    stats <- long |>
        group_by(sample_type, q_plot) |>
        summarise(med_div = median(diversity, na.rm = TRUE),
                  q25_div = quantile(diversity, 0.25, na.rm = TRUE),
                  q75_div = quantile(diversity, 0.75, na.rm = TRUE), .groups = "drop")
    ggplot() +
        # individual flowers: alpha / linewidth control how visible the traces are
        geom_line(data = long,
                  aes(x = q_plot, y = diversity, group = sample_id, colour = sample_type),
                  alpha = 0.10, linewidth = 0.3) +
        # interquartile range ribbon around the group median
        geom_ribbon(data = stats,
                    aes(x = q_plot, ymin = q25_div, ymax = q75_div, fill = sample_type),
                    alpha = 0.25, colour = NA) +
        # group median: straight segments between the orders, no smoothing
        geom_line(data = stats, aes(x = q_plot, y = med_div, colour = sample_type),
                  linewidth = 1.1) +
        geom_point(data = stats, aes(x = q_plot, y = med_div, colour = sample_type),
                   size = 2.5) +
        scale_x_continuous(breaks = HILL_Q_POSITIONS,
                           labels = c("0", "0.5", "1", "1.5", "2", "3", "5", "∞")) +
        scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
        scale_fill_manual(values   = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
        # log10 y axis; replace with scale_y_continuous() for a linear axis
        scale_y_log10() +
        coord_cartesian(ylim = ylim) +
        theme_bw() +
        labs(title = title, x = "Hill order (q)", y = y_label, colour = NULL, fill = NULL)
}

# Shared y range per marker (both depths)
ylim_hill_16S  <- range(c(hill_long(alpha_16S_2k)$diversity,  hill_long(alpha_16S_10k)$diversity),  na.rm = TRUE)
ylim_hill_ITS1 <- range(c(hill_long(alpha_ITS1_5k)$diversity, hill_long(alpha_ITS1_15k)$diversity), na.rm = TRUE)

p_hill <- (
    (plot_hill(alpha_16S_2k,   "16S, 2,000 reads",   "Effective number of ASVs",   ylim_hill_16S) |
     plot_hill(alpha_ITS1_5k,  "ITS1, 5,000 reads",  "Effective number of OTUs",   ylim_hill_ITS1)) /
    (plot_hill(alpha_16S_10k,  "16S, 10,000 reads",  "Effective number of ASVs",   ylim_hill_16S) |
     plot_hill(alpha_ITS1_15k, "ITS1, 15,000 reads", "Effective number of OTUs",   ylim_hill_ITS1))
) +
    plot_layout(guides = "collect") +
    plot_annotation(subtitle = paste0(
        "Thin lines: flowers that reached the depth (per flower, median of 999 rarefactions); ",
        "thick line: median across flowers; ribbon: interquartile range across flowers")) &
    theme(legend.position = "bottom")   # one shared legend under the grid

ggsave(file.path(RESULTS_DIR, "figures", paste0("07c_hill_profiles_16S_ITS1", SET_TAG, ".png")),
       plot = p_hill, width = 12, height = 9, dpi = 300)

################################################################################
#### 2. Forest plot: treatment effect (unbagged - bagged) with 95% CI ##########
#### Rows: "All farms" (global model, farm fixed + tree random) and the 7 farms #
#### (per-farm models, tree random). Primary depth only (16S 2k, ITS1 5k); ####
#### the 10k / 15k estimates are in the tables. ################################
#### Observed: NB GLMM -> ratio unbagged / bagged on a log axis (reference 1). ##
#### Other metrics: difference in metric units (reference 0). ##################
################################################################################

# Facet labels (column strips); the Observed label is also used to pick its log axis
METRIC_FOREST_LABELS <- c(
    Observed_median     = "Observed (ratio)",
    Shannon_median      = "Shannon H (difference)",
    BergerParker_median = "Berger-Parker (difference)",
    FaithPD_median      = "Faith's PD (difference)"
)
SIG_COLS  <- c(`TRUE` = "#B2182B", `FALSE` = "grey55")   # q < 0.05 vs not
Y_LEVELS  <- rev(c("All farms", FARM_LEVELS))            # "All farms" at the top

forest_df <- bind_rows(
    effects_global |> mutate(farm = "All farms"),
    lmm_16S_2k   |> mutate(marker = "16S",  depth = "2k"),
    lmm_ITS1_5k  |> mutate(marker = "ITS1", depth = "5k")
) |>
    # primary depth only
    filter(depth %in% c("2k", "5k")) |>
    select(marker, depth, farm, metric, estimate, conf_low, conf_high, q) |>
    mutate(
        # NB GLMM effect is on the log scale: exp() gives the ratio unbagged / bagged
        across(c(estimate, conf_low, conf_high),
               \(x) if_else(metric == "Observed_median", exp(x), x)),
        significant = factor(!is.na(q) & q < 0.05, levels = c(TRUE, FALSE)),
        # numeric y: row position ("All farms" at the top)
        y      = match(farm, Y_LEVELS),
        metric = factor(metric, levels = names(METRIC_FOREST_LABELS), labels = METRIC_FOREST_LABELS),
        marker = factor(marker, levels = c("16S", "ITS1"))
    )

# Reference line per panel: 1 for the ratio, 0 for the differences.
# Only marker x metric panels with data (no ITS1 x Faith's PD), so that panel stays blank.
forest_ref <- forest_df |>
    distinct(marker, metric) |>
    mutate(xint = if_else(metric == "Observed (ratio)", 1, 0))

p_forest <- forest_df |>
    ggplot(aes(y = y)) +
    geom_vline(data = forest_ref, aes(xintercept = xint), linetype = "dashed", colour = "grey40") +
    # line separating the global estimate from the per-farm estimates
    geom_hline(data = forest_ref, aes(yintercept = length(Y_LEVELS) - 0.5),
               colour = "grey75", linewidth = 0.4) +
    geom_linerange(aes(xmin = conf_low, xmax = conf_high, colour = significant), linewidth = 0.6) +
    geom_point(aes(x = estimate, colour = significant), size = 2.4) +
    # rows = marker, columns = metric; free x so each metric has its own range
    # render_empty = FALSE leaves the ITS1 x Faith's PD panel blank (no ITS1 phylogeny)
    ggh4x::facet_grid2(marker ~ metric, scales = "free_x", independent = "x",
                       render_empty = FALSE) +
    # Observed column on a log axis (ratio); others linear
    ggh4x::facetted_pos_scales(x = list(metric == "Observed (ratio)" ~ scale_x_log10())) +
    scale_y_continuous(breaks = seq_along(Y_LEVELS), labels = Y_LEVELS,
                       expand = expansion(add = 0.5)) +
    scale_colour_manual(values = SIG_COLS, labels = c(`TRUE` = "q < 0.05", `FALSE` = "q ≥ 0.05"),
                        drop = FALSE) +
    theme_bw() +
    theme(
        legend.position  = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
    ) +
    labs(
        title    = "Effect of insect access on alpha diversity (unbagged vs bagged flowers)",
        subtitle = paste0("16S 2,000 reads, ITS1 5,000 reads. Estimate and 95% CI. ",
                          "All farms: metric ~ treatment + farm + (1 | tree); ",
                          "per farm: metric ~ treatment + (1 | tree).\n",
                          "Observed: NB GLMM, ratio unbagged / bagged; others: LMM, difference. ",
                          "BH within marker."),
        x = "Unbagged vs bagged", y = NULL, colour = NULL
    )

ggsave(file.path(RESULTS_DIR, "figures", paste0("07c_forest_treatment_16S_ITS1", SET_TAG, ".png")),
       plot = p_forest, width = 13, height = 6.5, dpi = 300)

################################################################################
#### Helpers for the boxplot panels ##########################################
################################################################################

# One per-farm boxplot (bagged vs unbagged) for a single metric.
# Significance stars (q < 0.05) are read from the per-farm table.
plot_lmm_one <- function(alpha_sum, lmm_df, metric_col, y_label) {

    # Star y-position: 5% above per-farm max
    star_pos <- alpha_sum |>
        group_by(farm_id) |>
        summarise(star_y = max(.data[[metric_col]], na.rm = TRUE) * 1.05,
                  .groups = "drop")

    stars <- lmm_df |>
        filter(metric == metric_col, stars != "") |>
        left_join(star_pos, by = c("farm" = "farm_id")) |>
        mutate(farm_id = factor(farm, levels = FARM_LEVELS))

    alpha_sum |>
        mutate(farm_id = factor(farm_id, levels = FARM_LEVELS)) |>
        ggplot(aes(x = farm_id, y = .data[[metric_col]], fill = sample_type)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +
        geom_jitter(aes(colour = sample_type),
                    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
                    alpha = 0.4, size = 0.9, stroke = 0) +
        (if (nrow(stars) > 0)
            geom_text(data = stars, aes(x = farm_id, y = star_y, label = stars),
                      size = 4, inherit.aes = FALSE, vjust = 0)
        else NULL) +
        scale_fill_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
        scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 8)) +
        labs(x = "", y = y_label, fill = NULL, colour = NULL)
}

# One per-farm pooled boxplot (both treatments) with management type strips below.
# subtitle: farm LRT q from the global model of this metric.
plot_pooled_one <- function(alpha_sum, metric_col, y_label, subtitle) {
    alpha_sum |>
        mutate(farm_id = factor(farm_id, levels = FARM_LEVELS)) |>
        ggplot(aes(x = farm_id, y = .data[[metric_col]])) +
        geom_boxplot(outlier.shape = NA, alpha = 0.65, fill = POOL_COL,
                     colour = "grey30", width = 0.55) +
        geom_jitter(width = 0.15, alpha = 0.35, size = 0.9, stroke = 0, colour = POOL_COL) +
        # space = "free_x": panels proportional to number of farms per management type
        ggh4x::facet_nested(~ management_type, scales = "free_x", space = "free_x",
                            switch = "x", labeller = MGMT_LABELLER) +
        theme_bw() +
        theme(
            strip.placement  = "outside",
            panel.border     = element_blank(),
            panel.spacing.x  = unit(0.2, "lines"),
            strip.background = element_blank(),
            strip.text.x     = element_text(size = 9, margin = margin(t = 5, b = 5)),
            axis.text.x      = element_text(size = 8),
            plot.subtitle    = element_text(size = 9)
        ) +
        labs(x = "", y = y_label, subtitle = subtitle)
}

# Farm LRT q (global model) as a subtitle string
farm_q_label <- function(marker_name, depth_name, metric_col) {
    effects_global |>
        filter(marker == marker_name, depth == depth_name, metric == metric_col) |>
        pull(farm_q) |>
        (\(q) sprintf("farm LRT q = %.3f", q))()
}

################################################################################
#### 3. Panel 1: per-farm boxplots, bagged vs unbagged — primary depth ########
#### Column 4 of the ITS1 row is plot_spacer() (no Faith's PD for ITS1). ######
################################################################################

p_panel1 <- (
    (plot_lmm_one(alpha_16S_2k, lmm_16S_2k, "Observed_median",     "Observed\nASVs") |
     plot_lmm_one(alpha_16S_2k, lmm_16S_2k, "Shannon_median",      "Shannon H") |
     plot_lmm_one(alpha_16S_2k, lmm_16S_2k, "BergerParker_median", "Berger-Parker") |
     plot_lmm_one(alpha_16S_2k, lmm_16S_2k, "FaithPD_median",      "Faith's PD")) /
    (plot_lmm_one(alpha_ITS1_5k, lmm_ITS1_5k, "Observed_median",     "Observed\nOTUs") |
     plot_lmm_one(alpha_ITS1_5k, lmm_ITS1_5k, "Shannon_median",      "Shannon H") |
     plot_lmm_one(alpha_ITS1_5k, lmm_ITS1_5k, "BergerParker_median", "Berger-Parker") |
     plot_spacer())
) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

p_panel1 <- p_panel1 +
    plot_annotation(
        title    = "Alpha diversity by farm: bagged vs unbagged flowers",
        subtitle = paste0(
            "16S ASVs, 2,000 reads (top) / ITS1 OTU97, 5,000 reads (bottom) — ",
            "per-farm metric ~ treatment + (1 | tree) (Observed NB GLMM, others LMM); BH; * q < 0.05"
        ),
        caption  = "Farms ordered left to right: inside forest to full sun (deforestation gradient)"
    )

ggsave(file.path(RESULTS_DIR, "figures", paste0("07c_panel1_lmm_16S_ITS1", SET_TAG, ".png")),
       plot = p_panel1, width = 16, height = 10, dpi = 300)

################################################################################
#### 4. Panel 2: per-farm boxplots, treatments pooled — primary depth #########
#### Pooled for display of farm differences only; the treatment effect is in ##
#### the forest plot. Farm LRT q (global model) above each panel. #############
################################################################################

p_panel2 <- (
    (plot_pooled_one(alpha_16S_2k, "Observed_median",     "Observed ASVs", farm_q_label("16S", "2k", "Observed_median")) |
     plot_pooled_one(alpha_16S_2k, "Shannon_median",      "Shannon H",     farm_q_label("16S", "2k", "Shannon_median")) |
     plot_pooled_one(alpha_16S_2k, "BergerParker_median", "Berger-Parker", farm_q_label("16S", "2k", "BergerParker_median")) |
     plot_pooled_one(alpha_16S_2k, "FaithPD_median",      "Faith's PD",    farm_q_label("16S", "2k", "FaithPD_median"))) /
    (plot_pooled_one(alpha_ITS1_5k, "Observed_median",     "Observed OTUs", farm_q_label("ITS1", "5k", "Observed_median")) |
     plot_pooled_one(alpha_ITS1_5k, "Shannon_median",      "Shannon H",     farm_q_label("ITS1", "5k", "Shannon_median")) |
     plot_pooled_one(alpha_ITS1_5k, "BergerParker_median", "Berger-Parker", farm_q_label("ITS1", "5k", "BergerParker_median")) |
     plot_spacer())
) +
    plot_annotation(
        title    = "Alpha diversity by farm and management type (bagged and unbagged flowers pooled for display)",
        subtitle = paste0(
            "16S ASVs, 2,000 reads (top) / ITS1 OTU97, 5,000 reads (bottom) — farm effect: ",
            "likelihood ratio test in metric ~ treatment + farm + (1 | tree), BH within marker"
        ),
        caption  = "Farms ordered left to right: inside forest to full sun; management type in strip labels"
    )

ggsave(file.path(RESULTS_DIR, "figures", paste0("07c_panel2_pooled_16S_ITS1", SET_TAG, ".png")),
       plot = p_panel2, width = 16, height = 10, dpi = 300)

message("Done. Figures saved to ", file.path(RESULTS_DIR, "figures"))
