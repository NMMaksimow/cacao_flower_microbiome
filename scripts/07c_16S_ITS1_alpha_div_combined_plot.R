# ---- Libraries ----
library(here)
library(tidyverse)
library(ggh4x)
library(patchwork)

# ---- Parameters ----
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

# Treatment colours — used in panel 1 (bagged vs unbagged shown separately)
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

# ---- Load pre-computed data (produced by 07a / 07b) ----
alpha_16S_10k  <- readRDS(here("results", "rds", "alpha_16S_10k_999iters_summary.rds"))
alpha_ITS1_15k <- readRDS(here("results", "rds", "alpha_ITS1_otu97_15k_999iters_summary.rds"))

# Per-farm LMM results: cols farm, metric, estimate, se, t_value, p_raw, p_adj, stars
lmm_16S_10k  <- read_csv(here("results", "tables", "07a_10k_lmm_per_farm.csv"),
                          show_col_types = FALSE)
lmm_ITS1_15k <- read_csv(here("results", "tables", "07b_15k_lmm_per_farm.csv"),
                          show_col_types = FALSE)

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

################################################################################
#### Helpers ###################################################################
################################################################################

# One per-farm LMM boxplot (bagged vs unbagged) for a single metric.
# Significance stars are read from the pre-computed lmm_df CSV.
# No management nesting — farm IDs on x-axis only (same format as 07a/07b LMM figures).
# Parameters:
#   alpha_sum  — summary tibble (alpha_16S_10k or alpha_ITS1_15k)
#   lmm_df     — per-farm LMM results CSV (lmm_16S_10k or lmm_ITS1_15k)
#   metric_col — column name string, e.g. "Shannon_median"
#   y_label    — y-axis label string
plot_lmm_one <- function(alpha_sum, lmm_df, metric_col, y_label) {

    # Star y-position: 5% above per-farm max (matching 07a/07b star_pos logic)
    star_pos <- alpha_sum |>
        group_by(farm_id) |>
        summarise(star_y = max(.data[[metric_col]], na.rm = TRUE) * 1.05,
                  .groups = "drop")

    stars <- lmm_df |>
        filter(metric == metric_col, stars != "") |>
        left_join(star_pos, by = c("farm" = "farm_id")) |>
        mutate(farm_id = factor(farm, levels = FARM_LEVELS))

    p <- alpha_sum |>
        mutate(farm_id = factor(farm_id, levels = FARM_LEVELS)) |>
        ggplot(aes(x = farm_id, y = .data[[metric_col]], fill = sample_type)) +
        geom_boxplot(
            outlier.shape = NA,
            alpha         = 0.7,
            position      = position_dodge(width = 0.75)
        ) +
        geom_jitter(
            aes(colour = sample_type),
            position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
            alpha = 0.4, size = 0.9, stroke = 0
        ) +
        (if (nrow(stars) > 0)
            geom_text(data = stars, aes(x = farm_id, y = star_y, label = stars),
                      size = 4, inherit.aes = FALSE, vjust = 0)
        else NULL) +
        scale_fill_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
        scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 8)) +
        labs(x = "", y = y_label, fill = NULL, colour = NULL)
    p
}

# One per-farm pooled boxplot for a single metric.
# Both treatments combined — no fill split, single POOL_COL.
# Management type shown via facet_nested strip below x-axis
# (same theme as plot_alpha_farm_management() in 07a/07b).
# Parameters:
#   alpha_sum  — summary tibble (alpha_16S_10k or alpha_ITS1_15k)
#   metric_col — column name string
#   y_label    — y-axis label string
plot_pooled_one <- function(alpha_sum, metric_col, y_label) {
    alpha_sum |>
        mutate(farm_id = factor(farm_id, levels = FARM_LEVELS)) |>
        ggplot(aes(x = farm_id, y = .data[[metric_col]])) +
        geom_boxplot(
            outlier.shape = NA,
            alpha         = 0.65,
            fill          = POOL_COL,  # fixed colour — no treatment split
            colour        = "grey30",
            width         = 0.55
        ) +
        geom_jitter(
            width  = 0.15,
            alpha  = 0.35,
            size   = 0.9,
            stroke = 0,
            colour = POOL_COL
        ) +
        # space = "free_x": panels proportional to number of farms per management type
        ggh4x::facet_nested(
            ~ management_type,
            scales   = "free_x",
            space    = "free_x",
            switch   = "x",
            labeller = MGMT_LABELLER
        ) +
        theme_bw() +
        theme(
            strip.placement  = "outside",
            panel.border     = element_blank(),
            panel.spacing.x  = unit(0.2, "lines"),
            strip.background = element_blank(),
            strip.text.x     = element_text(size = 9, margin = margin(t = 5, b = 5)),
            axis.text.x      = element_text(size = 8)
        ) +
        labs(x = "", y = y_label)
}

################################################################################
#### Panel 1: per-farm LMM — 16S 10k (top) / ITS1 15k (bottom) ###############
#### Shows null bagged vs unbagged effect; justifies pooling in panel 2.      ##
#### Column 5 of ITS1 row is plot_spacer() so metric columns align vertically. ##
################################################################################

# Top row: 16S, 10k — all 5 metrics
p1_16s_obs <- plot_lmm_one(alpha_16S_10k, lmm_16S_10k, "Observed_median",     "Observed\nASVs")
p1_16s_sha <- plot_lmm_one(alpha_16S_10k, lmm_16S_10k, "Shannon_median",      "Shannon H")
p1_16s_inv <- plot_lmm_one(alpha_16S_10k, lmm_16S_10k, "InvSimpson_median",   "InvSimpson")
p1_16s_bp  <- plot_lmm_one(alpha_16S_10k, lmm_16S_10k, "BergerParker_median", "Berger-Parker")
p1_16s_fpd <- plot_lmm_one(alpha_16S_10k, lmm_16S_10k, "FaithPD_median",      "Faith's PD")

# Bottom row: ITS1, 15k — 4 metrics; plot_spacer() occupies column 5 (Faith's PD slot)
p1_its_obs <- plot_lmm_one(alpha_ITS1_15k, lmm_ITS1_15k, "Observed_median",     "Observed\nOTUs")
p1_its_sha <- plot_lmm_one(alpha_ITS1_15k, lmm_ITS1_15k, "Shannon_median",      "Shannon H")
p1_its_inv <- plot_lmm_one(alpha_ITS1_15k, lmm_ITS1_15k, "InvSimpson_median",   "InvSimpson")
p1_its_bp  <- plot_lmm_one(alpha_ITS1_15k, lmm_ITS1_15k, "BergerParker_median", "Berger-Parker")

p_panel1 <- (
    (p1_16s_obs | p1_16s_sha | p1_16s_inv | p1_16s_bp | p1_16s_fpd) /
    (p1_its_obs | p1_its_sha | p1_its_inv | p1_its_bp | plot_spacer())
) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

p_panel1 <- p_panel1 +
    plot_annotation(
        title    = "Alpha diversity by farm: bagged vs unbagged flowers",
        subtitle = paste0(
            "16S ASVs, 10k reads (top) / ITS1 OTU97, 15k reads (bottom) — ",
            "per-farm lmer(metric ~ sample_type + (1 | tree_id)); BH correction; * p_adj < 0.05"
        ),
        caption  = "Farms ordered left to right: inside forest to full sun (deforestation gradient)"
    )

ggsave(here("results", "figures", "07c_panel1_lmm_16S_ITS1.png"),
       plot = p_panel1, width = 20, height = 10, dpi = 300)

################################################################################
#### Panel 2: farm × management — bagged and unbagged pooled ##################
#### Null treatment effect (panel 1) justifies combining both groups.         ##
#### InvSimpson excluded — not a main-text metric.                            ##
#### Column 4 of ITS1 row is plot_spacer() (no Faith's PD for ITS1).         ##
################################################################################

# Top row: 16S, 10k — 4 main-text metrics (Observed, Shannon, Berger-Parker, Faith's PD)
p2_16s_obs <- plot_pooled_one(alpha_16S_10k, "Observed_median",     "Observed ASVs")
p2_16s_sha <- plot_pooled_one(alpha_16S_10k, "Shannon_median",      "Shannon H")
p2_16s_bp  <- plot_pooled_one(alpha_16S_10k, "BergerParker_median", "Berger-Parker")
p2_16s_fpd <- plot_pooled_one(alpha_16S_10k, "FaithPD_median",      "Faith's PD")

# Bottom row: ITS1, 15k — 3 metrics; plot_spacer() occupies column 4 (Faith's PD slot)
p2_its_obs <- plot_pooled_one(alpha_ITS1_15k, "Observed_median",     "Observed OTUs")
p2_its_sha <- plot_pooled_one(alpha_ITS1_15k, "Shannon_median",      "Shannon H")
p2_its_bp  <- plot_pooled_one(alpha_ITS1_15k, "BergerParker_median", "Berger-Parker")

p_panel2 <- (
    (p2_16s_obs | p2_16s_sha | p2_16s_bp | p2_16s_fpd) /
    (p2_its_obs | p2_its_sha | p2_its_bp | plot_spacer())
) +
    plot_annotation(
        title    = "Alpha diversity by farm and management type (bagged and unbagged flowers pooled)",
        subtitle = paste0(
            "16S ASVs, 10k reads (top) / ITS1 OTU97, 15k reads (bottom) — ",
            "treatments pooled (no significant insect-exclusion effect detected, see panel 1)"
        ),
        caption  = "Farms ordered left to right: inside forest to full sun; management type in strip labels"
    )

ggsave(here("results", "figures", "07c_panel2_pooled_16S_ITS1.png"),
       plot = p_panel2, width = 16, height = 10, dpi = 300)

message("Done. Figures saved to results/figures/")
