# ---- Libraries ----
library(here)
library(tidyverse)
library(lme4)
library(lmerTest)    # p-values for lmer()
library(glmmTMB)
library(DHARMa)
library(performance)
library(ggeffects)
library(ggh4x)
library(patchwork)
library(cowplot)

# ---- Parameters ----
# N_ITER and RAREFY_DEPTH* must match values used in 07b_ITS1_alpha_div_rarefaction.R
N_ITER           <- 999
RAREFY_DEPTH     <- 5000
RAREFY_DEPTH_15k <- 15000

# Farm order: inside_forest → full_sun (deforestation gradient)
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
MGMT_LEVELS <- c("inside_forest", "near_forest", "agroforest", "full_sun")

# Colour palette for bagged/unbagged (adjust hex if needed)
SAMPLE_TYPE_COLOURS <- c(bagged_flower = "#4C72B0", unbagged_flower = "#DD8452")
SAMPLE_TYPE_LABELS  <- c(bagged_flower = "bagged", unbagged_flower = "unbagged")

# ---- Input: HPC run (no Chao1) ----
alpha_iters_5k    <- readRDS(here("results", "rds", "alpha_ITS1_otu97_5k_999iters.rds"))
alpha_summary_5k  <- readRDS(here("results", "rds", "alpha_ITS1_otu97_5k_999iters_summary.rds"))
alpha_iters_15k   <- readRDS(here("results", "rds", "alpha_ITS1_otu97_15k_999iters.rds"))
alpha_summary_15k <- readRDS(here("results", "rds", "alpha_ITS1_otu97_15k_999iters_summary.rds"))

# ---- Local-run RDS (contain Chao1; loaded for reference, not used in plots or models) ----
alpha_iters_5k_local    <- readRDS(here("results", "rds", "alpha_ITS1_otu97_5k_999iters_local_run.rds"))
alpha_summary_5k_local  <- readRDS(here("results", "rds", "alpha_ITS1_otu97_5k_999iters_summary_local_run.rds"))
alpha_iters_15k_local   <- readRDS(here("results", "rds", "alpha_ITS1_otu97_15k_999iters_local_run.rds"))
alpha_summary_15k_local <- readRDS(here("results", "rds", "alpha_ITS1_otu97_15k_999iters_summary_local_run.rds"))

################################################################################
#### Shared helpers ############################################################
################################################################################

# Farm × management nested boxplot for one metric.
# metric_col: column name from alpha_sum (e.g. "Shannon_median")
# y_label:    y-axis label string
# alpha_sum:  summary tibble — either alpha_summary_5k or alpha_summary_15k
plot_alpha_farm_management <- function(metric_col, y_label, alpha_sum, ylim = NULL) {

    p <- alpha_sum |>
        ggplot(aes(x = farm_id, y = .data[[metric_col]], fill = sample_type)) +
        geom_boxplot(
            outlier.shape = NA,
            alpha         = 0.7,
            position      = position_dodge(width = 0.75)
        ) +
        geom_jitter(
            position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
            aes(colour = sample_type),
            alpha    = 0.4,
            size     = 0.8,
            stroke   = 0.2
        ) +
        ggh4x::facet_nested(
            ~ management_type,
            scales   = "free_x",
            space    = "free_x",
            switch   = "x",
            labeller = as_labeller(c(
                inside_forest = "inside\nforest",  # 1 farm — farm-level variance unidentifiable
                near_forest   = "near\nforest",
                agroforest    = "agroforest",
                full_sun      = "full\nsun"
            ))
        ) +
        scale_fill_manual(values   = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
        scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
        theme_bw() +
        theme(
            strip.placement  = "outside",
            panel.border     = element_blank(),
            panel.spacing.x  = unit(0.2, "lines"),
            strip.background = element_blank(),
            strip.text.x     = element_text(size = 9, margin = margin(t = 5, b = 5))
        ) +
        labs(x = "", y = y_label, fill = NULL)
    if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
    p
}

# Marginal means plot from ggpredict (one model × one term)
plot_marginal <- function(model, term, y_label) {
    ggeffects::ggpredict(model, terms = term) |>
        ggplot(aes(x = x, y = predicted)) +
        geom_point(size = 3) +
        geom_errorbar(
            aes(ymin = conf.low, ymax = conf.high),
            width = 0.25
        ) +
        theme_bw() +
        labs(x = "", y = y_label)
}

# x-axis positions for Hill profile
# Inf plotted at position 6 to avoid breaking the continuous scale
HILL_Q_POSITIONS <- c(`0` = 0, `0.5` = 0.5, `1` = 1, `1.5` = 1.5,
                      `2` = 2, `3` = 3, `5` = 5, `Inf` = 6)

################################################################################
#### 5. Convergence diagnostic — 5k ############################################
#### Per-iteration group median + cumulative running mean ######################
################################################################################

conv_long_5k <- alpha_iters_5k |>
    pivot_longer(
        cols      = c(Observed, Shannon, InvSimpson, BergerParker),
        names_to  = "metric",
        values_to = "value"
    ) |>
    # canonical left-to-right order: least sensitive → most sensitive to rarefaction depth
    mutate(metric = factor(metric,
        levels = c("Observed", "Shannon", "InvSimpson", "BergerParker"),
        labels = c("Observed", "Shannon", "InvSimpson", "Berger-Parker"))) |>
    group_by(iteration, sample_type, metric) |>
    summarise(iter_median = median(value, na.rm = TRUE), .groups = "drop") |>
    group_by(sample_type, metric) |>
    arrange(iteration) |>
    mutate(running_mean = cummean(iter_median)) |>
    ungroup()

# Top row: raw per-iteration values (shows natural noise level)
p_conv_raw_5k <- conv_long_5k |>
    ggplot(aes(x = iteration, y = iter_median, colour = sample_type)) +
    geom_line(alpha = 0.8, linewidth = 0.5) +
    facet_wrap(~ metric, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    labs(
        title  = "Per-iteration group median (5k)",
        x      = NULL,
        y      = "Median across samples",
        colour = NULL
    )

# Bottom row: cumulative running mean — flattens well before iteration 100
p_conv_cum_5k <- conv_long_5k |>
    ggplot(aes(x = iteration, y = running_mean, colour = sample_type)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    labs(
        title  = "Cumulative running mean (convergence diagnostic)",
        x      = "Iteration",
        y      = "Cumulative mean",
        colour = NULL
    )

p_conv_5k <- p_conv_raw_5k / p_conv_cum_5k +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

################################################################################
#### 5. Convergence diagnostic — 15k ###########################################
################################################################################

conv_long_15k <- alpha_iters_15k |>
    pivot_longer(
        cols      = c(Observed, Shannon, InvSimpson, BergerParker),
        names_to  = "metric",
        values_to = "value"
    ) |>
    mutate(metric = factor(metric,
        levels = c("Observed", "Shannon", "InvSimpson", "BergerParker"),
        labels = c("Observed", "Shannon", "InvSimpson", "Berger-Parker"))) |>
    group_by(iteration, sample_type, metric) |>
    summarise(iter_median = median(value, na.rm = TRUE), .groups = "drop") |>
    group_by(sample_type, metric) |>
    arrange(iteration) |>
    mutate(running_mean = cummean(iter_median)) |>
    ungroup()

p_conv_raw_15k <- conv_long_15k |>
    ggplot(aes(x = iteration, y = iter_median, colour = sample_type)) +
    geom_line(alpha = 0.8, linewidth = 0.5) +
    facet_wrap(~ metric, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    labs(
        title  = "Per-iteration group median (15k)",
        x      = NULL,
        y      = "Median across samples",
        colour = NULL
    )

p_conv_cum_15k <- conv_long_15k |>
    ggplot(aes(x = iteration, y = running_mean, colour = sample_type)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    labs(
        title  = "Cumulative running mean (convergence diagnostic)",
        x      = "Iteration",
        y      = "Cumulative mean",
        colour = NULL
    )

p_conv_15k <- p_conv_raw_15k / p_conv_cum_15k +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

################################################################################
#### 6a. Boxplots: bagged vs unbagged — 5k ####################################
################################################################################

p_sample_type_5k <- alpha_summary_5k |>
    pivot_longer(
        cols      = c(Observed_median, Shannon_median,
                      InvSimpson_median, BergerParker_median),
        names_to  = "metric",
        values_to = "value"
    ) |>
    mutate(metric = factor(
        metric,
        levels = c("Observed_median", "Shannon_median",
                   "InvSimpson_median", "BergerParker_median"),
        labels = c("Observed", "Shannon H", "InvSimpson", "Berger-Parker")
    )) |>
    ggplot(aes(x = sample_type, y = value, fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
    # ncol = 4: Observed, Shannon H, InvSimpson, Berger-Parker
    facet_wrap(~ metric, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = SAMPLE_TYPE_COLOURS) +
    scale_x_discrete(labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    theme(
        legend.position  = "none",
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold")
    ) +
    labs(
        title    = "Fungal alpha diversity (ITS1 OTU97) by flower treatment — 5k rarefaction",
        subtitle = sprintf("Median of %d rarefactions at %s reads per sample",
                           N_ITER, scales::comma(RAREFY_DEPTH)),
        x = "", y = ""
    )

################################################################################
#### 6a. Boxplots: bagged vs unbagged — 15k ###################################
################################################################################

p_sample_type_15k <- alpha_summary_15k |>
    pivot_longer(
        cols      = c(Observed_median, Shannon_median,
                      InvSimpson_median, BergerParker_median),
        names_to  = "metric",
        values_to = "value"
    ) |>
    mutate(metric = factor(
        metric,
        levels = c("Observed_median", "Shannon_median",
                   "InvSimpson_median", "BergerParker_median"),
        labels = c("Observed", "Shannon H", "InvSimpson", "Berger-Parker")
    )) |>
    ggplot(aes(x = sample_type, y = value, fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
    facet_wrap(~ metric, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = SAMPLE_TYPE_COLOURS) +
    scale_x_discrete(labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    theme(
        legend.position  = "none",
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold")
    ) +
    labs(
        title    = "Fungal alpha diversity (ITS1 OTU97) by flower treatment — 15k rarefaction",
        subtitle = sprintf("Median of %d rarefactions at %s reads (note: fewer samples retained)",
                           N_ITER, scales::comma(RAREFY_DEPTH_15k)),
        x = "", y = ""
    )

################################################################################
#### 6b. Farm × management nested boxplot — 5k #################################
################################################################################

p_farm_mgmt_5k <- (
    plot_alpha_farm_management("Observed_median",     "Observed",      alpha_summary_5k) |
    plot_alpha_farm_management("Shannon_median",      "Shannon H",     alpha_summary_5k) |
    plot_alpha_farm_management("InvSimpson_median",   "InvSimpson",    alpha_summary_5k) |
    plot_alpha_farm_management("BergerParker_median", "Berger-Parker", alpha_summary_5k)
) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

################################################################################
#### 6b. Farm × management nested boxplot — 15k ################################
################################################################################

p_farm_mgmt_15k <- (
    plot_alpha_farm_management("Observed_median",     "Observed",      alpha_summary_15k) |
    plot_alpha_farm_management("Shannon_median",      "Shannon H",     alpha_summary_15k) |
    plot_alpha_farm_management("InvSimpson_median",   "InvSimpson",    alpha_summary_15k) |
    plot_alpha_farm_management("BergerParker_median", "Berger-Parker", alpha_summary_15k)
) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

################################################################################
#### 6c. Hill diversity profile (q = 0 to ∞) — 5k #############################
#### Each curve = one sample; thick line = group mean; ribbon = ±1 SD ##########
################################################################################

hill_long_5k <- alpha_summary_5k |>
    mutate(Hill_qInf_median = 1 / BergerParker_median) |>
    select(
        sample_id, sample_type,
        `0`   = Observed_median,
        `0.5` = Hill_q0_5_median,
        `1`   = Hill_q1_median,
        `1.5` = Hill_q1_5_median,
        `2`   = InvSimpson_median,
        `3`   = Hill_q3_median,
        `5`   = Hill_q5_median,
        `Inf` = Hill_qInf_median
    ) |>
    pivot_longer(
        cols      = -c(sample_id, sample_type),
        names_to  = "q",
        values_to = "diversity"
    ) |>
    mutate(
        q      = as.numeric(q),
        # Replace Inf with plot position 6 so continuous x-axis works
        q_plot = HILL_Q_POSITIONS[as.character(q)]
    )

hill_group_stats_5k <- hill_long_5k |>
    group_by(sample_type, q, q_plot) |>
    summarise(
        mean_div = mean(diversity, na.rm = TRUE),
        sd_div   = sd(diversity, na.rm = TRUE),
        .groups  = "drop"
    )

p_hill_5k <- ggplot() +
    # individual sample traces (one connected line per sample across all q values)
    geom_line(
        data      = hill_long_5k,
        aes(x = q_plot, y = diversity, group = sample_id, colour = sample_type),
        alpha     = 0.10,
        linewidth = 0.3
    ) +
    # ±1 SD ribbon around group mean
    geom_ribbon(
        data  = hill_group_stats_5k,
        aes(x    = q_plot,
            ymin = mean_div - sd_div,
            ymax = mean_div + sd_div,
            fill = sample_type),
        alpha  = 0.25,
        colour = NA
    ) +
    # straight segments connecting group mean points — honest representation, no smoothing
    geom_line(
        data      = hill_group_stats_5k,
        aes(x = q_plot, y = mean_div, colour = sample_type),
        linewidth = 1.2
    ) +
    # group mean dots — one per q value per group
    geom_point(
        data = hill_group_stats_5k,
        aes(x = q_plot, y = mean_div, colour = sample_type),
        size = 3.5
    ) +
    scale_x_continuous(
        # 8 breaks at actual plot positions; label position 6 as ∞
        breaks = c(0, 0.5, 1, 1.5, 2, 3, 5, 6),
        labels = c("0", "0.5", "1", "1.5", "2", "3", "5", "∞")
    ) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    scale_fill_manual(values   = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(
        title    = "Hill diversity profile (q = 0 to ∞) — 5k rarefaction",
        subtitle = "Thick line = group mean; ribbon = ±1 SD; thin lines = individual samples",
        x        = "Hill order (q)",
        y        = "Effective number of OTUs",
        colour   = NULL,
        fill     = NULL
    )

################################################################################
#### 6c. Hill diversity profile (q = 0 to ∞) — 15k ############################
################################################################################

hill_long_15k <- alpha_summary_15k |>
    mutate(Hill_qInf_median = 1 / BergerParker_median) |>
    select(
        sample_id, sample_type,
        `0`   = Observed_median,
        `0.5` = Hill_q0_5_median,
        `1`   = Hill_q1_median,
        `1.5` = Hill_q1_5_median,
        `2`   = InvSimpson_median,
        `3`   = Hill_q3_median,
        `5`   = Hill_q5_median,
        `Inf` = Hill_qInf_median
    ) |>
    pivot_longer(
        cols      = -c(sample_id, sample_type),
        names_to  = "q",
        values_to = "diversity"
    ) |>
    mutate(
        q      = as.numeric(q),
        q_plot = HILL_Q_POSITIONS[as.character(q)]
    )

hill_group_stats_15k <- hill_long_15k |>
    group_by(sample_type, q, q_plot) |>
    summarise(
        mean_div = mean(diversity, na.rm = TRUE),
        sd_div   = sd(diversity, na.rm = TRUE),
        .groups  = "drop"
    )

p_hill_15k <- ggplot() +
    geom_line(
        data      = hill_long_15k,
        aes(x = q_plot, y = diversity, group = sample_id, colour = sample_type),
        alpha     = 0.10,
        linewidth = 0.3
    ) +
    geom_ribbon(
        data  = hill_group_stats_15k,
        aes(x    = q_plot,
            ymin = mean_div - sd_div,
            ymax = mean_div + sd_div,
            fill = sample_type),
        alpha  = 0.25,
        colour = NA
    ) +
    geom_line(
        data      = hill_group_stats_15k,
        aes(x = q_plot, y = mean_div, colour = sample_type),
        linewidth = 1.2
    ) +
    geom_point(
        data = hill_group_stats_15k,
        aes(x = q_plot, y = mean_div, colour = sample_type),
        size = 3.5
    ) +
    scale_x_continuous(
        breaks = c(0, 0.5, 1, 1.5, 2, 3, 5, 6),
        labels = c("0", "0.5", "1", "1.5", "2", "3", "5", "∞")
    ) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    scale_fill_manual(values   = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(
        title    = "Hill diversity profile (q = 0 to ∞) — 15k rarefaction",
        subtitle = "Thick line = group mean; ribbon = ±1 SD; thin lines = individual samples",
        x        = "Hill order (q)",
        y        = "Effective number of OTUs",
        colour   = NULL,
        fill     = NULL
    )

################################################################################
#### 7a. Statistical models: sample_type (main treatment effect) — 5k ##########
################################################################################

# Full nested random structure (1 | farm_id / tree_id).
# If farm_id variance = 0 (singular fit), fall back to (1 | tree_id).
lmm_shannon_treatment_5k <- lmer(
    Shannon_median ~ sample_type + (1 | farm_id / tree_id),
    data = alpha_summary_5k
)
summary(lmm_shannon_treatment_5k)
performance::check_singularity(lmm_shannon_treatment_5k)
performance::check_normality(lmm_shannon_treatment_5k)

# Simplified model used for marginal means if full model is singular
lmm_shannon_treatment_tree_5k <- lmer(
    Shannon_median ~ sample_type + (1 | tree_id),
    data = alpha_summary_5k
)
summary(lmm_shannon_treatment_tree_5k)

lmm_invsimpson_treatment_5k <- lmer(
    InvSimpson_median ~ sample_type + (1 | farm_id / tree_id),
    data = alpha_summary_5k
)
summary(lmm_invsimpson_treatment_5k)
performance::check_singularity(lmm_invsimpson_treatment_5k)

lmm_bergerparker_treatment_5k <- lmer(
    BergerParker_median ~ sample_type + (1 | farm_id / tree_id),
    data = alpha_summary_5k
)
summary(lmm_bergerparker_treatment_5k)
performance::check_singularity(lmm_bergerparker_treatment_5k)

# Observed OTUs: negative binomial GLMM (overdispersed integer counts)
glmm_observed_treatment_5k <- glmmTMB::glmmTMB(
    round(Observed_median) ~ sample_type + (1 | tree_id),
    family = nbinom2,
    data   = alpha_summary_5k
)
summary(glmm_observed_treatment_5k)

glmm_observed_treatment_5k |>
    DHARMa::simulateResiduals(n = 1000) |>
    (\(sim) {
        plot(sim)
        DHARMa::testUniformity(sim)
        DHARMa::testDispersion(sim)
    })()

################################################################################
#### 7a. Statistical models: sample_type (main treatment effect) — 15k #########
################################################################################

lmm_shannon_treatment_15k <- lmer(
    Shannon_median ~ sample_type + (1 | farm_id / tree_id),
    data = alpha_summary_15k
)
summary(lmm_shannon_treatment_15k)
performance::check_singularity(lmm_shannon_treatment_15k)
performance::check_normality(lmm_shannon_treatment_15k)

lmm_shannon_treatment_tree_15k <- lmer(
    Shannon_median ~ sample_type + (1 | tree_id),
    data = alpha_summary_15k
)
summary(lmm_shannon_treatment_tree_15k)

lmm_invsimpson_treatment_15k <- lmer(
    InvSimpson_median ~ sample_type + (1 | farm_id / tree_id),
    data = alpha_summary_15k
)
summary(lmm_invsimpson_treatment_15k)
performance::check_singularity(lmm_invsimpson_treatment_15k)

lmm_bergerparker_treatment_15k <- lmer(
    BergerParker_median ~ sample_type + (1 | farm_id / tree_id),
    data = alpha_summary_15k
)
summary(lmm_bergerparker_treatment_15k)
performance::check_singularity(lmm_bergerparker_treatment_15k)

glmm_observed_treatment_15k <- glmmTMB::glmmTMB(
    round(Observed_median) ~ sample_type + (1 | tree_id),
    family = nbinom2,
    data   = alpha_summary_15k
)
summary(glmm_observed_treatment_15k)

glmm_observed_treatment_15k |>
    DHARMa::simulateResiduals(n = 1000) |>
    (\(sim) {
        plot(sim)
        DHARMa::testUniformity(sim)
        DHARMa::testDispersion(sim)
    })()

################################################################################
#### 7c. Statistical models: CLSM subsets — 5k #################################
################################################################################

# pi_clsm 3-class: none (0), weak (1-39), strong (>= 40)
# Threshold >= 40 is the nearest integer to 43 pollen grains/stigma,
# the level corresponding to 50% of maximum fruit set (Falque et al. 1995).
# Binary is_pollination_clsm was removed: any detectable pollen conflates
# accidental/wind/self-deposition with insect-mediated pollination.
# NB: strong class n ~15 after 5k rarefaction — CIs are wide.
alpha_summary_5k |>
    filter(sample_type == "unbagged_flower", !is.na(pi_clsm)) |>
    mutate(pi_class = factor(
        case_when(
            pi_clsm == 0                ~ "none",
            pi_clsm > 0 & pi_clsm < 40 ~ "weak",
            pi_clsm >= 40               ~ "strong"
        ),
        levels = c("none", "weak", "strong")
    )) |>
    lmer(Shannon_median ~ pi_class + (1 | farm_id / tree_id), data = _) |>
    summary()

alpha_summary_5k |>
    filter(sample_type == "unbagged_flower", !is.na(pi_clsm)) |>
    mutate(pi_class = factor(
        case_when(
            pi_clsm == 0                ~ "none",
            pi_clsm > 0 & pi_clsm < 40 ~ "weak",
            pi_clsm >= 40               ~ "strong"
        ),
        levels = c("none", "weak", "strong")
    )) |>
    lmer(BergerParker_median ~ pi_class + (1 | farm_id / tree_id), data = _) |>
    summary()

# is_hyphae_present: all flowers with CLSM data; sample_type as covariate
alpha_summary_5k |>
    filter(!is.na(is_hyphae_present)) |>
    lmer(Shannon_median ~ is_hyphae_present + sample_type + (1 | farm_id / tree_id),
         data = _) |>
    summary()

alpha_summary_5k |>
    filter(!is.na(is_hyphae_present)) |>
    lmer(InvSimpson_median ~ is_hyphae_present + sample_type + (1 | farm_id / tree_id),
         data = _) |>
    summary()

alpha_summary_5k |>
    filter(!is.na(is_hyphae_present)) |>
    lmer(BergerParker_median ~ is_hyphae_present + sample_type + (1 | farm_id / tree_id),
         data = _) |>
    summary()

################################################################################
#### 7c. Statistical models: CLSM subsets — 15k ################################
################################################################################

# pi_clsm 3-class: none (0), weak (1-39), strong (>= 40)
# Threshold >= 40 is the nearest integer to 43 pollen grains/stigma,
# the level corresponding to 50% of maximum fruit set (Falque et al. 1995).
# Binary is_pollination_clsm was removed: any detectable pollen conflates
# accidental/wind/self-deposition with insect-mediated pollination.
# NB: strong class n ~15 after 15k rarefaction — CIs are wide.
alpha_summary_15k |>
    filter(sample_type == "unbagged_flower", !is.na(pi_clsm)) |>
    mutate(pi_class = factor(
        case_when(
            pi_clsm == 0                ~ "none",
            pi_clsm > 0 & pi_clsm < 40 ~ "weak",
            pi_clsm >= 40               ~ "strong"
        ),
        levels = c("none", "weak", "strong")
    )) |>
    lmer(Shannon_median ~ pi_class + (1 | farm_id / tree_id), data = _) |>
    summary()

alpha_summary_15k |>
    filter(sample_type == "unbagged_flower", !is.na(pi_clsm)) |>
    mutate(pi_class = factor(
        case_when(
            pi_clsm == 0                ~ "none",
            pi_clsm > 0 & pi_clsm < 40 ~ "weak",
            pi_clsm >= 40               ~ "strong"
        ),
        levels = c("none", "weak", "strong")
    )) |>
    lmer(BergerParker_median ~ pi_class + (1 | farm_id / tree_id), data = _) |>
    summary()

alpha_summary_15k |>
    filter(!is.na(is_hyphae_present)) |>
    lmer(Shannon_median ~ is_hyphae_present + sample_type + (1 | farm_id / tree_id),
         data = _) |>
    summary()

alpha_summary_15k |>
    filter(!is.na(is_hyphae_present)) |>
    lmer(InvSimpson_median ~ is_hyphae_present + sample_type + (1 | farm_id / tree_id),
         data = _) |>
    summary()

alpha_summary_15k |>
    filter(!is.na(is_hyphae_present)) |>
    lmer(BergerParker_median ~ is_hyphae_present + sample_type + (1 | farm_id / tree_id),
         data = _) |>
    summary()

################################################################################
#### 8a. Marginal means: sample_type (bagged vs unbagged) — 5k #################
#### Model-predicted means ± 95% CI, marginalised over tree and farm effects ####
################################################################################

p_obs_treat_5k <- plot_marginal(glmm_observed_treatment_5k,    "sample_type", "Observed OTUs")
p_sha_treat_5k <- plot_marginal(lmm_shannon_treatment_tree_5k, "sample_type", "Shannon H")
p_inv_treat_5k <- plot_marginal(lmm_invsimpson_treatment_5k,   "sample_type", "InvSimpson")
p_bp_treat_5k  <- plot_marginal(lmm_bergerparker_treatment_5k, "sample_type", "Berger-Parker")

p_treat_5k <- (p_obs_treat_5k | p_sha_treat_5k | p_inv_treat_5k | p_bp_treat_5k) +
    plot_annotation(
        title    = "Marginal means ± 95% CI: sample_type (bagged vs unbagged) — 5k rarefaction",
        subtitle = "Model-predicted means marginalised over tree and farm random effects"
    )

################################################################################
#### 8a. Marginal means: sample_type — 15k #####################################
################################################################################

p_obs_treat_15k <- plot_marginal(glmm_observed_treatment_15k,    "sample_type", "Observed OTUs")
p_sha_treat_15k <- plot_marginal(lmm_shannon_treatment_tree_15k, "sample_type", "Shannon H")
p_inv_treat_15k <- plot_marginal(lmm_invsimpson_treatment_15k,   "sample_type", "InvSimpson")
p_bp_treat_15k  <- plot_marginal(lmm_bergerparker_treatment_15k, "sample_type", "Berger-Parker")

p_treat_15k <- (p_obs_treat_15k | p_sha_treat_15k | p_inv_treat_15k | p_bp_treat_15k) +
    plot_annotation(
        title    = "Marginal means ± 95% CI: sample_type (bagged vs unbagged) — 15k rarefaction",
        subtitle = "Model-predicted means marginalised over tree and farm random effects"
    )

################################################################################
#### 9_icc. Within-tree ICC — are flowers from the same tree more similar? #####
#### Intercept-only LMM with (1 | tree_id) across all farms combined. ##########
#### ICC_adjusted = var(tree) / (var(tree) + var(residual)) ####################
#### Interpretation: < 0.1 negligible, 0.1–0.3 modest, > 0.3 meaningful. ######
#### High ICC → tree-level aggregation is justified; low ICC → aggregation ######
#### loses information without protecting against pseudoreplication. ############
################################################################################

METRICS_ICC_ITS1 <- c("Observed_median", "Shannon_median", "InvSimpson_median",
                      "BergerParker_median")

# Run at 5k depth (primary depth); ICC is a structural property of the sampling design
icc_df_5k <- expand.grid(
    metric      = METRICS_ICC_ITS1,
    sample_type = c("bagged_flower", "unbagged_flower"),
    stringsAsFactors = FALSE
) |>
    as_tibble() |>
    mutate(icc_adjusted = map2_dbl(metric, sample_type, \(m, st) {
        fit <- alpha_summary_5k |>
            filter(sample_type == st) |>
            lmer(as.formula(paste0(m, " ~ 1 + (1 | tree_id)")), data = _)
        performance::icc(fit)$ICC_adjusted
    }))

print(icc_df_5k)
write_csv(icc_df_5k, here("results", "tables", "07b_5k_icc_within_tree.csv"))

################################################################################
#### 9b. Per-farm LMM: bagged vs unbagged at flower level (1 | tree_id) ########
#### Replaces paired Wilcoxon on tree medians (removed based on ICC results). ##
#### ICC (section 9_icc) showed moderate within-tree clustering for ITS1;     ##
#### LMM models this directly via random intercept without aggregating.        ##
################################################################################

# 4 metrics for ITS1 (no Faith's PD — ITS too variable for phylogenetic inference)
METRICS_LMM_FARM <- c("Observed_median", "Shannon_median", "InvSimpson_median",
                      "BergerParker_median")
METRIC_LABELS_LMM <- c(
    Observed_median     = "Observed",
    Shannon_median      = "Shannon H",
    InvSimpson_median   = "InvSimpson",
    BergerParker_median = "Berger-Parker"
)

# 7 farms × 4 metrics = 28 tests for BH correction
n_tests_lmm_5k <- length(FARM_LEVELS) * length(METRICS_LMM_FARM)

# Per-farm LMM: sample_type fixed effect, tree_id random intercept.
# Reference level: bagged_flower; positive estimate = higher in unbagged_flower.
# lmerTest loaded at top → summary()$coefficients includes Pr(>|t|) (Satterthwaite).
# tryCatch handles degenerate farms (singular fit, missing data) → NA rows.
lmm_farm_df_5k <- expand.grid(
    farm   = FARM_LEVELS,
    metric = METRICS_LMM_FARM,
    stringsAsFactors = FALSE
) |>
    as_tibble() |>
    mutate(fit_result = map2(farm, metric, \(f, m) {
        dat <- alpha_summary_5k |>
            filter(farm_id == f) |>
            mutate(sample_type = factor(sample_type,
                                        levels = c("bagged_flower", "unbagged_flower")))
        fit <- tryCatch(
            lmer(as.formula(paste0(m, " ~ sample_type + (1 | tree_id)")), data = dat),
            error = \(e) NULL
        )
        if (is.null(fit))
            return(tibble(estimate = NA_real_, se = NA_real_, t_value = NA_real_, p_raw = NA_real_))
        coef_tbl <- summary(fit)$coefficients
        row_name <- "sample_typeunbagged_flower"
        if (!row_name %in% rownames(coef_tbl))
            return(tibble(estimate = NA_real_, se = NA_real_, t_value = NA_real_, p_raw = NA_real_))
        tibble(
            estimate = coef_tbl[row_name, "Estimate"],
            se       = coef_tbl[row_name, "Std. Error"],
            t_value  = coef_tbl[row_name, "t value"],
            p_raw    = coef_tbl[row_name, "Pr(>|t|)"]
        )
    })) |>
    unnest(fit_result) |>
    mutate(
        p_adj = p.adjust(p_raw, method = "BH"),
        stars = case_when(
            p_adj < 0.001 ~ "***",
            p_adj < 0.01  ~ "**",
            p_adj < 0.05  ~ "*",
            TRUE          ~ ""
        )
    )

write_csv(lmm_farm_df_5k, here("results", "tables", "07b_5k_lmm_per_farm.csv"))

# Star y-position: 5% above max flower-level value per farm × metric
star_pos_lmm_5k <- alpha_summary_5k |>
    pivot_longer(cols = all_of(METRICS_LMM_FARM), names_to = "metric", values_to = "value") |>
    group_by(farm_id, metric) |>
    summarise(star_y = max(value, na.rm = TRUE) * 1.05, .groups = "drop")

lmm_stars_5k <- lmm_farm_df_5k |>
    filter(stars != "") |>
    left_join(star_pos_lmm_5k, by = c("farm" = "farm_id", "metric"))

# Boxplot of raw flower-level data with LMM significance stars
p_lmm_farm_5k <- alpha_summary_5k |>
    mutate(farm_id = factor(farm_id, levels = FARM_LEVELS)) |>
    pivot_longer(cols = all_of(METRICS_LMM_FARM), names_to = "metric", values_to = "value") |>
    mutate(metric = factor(metric, levels = METRICS_LMM_FARM, labels = METRIC_LABELS_LMM)) |>
    ggplot(aes(x = farm_id, y = value, fill = sample_type)) +
    geom_boxplot(
        outlier.shape = NA,           # individual flowers shown via jitter below
        alpha         = 0.7,
        position      = position_dodge(width = 0.75)
    ) +
    geom_jitter(
        aes(colour = sample_type),
        position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
        alpha    = 0.4,
        size     = 0.9,
        stroke   = 0
    ) +
    # significance stars above each farm pair where p_adj < 0.05 after BH
    # guard: empty lmm_stars_5k skips layer to avoid NULL PANEL in ggh4x::FacetWrap2
    (if (nrow(lmm_stars_5k) > 0)
        geom_text(
            data = lmm_stars_5k |>
                mutate(
                    farm_id = factor(farm, levels = FARM_LEVELS),
                    metric  = factor(metric, levels = METRICS_LMM_FARM,
                                     labels = METRIC_LABELS_LMM)
                ),
            aes(x = farm_id, y = star_y, label = stars),
            size        = 4,
            inherit.aes = FALSE,
            vjust       = 0
        )
    else NULL) +
    # ncol = 4: Observed, Shannon H, InvSimpson, Berger-Parker
    # facet_wrap2 (ggh4x) required for facetted_pos_scales compatibility
    facet_wrap2(~ metric, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    scale_x_discrete() +
    theme_bw() +
    theme(
        legend.position  = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold"),
        axis.text.x      = element_text(size = 8)
    ) +
    labs(
        title    = "Fungal alpha diversity — per-farm LMM (bagged vs unbagged) — 5k rarefaction",
        subtitle = sprintf(
            "Flower-level data; lmer(metric ~ sample_type + (1|tree_id)) per farm; BH correction (%d tests); * p_adj < 0.05",
            n_tests_lmm_5k
        ),
        x = "", y = "", fill = NULL, colour = NULL
    )

# --- 15k depth: same LMM, secondary rarefaction depth ---

n_tests_lmm_15k <- length(FARM_LEVELS) * length(METRICS_LMM_FARM)

lmm_farm_df_15k <- expand.grid(
    farm   = FARM_LEVELS,
    metric = METRICS_LMM_FARM,
    stringsAsFactors = FALSE
) |>
    as_tibble() |>
    mutate(fit_result = map2(farm, metric, \(f, m) {
        dat <- alpha_summary_15k |>
            filter(farm_id == f) |>
            mutate(sample_type = factor(sample_type,
                                        levels = c("bagged_flower", "unbagged_flower")))
        fit <- tryCatch(
            lmer(as.formula(paste0(m, " ~ sample_type + (1 | tree_id)")), data = dat),
            error = \(e) NULL
        )
        if (is.null(fit))
            return(tibble(estimate = NA_real_, se = NA_real_, t_value = NA_real_, p_raw = NA_real_))
        coef_tbl <- summary(fit)$coefficients
        row_name <- "sample_typeunbagged_flower"
        if (!row_name %in% rownames(coef_tbl))
            return(tibble(estimate = NA_real_, se = NA_real_, t_value = NA_real_, p_raw = NA_real_))
        tibble(
            estimate = coef_tbl[row_name, "Estimate"],
            se       = coef_tbl[row_name, "Std. Error"],
            t_value  = coef_tbl[row_name, "t value"],
            p_raw    = coef_tbl[row_name, "Pr(>|t|)"]
        )
    })) |>
    unnest(fit_result) |>
    mutate(
        p_adj = p.adjust(p_raw, method = "BH"),
        stars = case_when(
            p_adj < 0.001 ~ "***",
            p_adj < 0.01  ~ "**",
            p_adj < 0.05  ~ "*",
            TRUE          ~ ""
        )
    )

write_csv(lmm_farm_df_15k, here("results", "tables", "07b_15k_lmm_per_farm.csv"))

star_pos_lmm_15k <- alpha_summary_15k |>
    pivot_longer(cols = all_of(METRICS_LMM_FARM), names_to = "metric", values_to = "value") |>
    group_by(farm_id, metric) |>
    summarise(star_y = max(value, na.rm = TRUE) * 1.05, .groups = "drop")

lmm_stars_15k <- lmm_farm_df_15k |>
    filter(stars != "") |>
    left_join(star_pos_lmm_15k, by = c("farm" = "farm_id", "metric"))

p_lmm_farm_15k <- alpha_summary_15k |>
    mutate(farm_id = factor(farm_id, levels = FARM_LEVELS)) |>
    pivot_longer(cols = all_of(METRICS_LMM_FARM), names_to = "metric", values_to = "value") |>
    mutate(metric = factor(metric, levels = METRICS_LMM_FARM, labels = METRIC_LABELS_LMM)) |>
    ggplot(aes(x = farm_id, y = value, fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7,
                 position = position_dodge(width = 0.75)) +
    geom_jitter(aes(colour = sample_type),
                position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
                alpha = 0.4, size = 0.9, stroke = 0) +
    # guard: empty lmm_stars_15k skips layer to avoid NULL PANEL in ggh4x::FacetWrap2
    (if (nrow(lmm_stars_15k) > 0)
        geom_text(
            data = lmm_stars_15k |>
                mutate(farm_id = factor(farm, levels = FARM_LEVELS),
                       metric  = factor(metric, levels = METRICS_LMM_FARM,
                                        labels = METRIC_LABELS_LMM)),
            aes(x = farm_id, y = star_y, label = stars),
            size = 4, inherit.aes = FALSE, vjust = 0
        )
    else NULL) +
    # facet_wrap2 (ggh4x) required for facetted_pos_scales compatibility
    facet_wrap2(~ metric, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    scale_x_discrete() +
    theme_bw() +
    theme(legend.position  = "bottom",
          strip.background = element_rect(fill = "grey90"),
          strip.text       = element_text(face = "bold"),
          axis.text.x      = element_text(size = 8)) +
    labs(
        title    = "Fungal alpha diversity — per-farm LMM (bagged vs unbagged) — 15k rarefaction",
        subtitle = sprintf(
            "Flower-level data; lmer(metric ~ sample_type + (1|tree_id)) per farm; BH correction (%d tests); * p_adj < 0.05",
            n_tests_lmm_15k
        ),
        x = "", y = "", fill = NULL, colour = NULL
    )

################################################################################
#### 10. Combined depth-comparison figures — ITS1 #############################
#### Left / top  = 5k (lower depth, more samples retained) ####################
#### Right / bottom = 15k (higher depth, fewer samples retained) ###############
################################################################################

################################################################################
#### Shared Y-axis limits for depth-comparison combined figures ################
#### Compute per-metric global range from both depths; apply before ggsave. ####
################################################################################

# ---- Hill profile (single Y axis, no facets) --------------------------------
ylim_hill_ITS1 <- range(c(hill_long_5k$diversity, hill_long_15k$diversity), na.rm = TRUE)
p_hill_5k  <- p_hill_5k  + coord_cartesian(ylim = ylim_hill_ITS1)
p_hill_15k <- p_hill_15k + coord_cartesian(ylim = ylim_hill_ITS1)

# ---- Convergence diagnostic -------------------------------------------------
# Separate limits for iter_median (raw) and running_mean (cumulative rows)
CONV_METRIC_ORDER_ITS1 <- c("Observed", "Shannon", "InvSimpson", "Berger-Parker")

conv_iter_lims_ITS1 <- bind_rows(
    conv_long_5k  |> select(metric, value = iter_median),
    conv_long_15k |> select(metric, value = iter_median)
) |>
    group_by(metric) |>
    summarise(lo = min(value, na.rm = TRUE), hi = max(value, na.rm = TRUE), .groups = "drop") |>
    arrange(match(metric, CONV_METRIC_ORDER_ITS1))

conv_cum_lims_ITS1 <- bind_rows(
    conv_long_5k  |> select(metric, value = running_mean),
    conv_long_15k |> select(metric, value = running_mean)
) |>
    group_by(metric) |>
    summarise(lo = min(value, na.rm = TRUE), hi = max(value, na.rm = TRUE), .groups = "drop") |>
    arrange(match(metric, CONV_METRIC_ORDER_ITS1))

y_scales_conv_iter_ITS1 <- setNames(
    lapply(seq_len(nrow(conv_iter_lims_ITS1)), \(i)
        scale_y_continuous(limits = c(conv_iter_lims_ITS1$lo[i], conv_iter_lims_ITS1$hi[i]))),
    CONV_METRIC_ORDER_ITS1
)
y_scales_conv_cum_ITS1 <- setNames(
    lapply(seq_len(nrow(conv_cum_lims_ITS1)), \(i)
        scale_y_continuous(limits = c(conv_cum_lims_ITS1$lo[i], conv_cum_lims_ITS1$hi[i]))),
    CONV_METRIC_ORDER_ITS1
)

p_conv_5k <- (
    (p_conv_raw_5k  + ggh4x::facetted_pos_scales(y = y_scales_conv_iter_ITS1)) /
    (p_conv_cum_5k  + ggh4x::facetted_pos_scales(y = y_scales_conv_cum_ITS1))
) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

p_conv_15k <- (
    (p_conv_raw_15k + ggh4x::facetted_pos_scales(y = y_scales_conv_iter_ITS1)) /
    (p_conv_cum_15k + ggh4x::facetted_pos_scales(y = y_scales_conv_cum_ITS1))
) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# ---- Farm management boxplots -----------------------------------------------
FARM_METRICS_ITS1 <- c("Observed_median", "Shannon_median", "InvSimpson_median",
                       "BergerParker_median")

farm_lims_ITS1 <- bind_rows(alpha_summary_5k, alpha_summary_15k) |>
    pivot_longer(cols = all_of(FARM_METRICS_ITS1), names_to = "metric", values_to = "value") |>
    group_by(metric) |>
    summarise(lo = min(value, na.rm = TRUE), hi = max(value, na.rm = TRUE), .groups = "drop") |>
    arrange(match(metric, FARM_METRICS_ITS1))

p_farm_mgmt_5k <- (
    plot_alpha_farm_management("Observed_median",     "Observed",      alpha_summary_5k,
                               ylim = c(farm_lims_ITS1$lo[1], farm_lims_ITS1$hi[1])) |
    plot_alpha_farm_management("Shannon_median",      "Shannon H",     alpha_summary_5k,
                               ylim = c(farm_lims_ITS1$lo[2], farm_lims_ITS1$hi[2])) |
    plot_alpha_farm_management("InvSimpson_median",   "InvSimpson",    alpha_summary_5k,
                               ylim = c(farm_lims_ITS1$lo[3], farm_lims_ITS1$hi[3])) |
    plot_alpha_farm_management("BergerParker_median", "Berger-Parker", alpha_summary_5k,
                               ylim = c(farm_lims_ITS1$lo[4], farm_lims_ITS1$hi[4]))
) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

p_farm_mgmt_15k <- (
    plot_alpha_farm_management("Observed_median",     "Observed",      alpha_summary_15k,
                               ylim = c(farm_lims_ITS1$lo[1], farm_lims_ITS1$hi[1])) |
    plot_alpha_farm_management("Shannon_median",      "Shannon H",     alpha_summary_15k,
                               ylim = c(farm_lims_ITS1$lo[2], farm_lims_ITS1$hi[2])) |
    plot_alpha_farm_management("InvSimpson_median",   "InvSimpson",    alpha_summary_15k,
                               ylim = c(farm_lims_ITS1$lo[3], farm_lims_ITS1$hi[3])) |
    plot_alpha_farm_management("BergerParker_median", "Berger-Parker", alpha_summary_15k,
                               ylim = c(farm_lims_ITS1$lo[4], farm_lims_ITS1$hi[4]))
) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# ---- Sample type boxplots ---------------------------------------------------
STYPE_METRICS_ITS1 <- c("Observed_median", "Shannon_median", "InvSimpson_median",
                        "BergerParker_median")
STYPE_LABELS_ITS1  <- c("Observed", "Shannon H", "InvSimpson", "Berger-Parker")

stype_lims_ITS1 <- bind_rows(alpha_summary_5k, alpha_summary_15k) |>
    pivot_longer(cols = all_of(STYPE_METRICS_ITS1), names_to = "metric", values_to = "value") |>
    group_by(metric) |>
    summarise(lo = min(value, na.rm = TRUE), hi = max(value, na.rm = TRUE), .groups = "drop") |>
    arrange(match(metric, STYPE_METRICS_ITS1))

y_scales_stype_ITS1 <- setNames(
    lapply(seq_len(nrow(stype_lims_ITS1)), \(i)
        scale_y_continuous(limits = c(stype_lims_ITS1$lo[i], stype_lims_ITS1$hi[i]))),
    STYPE_LABELS_ITS1
)

p_sample_type_5k  <- p_sample_type_5k  + ggh4x::facetted_pos_scales(y = y_scales_stype_ITS1)
p_sample_type_15k <- p_sample_type_15k + ggh4x::facetted_pos_scales(y = y_scales_stype_ITS1)

# ---- LMM farm boxplots (1.06× buffer to keep significance stars visible) ----
lmm_lims_ITS1 <- bind_rows(alpha_summary_5k, alpha_summary_15k) |>
    pivot_longer(cols = all_of(METRICS_LMM_FARM), names_to = "metric", values_to = "value") |>
    group_by(metric) |>
    summarise(lo = min(value, na.rm = TRUE),
              hi = max(value, na.rm = TRUE) * 1.06, .groups = "drop") |>
    arrange(match(metric, METRICS_LMM_FARM))

y_scales_lmm_ITS1 <- setNames(
    lapply(seq_len(nrow(lmm_lims_ITS1)), \(i)
        scale_y_continuous(limits = c(lmm_lims_ITS1$lo[i], lmm_lims_ITS1$hi[i]))),
    METRIC_LABELS_LMM
)

p_lmm_farm_5k  <- p_lmm_farm_5k  + ggh4x::facetted_pos_scales(y = y_scales_lmm_ITS1)
p_lmm_farm_15k <- p_lmm_farm_15k + ggh4x::facetted_pos_scales(y = y_scales_lmm_ITS1)

################################################################################
#### Depth-comparison combined figures #########################################
################################################################################

# Convergence diagnostic: side by side
p_conv_comb_ITS1 <- (p_conv_5k | p_conv_15k) +
    plot_annotation(title = "Rarefaction convergence diagnostic: 5k (left) vs 15k (right) — ITS1 OTU97")

ggsave(here("results", "figures", "07b_conv_depth_comparison.png"),
       plot = p_conv_comb_ITS1, width = 20, height = 8, dpi = 300)

# Hill diversity profile: side by side
p_hill_comb_ITS1 <- (p_hill_5k | p_hill_15k) +
    plot_annotation(title = "Hill diversity profile: 5k (left) vs 15k (right) — ITS1 OTU97")

ggsave(here("results", "figures", "07b_hill_depth_comparison.png"),
       plot = p_hill_comb_ITS1, width = 18, height = 5, dpi = 300)

# Farm management boxplots: stacked
p_farm_comb_ITS1 <- (p_farm_mgmt_5k / p_farm_mgmt_15k) +
    plot_annotation(title = "Alpha diversity by farm — 5k (top) vs 15k (bottom) — ITS1 OTU97")

ggsave(here("results", "figures", "07b_farm_mgmt_depth_comparison.png"),
       plot = p_farm_comb_ITS1, width = 14, height = 12, dpi = 300)

# Marginal means — sample_type: stacked
p_treat_comb_ITS1 <- (p_treat_5k / p_treat_15k) +
    plot_annotation(title = "Marginal means (sample_type): 5k (top) vs 15k (bottom) — ITS1 OTU97")

ggsave(here("results", "figures", "07b_marginal_treatment_depth_comparison.png"),
       plot = p_treat_comb_ITS1, width = 14, height = 10, dpi = 300)

# Overall sample_type boxplots (all farms pooled): stacked
p_sample_type_comb_ITS1 <- (p_sample_type_5k / p_sample_type_15k) +
    plot_annotation(title = "Alpha diversity (all farms): 5k (top) vs 15k (bottom) — ITS1 OTU97")

ggsave(here("results", "figures", "07b_sample_type_depth_comparison.png"),
       plot = p_sample_type_comb_ITS1, width = 16, height = 10, dpi = 300)

# LMM per-farm (bagged vs unbagged): stacked
p_lmm_farm_comb_ITS1 <- (p_lmm_farm_5k / p_lmm_farm_15k) +
    plot_annotation(title = "Per-farm LMM alpha diversity (bagged vs unbagged): 5k (top) vs 15k (bottom) — ITS1 OTU97")

ggsave(here("results", "figures", "07b_lmm_farm_depth_comparison.png"),
       plot = p_lmm_farm_comb_ITS1, width = 14, height = 10, dpi = 300)
