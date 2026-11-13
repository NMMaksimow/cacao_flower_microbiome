# ---- Libraries ----
library(here)
library(tidyverse)
library(lme4)
library(lmerTest)    # p-values for lmer()
library(glmmTMB)
library(DHARMa)
library(performance)
library(ggh4x)
library(patchwork)
library(cowplot)

# Data version and sample set: RESULTS_DIR, BIOSAMPLE_SUFFIX, SET_TAG
source(here("scripts", "utils_paths.R"))
# Model helpers shared with 07b: fit_alpha_model(), alpha_treatment_farm(),
# alpha_per_farm(), alpha_microscopy()
source(here("scripts", "utils_alpha.R"))

# ---- Parameters ----
# N_ITER and RAREFY_DEPTH* must match values used in 07a_16S_alpha_div_rarefaction.R
N_ITER           <- 999
RAREFY_DEPTH     <- 2000
RAREFY_DEPTH_10k <- 10000

# Farm order: inside_forest → full_sun (deforestation gradient)
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
MGMT_LEVELS <- c("inside_forest", "near_forest", "agroforest", "full_sun")

# Colour palette for bagged/unbagged (adjust hex if needed)
SAMPLE_TYPE_COLOURS <- c(bagged_flower = "#4C72B0", unbagged_flower = "#DD8452")
SAMPLE_TYPE_LABELS  <- c(bagged_flower = "bagged", unbagged_flower = "unbagged")

# ---- Input: load all four RDS outputs from rarefaction script ----
alpha_iters_2k    <- readRDS(file.path(RESULTS_DIR, "rds", paste0("alpha_16S_2k_999iters", BIOSAMPLE_SUFFIX, ".rds")))
alpha_summary_2k  <- readRDS(file.path(RESULTS_DIR, "rds", paste0("alpha_16S_2k_999iters_summary", BIOSAMPLE_SUFFIX, ".rds")))
alpha_iters_10k   <- readRDS(file.path(RESULTS_DIR, "rds", paste0("alpha_16S_10k_999iters", BIOSAMPLE_SUFFIX, ".rds")))
alpha_summary_10k <- readRDS(file.path(RESULTS_DIR, "rds", paste0("alpha_16S_10k_999iters_summary", BIOSAMPLE_SUFFIX, ".rds")))

################################################################################
#### 5. Convergence diagnostic — 2k ############################################
#### Per-iteration group median + cumulative running mean ######################
################################################################################

conv_long_2k <- alpha_iters_2k |>
    pivot_longer(
        cols      = c(Observed, Shannon, BergerParker),
        names_to  = "metric",
        values_to = "value"
    ) |>
    # canonical left-to-right order: least sensitive → most sensitive to rarefaction depth
    mutate(metric = factor(metric,
        levels = c("Observed", "Shannon", "BergerParker"),
        labels = c("Observed", "Shannon", "Berger-Parker"))) |>
    group_by(iteration, sample_type, metric) |>
    summarise(iter_median = median(value, na.rm = TRUE), .groups = "drop") |>
    group_by(sample_type, metric) |>
    arrange(iteration) |>
    mutate(running_mean = cummean(iter_median)) |>
    ungroup()

# Top row: raw per-iteration values (shows natural noise level)
p_conv_raw_2k <- conv_long_2k |>
    ggplot(aes(x = iteration, y = iter_median, colour = sample_type)) +
    geom_line(alpha = 0.8, linewidth = 0.5) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    labs(
        title  = "Per-iteration group median (2k)",
        x      = NULL,
        y      = "Median across samples",
        colour = NULL
    )

# Bottom row: cumulative running mean — flattens well before iteration 100
p_conv_cum_2k <- conv_long_2k |>
    ggplot(aes(x = iteration, y = running_mean, colour = sample_type)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    labs(
        title  = "Cumulative running mean (convergence diagnostic)",
        x      = "Iteration",
        y      = "Cumulative mean",
        colour = NULL
    )

p_conv_2k <- p_conv_raw_2k / p_conv_cum_2k +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

################################################################################
#### 5. Convergence diagnostic — 10k ###########################################
################################################################################

conv_long_10k <- alpha_iters_10k |>
    pivot_longer(
        cols      = c(Observed, Shannon, BergerParker),
        names_to  = "metric",
        values_to = "value"
    ) |>
    mutate(metric = factor(metric,
        levels = c("Observed", "Shannon", "BergerParker"),
        labels = c("Observed", "Shannon", "Berger-Parker"))) |>
    group_by(iteration, sample_type, metric) |>
    summarise(iter_median = median(value, na.rm = TRUE), .groups = "drop") |>
    group_by(sample_type, metric) |>
    arrange(iteration) |>
    mutate(running_mean = cummean(iter_median)) |>
    ungroup()

p_conv_raw_10k <- conv_long_10k |>
    ggplot(aes(x = iteration, y = iter_median, colour = sample_type)) +
    geom_line(alpha = 0.8, linewidth = 0.5) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    labs(
        title  = "Per-iteration group median (10k)",
        x      = NULL,
        y      = "Median across samples",
        colour = NULL
    )

p_conv_cum_10k <- conv_long_10k |>
    ggplot(aes(x = iteration, y = running_mean, colour = sample_type)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = SAMPLE_TYPE_COLOURS, labels = SAMPLE_TYPE_LABELS) +
    theme_bw() +
    labs(
        title  = "Cumulative running mean (convergence diagnostic)",
        x      = "Iteration",
        y      = "Cumulative mean",
        colour = NULL
    )

p_conv_10k <- p_conv_raw_10k / p_conv_cum_10k +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

################################################################################
#### 6a. Boxplots: bagged vs unbagged — 2k ####################################
################################################################################

p_sample_type_2k <- alpha_summary_2k |>
    pivot_longer(
        cols      = c(Observed_median, Shannon_median,
                      BergerParker_median, FaithPD_median),
        names_to  = "metric",
        values_to = "value"
    ) |>
    mutate(metric = factor(
        metric,
        levels = c("Observed_median", "Shannon_median",
                   "BergerParker_median", "FaithPD_median"),
        labels = c("Observed", "Shannon H", "Berger-Parker", "Faith's PD")
    )) |>
    ggplot(aes(x = sample_type, y = value, fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
    # ncol = 4: Observed, Shannon H, Berger-Parker, Faith's PD
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
        title    = "Bacterial alpha diversity (16S ASV) by flower treatment — 2k rarefaction",
        subtitle = sprintf("Median of %d rarefactions at %s reads per sample",
                           N_ITER, scales::comma(RAREFY_DEPTH)),
        x = "", y = ""
    )

################################################################################
#### 6a. Boxplots: bagged vs unbagged — 10k ###################################
################################################################################

p_sample_type_10k <- alpha_summary_10k |>
    pivot_longer(
        cols      = c(Observed_median, Shannon_median,
                      BergerParker_median, FaithPD_median),
        names_to  = "metric",
        values_to = "value"
    ) |>
    mutate(metric = factor(
        metric,
        levels = c("Observed_median", "Shannon_median",
                   "BergerParker_median", "FaithPD_median"),
        labels = c("Observed", "Shannon H", "Berger-Parker", "Faith's PD")
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
        title    = "Bacterial alpha diversity (16S ASV) by flower treatment — 10k rarefaction",
        subtitle = sprintf("Median of %d rarefactions at %s reads (note: fewer samples retained)",
                           N_ITER, scales::comma(RAREFY_DEPTH_10k)),
        x = "", y = ""
    )

################################################################################
#### 7a. Treatment and farm effects (global models) ###########################
#### metric ~ sample_type + farm_id + (1 | tree_id); Observed as NB GLMM. ######
#### Farm: likelihood ratio test (ML) against the model without farm_id. ######
#### BH within depth over the main-text metrics. ###############################
#### 2k = primary depth; 10k = depth sensitivity check. ########################
################################################################################

# Main-text metrics for 16S (InvSimpson dropped; it stays only as Hill q = 2)
METRICS_16S <- c("Observed_median", "Shannon_median", "BergerParker_median", "FaithPD_median")

treat_farm_16S <- alpha_treatment_farm(
    alpha_by_depth = list(`2k` = alpha_summary_2k, `10k` = alpha_summary_10k),
    metrics        = METRICS_16S,
    marker         = "16S"
)

print(treat_farm_16S$table, width = Inf)
write_csv(treat_farm_16S$table,
          file.path(RESULTS_DIR, "tables", paste0("07a_treatment_farm_effects", SET_TAG, ".csv")))

# Diagnostics: singular fits (zero tree variance) and residuals.
# LMM: normality of residuals; NB GLMM: DHARMa simulated residuals (uniformity, dispersion).
walk2(treat_farm_16S$models$fit, paste(treat_farm_16S$models$depth, treat_farm_16S$models$metric),
      \(fit, label) {
          cat("\n----", label, "----\n")
          if (inherits(fit, "glmmTMB")) {
              sim <- DHARMa::simulateResiduals(fit, n = 1000)
              print(DHARMa::testUniformity(sim, plot = FALSE))
              print(DHARMa::testDispersion(sim, plot = FALSE))
          } else {
              print(performance::check_singularity(fit))
              print(performance::check_normality(fit))
          }
      })

################################################################################
#### 7b. Sensitivity of the treatment effect to the model form (primary 2k) ####
#### Residuals of the Berger-Parker and Faith's PD LMMs are not normal: #######
####   Berger-Parker is a proportion -> beta GLMM (logit link); ################
####   Faith's PD is positive and right-skewed -> LMM on log(Faith's PD). ######
################################################################################

sens_data_2k <- alpha_summary_2k |>
    mutate(sample_type = factor(sample_type, levels = SAMPLE_TYPE_LEVELS))

sensitivity_16S <- bind_rows(
    primary_rows(treat_farm_16S$table, "2k", c("BergerParker_median", "FaithPD_median")),
    sensitivity_effect(BergerParker_median ~ sample_type + farm_id + (1 | tree_id),
                       sens_data_2k, family = glmmTMB::beta_family()) |>
        mutate(metric = "BergerParker_median", analysis = "beta GLMM (logit link)",
               scale = "log odds ratio"),
    sensitivity_effect(log(FaithPD_median) ~ sample_type + farm_id + (1 | tree_id),
                       sens_data_2k) |>
        mutate(metric = "FaithPD_median", analysis = "LMM on log(Faith's PD)",
               scale = "log ratio")
) |>
    mutate(marker = "16S", depth = "2k", .before = 1) |>
    arrange(metric)

print(sensitivity_16S, width = Inf)
write_csv(sensitivity_16S,
          file.path(RESULTS_DIR, "tables", paste0("07a_treatment_sensitivity", SET_TAG, ".csv")))

################################################################################
#### 7c. Microscopy features (primary depth 2k) ################################
#### Same tests as 10g: unbagged flowers (fungal colonization score, pollen, ##
#### callose plug, germination) and bagged flowers (pollen). ##################
#### metric ~ feature + farm_id + (1 | tree_id); BH over all tests. ###########
################################################################################

microscopy_16S <- alpha_microscopy(alpha_summary_2k, METRICS_16S) |>
    mutate(marker = "16S", depth = "2k", .before = 1)

print(microscopy_16S, width = Inf)
write_csv(microscopy_16S,
          file.path(RESULTS_DIR, "tables", paste0("07a_microscopy_effects", SET_TAG, ".csv")))

################################################################################
#### 9_icc. Within-tree ICC — are flowers from the same tree more similar? #####
#### Intercept-only LMM with (1 | tree_id) across all farms combined. ##########
#### ICC_adjusted = var(tree) / (var(tree) + var(residual)) ####################
#### Interpretation: < 0.1 negligible, 0.1–0.3 modest, > 0.3 meaningful. ######
#### High ICC → tree-level aggregation is justified; low ICC → aggregation ######
#### loses information without protecting against pseudoreplication. ############
################################################################################

METRICS_ICC_16S <- METRICS_16S

# Run at 2k depth (primary depth); ICC is a structural property of the sampling design
icc_df_2k <- expand.grid(
    metric      = METRICS_ICC_16S,
    sample_type = c("bagged_flower", "unbagged_flower"),
    stringsAsFactors = FALSE
) |>
    as_tibble() |>
    mutate(icc_adjusted = map2_dbl(metric, sample_type, \(m, st) {
        fit <- alpha_summary_2k |>
            filter(sample_type == st) |>
            lmer(as.formula(paste0(m, " ~ 1 + (1 | tree_id)")), data = _)
        # Zero tree variance (singular fit): icc() returns NA, not a table; ICC is 0
        icc_res <- performance::icc(fit)
        if (is.list(icc_res)) icc_res$ICC_adjusted else 0
    }))

print(icc_df_2k)
write_csv(icc_df_2k, file.path(RESULTS_DIR, "tables", paste0("07a_2k_icc_within_tree", SET_TAG, ".csv")))

################################################################################
#### 9b. Per-farm LMM: bagged vs unbagged at flower level (1 | tree_id) ########
#### Replaces paired Wilcoxon on tree medians (removed based on ICC results). ##
#### ICC (section 9_icc) showed low within-tree clustering, so flower-level ####
#### data are used directly; tree_id as random effect handles residual         ##
#### within-tree correlation without discarding observations.                  ##
################################################################################

# 4 main-text metrics for 16S (includes Faith's PD)
METRICS_LMM_FARM <- METRICS_16S
METRIC_LABELS_LMM <- c(
    Observed_median     = "Observed",
    Shannon_median      = "Shannon H",
    BergerParker_median = "Berger-Parker",
    FaithPD_median      = "Faith's PD"
)

# 7 farms × 4 metrics = 28 tests for BH correction
n_tests_lmm_2k <- length(FARM_LEVELS) * length(METRICS_LMM_FARM)

# Per-farm model (alpha_per_farm() in utils_alpha.R): metric ~ sample_type + (1 | tree_id)
# within each farm; Observed as NB GLMM (log scale), the others as LMM.
# Reference level: bagged_flower; positive estimate = higher in unbagged_flower.
# Degenerate farms (fit error) -> NA rows. Columns: estimate, se, conf_low, conf_high, p, q, stars.
lmm_farm_df_2k <- alpha_per_farm(alpha_summary_2k, METRICS_LMM_FARM, FARM_LEVELS)

write_csv(lmm_farm_df_2k, file.path(RESULTS_DIR, "tables", paste0("07a_2k_lmm_per_farm", SET_TAG, ".csv")))

# Star y-position: 5% above max flower-level value per farm × metric
star_pos_lmm_2k <- alpha_summary_2k |>
    pivot_longer(cols = all_of(METRICS_LMM_FARM), names_to = "metric", values_to = "value") |>
    group_by(farm_id, metric) |>
    summarise(star_y = max(value, na.rm = TRUE) * 1.05, .groups = "drop")

lmm_stars_2k <- lmm_farm_df_2k |>
    filter(stars != "") |>
    left_join(star_pos_lmm_2k, by = c("farm" = "farm_id", "metric"))

# Boxplot of raw flower-level data with LMM significance stars
p_lmm_farm_2k <- alpha_summary_2k |>
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
    # significance stars above each farm pair where q < 0.05 after BH
    # guard: empty lmm_stars_2k skips layer to avoid NULL PANEL in ggh4x::FacetWrap2
    (if (nrow(lmm_stars_2k) > 0)
        geom_text(
            data = lmm_stars_2k |>
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
    # ncol = 4: Observed, Shannon H, Berger-Parker, Faith's PD
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
        title    = "Bacterial alpha diversity — per-farm LMM (bagged vs unbagged) — 2k rarefaction",
        subtitle = sprintf(
            "Flower-level data; metric ~ sample_type + (1|tree_id) per farm (Observed NB GLMM, others LMM); BH (%d tests); * q < 0.05",
            n_tests_lmm_2k
        ),
        x = "", y = "", fill = NULL, colour = NULL
    )

# --- 10k depth: same LMM, secondary rarefaction depth ---

n_tests_lmm_10k <- length(FARM_LEVELS) * length(METRICS_LMM_FARM)

lmm_farm_df_10k <- alpha_per_farm(alpha_summary_10k, METRICS_LMM_FARM, FARM_LEVELS)

write_csv(lmm_farm_df_10k, file.path(RESULTS_DIR, "tables", paste0("07a_10k_lmm_per_farm", SET_TAG, ".csv")))

star_pos_lmm_10k <- alpha_summary_10k |>
    pivot_longer(cols = all_of(METRICS_LMM_FARM), names_to = "metric", values_to = "value") |>
    group_by(farm_id, metric) |>
    summarise(star_y = max(value, na.rm = TRUE) * 1.05, .groups = "drop")

lmm_stars_10k <- lmm_farm_df_10k |>
    filter(stars != "") |>
    left_join(star_pos_lmm_10k, by = c("farm" = "farm_id", "metric"))

p_lmm_farm_10k <- alpha_summary_10k |>
    mutate(farm_id = factor(farm_id, levels = FARM_LEVELS)) |>
    pivot_longer(cols = all_of(METRICS_LMM_FARM), names_to = "metric", values_to = "value") |>
    mutate(metric = factor(metric, levels = METRICS_LMM_FARM, labels = METRIC_LABELS_LMM)) |>
    ggplot(aes(x = farm_id, y = value, fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7,
                 position = position_dodge(width = 0.75)) +
    geom_jitter(aes(colour = sample_type),
                position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
                alpha = 0.4, size = 0.9, stroke = 0) +
    # guard: empty lmm_stars_10k skips layer to avoid NULL PANEL in ggh4x::FacetWrap2
    (if (nrow(lmm_stars_10k) > 0)
        geom_text(
            data = lmm_stars_10k |>
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
          axis.text.x      = element_text(size = 8, lineheight = 0.85)) +
    labs(
        title    = "Bacterial alpha diversity — per-farm LMM (bagged vs unbagged) — 10k rarefaction",
        subtitle = sprintf(
            "Flower-level data; metric ~ sample_type + (1|tree_id) per farm (Observed NB GLMM, others LMM); BH (%d tests); * q < 0.05",
            n_tests_lmm_10k
        ),
        x = "", y = "", fill = NULL, colour = NULL
    )

################################################################################
#### 10. Combined depth-comparison figures — 16S ##############################
#### Left / top  = 2k (lower depth, more samples retained) ####################
#### Right / bottom = 10k (higher depth, fewer samples retained) ###############
################################################################################

################################################################################
#### Shared Y-axis limits for depth-comparison combined figures ################
#### Compute per-metric global range from both depths; apply before ggsave. ####
################################################################################

# ---- Convergence diagnostic -------------------------------------------------
# Separate limits for iter_median (raw) and running_mean (cumulative rows)
CONV_METRIC_ORDER_16S <- c("Observed", "Shannon", "Berger-Parker")

conv_iter_lims_16S <- bind_rows(
    conv_long_2k  |> select(metric, value = iter_median),
    conv_long_10k |> select(metric, value = iter_median)
) |>
    group_by(metric) |>
    summarise(lo = min(value, na.rm = TRUE), hi = max(value, na.rm = TRUE), .groups = "drop") |>
    arrange(match(metric, CONV_METRIC_ORDER_16S))

conv_cum_lims_16S <- bind_rows(
    conv_long_2k  |> select(metric, value = running_mean),
    conv_long_10k |> select(metric, value = running_mean)
) |>
    group_by(metric) |>
    summarise(lo = min(value, na.rm = TRUE), hi = max(value, na.rm = TRUE), .groups = "drop") |>
    arrange(match(metric, CONV_METRIC_ORDER_16S))

y_scales_conv_iter_16S <- setNames(
    lapply(seq_len(nrow(conv_iter_lims_16S)), \(i)
        scale_y_continuous(limits = c(conv_iter_lims_16S$lo[i], conv_iter_lims_16S$hi[i]))),
    CONV_METRIC_ORDER_16S
)
y_scales_conv_cum_16S <- setNames(
    lapply(seq_len(nrow(conv_cum_lims_16S)), \(i)
        scale_y_continuous(limits = c(conv_cum_lims_16S$lo[i], conv_cum_lims_16S$hi[i]))),
    CONV_METRIC_ORDER_16S
)

p_conv_2k <- (
    (p_conv_raw_2k  + ggh4x::facetted_pos_scales(y = y_scales_conv_iter_16S)) /
    (p_conv_cum_2k  + ggh4x::facetted_pos_scales(y = y_scales_conv_cum_16S))
) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

p_conv_10k <- (
    (p_conv_raw_10k + ggh4x::facetted_pos_scales(y = y_scales_conv_iter_16S)) /
    (p_conv_cum_10k + ggh4x::facetted_pos_scales(y = y_scales_conv_cum_16S))
) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# ---- Sample type boxplots ---------------------------------------------------
STYPE_METRICS_16S <- METRICS_16S
STYPE_LABELS_16S  <- c("Observed", "Shannon H", "Berger-Parker", "Faith's PD")

stype_lims_16S <- bind_rows(alpha_summary_2k, alpha_summary_10k) |>
    pivot_longer(cols = all_of(STYPE_METRICS_16S), names_to = "metric", values_to = "value") |>
    group_by(metric) |>
    summarise(lo = min(value, na.rm = TRUE), hi = max(value, na.rm = TRUE), .groups = "drop") |>
    arrange(match(metric, STYPE_METRICS_16S))

y_scales_stype_16S <- setNames(
    lapply(seq_len(nrow(stype_lims_16S)), \(i)
        scale_y_continuous(limits = c(stype_lims_16S$lo[i], stype_lims_16S$hi[i]))),
    STYPE_LABELS_16S
)

p_sample_type_2k  <- p_sample_type_2k  + ggh4x::facetted_pos_scales(y = y_scales_stype_16S)
p_sample_type_10k <- p_sample_type_10k + ggh4x::facetted_pos_scales(y = y_scales_stype_16S)

# ---- LMM farm boxplots (1.06× buffer to keep significance stars visible) ----
lmm_lims_16S <- bind_rows(alpha_summary_2k, alpha_summary_10k) |>
    pivot_longer(cols = all_of(METRICS_LMM_FARM), names_to = "metric", values_to = "value") |>
    group_by(metric) |>
    summarise(lo = min(value, na.rm = TRUE),
              hi = max(value, na.rm = TRUE) * 1.06, .groups = "drop") |>
    arrange(match(metric, METRICS_LMM_FARM))

y_scales_lmm_16S <- setNames(
    lapply(seq_len(nrow(lmm_lims_16S)), \(i)
        scale_y_continuous(limits = c(lmm_lims_16S$lo[i], lmm_lims_16S$hi[i]))),
    METRIC_LABELS_LMM
)

p_lmm_farm_2k  <- p_lmm_farm_2k  + ggh4x::facetted_pos_scales(y = y_scales_lmm_16S)
p_lmm_farm_10k <- p_lmm_farm_10k + ggh4x::facetted_pos_scales(y = y_scales_lmm_16S)

################################################################################
#### Depth-comparison combined figures #########################################
################################################################################

# Convergence diagnostic: side by side
# Each sub-panel title already states the depth ("Per-iteration group median (2k)" etc.)
p_conv_comb_16S <- (p_conv_2k | p_conv_10k) +
    plot_annotation(title = "Rarefaction convergence diagnostic: 2k (left) vs 10k (right) — 16S")

ggsave(file.path(RESULTS_DIR, "figures", paste0("07a_conv_depth_comparison", SET_TAG, ".png")),
       plot = p_conv_comb_16S, width = 16, height = 8, dpi = 300)

# Hill diversity profile: moved to 07c (16S | ITS1, both depths)

# Treatment effects: forest plot in 07c (from 07a_treatment_farm_effects and the per-farm tables)

# Overall sample_type boxplots (all farms pooled): stacked
p_sample_type_comb_16S <- (p_sample_type_2k / p_sample_type_10k) +
    plot_annotation(title = "Alpha diversity (all farms): 2k (top) vs 10k (bottom) — 16S")

ggsave(file.path(RESULTS_DIR, "figures", paste0("07a_sample_type_depth_comparison", SET_TAG, ".png")),
       plot = p_sample_type_comb_16S, width = 13, height = 10, dpi = 300)

# LMM per-farm (bagged vs unbagged): stacked
p_lmm_farm_comb_16S <- (p_lmm_farm_2k / p_lmm_farm_10k) +
    plot_annotation(title = "Per-farm LMM alpha diversity (bagged vs unbagged): 2k (top) vs 10k (bottom) — 16S")

ggsave(file.path(RESULTS_DIR, "figures", paste0("07a_lmm_farm_depth_comparison", SET_TAG, ".png")),
       plot = p_lmm_farm_comb_16S, width = 15, height = 10, dpi = 300)
