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
# Model helpers shared with 07a: fit_alpha_model(), alpha_treatment_farm(),
# alpha_per_farm(), alpha_microscopy()
source(here("scripts", "utils_alpha.R"))

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
alpha_iters_5k    <- readRDS(file.path(RESULTS_DIR, "rds", paste0("alpha_ITS1_otu97_5k_999iters", BIOSAMPLE_SUFFIX, ".rds")))
alpha_summary_5k  <- readRDS(file.path(RESULTS_DIR, "rds", paste0("alpha_ITS1_otu97_5k_999iters_summary", BIOSAMPLE_SUFFIX, ".rds")))
alpha_iters_15k   <- readRDS(file.path(RESULTS_DIR, "rds", paste0("alpha_ITS1_otu97_15k_999iters", BIOSAMPLE_SUFFIX, ".rds")))
alpha_summary_15k <- readRDS(file.path(RESULTS_DIR, "rds", paste0("alpha_ITS1_otu97_15k_999iters_summary", BIOSAMPLE_SUFFIX, ".rds")))

################################################################################
#### 5. Convergence diagnostic — 5k ############################################
#### Per-iteration group median + cumulative running mean ######################
################################################################################

conv_long_5k <- alpha_iters_5k |>
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
p_conv_raw_5k <- conv_long_5k |>
    ggplot(aes(x = iteration, y = iter_median, colour = sample_type)) +
    geom_line(alpha = 0.8, linewidth = 0.5) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
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
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
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

p_conv_raw_15k <- conv_long_15k |>
    ggplot(aes(x = iteration, y = iter_median, colour = sample_type)) +
    geom_line(alpha = 0.8, linewidth = 0.5) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
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
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
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
                      BergerParker_median),
        names_to  = "metric",
        values_to = "value"
    ) |>
    mutate(metric = factor(
        metric,
        levels = c("Observed_median", "Shannon_median",
                   "BergerParker_median"),
        labels = c("Observed", "Shannon H", "Berger-Parker")
    )) |>
    ggplot(aes(x = sample_type, y = value, fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
    # ncol = 3: Observed, Shannon H, Berger-Parker
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
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
                      BergerParker_median),
        names_to  = "metric",
        values_to = "value"
    ) |>
    mutate(metric = factor(
        metric,
        levels = c("Observed_median", "Shannon_median",
                   "BergerParker_median"),
        labels = c("Observed", "Shannon H", "Berger-Parker")
    )) |>
    ggplot(aes(x = sample_type, y = value, fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
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
#### 7a. Treatment and farm effects (global models) ###########################
#### metric ~ sample_type + farm_id + (1 | tree_id); Observed as NB GLMM. ######
#### Farm: likelihood ratio test (ML) against the model without farm_id. ######
#### BH within depth over the main-text metrics. ###############################
#### 5k = primary depth; 15k = depth sensitivity check. ########################
################################################################################

# Main-text metrics for ITS1 (no Faith PD; InvSimpson dropped, it stays only as Hill q = 2)
METRICS_ITS1 <- c("Observed_median", "Shannon_median", "BergerParker_median")

treat_farm_ITS1 <- alpha_treatment_farm(
    alpha_by_depth = list(`5k` = alpha_summary_5k, `15k` = alpha_summary_15k),
    metrics        = METRICS_ITS1,
    marker         = "ITS1"
)

print(treat_farm_ITS1$table, width = Inf)
write_csv(treat_farm_ITS1$table,
          file.path(RESULTS_DIR, "tables", paste0("07b_treatment_farm_effects", SET_TAG, ".csv")))

# Diagnostics: singular fits (zero tree variance) and residuals.
# LMM: normality of residuals; NB GLMM: DHARMa simulated residuals (uniformity, dispersion).
walk2(treat_farm_ITS1$models$fit, paste(treat_farm_ITS1$models$depth, treat_farm_ITS1$models$metric),
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
#### 7b. Sensitivity of the treatment effect (primary depth 5k) ################
#### (1) Model form: Berger-Parker is a proportion -> beta GLMM (logit link). ##
#### (2) Depth and PCR yield: ITS1 reads and PCR product concentration are ####
####     lower in bagged flowers (04 check). Add them as covariates. One-way ###
####     test: an effect that stays is strong evidence; one that goes proves ###
####     nothing (both may be consequences of the treatment). #################
################################################################################

# Reads per flower before rarefaction (post-decontam fungal OTU97 object)
reads_ITS1 <- readRDS(file.path(RESULTS_DIR, "rds", paste0("ps_ITS1_otu97_fungi_biosamples", BIOSAMPLE_SUFFIX, ".rds"))) |>
    phyloseq::sample_sums()

sens_data_5k <- alpha_summary_5k |>
    mutate(
        sample_type = factor(sample_type, levels = SAMPLE_TYPE_LEVELS),
        log_reads   = log(reads_ITS1[sample_id]),        # reads before rarefaction
        log_pcr     = log(dna_concentration_hs_its1)     # pooled ITS1 PCR product, ng/ul
    )

# Response as written in the formula (Observed: integer counts for the NB GLMM)
SENS_RESPONSE_ITS1 <- c(Observed_median     = "round(Observed_median)",
                        Shannon_median      = "Shannon_median",
                        BergerParker_median = "BergerParker_median")
COVARIATE_LABELS <- c(log_reads = "log(reads before rarefaction)",
                      log_pcr   = "log(PCR product concentration)")

covariate_rows_ITS1 <- expand_grid(metric = METRICS_ITS1,
                                   covariates = c("log_reads", "log_reads + log_pcr")) |>
    mutate(
        effect = map2(metric, covariates, \(m, cv) sensitivity_effect(
            as.formula(paste(SENS_RESPONSE_ITS1[[m]], "~ sample_type +", cv, "+ farm_id + (1 | tree_id)")),
            sens_data_5k,
            family = if (m == "Observed_median") glmmTMB::nbinom2() else NULL)),
        analysis = paste("primary model +", str_replace_all(covariates, COVARIATE_LABELS)),
        scale    = if_else(metric == "Observed_median", "log ratio", "difference")
    ) |>
    unnest(effect) |>
    select(-covariates)

sensitivity_ITS1 <- bind_rows(
    primary_rows(treat_farm_ITS1$table, "5k", METRICS_ITS1),
    sensitivity_effect(BergerParker_median ~ sample_type + farm_id + (1 | tree_id),
                       sens_data_5k, family = glmmTMB::beta_family()) |>
        mutate(metric = "BergerParker_median", analysis = "beta GLMM (logit link)",
               scale = "log odds ratio"),
    covariate_rows_ITS1
) |>
    mutate(marker = "ITS1", depth = "5k", .before = 1) |>
    arrange(match(metric, METRICS_ITS1))

print(sensitivity_ITS1, width = Inf)
write_csv(sensitivity_ITS1,
          file.path(RESULTS_DIR, "tables", paste0("07b_treatment_sensitivity", SET_TAG, ".csv")))

################################################################################
#### 7c. Microscopy features (primary depth 5k) ################################
#### Same tests as 10g: unbagged flowers (fungal colonization score, pollen, ##
#### callose plug, germination) and bagged flowers (pollen). ##################
#### metric ~ feature + farm_id + (1 | tree_id); BH over all tests. ###########
################################################################################

microscopy_ITS1 <- alpha_microscopy(alpha_summary_5k, METRICS_ITS1) |>
    mutate(marker = "ITS1", depth = "5k", .before = 1)

print(microscopy_ITS1, width = Inf)
write_csv(microscopy_ITS1,
          file.path(RESULTS_DIR, "tables", paste0("07b_microscopy_effects", SET_TAG, ".csv")))

################################################################################
#### 9_icc. Within-tree ICC — are flowers from the same tree more similar? #####
#### Intercept-only LMM with (1 | tree_id) across all farms combined. ##########
#### ICC_adjusted = var(tree) / (var(tree) + var(residual)) ####################
#### Interpretation: < 0.1 negligible, 0.1–0.3 modest, > 0.3 meaningful. ######
#### High ICC → tree-level aggregation is justified; low ICC → aggregation ######
#### loses information without protecting against pseudoreplication. ############
################################################################################

METRICS_ICC_ITS1 <- METRICS_ITS1

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
        # Zero tree variance (singular fit): icc() returns NA, not a table; ICC is 0
        icc_res <- performance::icc(fit)
        if (is.list(icc_res)) icc_res$ICC_adjusted else 0
    }))

print(icc_df_5k)
write_csv(icc_df_5k, file.path(RESULTS_DIR, "tables", paste0("07b_5k_icc_within_tree", SET_TAG, ".csv")))

################################################################################
#### 9b. Per-farm LMM: bagged vs unbagged at flower level (1 | tree_id) ########
#### Replaces paired Wilcoxon on tree medians (removed based on ICC results). ##
#### ICC (section 9_icc) showed moderate within-tree clustering for ITS1;     ##
#### LMM models this directly via random intercept without aggregating.        ##
################################################################################

# 3 main-text metrics for ITS1 (no Faith's PD — ITS too variable for phylogenetic inference)
METRICS_LMM_FARM <- METRICS_ITS1
METRIC_LABELS_LMM <- c(
    Observed_median     = "Observed",
    Shannon_median      = "Shannon H",
    BergerParker_median = "Berger-Parker"
)

# 7 farms × 3 metrics = 21 tests for BH correction
n_tests_lmm_5k <- length(FARM_LEVELS) * length(METRICS_LMM_FARM)

# Per-farm model (alpha_per_farm() in utils_alpha.R): metric ~ sample_type + (1 | tree_id)
# within each farm; Observed as NB GLMM (log scale), the others as LMM.
# Reference level: bagged_flower; positive estimate = higher in unbagged_flower.
# Degenerate farms (fit error) -> NA rows. Columns: estimate, se, conf_low, conf_high, p, q, stars.
lmm_farm_df_5k <- alpha_per_farm(alpha_summary_5k, METRICS_LMM_FARM, FARM_LEVELS)

write_csv(lmm_farm_df_5k, file.path(RESULTS_DIR, "tables", paste0("07b_5k_lmm_per_farm", SET_TAG, ".csv")))

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
    # significance stars above each farm pair where q < 0.05 after BH
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
    # ncol = 3: Observed, Shannon H, Berger-Parker
    # facet_wrap2 (ggh4x) required for facetted_pos_scales compatibility
    facet_wrap2(~ metric, scales = "free_y", ncol = 3) +
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
            "Flower-level data; metric ~ sample_type + (1|tree_id) per farm (Observed NB GLMM, others LMM); BH (%d tests); * q < 0.05",
            n_tests_lmm_5k
        ),
        x = "", y = "", fill = NULL, colour = NULL
    )

# --- 15k depth: same LMM, secondary rarefaction depth ---

n_tests_lmm_15k <- length(FARM_LEVELS) * length(METRICS_LMM_FARM)

lmm_farm_df_15k <- alpha_per_farm(alpha_summary_15k, METRICS_LMM_FARM, FARM_LEVELS)

write_csv(lmm_farm_df_15k, file.path(RESULTS_DIR, "tables", paste0("07b_15k_lmm_per_farm", SET_TAG, ".csv")))

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
    facet_wrap2(~ metric, scales = "free_y", ncol = 3) +
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
            "Flower-level data; metric ~ sample_type + (1|tree_id) per farm (Observed NB GLMM, others LMM); BH (%d tests); * q < 0.05",
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

# ---- Convergence diagnostic -------------------------------------------------
# Separate limits for iter_median (raw) and running_mean (cumulative rows)
CONV_METRIC_ORDER_ITS1 <- c("Observed", "Shannon", "Berger-Parker")

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

# ---- Sample type boxplots ---------------------------------------------------
STYPE_METRICS_ITS1 <- METRICS_ITS1
STYPE_LABELS_ITS1  <- c("Observed", "Shannon H", "Berger-Parker")

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

ggsave(file.path(RESULTS_DIR, "figures", paste0("07b_conv_depth_comparison", SET_TAG, ".png")),
       plot = p_conv_comb_ITS1, width = 16, height = 8, dpi = 300)

# Hill diversity profile: moved to 07c (16S | ITS1, both depths)

# Treatment effects: forest plot in 07c (from 07b_treatment_farm_effects and the per-farm tables)

# Overall sample_type boxplots (all farms pooled): stacked
p_sample_type_comb_ITS1 <- (p_sample_type_5k / p_sample_type_15k) +
    plot_annotation(title = "Alpha diversity (all farms): 5k (top) vs 15k (bottom) — ITS1 OTU97")

ggsave(file.path(RESULTS_DIR, "figures", paste0("07b_sample_type_depth_comparison", SET_TAG, ".png")),
       plot = p_sample_type_comb_ITS1, width = 10, height = 10, dpi = 300)

# LMM per-farm (bagged vs unbagged): stacked
p_lmm_farm_comb_ITS1 <- (p_lmm_farm_5k / p_lmm_farm_15k) +
    plot_annotation(title = "Per-farm LMM alpha diversity (bagged vs unbagged): 5k (top) vs 15k (bottom) — ITS1 OTU97")

ggsave(file.path(RESULTS_DIR, "figures", paste0("07b_lmm_farm_depth_comparison", SET_TAG, ".png")),
       plot = p_lmm_farm_comb_ITS1, width = 12, height = 10, dpi = 300)
