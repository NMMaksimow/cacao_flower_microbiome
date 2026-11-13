# Requires tidyverse (dplyr, purrr, tidyr, tibble) to be loaded before sourcing.
# Shared model helpers for alpha diversity (07a 16S, 07b ITS1).
#
# One model family for all metrics:
#   Observed (integer richness, overdispersed) -> negative binomial GLMM (glmmTMB, nbinom2),
#                                                 effect on the log scale (exp = ratio)
#   Shannon, Berger-Parker, Faith PD           -> LMM (lmerTest::lmer), effect as a difference
# Random intercept for tree: flowers of one tree are not independent.
# Effects: estimate, SE, Wald 95% CI and p from summary() of the fitted model
#   (lmerTest: Satterthwaite df, t quantile; glmmTMB: z).

library(lmerTest)   # lmer() with Satterthwaite p-values
library(glmmTMB)

# Treatment reference level: positive estimate = higher in unbagged flowers
SAMPLE_TYPE_LEVELS <- c("bagged_flower", "unbagged_flower")

# ---- Fit one model ----------------------------------------------------------
# metric: column name, e.g. "Shannon_median"; rhs: right-hand side of the formula
# REML = FALSE for likelihood ratio tests of fixed effects (ignored for the NB GLMM,
# which glmmTMB fits by ML).
fit_alpha_model <- function(metric, rhs, data, REML = TRUE) {
    if (metric == "Observed_median") {
        glmmTMB::glmmTMB(as.formula(paste("round(Observed_median) ~", rhs)),
                         family = glmmTMB::nbinom2, data = data)
    } else {
        lmerTest::lmer(as.formula(paste(metric, "~", rhs)), data = data, REML = REML)
    }
}

# ---- Effect of one term -------------------------------------------------------
# Returns estimate, SE, Wald 95% CI and p for the coefficient named `term`.
# NA row if the term is missing (e.g. a level absent in a small subset).
tidy_term <- function(fit, term) {
    is_nb    <- inherits(fit, "glmmTMB")
    coef_tbl <- if (is_nb) summary(fit)$coefficients$cond else summary(fit)$coefficients
    if (!term %in% rownames(coef_tbl))
        return(tibble(estimate = NA_real_, se = NA_real_, conf_low = NA_real_,
                      conf_high = NA_real_, p = NA_real_))
    est <- coef_tbl[term, "Estimate"]
    se  <- coef_tbl[term, "Std. Error"]
    # critical value: z for the GLMM, t with Satterthwaite df for the LMM
    crit <- if (is_nb) qnorm(0.975) else qt(0.975, df = coef_tbl[term, "df"])
    tibble(
        estimate  = est,
        se        = se,
        conf_low  = est - crit * se,
        conf_high = est + crit * se,
        p         = coef_tbl[term, ncol(coef_tbl)]   # last column: Pr(>|t|) or Pr(>|z|)
    )
}

# ---- Likelihood ratio test of one fixed term ---------------------------------
# Full vs reduced model, both fitted by ML. Returns chi-square, df, p.
lrt_term <- function(metric, rhs_full, rhs_reduced, data) {
    full    <- fit_alpha_model(metric, rhs_full,    data, REML = FALSE)
    reduced <- fit_alpha_model(metric, rhs_reduced, data, REML = FALSE)
    ll_full <- logLik(full)
    ll_red  <- logLik(reduced)
    chisq   <- 2 * (as.numeric(ll_full) - as.numeric(ll_red))
    df      <- attr(ll_full, "df") - attr(ll_red, "df")
    tibble(chisq = chisq, df = df, p = pchisq(chisq, df, lower.tail = FALSE))
}

# ---- Global treatment + farm table -------------------------------------------
# alpha_by_depth: named list of summary tibbles, e.g. list(`2k` = ..., `10k` = ...)
# Model: metric ~ sample_type + farm_id + (1 | tree_id)
# Farm: LRT against metric ~ sample_type + (1 | tree_id)
# BH within depth over the metrics, separately for treatment and farm.
# Returns list(table = tibble, models = tibble with the fitted REML models).
alpha_treatment_farm <- function(alpha_by_depth, metrics, marker) {
    models <- tidyr::expand_grid(depth = names(alpha_by_depth), metric = metrics) |>
        mutate(
            data = map(depth, \(d) alpha_by_depth[[d]] |>
                           mutate(sample_type = factor(sample_type, levels = SAMPLE_TYPE_LEVELS))),
            fit  = map2(metric, data, \(m, d)
                        fit_alpha_model(m, "sample_type + farm_id + (1 | tree_id)", d))
        )

    table <- models |>
        mutate(
            n_flowers = map_int(fit, nobs),
            model     = if_else(metric == "Observed_median", "NB GLMM (log link)", "LMM"),
            effect    = map(fit, \(f) tidy_term(f, "sample_typeunbagged_flower")),
            farm      = map2(metric, data, \(m, d)
                             lrt_term(m, "sample_type + farm_id + (1 | tree_id)",
                                      "sample_type + (1 | tree_id)", d) |>
                                 rename_with(\(x) paste0("farm_", x)))
        ) |>
        select(depth, metric, n_flowers, model, effect, farm) |>
        unnest(c(effect, farm)) |>
        group_by(depth) |>
        mutate(q = p.adjust(p, method = "BH"), farm_q = p.adjust(farm_p, method = "BH")) |>
        ungroup() |>
        mutate(
            marker = marker,
            # Observed: ratio unbagged / bagged (exp of the log-scale effect)
            ratio      = if_else(metric == "Observed_median", exp(estimate),  NA_real_),
            ratio_low  = if_else(metric == "Observed_median", exp(conf_low),  NA_real_),
            ratio_high = if_else(metric == "Observed_median", exp(conf_high), NA_real_)
        ) |>
        select(marker, depth, metric, n_flowers, model, estimate, se, conf_low, conf_high,
               p, q, ratio, ratio_low, ratio_high, farm_chisq, farm_df, farm_p, farm_q)

    list(table = table, models = models)
}

# ---- Sensitivity analyses of the treatment effect --------------------------------
# Refit the global model in another form and return the treatment effect, to compare
# with the primary model. formula: full model formula, written out by the caller;
# family: NULL = LMM (lmerTest), otherwise a glmmTMB family (e.g. glmmTMB::beta_family()).
sensitivity_effect <- function(formula, data, family = NULL) {
    fit <- if (is.null(family)) {
        lmerTest::lmer(formula, data = data)
    } else {
        glmmTMB::glmmTMB(formula, family = family, data = data)
    }
    tidy_term(fit, "sample_typeunbagged_flower") |>
        mutate(n_flowers = nobs(fit), .before = 1)
}

# Primary-model rows of alpha_treatment_farm()$table in the layout of the sensitivity table
primary_rows <- function(effects_table, depth_name, metrics) {
    effects_table |>
        filter(depth == depth_name, metric %in% metrics) |>
        transmute(metric,
                  analysis = paste("primary:", model),
                  scale    = if_else(metric == "Observed_median", "log ratio", "difference"),
                  n_flowers, estimate, se, conf_low, conf_high, p)
}

# ---- Per-farm treatment effects ------------------------------------------------
# Within each farm: metric ~ sample_type + (1 | tree_id); BH over farms x metrics.
# Degenerate fits (error) give NA rows. `stars` marks q < 0.05 for the boxplots.
alpha_per_farm <- function(alpha_sum, metrics, farm_levels) {
    tidyr::expand_grid(farm = farm_levels, metric = metrics) |>
        mutate(effect = map2(farm, metric, \(f, m) {
            dat <- alpha_sum |>
                filter(farm_id == f) |>
                mutate(sample_type = factor(sample_type, levels = SAMPLE_TYPE_LEVELS))
            fit <- tryCatch(fit_alpha_model(m, "sample_type + (1 | tree_id)", dat),
                            error = \(e) NULL)
            if (is.null(fit))
                return(tibble(n_flowers = nrow(dat), estimate = NA_real_, se = NA_real_,
                              conf_low = NA_real_, conf_high = NA_real_, p = NA_real_))
            tidy_term(fit, "sample_typeunbagged_flower") |>
                mutate(n_flowers = nobs(fit), .before = 1)
        })) |>
        unnest(effect) |>
        mutate(
            model = if_else(metric == "Observed_median", "NB GLMM (log link)", "LMM"),
            q     = p.adjust(p, method = "BH"),
            stars = case_when(
                q < 0.001 ~ "***",
                q < 0.01  ~ "**",
                q < 0.05  ~ "*",
                TRUE      ~ ""
            )
        )
}

# ---- Microscopy features ----------------------------------------------------------
# Same tests and transforms as 10g (scripts/10g_16S_ITS1_CLR_permanova_microscopy.R).
MICROSCOPY_TESTS <- tibble::tribble(
    ~flowers,           ~feature,                     ~transform,
    "unbagged_flower",  "fungal_colonization_score",  "numeric",
    "unbagged_flower",  "pi_clsm",                    "log1p",
    "unbagged_flower",  "callose_plug_index",         "log1p",
    "unbagged_flower",  "is_germination",             "factor",
    "bagged_flower",    "pi_clsm",                    "log1p"
)

transform_feature <- function(x, transform) {
    switch(transform,
           numeric = as.numeric(as.character(x)),   # ordered factor 0-5 -> 0..5
           log1p   = log1p(x),
           factor  = factor(x))
}

# One microscopy test: flowers of group `flowers` with a value of `feature`,
# feature transformed to x; returns n_flowers and the effect of x.
fit_microscopy_test <- function(alpha_sum, flowers, feature, transform, metric) {
    keep <- alpha_sum$sample_type == flowers & !is.na(alpha_sum[[feature]])
    dat  <- alpha_sum[keep, ]
    dat$x <- transform_feature(dat[[feature]], transform)
    fit  <- fit_alpha_model(metric, "x + farm_id + (1 | tree_id)", dat)
    coef_names <- rownames(if (inherits(fit, "glmmTMB"))
        summary(fit)$coefficients$cond else summary(fit)$coefficients)
    tidy_term(fit, grep("^x", coef_names, value = TRUE)[1]) |>
        mutate(n_flowers = nobs(fit), .before = 1)
}

# Within one flower group, flowers with a value: metric ~ x + farm_id + (1 | tree_id).
# The effect row is the coefficient of x (for a factor: its second level).
# BH over all tests (features x metrics) within the marker, as in 10g.
alpha_microscopy <- function(alpha_sum, metrics, tests = MICROSCOPY_TESTS) {
    tidyr::expand_grid(tests, metric = metrics) |>
        mutate(effect = pmap(list(flowers, feature, transform, metric),
                             \(fl, ft, tr, m) fit_microscopy_test(alpha_sum, fl, ft, tr, m))) |>
        unnest(effect) |>
        mutate(
            model = if_else(metric == "Observed_median", "NB GLMM (log link)", "LMM"),
            q     = p.adjust(p, method = "BH")
        )
}
