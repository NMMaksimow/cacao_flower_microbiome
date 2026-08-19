# ---- 0. Parameters -------------------------------------------------------
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
# Max zeta order: 10 is appropriate for 21 samples per group (3 flowers × 7 trees).
# Order 3 ≈ within-tree scale; orders 4–10 ≈ across-tree within-farm.
ZETA_ORDERS <- 1:10
# MC iterations: 1000 gives smooth curves for n=21; increase to 2000 if noisy.
SAM <- 1000L

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(phyloseq)
library(zetadiv)   # Zeta.decline.mc()

# ---- 2. Load data --------------------------------------------------------
# OTU97 phyloseq: consistent with all other ITS1 diversity analyses.
# No prevalence filter before zeta: presence/absence method; rare OTUs contribute
# to turnover and should not be removed.
ps <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds")) |>
  subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))

# Incidence matrix: samples × OTUs (0/1)
incidence <- ps |>
  otu_table() |>
  as("matrix") |>
  (\(m) if (taxa_are_rows(ps)) t(m) else m)() |>
  (\(m) { m[m > 0] <- 1L; as.data.frame(m) })()

meta <- as(sample_data(ps), "data.frame")

cat(sprintf("Loaded: %d samples × %d OTUs\n", nrow(incidence), ncol(incidence)))
stopifnot(all(rownames(incidence) == rownames(meta)))

# ---- 3. Per-farm zeta decline --------------------------------------------
# Same approach as 11a_16S_zeta_by_farm.R; results stored in nested list.

dir.create(here("results", "rds"),    showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)

zeta_results <- list()

for (farm in FARM_LEVELS) {
  zeta_results[[farm]] <- list()
  farm_rows <- which(meta$farm_id == farm)
  inc_farm  <- incidence[farm_rows, ]
  meta_farm <- meta[farm_rows, ]

  for (grp in c("bagged_flower", "unbagged_flower")) {
    grp_rows  <- which(meta_farm$sample_type == grp)
    n_samples <- length(grp_rows)
    cat(sprintf("[%s × %s] %d samples — computing orders 1:%d, sam=%d ...\n",
                farm, grp, n_samples, max(ZETA_ORDERS), SAM))

    if (n_samples < max(ZETA_ORDERS)) {
      warning(sprintf("Fewer samples (%d) than max order (%d); capping orders at %d",
                      n_samples, max(ZETA_ORDERS), n_samples))
      orders_use <- 1:n_samples
    } else {
      orders_use <- ZETA_ORDERS
    }

    sub <- inc_farm[grp_rows, ]
    # Trim OTU columns absent from this farm-group subset before passing to
    # Zeta.decline.mc() — avoids operating on a 21 × 5k sparse matrix.
    sub <- sub[, colSums(sub) > 0, drop = FALSE]
    set.seed(42)
    zeta_results[[farm]][[grp]] <- Zeta.decline.mc(
      data.spec = sub,
      orders    = orders_use,
      sam       = SAM,
      plot      = FALSE
    )
  }
}

saveRDS(zeta_results,
        here("results", "rds", "zeta_ITS1_by_farm.rds"),
        compress = "xz")

cat("Saved: results/rds/zeta_ITS1_by_farm.rds\n")

# ---- 4. Model fit table (power-law vs exponential) -----------------------
# Power-law fit:    log10(zeta) ~ log10(order)  — flat decline, deterministic
# Exponential fit:  log10(zeta) ~ order          — steep decline, stochastic

fit_models <- function(z) {
  ord <- z$zeta.order
  val <- z$zeta.val
  pl  <- lm(log10(val) ~ log10(ord))
  ep  <- lm(log10(val) ~ ord)
  data.frame(
    power_r2    = summary(pl)$r.squared,
    exp_r2      = summary(ep)$r.squared,
    aic_power   = AIC(pl),
    aic_exp     = AIC(ep),
    best_fit    = ifelse(AIC(pl) < AIC(ep), "power_law", "exponential")
  )
}

model_fits <- lapply(FARM_LEVELS, \(farm) {
  lapply(c("bagged_flower", "unbagged_flower"), \(grp) {
    fit_models(zeta_results[[farm]][[grp]]) |>
      mutate(farm = farm, group = grp, .before = 1)
  }) |> bind_rows()
}) |> bind_rows()

write_csv(model_fits,
          here("results", "tables", "11b_ITS1_zeta_by_farm_model_fits.csv"))

cat("Saved: results/tables/11b_ITS1_zeta_by_farm_model_fits.csv\n")
cat("Done.\n")
