################################################################################
####################### ITS1 Zeta-diversity decline ###########################
################################################################################
# Question: Do flower visitors homogenize or heterogenize the fungal microbiome?
# Prediction:
#   - unbagged (visitor-exposed) ~ power law  → homogenization
#   - unbagged (visitor-exposed) ~ exp law    → heterogenization
#   - parallel curves            → visitors have negligible effect
#
# Input:  results/rds/phyloseq_ITS1_fungi_biosamples.rds
# Output: results/rds/zeta_ITS1_*.rds
#         results/figures/zeta_ITS1_*.png

# ---- Packages ----
library(phyloseq)
library(zetadiv)
library(dplyr)
library(ggplot2)

# ---- Paths ----
project_dir <- here::here()
rds_dir     <- file.path(project_dir, "results", "rds")
fig_dir     <- file.path(project_dir, "results", "figures")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load phyloseq ----
ps <- readRDS(file.path(rds_dir, "phyloseq_ITS1_fungi_biosamples.rds"))

# ---- ASV incidence matrix (samples × ASVs, presence/absence) ----
incidence_ITS1 <-
  ps |>
  otu_table() |>
  as("matrix") |>
  t() |>
  (\(x) { x[x > 0] <- 1; as.data.frame(x) })()

metadata_ITS1 <-
  sample_data(ps) |>
  as("data.frame")

stopifnot(all(rownames(incidence_ITS1) == rownames(metadata_ITS1)))

################################################################################
# 1.  Zeta decline on ASV level (orders 1–10, MC sam = 1 000)
################################################################################

set.seed(42)
zeta_bagged_ITS1 <-
  incidence_ITS1[metadata_ITS1$sample_type == "bagged_flower", ] |>
  Zeta.decline.mc(orders = 1:10, sam = 1000, plot = FALSE)

set.seed(42)
zeta_unbagged_ITS1 <-
  incidence_ITS1[metadata_ITS1$sample_type == "unbagged_flower", ] |>
  Zeta.decline.mc(orders = 1:10, sam = 1000, plot = FALSE)

saveRDS(zeta_bagged_ITS1,   file.path(rds_dir, "zeta_bagged_ITS1_1k.rds"),   compress = "xz")
saveRDS(zeta_unbagged_ITS1, file.path(rds_dir, "zeta_unbagged_ITS1_1k.rds"), compress = "xz")

# ---- Model fit diagnostics (power-law vs exponential) ----
extract_zeta_fits <- function(zeta_result) {
  power_r2   <- summary(zeta_result$zeta.pl)$r.squared
  exp_r2     <- summary(zeta_result$zeta.exp)$r.squared
  power_coef <- coef(zeta_result$zeta.pl)
  exp_coef   <- coef(zeta_result$zeta.exp)
  list(
    power_law = list(
      r2          = power_r2,
      a_intercept = 10^(power_coef[["(Intercept)"]]),
      b_coef      = power_coef[["log10(c(orders))"]]
    ),
    exponential = list(
      r2          = exp_r2,
      a_intercept = 10^(exp_coef[["(Intercept)"]]),
      b_coef      = exp_coef[["c(orders)"]]
    ),
    aic = zeta_result$aic
  )
}

cat("=== Bagged flowers – ITS1 ASV ===\n")
print(extract_zeta_fits(zeta_bagged_ITS1))

cat("=== Unbagged flowers – ITS1 ASV ===\n")
print(extract_zeta_fits(zeta_unbagged_ITS1))

# ---- Plot ----
zeta_tbl_ITS1 <-
  dplyr::bind_rows(
    data.frame(order    = zeta_bagged_ITS1$zeta.order,
               zeta     = zeta_bagged_ITS1$zeta.val,
               zeta_sd  = zeta_bagged_ITS1$zeta.val.sd,
               group    = "bagged_flower"),
    data.frame(order    = zeta_unbagged_ITS1$zeta.order,
               zeta     = zeta_unbagged_ITS1$zeta.val,
               zeta_sd  = zeta_unbagged_ITS1$zeta.val.sd,
               group    = "unbagged_flower")
  )

p_ITS1_asv <-
  ggplot(zeta_tbl_ITS1, aes(x = order, y = zeta, colour = group)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = pmax(zeta - zeta_sd, .Machine$double.eps),
        ymax = zeta + zeta_sd),
    width = 0.15, alpha = 0.9
  ) +
  theme_bw() +
  scale_y_log10() +
  scale_x_continuous(breaks = 1:10) +
  labs(
    x      = "Zeta order (number of flowers)",
    y      = "Zeta diversity component (average shared ASVs)",
    title  = "ITS1 zeta-diversity decline – ASV level",
    colour = "Flower treatment"
  )

ggsave(file.path(fig_dir, "zeta_ITS1_asv_1k.png"),
       p_ITS1_asv, width = 7, height = 5, dpi = 300)

cat("Done. Results saved to results/rds/ and results/figures/\n")
