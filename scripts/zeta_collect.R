################################################################################
# zeta_collect.R — собрать per-order RDS файлы и построить финальные графики
#
# Запускать ПОСЛЕ того как все array tasks завершились:
#   Rscript scripts/zeta_collect.R
#
# или в интерактивной сессии R (локально или на кластере)
################################################################################

library(dplyr)
library(ggplot2)

# ---- Paths ----
project_dir <- here::here()    # локально
# project_dir <- "/nobackup/jcnx71/cacao_flower_microbiome"  # на кластере
rds_dir  <- file.path(project_dir, "results", "rds")
tmp_dir  <- file.path(rds_dir, "zeta_tmp")
fig_dir  <- file.path(project_dir, "results", "figures")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Функция: собрать все orders в таблицу ----
collect_zeta <- function(amplicon, group) {

  pattern <- sprintf("zeta_%s_%s_order*.rds", amplicon, group)
  files   <- sort(list.files(tmp_dir, pattern = pattern, full.names = TRUE))

  if (length(files) == 0) {
    stop(sprintf("No files found for %s %s in %s", amplicon, group, tmp_dir))
  }

  cat(sprintf("Found %d order files for %s %s\n", length(files), amplicon, group))

  rows <- lapply(files, function(f) {
    r <- readRDS(f)
    data.frame(
      order   = r$zeta.order,
      zeta    = r$zeta.val,
      zeta_sd = r$zeta.val.sd
    )
  })

  bind_rows(rows) |> arrange(order)
}

# ---- Функция: fit power-law and exponential, return diagnostics ----
fit_zeta_models <- function(tbl) {
  orders <- tbl$order
  vals   <- tbl$zeta

  pl  <- lm(log10(vals) ~ log10(orders))
  exp <- lm(log10(vals) ~ orders)

  list(
    power_law = list(
      r2          = summary(pl)$r.squared,
      a_intercept = 10^coef(pl)[["(Intercept)"]],
      b_coef      = coef(pl)[["log10(orders)"]]
    ),
    exponential = list(
      r2          = summary(exp)$r.squared,
      a_intercept = 10^coef(exp)[["(Intercept)"]],
      b_coef      = coef(exp)[["orders"]]
    ),
    aic = AIC(pl, exp)
  )
}

################################################################################
# 16S
################################################################################

tbl_16S <-
  bind_rows(
    collect_zeta("16S", "bagged")   |> mutate(group = "bagged_flower"),
    collect_zeta("16S", "unbagged") |> mutate(group = "unbagged_flower")
  )

saveRDS(tbl_16S, file.path(rds_dir, "zeta_16S_results.rds"), compress = "xz")

cat("=== 16S model fits ===\n")
cat("-- Bagged --\n");   print(fit_zeta_models(filter(tbl_16S, group == "bagged_flower")))
cat("-- Unbagged --\n"); print(fit_zeta_models(filter(tbl_16S, group == "unbagged_flower")))

p_16S <-
  ggplot(tbl_16S, aes(x = order, y = zeta, colour = group)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = pmax(zeta - zeta_sd, .Machine$double.eps),
        ymax = zeta + zeta_sd),
    width = 0.15, alpha = 0.9
  ) +
  theme_bw() +
  scale_y_log10() +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    x      = "Zeta order (number of flowers)",
    y      = "Zeta diversity (average shared ASVs)",
    title  = "16S zeta-diversity decline",
    colour = "Flower treatment"
  )

ggsave(file.path(fig_dir, "zeta_16S_10k.png"), p_16S, width = 7, height = 5, dpi = 300)

################################################################################
# ITS1
################################################################################

tbl_ITS1 <-
  bind_rows(
    collect_zeta("ITS1", "bagged")   |> mutate(group = "bagged_flower"),
    collect_zeta("ITS1", "unbagged") |> mutate(group = "unbagged_flower")
  )

saveRDS(tbl_ITS1, file.path(rds_dir, "zeta_ITS1_results.rds"), compress = "xz")

cat("=== ITS1 model fits ===\n")
cat("-- Bagged --\n");   print(fit_zeta_models(filter(tbl_ITS1, group == "bagged_flower")))
cat("-- Unbagged --\n"); print(fit_zeta_models(filter(tbl_ITS1, group == "unbagged_flower")))

p_ITS1 <-
  ggplot(tbl_ITS1, aes(x = order, y = zeta, colour = group)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = pmax(zeta - zeta_sd, .Machine$double.eps),
        ymax = zeta + zeta_sd),
    width = 0.15, alpha = 0.9
  ) +
  theme_bw() +
  scale_y_log10() +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    x      = "Zeta order (number of flowers)",
    y      = "Zeta diversity component (average shared ASVs)",
    title  = "ITS1 zeta-diversity decline",
    colour = "Flower treatment"
  )

ggsave(file.path(fig_dir, "zeta_ITS1_10k.png"), p_ITS1, width = 7, height = 5, dpi = 300)

cat("Done. Figures saved to results/figures/\n")
