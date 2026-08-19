# ---- 0. Parameters -------------------------------------------------------
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
MGMT        <- c(ib = "inside_forest", vr = "near_forest",  sa = "near_forest",
                 kk = "agroforest",   mt = "agroforest",   vi = "full_sun",
                 yb = "full_sun")
# Colours match Python notebooks and 09c
COLS <- c(bagged_flower = "#4C72B0", unbagged_flower = "#DD8452")

# Strip labels for the top row (16S): "farm | management"
# ITS1 row strips are suppressed (same column = same farm info)
FARM_LABELS <- setNames(paste0(FARM_LEVELS, "\n", MGMT[FARM_LEVELS]),
                        FARM_LEVELS)

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(patchwork)

# ---- 2. Load RDS results -------------------------------------------------
zeta_16s <- readRDS(here("results", "rds", "zeta_16S_by_farm.rds"))
zeta_its <- readRDS(here("results", "rds", "zeta_ITS1_by_farm.rds"))

# ---- 3. Flatten to tidy tibbles ------------------------------------------
# zeta_results[[farm]][[group]] stores a zetadiv object with:
#   $zeta.order — integer vector of orders
#   $zeta.val   — mean shared taxa at each order
#   $zeta.val.sd — standard deviation across MC replicates

flatten_zeta <- function(zeta_list, farm_levels) {
  lapply(farm_levels, \(farm) {
    lapply(c("bagged_flower", "unbagged_flower"), \(grp) {
      z <- zeta_list[[farm]][[grp]]
      tibble(
        farm_id  = farm,
        group    = grp,
        order    = z$zeta.order,
        zeta     = z$zeta.val,
        zeta_sd  = z$zeta.val.sd
      )
    }) |> bind_rows()
  }) |> bind_rows() |>
    mutate(
      farm_id = factor(farm_id, levels = farm_levels),
      group   = factor(group, levels = c("bagged_flower", "unbagged_flower"))
    )
}

tbl_16s <- flatten_zeta(zeta_16s, FARM_LEVELS)
tbl_its <- flatten_zeta(zeta_its, FARM_LEVELS)

# ---- 4. Model fits (power-law vs exponential) ----------------------------
# Power-law:   log10(zeta) ~ log10(order) — flat decline, deterministic filtering
# Exponential: log10(zeta) ~ order        — steep decline, stochastic assembly
# AIC winner characterises dominant assembly mechanism per farm × treatment.

fit_models <- function(z) {
  ord  <- z$zeta.order
  val  <- z$zeta.val
  # Drop orders where zeta = 0: log10(0) = -Inf breaks lm.
  keep <- val > 0
  ord  <- ord[keep]
  val  <- val[keep]
  if (length(val) < 3)
    return(data.frame(power_r2 = NA, exp_r2 = NA,
                      aic_power = NA, aic_exp = NA, best_fit = NA))
  pl <- lm(log10(val) ~ log10(ord))
  ep <- lm(log10(val) ~ ord)
  data.frame(
    power_r2  = summary(pl)$r.squared,
    exp_r2    = summary(ep)$r.squared,
    aic_power = AIC(pl),
    aic_exp   = AIC(ep),
    best_fit  = ifelse(AIC(pl) < AIC(ep), "power_law", "exponential")
  )
}

run_model_fits <- function(zeta_list, farm_levels) {
  do.call(rbind, lapply(farm_levels, \(farm) {
    do.call(rbind, lapply(c("bagged_flower", "unbagged_flower"), \(grp) {
      row <- fit_models(zeta_list[[farm]][[grp]])
      cbind(data.frame(farm = farm, group = grp), row)
    }))
  }))
}

dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)

fits_16s <- run_model_fits(zeta_16s, FARM_LEVELS)
write.csv(fits_16s, here("results", "tables", "11c_16S_zeta_model_fits.csv"),
          row.names = FALSE)

fits_its <- run_model_fits(zeta_its, FARM_LEVELS)
write.csv(fits_its, here("results", "tables", "11c_ITS1_zeta_model_fits.csv"),
          row.names = FALSE)

message("Model fit tables saved to results/tables/")

# ---- 5. Helper: build one facet_wrap row ---------------------------------
make_row <- function(tbl, y_label, show_strip = TRUE) {
  ggplot(tbl, aes(x = order, y = zeta, colour = group, fill = group)) +
    # ±1 SD ribbon; floor at 1e-3 to stay positive on log scale
    geom_ribbon(
      aes(ymin = pmax(zeta - zeta_sd, 1e-3), ymax = zeta + zeta_sd),
      alpha = 0.20,
      colour = NA
    ) +
    # Mean zeta decline line
    geom_line(linewidth = 0.85) +
    # Individual order points
    geom_point(size = 1.2) +
    facet_wrap(~ farm_id, nrow = 1,
               labeller = as_labeller(FARM_LABELS)) +
    # Log10 y: shared taxa span orders of magnitude across zeta orders
    scale_y_log10(
      # Minor breaks give grid lines at log10 sub-divisions
      breaks = scales::trans_breaks("log10", \(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    # Log x: power-law decline appears as a straight line on log-log (diagnostic).
    scale_x_log10(
      breaks = c(1, 2, 5, 10, 21),
      labels = c("1", "2", "5", "10", "21")
    ) +
    scale_colour_manual(
      values = COLS,
      labels = c("bagged", "unbagged"),
      name   = NULL
    ) +
    scale_fill_manual(
      values = COLS,
      labels = c("bagged", "unbagged"),
      name   = NULL
    ) +
    theme_bw(base_size = 8) +
    theme(
      # Strip text: bold for top row only; hide for bottom row
      strip.text       = if (show_strip)
        element_text(size = 7, face = "bold", margin = margin(b = 3))
      else
        element_blank(),
      strip.background = element_blank(),
      # Legend handled by patchwork collect at bottom
      legend.position  = "none",
      # Tighter axis labels to conserve vertical space
      axis.text        = element_text(size = 6.5),
      axis.title       = element_text(size = 7.5),
      # Panels separated by a thin border
      panel.border     = element_rect(colour = "grey70", linewidth = 0.4),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Zeta order",
      y = y_label
    )
}

p_16s <- make_row(tbl_16s,
                  y_label    = "16S bacteria (ASV)\nShared ASVs (log scale)",
                  show_strip = TRUE) +
  # Shared y limits with ITS1 row (10^0.5 to 10^3) for cross-kingdom comparison.
  # coord_cartesian clips without removing data from ribbon/scale fitting.
  coord_cartesian(ylim = c(10^0.5, 10^3))

p_its <- make_row(tbl_its,
                  y_label    = "ITS1 OTU97 fungi\nShared OTUs (log scale)",
                  show_strip = FALSE) +
  coord_cartesian(ylim = c(10^0.5, 10^3))

# ---- 6. Assemble with patchwork ------------------------------------------
# (p_16s / p_its): stacks two rows vertically
# guides = "collect": pulls all identical legends into one shared legend
# & theme(legend.position = "bottom"): applies to both assembled rows
p_final <- (p_16s / p_its) +
  plot_layout(guides = "collect") &
  theme(
    legend.position  = "bottom",
    # Legend key size: adjust width/height (in lines) to control swatch size
    legend.key.size  = unit(0.8, "lines"),
    legend.text      = element_text(size = 8),
    legend.margin    = margin(t = 2)
  )

# ---- 7. Save figure ------------------------------------------------------
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

# Width 20 in × 7 farms works well; height 7 in gives each row ~3.5 in.
# Increase height (e.g. 8) if y-axis labels or tick numbers feel crowded.
ggsave(
  here("results", "figures", "11c_16S_ITS1_zeta_by_farm.png"),
  plot   = p_final,
  width  = 20,
  height = 7,
  dpi    = 300
)

message("Done. Figure saved to results/figures/11c_16S_ITS1_zeta_by_farm.png")
