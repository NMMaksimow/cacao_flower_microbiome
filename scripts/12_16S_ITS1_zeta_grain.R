# ---- 0. Parameters -------------------------------------------------------
# Zeta orders capped at 7: farm grain has only 7 samples (one per farm),
# so all three grains share the same x-axis range for direct comparison.
ZETA_ORDERS <- 1:7
SAM         <- 1000L

# Grain colour palette; line type distinguishes bagged vs unbagged.
GRAIN_COLS <- c(flower = "#4dac26", tree = "#b8860b", farm = "#7b2d8b")
GRAIN_LTYS <- c(bagged_flower = "solid", unbagged_flower = "dotted")

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(phyloseq)
library(zetadiv)
library(dplyr)
library(tibble)
library(tidyverse)
library(patchwork)

# ---- 2. Load phyloseq objects --------------------------------------------
ps_16s <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds")) |>
  subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))

ps_its <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds")) |>
  subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))

# ---- 3. Aggregation helper -----------------------------------------------
# Union rule: taxon is present at grain G if it appears on any sample within G.
# incidence_df: samples × taxa data.frame (values 0/1)
# meta:         sample metadata with rownames matching incidence_df
# group_col:    column in meta defining the new grain units (tree_id / farm_id)
aggregate_grain <- function(incidence_df, meta, group_col) {
  grp         <- meta[[group_col]]
  unique_grps <- unique(grp)
  do.call(rbind, lapply(unique_grps, \(g) {
    rows <- which(grp == g)
    as.integer(colSums(incidence_df[rows, , drop = FALSE]) > 0)
  })) |>
    as.data.frame() |>
    setNames(colnames(incidence_df)) |>
    (\(df) { rownames(df) <- unique_grps; df })()
}

# ---- 4. Build incidence at 3 grains for one phyloseq object --------------
# Returns nested list: [[treatment]][[grain]] = trimmed incidence data.frame
build_grain_combos <- function(ps, marker_label) {
  lapply(c("bagged_flower", "unbagged_flower"), \(trt) {
    # prune_samples() takes a logical vector directly — avoids NSE scope issue
    # that causes subset_samples() to fail inside lapply lambdas.
    ps_t <- prune_samples(sample_data(ps)$sample_type == trt, ps)

    inc  <- otu_table(ps_t) |>
      as("matrix") |>
      (\(m) if (taxa_are_rows(ps_t)) t(m) else m)() |>
      (\(m) { m[m > 0] <- 1L; as.data.frame(m) })()

    meta <- as(sample_data(ps_t), "data.frame")
    stopifnot(all(rownames(inc) == rownames(meta)))

    inc_flower <- inc[, colSums(inc) > 0, drop = FALSE]
    inc_tree   <- aggregate_grain(inc, meta, "tree_id")
    inc_tree   <- inc_tree[, colSums(inc_tree) > 0, drop = FALSE]
    inc_farm   <- aggregate_grain(inc, meta, "farm_id")
    inc_farm   <- inc_farm[, colSums(inc_farm) > 0, drop = FALSE]

    cat(sprintf("[%s | %s]  flower: %d × %d  tree: %d × %d  farm: %d × %d\n",
                marker_label, trt,
                nrow(inc_flower), ncol(inc_flower),
                nrow(inc_tree),   ncol(inc_tree),
                nrow(inc_farm),   ncol(inc_farm)))

    list(flower = inc_flower, tree = inc_tree, farm = inc_farm)
  }) |>
    setNames(c("bagged_flower", "unbagged_flower"))
}

# ---- 5. Run Zeta.decline.mc for every grain × treatment combination ------
run_zeta_combos <- function(combo_list, marker_label) {
  # Reproducible seeds: one per (treatment, grain) combo
  seeds <- list(
    bagged_flower   = c(flower = 101L, tree = 201L, farm = 301L),
    unbagged_flower = c(flower = 102L, tree = 202L, farm = 302L)
  )
  lapply(names(combo_list), \(trt) {
    lapply(names(combo_list[[trt]]), \(grain) {
      n <- nrow(combo_list[[trt]][[grain]])
      cat(sprintf("  zeta [%s | %s | %s]  n = %d\n", marker_label, trt, grain, n))
      set.seed(seeds[[trt]][[grain]])
      z <- Zeta.decline.mc(
        data.spec = combo_list[[trt]][[grain]],
        orders    = ZETA_ORDERS,
        sam       = SAM,
        plot      = FALSE
      )
      tibble(
        treatment = trt,
        grain     = grain,
        order     = z$zeta.order,
        zeta      = z$zeta.val,
        zeta_sd   = z$zeta.val.sd
      )
    }) |> bind_rows()
  }) |> bind_rows()
}

# ---- 6. Compute zeta for both markers ------------------------------------
cat("=== 16S ===\n")
combos_16s <- build_grain_combos(ps_16s, "16S")
tbl_16s    <- run_zeta_combos(combos_16s, "16S") |>
  mutate(
    grain     = factor(grain,     levels = c("flower", "tree", "farm")),
    treatment = factor(treatment, levels = c("bagged_flower", "unbagged_flower"))
  )

cat("\n=== ITS1 ===\n")
combos_its <- build_grain_combos(ps_its, "ITS1")
tbl_its    <- run_zeta_combos(combos_its, "ITS1") |>
  mutate(
    grain     = factor(grain,     levels = c("flower", "tree", "farm")),
    treatment = factor(treatment, levels = c("bagged_flower", "unbagged_flower"))
  )

# ---- 7. Build one-panel plot for a single marker -------------------------
make_grain_panel <- function(tbl, y_label, show_strip = TRUE) {
  ggplot(
    tbl,
    aes(
      x      = order,
      y      = zeta,
      colour = grain,
      fill   = grain,
      # group keeps each grain × treatment combo on its own ribbon/line
      group  = interaction(grain, treatment)
    )
  ) +
    # ±1 SD ribbon per grain; no border, semi-transparent
    geom_ribbon(
      aes(ymin = pmax(zeta - zeta_sd, 1e-3), ymax = zeta + zeta_sd),
      alpha  = 0.15,
      colour = NA
    ) +
    # Solid/dashed lines distinguish bagged vs unbagged
    geom_line(aes(linetype = treatment), linewidth = 0.9) +
    geom_point(aes(shape = treatment), size = 1.5) +
    # Log x: power-law appears as a straight line (diagnostic)
    scale_x_log10(
      breaks = c(1, 2, 3, 5, 7),
      labels = c("1", "2", "3", "5", "7")
    ) +
    # Log y: shared taxa span orders of magnitude
    scale_y_log10(
      breaks = scales::trans_breaks("log10", \(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    # Floor y at ~3 shared taxa; ceiling auto-scales to data
    coord_cartesian(ylim = c(10^0.5, NA)) +
    # Grain colours
    scale_colour_manual(
      name   = "Grain",
      values = GRAIN_COLS,
      labels = c("Flower", "Tree", "Farm")
    ) +
    scale_fill_manual(
      values = GRAIN_COLS,
      labels = c("Flower", "Tree", "Farm"),
      guide  = "none"    # fill shares meaning with colour; suppress duplicate
    ) +
    # Treatment line types
    scale_linetype_manual(
      name   = "Treatment",
      values = GRAIN_LTYS,
      labels = c("bagged", "unbagged")
    ) +
    scale_shape_manual(
      name   = "Treatment",
      values = c(bagged_flower = 16L, unbagged_flower = 1L),
      labels = c("bagged", "unbagged")
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position  = "right",
      legend.key.size  = unit(1.4, "lines"),
      legend.text      = element_text(size = 9),
      legend.title     = element_text(size = 9, face = "bold"),
      axis.text        = element_text(size = 7),
      axis.title       = element_text(size = 8.5),
      panel.border     = element_rect(colour = "grey70", linewidth = 0.4),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Zeta order", y = y_label)
}

p_16s <- make_grain_panel(tbl_16s,
                          y_label = "16S bacteria (ASV)\nShared ASVs (log scale)")
p_its <- make_grain_panel(tbl_its,
                          y_label = "ITS1 OTU97 fungi\nShared OTUs (log scale)")

# ---- 8. Assemble and save ------------------------------------------------
# Stack panels vertically; legends are identical so patchwork collects them.
p_final <- (p_16s / p_its) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.key.size = unit(1.4, "lines"),
    legend.text     = element_text(size = 9),
    legend.title    = element_text(size = 9, face = "bold")
  )

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

# Width 10 in: no farm facets, figure is narrower than 11c
ggsave(
  here("results", "figures", "12_16S_ITS1_zeta_grain.png"),
  plot   = p_final,
  width  = 10,
  height = 7,
  dpi    = 300
)

message("Done. Figure saved to results/figures/12_16S_ITS1_zeta_grain.png")
