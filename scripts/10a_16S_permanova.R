# ---- 0. Parameters -------------------------------------------------------
N_PERM      <- 999
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
# Canonical metric order (top → bottom in figures): non-rarefied Aitchison first,
# then abundance-weighted rarefied, then binary rarefied.
DIST_LEVELS <- c("aitchison_rCLR", "aitchison_CLR",
                 "bray_curtis", "unifrac_w", "jaccard", "unifrac_u")

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(vegan)
library(permute)   # how() for blocked permutations; preserves term names in adonis2

# ---- 2. Load distance matrices -------------------------------------------
# Computed and saved by 08a_16S_pcoa_global.R.
# Re-run 08a if these files are missing or outdated.
dist_list <- readRDS(here("results", "rds", "16S_dist_list.rds"))
meta_list <- readRDS(here("results", "rds", "16S_meta_list.rds"))

# ---- 3. Per-farm PERMANOVA -----------------------------------------------
# adonis2(d ~ sample_type, by = "terms", permutations = how(blocks = tree_id))
# blocked permutations: shuffles stay within trees (3 bagged + 3 unbagged per tree).
# Tests whether bagged/unbagged differ in community composition within each farm.

cat("Running per-farm PERMANOVA...\n")

results_perm_farm <- lapply(FARM_LEVELS, \(farm) {
  lapply(names(dist_list), \(dname) {
    meta_d    <- meta_list[[dname]]
    farm_idx  <- which(!is.na(meta_d$farm_id) & meta_d$farm_id == farm)
    meta_farm <- meta_d[farm_idx, ]
    # Require >= 2 trees with both sample types for blocked permutations to be feasible.
    # Rarefaction may leave some farms with 0 or 1 complete trees → skip.
    trees_ok <- meta_farm |>
      group_by(tree_id) |>
      summarise(n_types = n_distinct(sample_type), .groups = "drop") |>
      filter(n_types == 2) |>
      nrow()
    if (trees_ok < 2) {
      message(sprintf("skip %s x %s: only %d trees with both sample types",
                      farm, dname, trees_ok))
      return(NULL)
    }
    d_mat  <- as.matrix(dist_list[[dname]])
    d_farm <- as.dist(d_mat[farm_idx, farm_idx])
    set.seed(42)
    res <- adonis2(d_farm ~ sample_type,
                   data         = meta_farm,
                   by           = "terms",
                   permutations = how(nperm = N_PERM, blocks = meta_farm$tree_id))
    as.data.frame(res) |>
      rownames_to_column("term") |>
      filter(term == "sample_type") |>
      mutate(farm = farm, distance = dname, model = "per_farm")
  }) |> bind_rows()
}) |> bind_rows()

# ---- 3b. Per-farm PERMDISP -----------------------------------------------
# Multivariate homoscedasticity test within each farm.
# Tests whether bagged/unbagged differ in *spread* around their centroid.
# If sig: PERMANOVA result for that farm may partly reflect dispersion, not
# centroid shift. Compare to PERMANOVA result when interpreting per-farm heatmap.

cat("Running per-farm PERMDISP...\n")

results_permdisp_farm <- lapply(FARM_LEVELS, \(farm) {
  lapply(names(dist_list), \(dname) {
    meta_d    <- meta_list[[dname]]
    farm_idx  <- which(!is.na(meta_d$farm_id) & meta_d$farm_id == farm)
    meta_farm <- meta_d[farm_idx, ]
    # Require >= 2 trees with both sample types for blocked permutations to be feasible.
    # Rarefaction may leave some farms with 0 or 1 complete trees → skip.
    trees_ok <- meta_farm |>
      group_by(tree_id) |>
      summarise(n_types = n_distinct(sample_type), .groups = "drop") |>
      filter(n_types == 2) |>
      nrow()
    if (trees_ok < 2) {
      message(sprintf("skip %s x %s: only %d trees with both sample types",
                      farm, dname, trees_ok))
      return(NULL)
    }
    d_mat  <- as.matrix(dist_list[[dname]])
    d_farm <- as.dist(d_mat[farm_idx, farm_idx])
    bd    <- betadisper(d_farm, meta_farm$sample_type)
    ptest <- permutest(bd, permutations = how(nperm = N_PERM,
                                              blocks = meta_farm$tree_id))
    data.frame(
      farm     = farm,
      distance = dname,
      F_value  = ptest$tab[["F"]][1],
      p_value  = ptest$tab[["Pr(>F)"]][1],
      df_group = ptest$tab[["Df"]][1],
      df_resid = ptest$tab[["Df"]][2]
    )
  }) |> bind_rows()
}) |> bind_rows()

# ---- 4. Combined PERMANOVA -----------------------------------------------
# Model A (additive): farm_id controls between-farm variance; tests overall
#   sample_type effect after accounting for farm.
# Model B (interaction): tests whether the sample_type effect is heterogeneous
#   across farms (the variation visible in the per-farm PCA panels).

cat("Running combined PERMANOVA...\n")

results_combined <- lapply(names(dist_list), \(dname) {
  d      <- dist_list[[dname]]
  meta_d <- meta_list[[dname]]

  set.seed(42)
  res_add <- adonis2(d ~ farm_id + sample_type,
                     data         = meta_d,
                     by           = "terms",
                     permutations = how(nperm = N_PERM, blocks = meta_d$tree_id))

  set.seed(42)
  res_int <- adonis2(d ~ farm_id * sample_type,
                     data         = meta_d,
                     by           = "terms",
                     permutations = how(nperm = N_PERM, blocks = meta_d$tree_id))

  bind_rows(
    as.data.frame(res_add) |> rownames_to_column("term") |>
      mutate(distance = dname, model = "farm_add_treatment"),
    as.data.frame(res_int) |> rownames_to_column("term") |>
      mutate(distance = dname, model = "farm_x_treatment")
  )
}) |> bind_rows()

# ---- 5. PERMDISP ---------------------------------------------------------
# Multivariate homoscedasticity test (betadisper).
# Significant PERMANOVA + non-significant PERMDISP → true location shift.
# Both significant → PERMANOVA result may partly reflect dispersion differences.

cat("Running PERMDISP...\n")

results_permdisp <- lapply(names(dist_list), \(dname) {
  meta_d <- meta_list[[dname]]
  bd    <- betadisper(dist_list[[dname]], meta_d$sample_type)
  ptest <- permutest(bd, permutations = N_PERM)
  data.frame(
    distance = dname,
    F_value  = ptest$tab[["F"]][1],
    p_value  = ptest$tab[["Pr(>F)"]][1],
    df_group = ptest$tab[["Df"]][1],
    df_resid = ptest$tab[["Df"]][2]
  )
}) |> bind_rows()

# ---- 6. Visualisation ----------------------------------------------------
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

# Figure A: per-farm significance heatmap (farm × distance; fill = PERMANOVA R²)
# Bold black border = PERMDISP p < 0.05 (unequal dispersion between bagged/unbagged)
heatmap_data <- results_perm_farm |>
  left_join(
    results_permdisp_farm |> dplyr::select(farm, distance, permdisp_p = p_value),
    by = c("farm", "distance")
  ) |>
  mutate(
    # farms left → right in established order; distances top → bottom
    farm         = factor(farm, levels = FARM_LEVELS),
    distance     = factor(distance, levels = rev(DIST_LEVELS)),
    permdisp_sig = permdisp_p < 0.05,
    p_label      = case_when(
      `Pr(>F)` < 0.001 ~ "<0.001 ***",
      `Pr(>F)` < 0.01  ~ paste0(sprintf("%.3f", `Pr(>F)`), " **"),
      `Pr(>F)` < 0.05  ~ paste0(sprintf("%.3f", `Pr(>F)`), " *"),
      TRUE             ~ sprintf("%.3f", `Pr(>F)`)
    )
  )

p_farm <- heatmap_data |>
  ggplot(aes(x = farm, y = distance)) +
  # Layer 1: tile fill by PERMANOVA R² (thin grey grid between tiles)
  geom_tile(aes(fill = R2), colour = "grey80", linewidth = 0.3) +
  # Layer 2: bold border on tiles where PERMDISP is significant (p < 0.05)
  geom_tile(data = filter(heatmap_data, permdisp_sig),
            fill = NA, colour = "black", linewidth = 1.2) +
  geom_text(aes(label = p_label), size = 2.5) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = expression(R^2), limits = c(0, NA)) +
  labs(
    title    = "PERMANOVA per farm — 16S ASV bacteria",
    subtitle = sprintf("adonis2(d ~ sample_type, blocks = tree_id, nperm = %d)\nFill = R²; label = PERMANOVA p-value; bold border = PERMDISP p < 0.05",
                       N_PERM),
    x = "Farm", y = "Distance matrix"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here("results", "figures", "10a_16S_perm_farm_heatmap.png"),
       plot = p_farm, width = 9, height = 6, dpi = 300)

# Figure B: combined model R² bar chart (term × distance; facet = model × distance)
p_comb <- results_combined |>
  filter(!term %in% c("Residual", "Total")) |>
  mutate(
    distance = factor(distance, levels = DIST_LEVELS),
    model    = recode(model,
                      "farm_add_treatment" = "farm + treatment",
                      "farm_x_treatment"   = "farm × treatment"),
    sig      = case_when(
      `Pr(>F)` < 0.001 ~ "***",
      `Pr(>F)` < 0.01  ~ "**",
      `Pr(>F)` < 0.05  ~ "*",
      TRUE             ~ ""
    )
  ) |>
  ggplot(aes(x = term, y = R2, fill = term)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sig), vjust = -0.3, size = 4) +
  facet_grid(model ~ distance) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "PERMANOVA combined models — 16S ASV bacteria",
    subtitle = sprintf("R² per term (blocked permutations, tree_id, %d permutations)", N_PERM),
    x = NULL, y = expression(R^2)
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here("results", "figures", "10a_16S_perm_combined_r2.png"),
       plot = p_comb, width = 15, height = 5, dpi = 300)

# ---- 7. Save outputs -----------------------------------------------------
write_csv(results_perm_farm,
          here("results", "tables", "10a_16S_permanova_per_farm.csv"))
write_csv(results_combined,
          here("results", "tables", "10a_16S_permanova_combined.csv"))
write_csv(results_permdisp,
          here("results", "tables", "10a_16S_permdisp.csv"))
write_csv(results_permdisp_farm,
          here("results", "tables", "10a_16S_permdisp_per_farm.csv"))

cat("\nPERMDISP summary:\n")
print(results_permdisp)

message("\nDone. CSVs and figures saved.")
