# ---- 0. Parameters -------------------------------------------------------
N_PERM      <- 999
ALPHA       <- 0.05   # single significance threshold used in all figures
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
# Canonical metric order (top → bottom in figures): primary metric (CLR + CZM) first,
# then rCLR (sensitivity), then abundance-weighted rarefied, then binary rarefied.
DIST_LEVELS <- c("aitchison_CLR", "aitchison_rCLR",
                 "bray_curtis", "unifrac_w", "jaccard", "unifrac_u")
# Euclidean metrics: group centroids are exact in PCoA space, so tree-level tests
# (section 5) are run for these only.
EUCLID_DISTS <- c("aitchison_CLR", "aitchison_rCLR")

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

# Rarefaction settings saved by 08a; used in figure captions
rarefy    <- attr(dist_list, "rarefy")
RARE_NOTE <- sprintf(
  "BC/UniFrac/Jaccard: rarefied to %d reads, mean of %d rarefactions (vegan::avgdist); n = %d",
  rarefy$depth, rarefy$iterations, nrow(meta_list[["bray_curtis"]]))
AIT_NOTE  <- sprintf("Aitchison (CLR, rCLR): non-rarefied, n = %d",
                     nrow(meta_list[["aitchison_CLR"]]))

# Distances between group centroids.
# wcmdscale() returns PCoA coordinates on all positive-eigenvalue axes; for a
# Euclidean distance these reproduce the original distances exactly, so the mean
# coordinates of a group are its true centroid.
group_centroid_dist <- function(d, group) {
  coords <- vegan::wcmdscale(d)
  sums   <- rowsum(coords, as.character(group))
  dist(sums / as.vector(table(as.character(group))[rownames(sums)]))
}

# ---- 3. Per-farm PERMANOVA -----------------------------------------------
# adonis2(d ~ tree_id + sample_type, by = "terms", permutations = how(blocks = tree_id))
# Bagging was applied to flowers within each tree, so sample_type labels are
# permuted within trees. tree_id in the formula removes between-tree variation
# from the residual (blocking alone only restricts the permutations).
# Tests whether bagged/unbagged differ in community composition within each farm.

cat("Running per-farm PERMANOVA...\n")

results_perm_farm <- lapply(FARM_LEVELS, \(farm) {
  lapply(names(dist_list), \(dname) {
    meta_d    <- meta_list[[dname]]
    farm_idx  <- which(!is.na(meta_d$farm_id) & meta_d$farm_id == farm)
    meta_farm <- droplevels(meta_d[farm_idx, ])
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
    res <- tryCatch(
      adonis2(d_farm ~ tree_id + sample_type,
              data         = meta_farm,
              by           = "terms",
              permutations = how(nperm = N_PERM, blocks = meta_farm$tree_id)),
      error = function(e) {
        message(sprintf("PERMANOVA error %s x %s (%d samples, %d complete trees): %s",
                        farm, dname, nrow(meta_farm), trees_ok, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res)) return(NULL)
    as.data.frame(res) |>
      rownames_to_column("term") |>
      filter(term == "sample_type") |>
      mutate(farm = farm, distance = dname, n_samples = nrow(meta_farm),
             level = "flower")
  }) |> bind_rows()
}) |> bind_rows() |>
  # BH correction across the 7 farms, separately for each distance metric
  group_by(distance) |>
  mutate(p_adj = p.adjust(`Pr(>F)`, method = "BH")) |>
  ungroup()

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
    meta_farm <- droplevels(meta_d[farm_idx, ])
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
    tryCatch({
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
    }, error = function(e) {
      message(sprintf("PERMDISP error %s x %s (%d samples, %d complete trees): %s",
                      farm, dname, nrow(meta_farm), trees_ok, conditionMessage(e)))
      NULL
    })
  }) |> bind_rows()
}) |> bind_rows() |>
  # BH correction across the 7 farms, separately for each distance metric
  group_by(distance) |>
  mutate(p_adj = p.adjust(p_value, method = "BH")) |>
  ungroup()

# ---- 4. Combined PERMANOVA, flower level (all metrics) -------------------
# adonis2(d ~ tree_id + sample_type, permutations within tree).
# Overall bagged/unbagged effect across all farms; farm is nested in tree, so
# tree_id also absorbs between-farm variation. Valid for every metric, including
# the rarefied ones where trees are no longer balanced.

cat("Running combined PERMANOVA (flower level)...\n")

results_flower <- lapply(names(dist_list), \(dname) {
  d      <- dist_list[[dname]]
  meta_d <- meta_list[[dname]]
  set.seed(42)
  adonis2(d ~ tree_id + sample_type,
          data         = meta_d,
          by           = "terms",
          permutations = how(nperm = N_PERM, blocks = meta_d$tree_id)) |>
    as.data.frame() |>
    rownames_to_column("term") |>
    filter(term == "sample_type") |>
    mutate(distance = dname, level = "flower")
}) |> bind_rows()

# ---- 5. Combined PERMANOVA, tree level (Euclidean metrics) ---------------
# The tree is the unit of replication: flowers of one tree are subsamples.
# Farm: 49 tree centroids ~ farm_id, trees permuted freely between farms.
# Treatment and farm x treatment: 98 tree x treatment centroids
#   ~ farm_id + tree_id + sample_type + farm_id:sample_type, permutations within tree.
#   The residual is the tree x treatment variation (df 42), i.e. how much the
#   bagging effect varies between trees of the same farm. Only the sample_type and
#   farm_id:sample_type rows are reported; within-tree permutations do not test
#   farm_id or tree_id.

cat("Running combined PERMANOVA (tree level)...\n")

results_tree <- lapply(EUCLID_DISTS, \(dname) {
  d      <- dist_list[[dname]]
  meta_d <- meta_list[[dname]]

  # Farm effect on tree centroids
  d_tree    <- group_centroid_dist(d, meta_d$tree_id)
  meta_tree <- meta_d |>
    distinct(tree_id, farm_id) |>
    (\(x) x[match(labels(d_tree), x$tree_id), ])()
  set.seed(42)
  res_farm <- adonis2(d_tree ~ farm_id,
                      data         = meta_tree,
                      by           = "terms",
                      permutations = how(nperm = N_PERM))

  # Treatment and farm x treatment on tree x treatment centroids
  meta_cell <- meta_d |>
    mutate(cell_id = paste(tree_id, sample_type, sep = "__"))
  d_cell    <- group_centroid_dist(d, meta_cell$cell_id)
  meta_cell <- meta_cell |>
    distinct(cell_id, tree_id, farm_id, sample_type) |>
    (\(x) x[match(labels(d_cell), x$cell_id), ])()
  set.seed(42)
  res_cell <- adonis2(d_cell ~ farm_id + tree_id + sample_type + farm_id:sample_type,
                      data         = meta_cell,
                      by           = "terms",
                      permutations = how(nperm = N_PERM, blocks = meta_cell$tree_id))

  bind_rows(
    as.data.frame(res_farm) |> rownames_to_column("term") |>
      filter(term == "farm_id") |>
      mutate(level = "tree"),
    as.data.frame(res_cell) |> rownames_to_column("term") |>
      filter(term %in% c("sample_type", "farm_id:sample_type")) |>
      mutate(level = "tree_x_treatment")
  ) |>
    mutate(distance = dname)
}) |> bind_rows()

results_combined <- bind_rows(results_tree, results_flower)

# ---- 6. Global PERMDISP --------------------------------------------------
# Dispersion of each flower around its own farm x treatment centroid, so that
# between-farm differences do not count as "spread". The distances to centroid
# are then compared between bagged and unbagged after accounting for farm,
# with the same within-tree permutations as the PERMANOVA.
# Significant PERMANOVA + non-significant PERMDISP → true location shift.
# Both significant → PERMANOVA result may partly reflect dispersion differences.

cat("Running global PERMDISP...\n")

results_permdisp <- lapply(names(dist_list), \(dname) {
  meta_d <- meta_list[[dname]]
  bd     <- betadisper(dist_list[[dname]],
                       droplevels(interaction(meta_d$farm_id, meta_d$sample_type)))
  d_disp <- dist(bd$distances)   # univariate distance-to-centroid as a distance matrix
  set.seed(42)
  res <- adonis2(d_disp ~ farm_id + sample_type,
                 data         = meta_d,
                 by           = "terms",
                 permutations = how(nperm = N_PERM, blocks = meta_d$tree_id))
  data.frame(
    distance         = dname,
    disp_bagged      = mean(bd$distances[meta_d$sample_type == "bagged_flower"]),
    disp_unbagged    = mean(bd$distances[meta_d$sample_type == "unbagged_flower"]),
    F_value          = res["sample_type", "F"],
    p_value          = res["sample_type", "Pr(>F)"],
    df_group         = res["sample_type", "Df"],
    df_resid         = res["Residual", "Df"]
  )
}) |> bind_rows()

# ---- 7. Visualisation ----------------------------------------------------
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

# Figure A: per-farm significance heatmap (farm × distance; fill = PERMANOVA R²)
# Label = BH-adjusted p; "*" = p_adj < ALPHA.
# Bold black border = PERMDISP p_adj < ALPHA (unequal dispersion between bagged/unbagged)
heatmap_data <- results_perm_farm |>
  left_join(
    results_permdisp_farm |> dplyr::select(farm, distance, permdisp_p_adj = p_adj),
    by = c("farm", "distance")
  ) |>
  mutate(
    # farms left → right in established order; distances top → bottom
    farm         = factor(farm, levels = FARM_LEVELS),
    distance     = factor(distance, levels = rev(DIST_LEVELS)),
    permdisp_sig = permdisp_p_adj < ALPHA,
    p_label      = paste0(sprintf("%.3f", p_adj), if_else(p_adj < ALPHA, " *", ""))
  )

p_farm <- heatmap_data |>
  ggplot(aes(x = farm, y = distance)) +
  # Layer 1: tile fill by PERMANOVA R² (thin grey grid between tiles)
  geom_tile(aes(fill = R2), colour = "grey80", linewidth = 0.3) +
  # Layer 2: bold border on tiles where PERMDISP is significant
  geom_tile(data = filter(heatmap_data, permdisp_sig),
            fill = NA, colour = "black", linewidth = 1.2) +
  # Layer 3: adjusted p-value text; size controls label font size
  geom_text(aes(label = p_label), size = 2.5) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = expression(R^2), limits = c(0, NA)) +
  labs(
    title    = "PERMANOVA per farm — 16S ASV bacteria",
    subtitle = sprintf("Fill = R² of sample_type; label = BH-adjusted p (7 farms per metric); * p_adj < %.2f; bold border = PERMDISP p_adj < %.2f",
                       ALPHA, ALPHA),
    caption  = paste(
      sprintf("adonis2(d ~ tree_id + sample_type, by = \"terms\"), %d permutations within tree: bagging was applied to flowers within each tree.", N_PERM),
      "PERMDISP: betadisper(d, sample_type), permutations within tree.",
      AIT_NOTE, RARE_NOTE, sep = "\n"),
    x = "Farm", y = "Distance matrix"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1),
        # caption: left-aligned, smaller font
        plot.caption = element_text(hjust = 0, size = 7))

ggsave(here("results", "figures", "10a_16S_perm_farm_heatmap.png"),
       plot = p_farm, width = 9, height = 6.5, dpi = 300)

# Figure B: combined model R² bar chart (term × distance; one panel per distance)
# Tree-level terms exist for Euclidean metrics only; flower-level treatment for all.
TERM_LABELS <- c(
  "tree.farm_id"                         = "farm\n(tree)",
  "tree_x_treatment.sample_type"         = "treatment\n(tree x trt)",
  "tree_x_treatment.farm_id:sample_type" = "farm x treatment\n(tree x trt)",
  "flower.sample_type"                   = "treatment\n(flower)"
)

p_comb <- results_combined |>
  mutate(
    distance   = factor(distance, levels = DIST_LEVELS),
    term_label = factor(TERM_LABELS[paste(level, term, sep = ".")],
                        levels = TERM_LABELS),
    sig        = if_else(`Pr(>F)` < ALPHA, "*", "")
  ) |>
  ggplot(aes(x = term_label, y = R2, fill = term_label)) +
  geom_col(width = 0.6) +
  # significance marker above each bar; vjust moves it up from the bar top
  geom_text(aes(label = sig), vjust = -0.3, size = 5) +
  # one panel per distance; free x drops terms that a metric does not have
  facet_wrap(~ distance, nrow = 1, scales = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  # 12% headroom above the tallest bar so the "*" is not clipped
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "PERMANOVA combined models — 16S ASV bacteria",
    subtitle = sprintf("R² per term; * p < %.2f (%d permutations)", ALPHA, N_PERM),
    caption  = paste(
      "farm (tree): 49 tree centroids ~ farm_id; trees permuted freely between farms (tree = unit of replication).",
      "treatment, farm x treatment (tree x trt): 98 tree x treatment centroids ~ farm_id + tree_id + sample_type + farm_id:sample_type;",
      "   permutations within tree; residual = tree x treatment variation. Tree-level tests for Euclidean metrics (CLR, rCLR) only.",
      "treatment (flower): d ~ tree_id + sample_type, permutations within tree; all metrics.",
      AIT_NOTE, RARE_NOTE, sep = "\n"),
    x = NULL, y = expression(R^2)
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x  = element_text(angle = 30, hjust = 1),
        plot.caption = element_text(hjust = 0, size = 7))

ggsave(here("results", "figures", "10a_16S_perm_combined_r2.png"),
       plot = p_comb, width = 15, height = 6, dpi = 300)

# ---- 8. Save outputs -----------------------------------------------------
write_csv(results_perm_farm,
          here("results", "tables", "10a_16S_permanova_per_farm.csv"))
write_csv(results_combined,
          here("results", "tables", "10a_16S_permanova_combined.csv"))
write_csv(results_permdisp,
          here("results", "tables", "10a_16S_permdisp.csv"))
write_csv(results_permdisp_farm,
          here("results", "tables", "10a_16S_permdisp_per_farm.csv"))

cat("\nCombined PERMANOVA:\n")
print(results_combined)
cat("\nPERMDISP summary:\n")
print(results_permdisp)

message("\nDone. CSVs and figures saved.")
