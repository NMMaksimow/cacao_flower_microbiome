## ============================================================================
## 10e — Pairwise farm comparison: PERMANOVA and PERMDISP, 16S and ITS1
##
## Which farms differ in community composition, pair by pair, and do farms
## differ in how variable their trees are?
##
## Unit of replication: the tree, as in the farm test of 10a/10b (section 5).
## Each tree is its centroid in Aitchison (CLR) space over all its flowers;
## every tree holds bagged and unbagged flowers, so the centroid averages over
## treatment. 49 trees, 7 per farm.
##   PERMANOVA:  for each of the 21 farm pairs, 14 tree centroids ~ farm_id,
##               trees permuted freely between the two farms; BH over the 21 pairs.
##   PERMDISP:   distance of each tree centroid to its farm centroid
##               (betadisper), global and pairwise permutation tests; BH.
##
## Input (08a/08b):
##   results/rds/16S_dist_list.rds,        16S_meta_list.rds
##   results/rds/ITS1_otu97_dist_list.rds, ITS1_otu97_meta_list.rds
## Output:
##   results/tables/10e_farm_pairwise_permanova.csv
##   results/tables/10e_farm_permdisp.csv
##   results/figures/10e_farm_pairwise_permanova.png
## ============================================================================

# ---- 0. Parameters -------------------------------------------------------
N_PERM      <- 999
ALPHA       <- 0.05     # BH-adjusted threshold for the stars in the figure
DIST_NAME   <- "aitchison_CLR"   # primary metric (CLR + CZM zero replacement)
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

# Marker label -> RDS prefix written by 08a/08b
MARKERS <- c("16S" = "16S", "ITS1" = "ITS1_otu97")

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(vegan)
library(permute)     # how() for the permutation design
library(patchwork)

# Data version and sample set: RESULTS_DIR, BIOSAMPLE_SUFFIX, SET_TAG
source(here("scripts", "utils_paths.R"))

# ---- 2. Tree centroids ---------------------------------------------------
# Same as group_centroid_dist() in 10a: for a Euclidean distance the PCoA
# coordinates on all positive axes reproduce the distances exactly, so the mean
# coordinates of a group are its true centroid.
group_centroid_dist <- function(d, group) {
  coords <- vegan::wcmdscale(d)
  sums   <- rowsum(coords, as.character(group))
  dist(sums / as.vector(table(as.character(group))[rownames(sums)]))
}

tree_level <- function(prefix) {
  d    <- readRDS(file.path(RESULTS_DIR, "rds", paste0(prefix, "_dist_list", BIOSAMPLE_SUFFIX, ".rds")))[[DIST_NAME]]
  meta <- readRDS(file.path(RESULTS_DIR, "rds", paste0(prefix, "_meta_list", BIOSAMPLE_SUFFIX, ".rds")))[[DIST_NAME]]
  d_tree <- group_centroid_dist(d, meta$tree_id)
  meta_tree <- meta |>
    distinct(tree_id, farm_id) |>
    (\(x) x[match(labels(d_tree), x$tree_id), ])() |>
    mutate(farm_id = factor(farm_id, levels = FARM_LEVELS))
  list(d = d_tree, meta = meta_tree)
}

trees <- map(MARKERS, tree_level)
walk2(trees, names(trees), \(t, m) cat(sprintf("%s: %d tree centroids\n", m, nrow(t$meta))))

# ---- 3. Pairwise PERMANOVA -----------------------------------------------
pairwise_permanova <- function(t) {
  combn(FARM_LEVELS, 2, simplify = FALSE) |>
    map(\(pair) {
      keep <- t$meta$farm_id %in% pair
      set.seed(42)
      res <- adonis2(as.dist(as.matrix(t$d)[keep, keep]) ~ farm_id,
                     data         = droplevels(t$meta[keep, ]),
                     permutations = how(nperm = N_PERM))
      tibble(farm_a = pair[1], farm_b = pair[2], n_trees = sum(keep),
             F = res$F[1], R2 = res$R2[1], p = res$`Pr(>F)`[1])
    }) |>
    list_rbind() |>
    mutate(q = p.adjust(p, method = "BH"))
}

permanova_pairs <- map(trees, pairwise_permanova) |> list_rbind(names_to = "marker")

# Global farm test on the same tree centroids, for the figure subtitle
permanova_global <- map(trees, \(t) {
  set.seed(42)
  res <- adonis2(t$d ~ farm_id, data = t$meta, permutations = how(nperm = N_PERM))
  tibble(R2 = res$R2[1], p = res$`Pr(>F)`[1])
}) |> list_rbind(names_to = "marker")

# ---- 4. PERMDISP between farms -------------------------------------------
farm_dispersion <- function(t) {
  bd <- betadisper(t$d, t$meta$farm_id)
  set.seed(42)
  pt <- permutest(bd, pairwise = TRUE, permutations = how(nperm = N_PERM))
  list(
    distances = tibble(tree_id = labels(t$d), farm_id = t$meta$farm_id, distance = bd$distances),
    tests = bind_rows(
      tibble(comparison = "global", F = pt$tab["Groups", "F"], p = pt$tab["Groups", "Pr(>F)"]),
      tibble(comparison = names(pt$pairwise$permuted), F = NA_real_, p = pt$pairwise$permuted) |>
        mutate(q = p.adjust(p, method = "BH"))
    )
  )
}

dispersion <- map(trees, farm_dispersion)
permdisp_tests <- map(dispersion, "tests") |> list_rbind(names_to = "marker")

# ---- 5. Figure -----------------------------------------------------------
# One row per marker: left, pairwise R2 between farms (lower triangle);
# right, spread of the tree centroids within each farm.

plot_pairs <- function(m) {
  df <- filter(permanova_pairs, marker == m) |>
    mutate(label = sprintf("%.2f%s", R2, if_else(q < ALPHA, "*", "")))
  g  <- filter(permanova_global, marker == m)
  ggplot(df, aes(x = farm_a, y = farm_b, fill = R2)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 3.2) +           # size: R2 label text
    scale_x_discrete(limits = FARM_LEVELS[-length(FARM_LEVELS)]) +
    scale_y_discrete(limits = rev(FARM_LEVELS[-1])) +     # lower triangle, ib first
    # Fill: sequential palette; limits shared by both markers so the colours compare
    scale_fill_distiller(palette = "YlOrRd", direction = 1, name = expression(R^2),
                         limits = c(0, max(permanova_pairs$R2))) +
    coord_equal() +
    labs(title    = sprintf("%s: pairwise PERMANOVA between farms", m),
         subtitle = sprintf("Global farm test: R2 %.2f, p %.3f; * BH q < %.2f (21 pairs)",
                            g$R2, g$p, ALPHA),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.subtitle = element_text(size = 9, colour = "grey40"))
}

plot_dispersion <- function(m) {
  g <- filter(permdisp_tests, marker == m, comparison == "global")
  ggplot(dispersion[[m]]$distances, aes(x = farm_id, y = distance)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", width = 0.6) +
    geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.7) +  # one point per tree
    labs(title    = sprintf("%s: spread of trees within farms (PERMDISP)", m),
         subtitle = sprintf("Global test: F %.2f, p %.3f", g$F, g$p),
         x = NULL, y = "Distance of tree centroid\nto farm centroid (Aitchison)") +
    theme_bw(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9, colour = "grey40"))
}

p_farms <- map(names(MARKERS), \(m) plot_pairs(m) | plot_dispersion(m)) |>
  wrap_plots(ncol = 1) +
  plot_annotation(
    caption = sprintf(paste0(
      "Unit: tree centroid in Aitchison (CLR) space, 7 trees per farm; ",
      "PERMANOVA %d permutations, trees permuted freely between the two farms. %s set."),
      N_PERM, sub("_", "", SET_TAG)),
    theme = theme(plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))
  )

# ---- 6. Save -------------------------------------------------------------
dir.create(file.path(RESULTS_DIR, "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RESULTS_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

write_csv(permanova_pairs,
          file.path(RESULTS_DIR, "tables", paste0("10e_farm_pairwise_permanova", SET_TAG, ".csv")))
write_csv(permdisp_tests,
          file.path(RESULTS_DIR, "tables", paste0("10e_farm_permdisp", SET_TAG, ".csv")))
ggsave(file.path(RESULTS_DIR, "figures", paste0("10e_farm_pairwise_permanova", SET_TAG, ".png")),
       plot = p_farms, width = 12, height = 10, dpi = 300)

permanova_pairs |>
  summarise(pairs_q05 = sum(q < ALPHA), median_R2 = median(R2), .by = marker) |>
  print()
cat("10e done.\n")
