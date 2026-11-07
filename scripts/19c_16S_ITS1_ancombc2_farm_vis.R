## ============================================================================
## 19c — ANCOM-BC2 between farms: heatmap of the taxa that differ
##
## Rows: taxa passing the global farm test (FDR and the pseudo-count
## sensitivity filter) in 19a/19b. Columns: farms. Colour: how far a farm sits
## from the mean of the seven farms, LFC on the natural-log scale (the scale
## ANCOM-BC2 reports). Centring on the farm mean removes the arbitrary
## reference farm (ib) of the model output.
##
## The pairwise test (mdFDR) passes one pair only (ITS1 Cladosporium, mt vs ib),
## so the figure is descriptive: it shows which farm a taxon departs in, not
## which pairs differ significantly.
##
## Input (19a/19b, HPC):
##   results/rds/19a_ancombc2_farm.rds   16S
##   results/rds/19b_ancombc2_farm.rds   ITS1
##   results/rds/ps_16S_bacteria_biosamples.rds, ps_ITS1_otu97_fungi_biosamples.rds  (labels)
## Output:
##   results/tables/19c_ancombc2_farm_lfc.csv
##   results/figures/19c_ancombc2_farm_heatmap.png
## ============================================================================

# ---- 0. Parameters -------------------------------------------------------
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
REF_FARM    <- "ib"   # reference level of farm_id in 19a/19b (lfc = 0 by definition)
HASH_CHARS  <- 6      # characters of the feature hash shown after the label

# Marker label -> ANCOM-BC2 output (19a/19b) and phyloseq object for labels
MARKERS <- list(
  "16S"  = list(ancom = "19a_ancombc2_farm", ps = "ps_16S_bacteria_biosamples"),
  "ITS1" = list(ancom = "19b_ancombc2_farm", ps = "ps_ITS1_otu97_fungi_biosamples")
)

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(phyloseq)
library(patchwork)

# Data version and sample set: RESULTS_DIR, BIOSAMPLE_SUFFIX, SET_TAG
source(here("scripts", "utils_paths.R"))
source(here("scripts", "utils_tax_label.R"))

# ---- 2. Farm deviations of the significant taxa --------------------------
farm_deviation <- function(m) {
  out <- readRDS(file.path(RESULTS_DIR, "rds", paste0(m$ancom, BIOSAMPLE_SUFFIX, ".rds")))
  ps  <- readRDS(file.path(RESULTS_DIR, "rds", paste0(m$ps, BIOSAMPLE_SUFFIX, ".rds")))

  keep <- out$res_global |> filter(diff_robust_abn) |> select(taxon, q_global = q_val)

  out$res |>
    filter(taxon %in% keep$taxon) |>
    select(taxon, starts_with("lfc_farm_id")) |>
    mutate(!!paste0("lfc_farm_id", REF_FARM) := 0) |>   # reference farm
    pivot_longer(-taxon, names_to = "farm_id", names_prefix = "lfc_farm_id", values_to = "lfc") |>
    mutate(deviation = lfc - mean(lfc), .by = taxon) |>  # distance from the mean of 7 farms
    left_join(keep, by = "taxon") |>
    mutate(label   = paste0(make_tax_label(ps, taxon), " [", substr(taxon, 1, HASH_CHARS), "]"),
           farm_id = factor(farm_id, levels = FARM_LEVELS))
}

lfc_long <- map(MARKERS, farm_deviation) |> list_rbind(names_to = "marker")

# ---- 3. Figure -----------------------------------------------------------
plot_marker <- function(mk) {
  df <- filter(lfc_long, marker == mk)

  # Row order: taxa with a similar farm profile next to each other
  # (hierarchical clustering on the deviations; the dendrogram is not drawn)
  wide      <- df |> select(label, farm_id, deviation) |>
    pivot_wider(names_from = farm_id, values_from = deviation) |> column_to_rownames("label")
  row_order <- rownames(wide)[hclust(dist(wide))$order]

  lim <- max(abs(df$deviation))   # symmetric colour limits, per panel
  ggplot(df, aes(x = farm_id, y = factor(label, levels = row_order), fill = deviation)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", deviation)), size = 2.6) +   # size: value labels
    # Diverging palette: red = more abundant than the farm mean, blue = less
    scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-lim, lim),
                         name = "LFC vs\nfarm mean") +
    labs(title = sprintf("%s: %d taxa differ between farms", mk, n_distinct(df$taxon)),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid  = element_blank(),
          axis.text.y = element_text(face = "italic", size = 8),   # size: taxon labels
          plot.title  = element_text(face = "bold", size = 11))
}

# Panel heights follow the number of taxa so the rows have equal height
n_rows  <- lfc_long |> distinct(marker, taxon) |> count(marker) |> pull(n, name = marker)
p_heat  <- wrap_plots(map(names(MARKERS), plot_marker), ncol = 1, heights = n_rows[names(MARKERS)]) +
  plot_annotation(
    caption = paste0(
      "ANCOM-BC2 global test, sample_type + farm_id, (1 | tree_id); taxa passing BH FDR < 0.05 ",
      "and the pseudo-count sensitivity filter.\nLFC on the natural-log scale, centred on the mean ",
      "of the seven farms. Pairwise test (mdFDR): one pair only, ITS1 Cladosporium mt vs ib. ",
      sub("_", "", SET_TAG), " set."),
    theme = theme(plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))
  )

# ---- 4. Save -------------------------------------------------------------
dir.create(file.path(RESULTS_DIR, "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RESULTS_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

lfc_long |>
  arrange(marker, label, farm_id) |>
  write_csv(file.path(RESULTS_DIR, "tables", paste0("19c_ancombc2_farm_lfc", SET_TAG, ".csv")))

ggsave(file.path(RESULTS_DIR, "figures", paste0("19c_ancombc2_farm_heatmap", SET_TAG, ".png")),
       plot = p_heat, width = 8, height = 2 + 0.28 * sum(n_rows), dpi = 300)

print(count(distinct(lfc_long, marker, taxon), marker))
cat("19c done.\n")
