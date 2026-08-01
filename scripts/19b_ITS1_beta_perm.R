# ---- 0. Parameters -------------------------------------------------------
PREV_THRESHOLD <- 0.10   # global prevalence filter: ≥10% of biosamples = ≥29 of 294
N_PERM         <- 999
FARM_LEVELS    <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
RAREFY_DEPTH   <- 15000L  # rarefaction depth for depth-sensitive metrics (BC, Jaccard)
# Canonical metric order (top → bottom in figures): non-rarefied Aitchison first,
# then rarefied abundance-weighted, then rarefied binary.
DIST_LEVELS    <- c("aitchison_rCLR", "aitchison_CLR", "bray_curtis", "jaccard")

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(phyloseq)
library(vegan)
library(permute)       # how() for blocked permutations; preserves term names in adonis2
library(zCompositions) # cmultRepl() for multiplicative zero replacement (CLR pipeline)

# ---- 2. Load & filter ----------------------------------------------------
ps <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds")) |>
  subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))

# Global prevalence filter: keep features present in ≥10% of all biosamples.
# Mirrors filter_counts(counts_raw, 0.10) from notebooks/01b_its1_pca_by_farm.ipynb.
ps_prev <- ps |>
  (\(x) prune_taxa(
    rowSums(otu_table(x) > 0) / nsamples(x) >= PREV_THRESHOLD, x
  ))()

cat(sprintf("After prevalence filter: %d OTUs × %d samples\n",
            ntaxa(ps_prev), nsamples(ps_prev)))

# ---- 3. Distance matrices ------------------------------------------------
otu_mat <- as(otu_table(ps_prev), "matrix")          # features x samples
if (!taxa_are_rows(ps_prev)) otu_mat <- t(otu_mat)
otu_t   <- t(otu_mat)                                 # samples x features

# Create meta before distances so we can use its rownames as the canonical
# sample ID source — sample_data() preserves the original Excel IDs reliably.
meta    <- as(sample_data(ps_prev), "data.frame")
rownames(otu_t) <- rownames(meta)

# Robust Aitchison distance: rCLR + Euclidean, WITHOUT matrix completion.
# vegan >= 2.6-2 implements rCLR natively (geometric mean from non-zero values only;
# zeros receive 0). This differs from full RPCA (deicode, Python notebooks):
# deicode applies OptSpace matrix completion before rCLR. Use vegan rCLR for
# PERMANOVA — matrix completion is not designed for distance-matrix-based testing.
dist_ait <- vegdist(otu_t, method = "robust.aitchison")
# vegdist(robust.aitchison) sets Labels = "1","2",... instead of sample IDs
# (internal rCLR C routine strips rownames). Force correct Labels.
attr(dist_ait, "Labels") <- rownames(meta)

# CLR with multiplicative replacement (Martín-Fernández et al. 2003).
# cmultRepl replaces zeros while preserving compositional totals (CZM method).
# decostand(method="clr") is safe here because cmultRepl has already removed all zeros.
otu_t_mr     <- zCompositions::cmultRepl(otu_t, label = 0, method = "CZM",
                                          output = "p-counts", z.delete = FALSE)
clr_mat      <- vegan::decostand(otu_t_mr, method = "clr")  # log(x) - mean(log(x)) per row
dist_ait_clr <- dist(clr_mat)
# dist() already sets Labels = rownames(clr_mat); no manual override needed

# Rarefaction for depth-sensitive metrics.
# Rationale: unbagged flowers have systematically more reads than bagged across all farms
# (confirmed in 04b depth boxplots). Without rarefaction, depth confounds the bagged vs
# unbagged comparison for Bray-Curtis and Jaccard.
# Aitchison metrics above are inherently depth-robust and stay non-rarefied.
# ITS1 uses 15000 reads, consistent with alpha diversity rarefaction depth; ITS1 libraries are deeper than 16S.
set.seed(42)
ps_rare <- suppressMessages(
  phyloseq::rarefy_even_depth(ps_prev, sample.size = RAREFY_DEPTH,
                               rngseed = 42L, replace = FALSE,
                               trimOTUs = TRUE, verbose = FALSE)
)
cat(sprintf("After rarefaction to %d reads: %d samples remaining (dropped %d)\n",
            RAREFY_DEPTH, nsamples(ps_rare),
            nsamples(ps_prev) - nsamples(ps_rare)))

meta_rare  <- as(sample_data(ps_rare), "data.frame")
otu_t_rare <- as(otu_table(ps_rare), "matrix") |>
  (\(m) if (taxa_are_rows(ps_rare)) t(m) else m)()
rownames(otu_t_rare) <- rownames(meta_rare)

# Bray-Curtis on rarefied counts (abundance-weighted)
dist_bc  <- vegdist(otu_t_rare, method = "bray")
attr(dist_bc,  "Labels") <- rownames(meta_rare)

# Jaccard on rarefied presence/absence (binary)
dist_jac <- vegdist(otu_t_rare, method = "jaccard", binary = TRUE)
attr(dist_jac, "Labels") <- rownames(meta_rare)

dist_list <- list(
  aitchison_rCLR = dist_ait,      # non-rarefied: robust rCLR (zeros as missing)
  aitchison_CLR  = dist_ait_clr,  # non-rarefied: CLR + multiplicative zero replacement
  bray_curtis    = dist_bc,        # rarefied 15k: abundance-weighted dissimilarity
  jaccard        = dist_jac        # rarefied 15k: binary presence/absence
)

# meta_list maps each distance key to the correct metadata frame for that distance.
# Rarefied distances use meta_rare (samples that survived rarefaction); non-rarefied
# distances use meta (all biological samples after prevalence filter).
meta_list <- list(
  aitchison_rCLR = meta,
  aitchison_CLR  = meta,
  bray_curtis    = meta_rare,
  jaccard        = meta_rare
)

saveRDS(dist_list, here("results", "rds", "19b_ITS1_dist_list.rds"))
saveRDS(meta_list, here("results", "rds", "19b_ITS1_meta_list.rds"))

# ---- 4. Per-farm PERMANOVA -----------------------------------------------
# adonis2(d ~ sample_type, by = "terms", permutations = how(blocks = tree_id))
# blocked permutations: shuffles stay within trees (3 bagged + 3 unbagged per tree).
# Tests whether bagged/unbagged differ in community composition within each farm.

cat("Running per-farm PERMANOVA...\n")

results_perm_farm <- lapply(FARM_LEVELS, \(farm) {
  lapply(names(dist_list), \(dname) {
    meta_d    <- meta_list[[dname]]
    farm_idx  <- which(!is.na(meta_d$farm_id) & meta_d$farm_id == farm)
    meta_farm <- meta_d[farm_idx, ]
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

# ---- 4b. Per-farm PERMDISP -----------------------------------------------
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

# ---- 5. Combined PERMANOVA -----------------------------------------------
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

# ---- 6. PERMDISP ---------------------------------------------------------
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

# ---- 7. Visualisation ----------------------------------------------------
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
    title    = "PERMANOVA per farm — ITS1 OTU97 fungi",
    subtitle = sprintf("adonis2(d ~ sample_type, blocks = tree_id, nperm = %d)\nFill = R²; label = PERMANOVA p-value; bold border = PERMDISP p < 0.05",
                       N_PERM),
    x = "Farm", y = "Distance matrix"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here("results", "figures", "19b_ITS1_perm_farm_heatmap.png"),
       plot = p_farm, width = 9, height = 4, dpi = 300)

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
    title    = "PERMANOVA combined models — ITS1 OTU97 fungi",
    subtitle = sprintf("R² per term (blocked permutations, tree_id, %d permutations)", N_PERM),
    x = NULL, y = expression(R^2)
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here("results", "figures", "19b_ITS1_perm_combined_r2.png"),
       plot = p_comb, width = 12, height = 5, dpi = 300)

# ---- 8. Save outputs -----------------------------------------------------
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

write_csv(results_perm_farm,
          here("results", "tables", "19b_ITS1_permanova_per_farm.csv"))
write_csv(results_combined,
          here("results", "tables", "19b_ITS1_permanova_combined.csv"))
write_csv(results_permdisp,
          here("results", "tables", "19b_ITS1_permdisp.csv"))
write_csv(results_permdisp_farm,
          here("results", "tables", "19b_ITS1_permdisp_per_farm.csv"))

cat("\nPERMDISP summary:\n")
print(results_permdisp)

message("\nDone. CSVs and figures saved.")
