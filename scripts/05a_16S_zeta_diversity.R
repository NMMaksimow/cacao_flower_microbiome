################################################################################
####################### 16S Zeta-diversity decline ############################
################################################################################
# Question: Do flower visitors homogenize or heterogenize the bacterial microbiome?
# Prediction:
#   - unbagged (visitor-exposed) ~ power law  → homogenization
#   - unbagged (visitor-exposed) ~ exp law    → heterogenization
#   - parallel curves            → visitors have negligible effect
#
# Input:  results/rds/phyloseq_16S_bacteria_biosamples.rds
# Output: results/rds/zeta_16S_*.rds
#         results/figures/zeta_16S_*.png

# ---- Packages ----
library(phyloseq)
library(zetadiv)
library(dplyr)
library(ggplot2)

# ---- Paths ----
project_dir <- here::here()           # project root (works with RStudio / Rscript)
rds_dir     <- file.path(project_dir, "results", "rds")
fig_dir     <- file.path(project_dir, "results", "figures")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load phyloseq ----
ps <- readRDS(file.path(rds_dir, "phyloseq_16S_bacteria_biosamples.rds"))

# ---- ASV incidence matrix (samples × ASVs, presence/absence) ----
incidence_16S <-
  ps |>
  otu_table() |>
  as("matrix") |>
  t() |>
  (\(x) { x[x > 0] <- 1; as.data.frame(x) })()

metadata_16S <-
  sample_data(ps) |>
  as("data.frame")

stopifnot(all(rownames(incidence_16S) == rownames(metadata_16S)))

################################################################################
# 1.  Zeta decline on ASV level (orders 1–10, MC sam = 1 000)
################################################################################

set.seed(42)
zeta_bagged_16S <-
  incidence_16S[metadata_16S$sample_type == "bagged_flower", ] |>
  Zeta.decline.mc(orders = 1:10, sam = 1000, plot = FALSE)

set.seed(42)
zeta_unbagged_16S <-
  incidence_16S[metadata_16S$sample_type == "unbagged_flower", ] |>
  Zeta.decline.mc(orders = 1:10, sam = 1000, plot = FALSE)

saveRDS(zeta_bagged_16S,   file.path(rds_dir, "zeta_bagged_16S_1k.rds"),   compress = "xz")
saveRDS(zeta_unbagged_16S, file.path(rds_dir, "zeta_unbagged_16S_1k.rds"), compress = "xz")

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

cat("=== Bagged flowers – 16S ASV ===\n")
print(extract_zeta_fits(zeta_bagged_16S))

cat("=== Unbagged flowers – 16S ASV ===\n")
print(extract_zeta_fits(zeta_unbagged_16S))

# ---- Plot ----
zeta_tbl_16S <-
  dplyr::bind_rows(
    data.frame(order    = zeta_bagged_16S$zeta.order,
               zeta     = zeta_bagged_16S$zeta.val,
               zeta_sd  = zeta_bagged_16S$zeta.val.sd,
               group    = "bagged_flower"),
    data.frame(order    = zeta_unbagged_16S$zeta.order,
               zeta     = zeta_unbagged_16S$zeta.val,
               zeta_sd  = zeta_unbagged_16S$zeta.val.sd,
               group    = "unbagged_flower")
  )

p_16S_asv <-
  ggplot(zeta_tbl_16S, aes(x = order, y = zeta, colour = group)) +
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
    y      = "Zeta diversity (average shared ASVs)",
    title  = "16S zeta-diversity decline – ASV level",
    colour = "Flower treatment"
  )

ggsave(file.path(fig_dir, "zeta_16S_asv_1k.png"),
       p_16S_asv, width = 7, height = 5, dpi = 300)

################################################################################
# 2.  Zeta decline on Genus level  (sensitivity / noise-reduction)
################################################################################

# Agglomerate to genus, assigning unique labels to NA genera so tax_glom()
# does not merge all unannotated ASVs into one phantom genus.
ps_genus <-
  ps |>
  (\(physeq) {
    tax_mat <- tax_table(physeq) |> as.matrix()
    na_genus <- is.na(tax_mat[, "Genus"])
    tax_mat[na_genus, "Genus"] <-
      paste0("Unclassified_g__ASV_", rownames(tax_mat)[na_genus])
    tax_table(physeq) <- tax_table(tax_mat)
    phyloseq::tax_glom(physeq, taxrank = "Genus", NArm = FALSE)
  })()

cat(sprintf(
  "ASV count: %d  →  Genus count: %d  (of which unclassified: %d)\n",
  ntaxa(ps),
  ntaxa(ps_genus),
  sum(grepl("^Unclassified_g__ASV_", tax_table(ps_genus)[, "Genus"]))
))

incidence_16S_genus <-
  ps_genus |>
  otu_table() |>
  as("matrix") |>
  t() |>
  (\(x) { x[x > 0] <- 1; as.data.frame(x) })()

set.seed(42)
zeta_bagged_16S_genus <-
  incidence_16S_genus[metadata_16S$sample_type == "bagged_flower", ] |>
  Zeta.decline.mc(orders = 1:10, sam = 1000, plot = FALSE)

set.seed(42)
zeta_unbagged_16S_genus <-
  incidence_16S_genus[metadata_16S$sample_type == "unbagged_flower", ] |>
  Zeta.decline.mc(orders = 1:10, sam = 1000, plot = FALSE)

saveRDS(zeta_bagged_16S_genus,   file.path(rds_dir, "zeta_bagged_16S_genus_1k.rds"),   compress = "xz")
saveRDS(zeta_unbagged_16S_genus, file.path(rds_dir, "zeta_unbagged_16S_genus_1k.rds"), compress = "xz")

cat("=== Bagged flowers – 16S Genus ===\n")
print(extract_zeta_fits(zeta_bagged_16S_genus))
cat("=== Unbagged flowers – 16S Genus ===\n")
print(extract_zeta_fits(zeta_unbagged_16S_genus))

# ---- Plot genus ----
zeta_tbl_16S_genus <-
  dplyr::bind_rows(
    data.frame(order    = zeta_bagged_16S_genus$zeta.order,
               zeta     = zeta_bagged_16S_genus$zeta.val,
               zeta_sd  = zeta_bagged_16S_genus$zeta.val.sd,
               group    = "bagged_flower"),
    data.frame(order    = zeta_unbagged_16S_genus$zeta.order,
               zeta     = zeta_unbagged_16S_genus$zeta.val,
               zeta_sd  = zeta_unbagged_16S_genus$zeta.val.sd,
               group    = "unbagged_flower")
  )

p_16S_genus <-
  ggplot(zeta_tbl_16S_genus, aes(x = order, y = zeta, colour = group)) +
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
    y      = "Zeta diversity (average shared genera)",
    title  = "16S zeta-diversity decline – Genus level",
    colour = "Flower treatment"
  )

ggsave(file.path(fig_dir, "zeta_16S_genus_1k.png"),
       p_16S_genus, width = 7, height = 5, dpi = 300)

cat("Done. Results saved to results/rds/ and results/figures/\n")
