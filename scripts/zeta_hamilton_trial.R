################################################################################
######################## Zeta diversity on Hamilton HPC #######################
################################################################################

# ---- Load packages ----
library(phyloseq)
library(zetadiv)
library(dplyr)
library(ggplot2)

# ---- Define project paths ----
project_dir <- "/nobackup/jcnx71/cacao_flower_microbiome"
rds_dir     <- file.path(project_dir, "results", "rds")
results_dir <- file.path(project_dir, "results", "zeta")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load phyloseq object ----
phyloseq_16S_bacteria_biosamples <-
  readRDS(file.path(rds_dir, "phyloseq_16S_bacteria_biosamples.rds"))

# ---- Convert to incidence matrix ----
incidence_df_16S <-
  phyloseq_16S_bacteria_biosamples |>
  otu_table() |>
  as("matrix") |>
  t() |>
  (\(x) {
    x[x > 0] <- 1
    as.data.frame(x)
  })()

# ---- Metadata ----
metadata_16S <-
  sample_data(phyloseq_16S_bacteria_biosamples) |>
  as("data.frame")

stopifnot(
  nrow(incidence_df_16S) == nrow(metadata_16S)
)

################################################################################
########################## Zeta Monte Carlo ###################################
################################################################################

cat("Starting zeta calculation...\n")
start_time <- Sys.time()

set.seed(42)
zeta_bagged_16S_10k <-
  incidence_df_16S[
    metadata_16S$sample_type == "bagged_flower",
  ] |>
  Zeta.decline.mc(
    orders = 1:100,
    sam = 10000,
    plot = FALSE
  )

saveRDS(
  zeta_bagged_16S_10k,
  file.path(results_dir, "zeta_bagged_16S_10k.rds"),
  compress = "xz"
)

set.seed(42)
zeta_unbagged_16S_10k <-
  incidence_df_16S[
    metadata_16S$sample_type == "unbagged_flower",
  ] |>
  Zeta.decline.mc(
    orders = 1:100,
    sam = 10000,
    plot = FALSE
  )

saveRDS(
  zeta_unbagged_16S_10k,
  file.path(results_dir, "zeta_unbagged_16S_10k.rds"),
  compress = "xz"
)

end_time <- Sys.time()
cat("Zeta finished in:", end_time - start_time, "\n")

################################################################################
########################## Plot results #######################################
################################################################################

zeta_plot_tbl_16S_10k <-
  bind_rows(
    data.frame(
      order = zeta_bagged_16S_10k$zeta.order,
      zeta = pmax(zeta_bagged_16S_10k$zeta.val, .Machine$double.eps),
      zeta_sd = zeta_bagged_16S_10k$zeta.val.sd,
      group = "bagged_flower"
    ),
    data.frame(
      order = zeta_unbagged_16S_10k$zeta.order,
      zeta = pmax(zeta_unbagged_16S_10k$zeta.val, .Machine$double.eps),
      zeta_sd = zeta_unbagged_16S_10k$zeta.val.sd,
      group = "unbagged_flower"
    )
  )

p <-
  ggplot(zeta_plot_tbl_16S_10k,
         aes(x = order, y = zeta, colour = group)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5) +
  geom_errorbar(
    aes(ymin = pmax(zeta - zeta_sd, .Machine$double.eps),
        ymax = zeta + zeta_sd),
    width = 0.15
  ) +
  theme_bw() +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  labs(
    x = "Zeta order",
    y = "Shared ASVs",
    title = "Zeta diversity decline (16S)",
    colour = "Treatment"
  )

ggsave(
  file.path(results_dir, "zeta_16S_10k_plot.png"),
  p,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Plot saved.\n")
