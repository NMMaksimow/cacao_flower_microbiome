################################################################################
# zeta_array.R — compute ONE zeta order (MC) for bagged and unbagged flowers
#
# Called by zeta_16S_array.sh / zeta_ITS1_array.sh via:
#   Rscript scripts/zeta_array.R --amplicon 16S --order 5 --sam 10000
#
# Output: results/rds/zeta_tmp/zeta_{amplicon}_{group}_order{i}.rds
################################################################################

# ---- Parse command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, args) {
  idx <- which(args == flag)
  if (length(idx) == 0) stop(paste("Missing argument:", flag))
  args[idx + 1]
}

amplicon <- get_arg("--amplicon", args)   # "16S" or "ITS1"
order_i  <- as.integer(get_arg("--order", args))
sam      <- as.integer(get_arg("--sam",   args))

cat(sprintf("Amplicon: %s | Order: %d | sam: %d\n", amplicon, order_i, sam))

# ---- Packages ----
library(phyloseq)
library(zetadiv)

# ---- Paths ----
project_dir <- "/nobackup/jcnx71/cacao_flower_microbiome"
rds_dir     <- file.path(project_dir, "results", "rds")
tmp_dir     <- file.path(project_dir, "results", "rds", "zeta_tmp")

dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load correct phyloseq object ----
if (amplicon == "16S") {
  ps <- readRDS(file.path(rds_dir, "phyloseq_16S_bacteria_biosamples.rds"))
} else if (amplicon == "ITS1") {
  ps <- readRDS(file.path(rds_dir, "phyloseq_ITS1_fungi_biosamples.rds"))
} else {
  stop(paste("Unknown amplicon:", amplicon))
}

# ---- Build incidence matrix ----
incidence_mat <-
  ps |>
  otu_table() |>
  as("matrix") |>
  t() |>
  (\(x) { x[x > 0] <- 1; as.data.frame(x) })()

metadata <-
  sample_data(ps) |>
  as("data.frame")

stopifnot(all(rownames(incidence_mat) == rownames(metadata)))

# ---- Compute zeta for this order — bagged flowers ----
cat(sprintf("Computing order %d | bagged_flower ...\n", order_i))
set.seed(order_i * 100 + 1)   # reproducible per-order seed

result_bagged <- Zeta.order.mc(
  data.spec = incidence_mat[metadata$sample_type == "bagged_flower", ],
  order     = order_i,
  sam       = sam
)

saveRDS(
  result_bagged,
  file.path(tmp_dir, sprintf("zeta_%s_bagged_order%03d.rds", amplicon, order_i)),
  compress = "xz"
)

# ---- Compute zeta for this order — unbagged flowers ----
cat(sprintf("Computing order %d | unbagged_flower ...\n", order_i))
set.seed(order_i * 100 + 2)

result_unbagged <- Zeta.order.mc(
  data.spec = incidence_mat[metadata$sample_type == "unbagged_flower", ],
  order     = order_i,
  sam       = sam
)

saveRDS(
  result_unbagged,
  file.path(tmp_dir, sprintf("zeta_%s_unbagged_order%03d.rds", amplicon, order_i)),
  compress = "xz"
)

cat(sprintf("Done: order %d saved.\n", order_i))
