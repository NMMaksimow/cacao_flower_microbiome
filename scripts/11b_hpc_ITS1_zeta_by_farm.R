# ---- 0. Parameters -------------------------------------------------------
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
# Orders 1:21 = full range up to the group size (21 samples per farm × treatment).
ZETA_ORDERS <- 1:21
# MC iterations: 1000 per order.
SAM <- 1000L

# ---- 1. Libraries (main thread only) -------------------------------------
library(phyloseq)
library(dplyr)
library(tibble)
library(parallel)

# ---- 2. Paths ------------------------------------------------------------
project_dir <- "/nobackup/jcnx71/cacao_flower_microbiome"
rds_dir     <- file.path(project_dir, "results", "rds")
tables_dir  <- file.path(project_dir, "results", "tables")
dir.create(rds_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 3. Load and pre-process data ----------------------------------------
# OTU97: consistent with all other ITS1 diversity analyses.
# Note: use ps_ITS1_otu97_fungi_biosamples.rds — NOT ps_ITS1_fungi_biosamples.rds
# (the older non-otu97 file used in zeta_array.R is a pre-existing bug there).
cat("Loading phyloseq...\n")
ps  <- readRDS(file.path(rds_dir, "ps_ITS1_otu97_fungi_biosamples.rds"))
ps  <- subset_samples(ps, sample_type %in% c("bagged_flower", "unbagged_flower"))

incidence_full <- otu_table(ps) |>
  as("matrix") |>
  (\(m) if (taxa_are_rows(ps)) t(m) else m)() |>
  (\(m) { m[m > 0] <- 1L; as.data.frame(m) })()

meta <- as(sample_data(ps), "data.frame")

cat(sprintf("Loaded: %d samples × %d OTUs\n", nrow(incidence_full), ncol(incidence_full)))
stopifnot(all(rownames(incidence_full) == rownames(meta)))

combo_df <- expand.grid(
  farm  = FARM_LEVELS,
  group = c("bagged_flower", "unbagged_flower"),
  stringsAsFactors = FALSE
)

combos <- lapply(seq_len(nrow(combo_df)), \(i) {
  farm      <- combo_df$farm[i]
  grp       <- combo_df$group[i]
  farm_rows <- which(meta$farm_id == farm)
  grp_rows  <- farm_rows[meta$sample_type[farm_rows] == grp]
  sub       <- incidence_full[grp_rows, ]
  # Trim OTU columns absent from this farm-group subset — critical for speed.
  sub       <- sub[, colSums(sub) > 0, drop = FALSE]
  cat(sprintf("  prepared [%s × %s]: %d samples × %d OTUs\n",
              farm, grp, nrow(sub), ncol(sub)))
  list(farm = farm, group = grp, data = sub, seed = i * 100L + 1L)
})

# ---- 4. PSOCK parallel zeta computation ----------------------------------
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat(sprintf("Starting PSOCK cluster: %d cores for %d farm × group combinations\n",
            N_CORES, nrow(combo_df)))

cl <- parallel::makePSOCKcluster(N_CORES)
on.exit(parallel::stopCluster(cl), add = TRUE)

parallel::clusterExport(cl, c("ZETA_ORDERS", "SAM"))
parallel::clusterEvalQ(cl, library(zetadiv))

results_flat <- parallel::parLapply(cl, combos, function(item) {
  set.seed(item$seed)
  list(
    farm  = item$farm,
    group = item$group,
    zeta  = Zeta.decline.mc(
      data.spec = item$data,
      orders    = ZETA_ORDERS,
      sam       = SAM,
      plot      = FALSE
    )
  )
})

# ---- 5. Reassemble nested list [[farm]][[group]] -------------------------
zeta_results <- list()
for (r in results_flat) {
  if (is.null(zeta_results[[r$farm]])) zeta_results[[r$farm]] <- list()
  zeta_results[[r$farm]][[r$group]] <- r$zeta
}

saveRDS(zeta_results,
        file.path(rds_dir, "zeta_ITS1_by_farm.rds"),
        compress = "xz")
cat("Saved: results/rds/zeta_ITS1_by_farm.rds\n")

cat("Done. Model fits are computed in 11c_16S_ITS1_zeta_by_farm.R.\n")
