## ============================================================================
## 08c — CLR matrices (CZM zero replacement) for the random forest notebooks
##
## Single source of truth for the CLR features: the Python random forest (notebooks
## 03, 03b, 03c) reads these matrices instead of computing its own CLR. The
## computation is clr_czm() of utils_clr.R, the same function 08a/08b use for the
## Aitchison CLR distance, so config "all10" reproduces the 08a/08b CLR exactly.
##
## No distances, no rarefaction: prevalence filter -> cmultRepl (CZM) -> CLR only.
##
## Configurations (prevalence filter before zero replacement):
##   all10       >= 10% of all flowers (bagged + unbagged)   main, as 08a/08b
##   all05       >= 5% of all flowers                        sensitivity (03c)
##   all20       >= 20% of all flowers                       sensitivity (03c)
##   unbagged05  >= 5% of unbagged flowers                   sensitivity (03c)
##   unbagged10  >= 10% of unbagged flowers                  sensitivity (03c)
## The unbagged configs test whether taxa of the unbagged-only targets (pollination,
## fungal colonisation) are lost because prevalence is computed over all flowers.
## CLR is always computed on all flowers; the filter only chooses the features.
##
## No 0% matrix: cmultRepl time grows with the number of features (16S all10, 209
## features: ~20 s; 16S at 0%, ~26 700 features: ~45 min estimated). 03c tests the
## 0% threshold on TSS features, which need no zero replacement.
##
## Input:
##   results/<version>/rds/ps_16S_bacteria_biosamples<suffix>.rds
##   results/<version>/rds/ps_ITS1_otu97_fungi_biosamples<suffix>.rds
## Output (samples x features, first column sample_id):
##   results/<version>/tables/clr/08c_<marker>_clr_czm_<config><SET_TAG>.csv.gz
##
## Run for both sample sets (CFM_SAMPLE_SET = main | with_deviations), e.g.
##   WSLENV=CFM_SAMPLE_SET CFM_SAMPLE_SET=with_deviations Rscript.exe scripts/08c_16S_ITS1_clr_export.R
## ============================================================================

# ---- 0. Parameters -------------------------------------------------------
# Marker label (file name) -> phyloseq object (01a/01b)
MARKERS <- c("16S"  = "ps_16S_bacteria_biosamples",
             "ITS1" = "ps_ITS1_otu97_fungi_biosamples")

# Prevalence configurations: threshold and the flowers prevalence is computed in
CONFIGS <- tibble::tribble(
  ~config,      ~thresh, ~ref,
  "all10",      0.10,    "all",
  "all05",      0.05,    "all",
  "all20",      0.20,    "all",
  "unbagged05", 0.05,    "unbagged",
  "unbagged10", 0.10,    "unbagged"
)

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(phyloseq)

source(here("scripts", "utils_paths.R"))   # RESULTS_DIR, BIOSAMPLE_SUFFIX, SET_TAG
source(here("scripts", "utils_clr.R"))     # prevalence_keep(), clr_czm()

OUT_DIR <- file.path(RESULTS_DIR, "tables", "clr")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Counts, as in 08a/08b --------------------------------------------
# samples x features count matrix of bagged + unbagged flowers, rownames = sample IDs
read_counts <- function(ps_name) {
  ps <- readRDS(file.path(RESULTS_DIR, "rds", paste0(ps_name, BIOSAMPLE_SUFFIX, ".rds"))) |>
    subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  rownames(otu) <- rownames(as(sample_data(ps), "data.frame"))
  list(otu = otu, sample_type = as(sample_data(ps), "data.frame")$sample_type)
}

# ---- 3. Filter -> CZM -> CLR -> csv.gz -------------------------------------
for (marker in names(MARKERS)) {
  dat <- read_counts(MARKERS[[marker]])
  cat(sprintf("%s: %d flowers x %d features\n", marker, nrow(dat$otu), ncol(dat$otu)))

  CONFIGS |>
    pwalk(\(config, thresh, ref) {
      ref_rows <- if (ref == "all") NULL else which(dat$sample_type == "unbagged_flower")
      keep <- prevalence_keep(dat$otu, thresh, ref_rows)
      t0   <- Sys.time()
      # cmultRepl warns for every row with > 80% zeros (most rows at low thresholds);
      # the warnings carry no decision here, z.delete = FALSE keeps everything.
      clr  <- suppressWarnings(clr_czm(dat$otu[, keep, drop = FALSE]))
      out  <- file.path(OUT_DIR, sprintf("08c_%s_clr_czm_%s%s.csv.gz", marker, config, SET_TAG))
      clr |>
        as.data.frame() |>
        rownames_to_column("sample_id") |>
        write_csv(out)
      cat(sprintf("  %-10s %5d features  %5.1f s  -> %s\n", config, sum(keep),
                  as.numeric(difftime(Sys.time(), t0, units = "secs")), basename(out)))
    })
}
