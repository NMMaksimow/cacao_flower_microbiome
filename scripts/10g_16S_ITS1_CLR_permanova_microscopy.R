## ============================================================================
## 10g — Microscopy phenotypes and community composition: PERMANOVA and PERMDISP
##
## Do flowers with more pollen, more callose, germinated pollen or a heavier
## fungal colonisation of the stigma carry a different microbial community?
##
## Unbagged flowers only for the main tests: in bagged flowers the fungal score
## is 0 in 77 of 81 and callose is near absent, so these phenotypes are almost
## nested in treatment. Bagged flowers: pi_clsm only (pollen reached some bagged
## stigmas; 29 of 90 have pollen).
##
## Model, Aitchison CLR distance, flowers with a value for the feature:
##   PERMANOVA  d ~ farm_id + feature, by = "terms", permutations within farm
##   PERMDISP   distance of each flower to its farm centroid ~ farm_id + feature,
##              same permutations (does the spread around the farm centroid change
##              with the feature?)
## BH over the five tests of a marker, separately for PERMANOVA and PERMDISP.
## Flowers of one tree are not independent; permutations within farm do not
## account for that, so p-values lean liberal.
##
## Features (transform):
##   fungal_colonization_score  0-5 as a number (linear trend), unbagged
##   pi_clsm                    log1p (pollen count; median 2, max 349), unbagged / bagged
##   callose_plug_index         log1p (median 19, max 443), unbagged
##   is_germination             yes / no, unbagged (secondary)
##
## Input:
##   results/rds/{16S,ITS1_otu97}_dist_list.rds, _meta_list.rds (08a/08b)
## Output:
##   results/tables/10g_permanova_microscopy.csv
## ============================================================================

# ---- 0. Parameters -------------------------------------------------------
N_PERM    <- 999
DIST_NAME <- "aitchison_CLR"   # primary metric

# Marker label -> distance RDS prefix (08a/08b)
MARKERS <- c("16S" = "16S", "ITS1" = "ITS1_otu97")

# Tests: which flowers, which feature, how it enters the model
TESTS <- tibble::tribble(
  ~flowers,           ~feature,                     ~transform,
  "unbagged_flower",  "fungal_colonization_score",  "numeric",
  "unbagged_flower",  "pi_clsm",                    "log1p",
  "unbagged_flower",  "callose_plug_index",         "log1p",
  "unbagged_flower",  "is_germination",             "factor",
  "bagged_flower",    "pi_clsm",                    "log1p"
)

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(vegan)
library(permute)   # how() for permutations within farm

# Data version and sample set: RESULTS_DIR, BIOSAMPLE_SUFFIX, SET_TAG
source(here("scripts", "utils_paths.R"))

# ---- 2. One test ---------------------------------------------------------
transform_feature <- function(x, transform) {
  switch(transform,
         numeric = as.numeric(as.character(x)),   # ordered factor 0-5 -> 0..5
         log1p   = log1p(x),
         factor  = factor(x))
}

run_test <- function(d, meta, flowers, feature, transform) {
  keep  <- meta$sample_type == flowers & !is.na(meta[[feature]])
  m     <- meta[keep, ] |> mutate(x = transform_feature(.data[[feature]], transform))
  d_sub <- as.dist(as.matrix(d)[keep, keep])
  perm  <- how(nperm = N_PERM, blocks = m$farm_id)

  set.seed(42)
  res <- adonis2(d_sub ~ farm_id + x, data = m, by = "terms", permutations = perm)

  # Dispersion: distance to the farm centroid, then the same model on that distance
  bd <- betadisper(d_sub, m$farm_id)
  set.seed(42)
  res_disp <- adonis2(dist(bd$distances) ~ farm_id + x, data = m, by = "terms", permutations = perm)

  tibble(n = nrow(m),
         R2 = res["x", "R2"], F = res["x", "F"], p = res["x", "Pr(>F)"],
         farm_R2 = res["farm_id", "R2"],
         disp_F = res_disp["x", "F"], disp_p = res_disp["x", "Pr(>F)"])
}

# ---- 3. All tests --------------------------------------------------------
results <- imap(MARKERS, \(prefix, marker) {
  d    <- readRDS(file.path(RESULTS_DIR, "rds", paste0(prefix, "_dist_list", BIOSAMPLE_SUFFIX, ".rds")))[[DIST_NAME]]
  meta <- readRDS(file.path(RESULTS_DIR, "rds", paste0(prefix, "_meta_list", BIOSAMPLE_SUFFIX, ".rds")))[[DIST_NAME]]
  stopifnot(identical(labels(d), rownames(meta)))

  TESTS |>
    mutate(res = pmap(list(flowers, feature, transform), \(f, v, t) run_test(d, meta, f, v, t))) |>
    unnest(res) |>
    mutate(q = p.adjust(p, method = "BH"), disp_q = p.adjust(disp_p, method = "BH"))
}) |> list_rbind(names_to = "marker")

# ---- 4. Save and print ---------------------------------------------------
dir.create(file.path(RESULTS_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
write_csv(results, file.path(RESULTS_DIR, "tables", paste0("10g_permanova_microscopy", SET_TAG, ".csv")))

results |>
  mutate(across(c(R2, F, farm_R2, disp_F), \(x) round(x, 3)),
         across(c(q, disp_q), \(x) round(x, 3))) |>
  select(marker, flowers, feature, n, R2, F, p, q, disp_p, disp_q) |>
  print(n = Inf, width = Inf)
cat("10g done.\n")
