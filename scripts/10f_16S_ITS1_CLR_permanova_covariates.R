## ============================================================================
## 10f — Does the treatment effect survive technical covariates? 16S and ITS1
##
## Within a farm, bagged and unbagged flowers were extracted (and, for ITS1,
## amplified) in separate batches, and ITS1 DNA yield and read depth are higher
## in unbagged flowers (step 04). This step asks whether the treatment term of
## 10a/10b stays when these technical variables are fitted first.
##
## Reading: the test works one way. A treatment effect that stays after the
## covariates is strong evidence; one that goes proves nothing, because DNA
## yield and depth may themselves be consequences of insect visits.
##
## Model (as the flower-level PERMANOVA of 10a/10b), Aitchison CLR distance:
##   d ~ tree_id + <covariates> + sample_type
##   by = "terms" (sequential), permutations within tree. Treatment comes last, so
##   it is tested after tree and covariates. Farm needs no own term: each tree
##   belongs to one farm, so tree_id already holds the farm differences.
## Covariates, alone and together:
##   log_depth  log of the read count of the flower in the biosample object
##   log_dna    log of the PCR product concentration (dna_concentration_hs_*)
##   pcr_redo   at least one PCR replicate repeated on another day (01c), yes / no
##
## Input:
##   results/rds/{16S,ITS1_otu97}_dist_list.rds, _meta_list.rds (08a/08b)
##   results/rds/ps_16S_bacteria_biosamples.rds, ps_ITS1_otu97_fungi_biosamples.rds (depth)
## Output:
##   results/tables/10f_permanova_covariates.csv
## ============================================================================

# ---- 0. Parameters -------------------------------------------------------
N_PERM    <- 999
DIST_NAME <- "aitchison_CLR"   # primary metric

# Marker label -> distance RDS prefix (08a/08b), phyloseq object, metadata columns
MARKERS <- list(
  "16S"  = list(prefix = "16S",        ps = "ps_16S_bacteria_biosamples",
                dna = "dna_concentration_hs_16s",  redo = "pcr_redo_16s"),
  "ITS1" = list(prefix = "ITS1_otu97", ps = "ps_ITS1_otu97_fungi_biosamples",
                dna = "dna_concentration_hs_its1", redo = "pcr_redo_its1")
)

# Covariate sets fitted between tree_id and sample_type
COVARIATE_SETS <- list(
  none     = character(0),
  depth    = "log_depth",
  dna      = "log_dna",
  pcr_redo = "pcr_redo",
  all      = c("log_depth", "log_dna", "pcr_redo")
)

# ---- 1. Libraries --------------------------------------------------------
library(here)
library(tidyverse)
library(phyloseq)
library(vegan)
library(permute)   # how() for permutations within tree

# Data version and sample set: RESULTS_DIR, BIOSAMPLE_SUFFIX, SET_TAG
source(here("scripts", "utils_paths.R"))

# ---- 2. Data -------------------------------------------------------------
load_marker <- function(m) {
  d    <- readRDS(file.path(RESULTS_DIR, "rds", paste0(m$prefix, "_dist_list", BIOSAMPLE_SUFFIX, ".rds")))[[DIST_NAME]]
  meta <- readRDS(file.path(RESULTS_DIR, "rds", paste0(m$prefix, "_meta_list", BIOSAMPLE_SUFFIX, ".rds")))[[DIST_NAME]]
  depth <- readRDS(file.path(RESULTS_DIR, "rds", paste0(m$ps, BIOSAMPLE_SUFFIX, ".rds"))) |> sample_sums()

  stopifnot(identical(labels(d), rownames(meta)), all(rownames(meta) %in% names(depth)),
            all(meta[[m$dna]] > 0))   # log needs positive concentrations

  meta <- meta |>
    mutate(log_depth = log(depth[rownames(meta)]),
           log_dna   = log(.data[[m$dna]]),
           pcr_redo  = factor(if_else(.data[[m$redo]] > 0, "yes", "no")))
  list(d = d, meta = meta)
}

data_by_marker <- map(MARKERS, load_marker)

# ---- 3. PERMANOVA with each covariate set --------------------------------
fit_set <- function(dat, covariates) {
  rhs <- paste(c("tree_id", covariates, "sample_type"), collapse = " + ")
  set.seed(42)
  res <- adonis2(as.formula(paste("dat$d ~", rhs)),
                 data         = dat$meta,
                 by           = "terms",
                 permutations = how(nperm = N_PERM, blocks = dat$meta$tree_id))
  as.data.frame(res) |>
    rownames_to_column("term") |>
    filter(term %in% c(covariates, "sample_type")) |>
    transmute(term, df = Df, R2, F, p = `Pr(>F)`, formula = rhs)
}

results <- map(data_by_marker, \(dat) {
  map(COVARIATE_SETS, \(cv) fit_set(dat, cv)) |> list_rbind(names_to = "covariates")
}) |> list_rbind(names_to = "marker")

# ---- 4. Save and print ---------------------------------------------------
dir.create(file.path(RESULTS_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
write_csv(results, file.path(RESULTS_DIR, "tables", paste0("10f_permanova_covariates", SET_TAG, ".csv")))

results |>
  as_tibble() |>
  mutate(across(c(R2, F), \(x) round(x, 4))) |>
  select(marker, covariates, term, R2, F, p) |>
  print(n = Inf)
cat("10f done.\n")
