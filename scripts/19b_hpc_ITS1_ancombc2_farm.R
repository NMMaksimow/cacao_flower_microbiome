## ============================================================================
## 19b — ANCOM-BC2: ITS1 OTU97 fungi — differences between farms (HPC)
##
## Which taxa differ between farms, over all farms (global test) and for each
## pair of farms (pairwise test, 21 pairs). Same feature set as steps 16-18.
##
## Two models, each with its own global and pairwise test:
##   farm:       sample_type + farm_id, (1 | tree_id)
##   management: sample_type + management_type, (1 | farm_id/tree_id)
##               the replicate is the farm (7), so descriptive at most
##
## A taxon counts as different only if it passes FDR and the pseudo-count
## sensitivity filter (passed_ss), as in steps 16-18.
##
## The run is slow (mixed model per taxon, then mdFDR permutations), so each
## model saves its complete output as soon as it finishes.
##
## Input:
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output:
##   results/rds/19b_ancombc2_{farm,management}.rds            complete ancombc2 output
##   results/tables/19b_ancombc2_{farm,management}_global.csv  global test
##   results/tables/19b_ancombc2_{farm,management}_pairwise.csv pairwise test
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

FDR_ALPHA   <- 0.05   # BH-corrected p-value threshold
PREV_THRESH <- 0.10   # min fraction of samples a taxon must be present in
MIN_LIB     <- 500    # min library size; samples below are dropped
N_CORES     <- 16     # ancombc2 workers; keep equal to -c in the job script
SEED        <- 123

FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

# ── 1. Libraries & data ───────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

# Data version and sample set: RESULTS_DIR, BIOSAMPLE_SUFFIX, SET_TAG
source(here("scripts", "utils_paths.R"))

# Taxonomy labels, shared with the other ANCOM-BC2 scripts
source(here("scripts", "utils_tax_label.R"))

cat(sprintf("ANCOMBC %s\n", packageVersion("ANCOMBC")))

ps_raw <- readRDS(file.path(RESULTS_DIR, "rds", paste0("ps_ITS1_otu97_fungi_biosamples", BIOSAMPLE_SUFFIX, ".rds")))

cat(sprintf("Loaded: %d taxa x %d samples\n", ntaxa(ps_raw), nsamples(ps_raw)))

# ── 2. Feature set (as in 16-18) ──────────────────────────────────────────────

filter_ps <- function(physeq, prev = PREV_THRESH, min_lib = MIN_LIB) {
        physeq <- prune_samples(sample_sums(physeq) >= min_lib, physeq)
        keep   <- rowSums(otu_table(physeq) > 0) / nsamples(physeq) >= prev
        prune_taxa(keep, physeq)
}

ps_f <- filter_ps(ps_raw)

sample_data(ps_f)$farm_id     <- factor(sample_data(ps_f)$farm_id, levels = FARM_LEVELS)
sample_data(ps_f)$sample_type <- factor(sample_data(ps_f)$sample_type,
                                        levels = c("bagged_flower", "unbagged_flower"))
sample_data(ps_f)$management_type <- factor(sample_data(ps_f)$management_type)

cat(sprintf("After prevalence/library filter: %d taxa x %d samples\n", ntaxa(ps_f), nsamples(ps_f)))

# ── 3. Models ─────────────────────────────────────────────────────────────────

# The tested group must come LAST in fix_formula. With a random effect the global
# test (.ancombc_global_LRT, ANCOMBC 2.14.0) builds the reduced model by deleting
# "+ group" from the formula; a group in first place has no "+" before it, so the
# reduced model equals the full one and every taxon gets W = 0, p = NA, no error.
models <- list(
        farm = list(
                group        = "farm_id",
                fix_formula  = "sample_type + farm_id",
                rand_formula = "(1|tree_id)"
        ),
        management = list(
                group        = "management_type",
                fix_formula  = "sample_type + management_type",
                rand_formula = "(1|farm_id/tree_id)"
        )
)

# ── 4. Run and save each model ────────────────────────────────────────────────

for (model_name in names(models)) {
        model <- models[[model_name]]
        cat(sprintf("\n── %s: %s, %s ──\n", model_name, model$fix_formula, model$rand_formula))
        started <- Sys.time()

        out <- ancombc2(
                data         = ps_f,
                fix_formula  = model$fix_formula,
                rand_formula = model$rand_formula,
                group        = model$group,
                p_adj_method = "BH",
                alpha        = FDR_ALPHA,
                global       = TRUE,    # one test per taxon over all groups
                pairwise     = TRUE,    # every pair of groups, mdFDR controlled
                n_cl         = N_CORES,
                verbose      = FALSE
        )

        # Stop rather than save an empty global test (see the note in section 3)
        stopifnot(any(out$res_global$W > 0, na.rm = TRUE))

        saveRDS(out, file.path(RESULTS_DIR, "rds", paste0("19b_ancombc2_", model_name, BIOSAMPLE_SUFFIX, ".rds")))

        out$res_global |>
                dplyr::mutate(tax_label = make_tax_label(ps_raw, taxon), .after = taxon) |>
                dplyr::arrange(q_val) |>
                readr::write_csv(file.path(RESULTS_DIR, "tables",
                                           paste0("19b_ancombc2_", model_name, "_global", SET_TAG, ".csv")))

        out$res_pair |>
                dplyr::mutate(tax_label = make_tax_label(ps_raw, taxon), .after = taxon) |>
                readr::write_csv(file.path(RESULTS_DIR, "tables",
                                           paste0("19b_ancombc2_", model_name, "_pairwise", SET_TAG, ".csv")))

        cat(sprintf("  global: %d taxa pass FDR, %d also pass the sensitivity filter (of %d)\n",
                    sum(out$res_global$diff_abn, na.rm = TRUE),
                    sum(out$res_global$diff_robust_abn, na.rm = TRUE),
                    nrow(out$res_global)))
        cat(sprintf("  done in %.1f min\n", as.numeric(difftime(Sys.time(), started, units = "mins"))))
}

cat("\n19b done.\n")
