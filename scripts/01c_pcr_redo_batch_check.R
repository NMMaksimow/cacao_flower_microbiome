## ============================================================================
## 01c_pcr_redo_batch_check — is the PCR redo confounded with farm or treatment?
##
## Failed PCR replicates were repeated on another day and pooled with the rest
## of the sample (pcr_redo_16s, pcr_redo_its1: number of redone replicates; mocks
## are single PCRs). If redone samples cluster in one farm or treatment, a PCR-day
## effect could look like a farm or treatment effect. Tested here on the metadata
## (redo = at least one redone replicate); the effect on the communities is
## tested as a PERMANOVA term later (10a/10b).
##
## Expected cell counts are far below 5, so Fisher's exact test is used instead
## of chi-square (exact for 2 x 2; Monte Carlo p-value for larger tables).
## Cramer's V gives the strength of the association: with ~10 redo samples the
## test has little power, so V is read next to p.
##
## Input:
##   data/sample_metadata_with_microscopy.txt
## Output:
##   results/<version>/tables/01c_pcr_redo_counts.csv
##   results/<version>/tables/01c_pcr_redo_association.csv
## ============================================================================

library(here)
library(tidyverse)
library(DescTools)  # CramerV()

# Data version switch: CFM_VERSION, EXPORT_DIR, RESULTS_DIR
source(here("scripts", "utils_paths.R"))
# read_sample_metadata(), EXCLUDE_DEVIATIONS
source(here("scripts", "utils_read_metadata.R"))

# ── Parameters ────────────────────────────────────────────────────────────────
N_SIM <- 1e5   # Monte Carlo replicates for Fisher p-values of tables > 2 x 2
SEED  <- 42

# Metadata column of each marker
REDO_COLUMNS <- c(`16S` = "pcr_redo_16s", ITS1 = "pcr_redo_its1")

# Factors tested against redo
FACTORS <- c("farm_id", "sample_type", "farm_treatment")

# ── Flowers: 288 main set and 294 with deviations ─────────────────────────────
flowers <- read_sample_metadata() |>
        filter(sample_type %in% c("bagged_flower", "unbagged_flower")) |>
        mutate(farm_treatment = paste(farm_id, sample_type, sep = ":"))

sets <- list(
        main_288            = filter(flowers, !sampling_deviation %in% EXCLUDE_DEVIATIONS),
        with_deviations_294 = flowers
)

# ── Counts per farm x treatment ───────────────────────────────────────────────
redo_counts <- sets |>
        imap(\(d, set) d |>
                     group_by(farm_id, sample_type) |>
                     summarise(n_flowers      = n(),
                               n_redo_16S     = sum(pcr_redo_16s >= 1),
                               n_redo_ITS1    = sum(pcr_redo_its1 >= 1),
                               .groups = "drop") |>
                     mutate(set = set, .before = 1)) |>
        bind_rows()

# ── Association tests ─────────────────────────────────────────────────────────
# One factor against redo (yes/no) in one set of flowers
test_redo <- function(d, redo_column, factor) {
        redo <- d[[redo_column]] >= 1
        tab  <- table(d[[factor]], redo)
        set.seed(SEED)
        fisher <- if (all(dim(tab) == 2)) {
                fisher.test(tab)
        } else {
                fisher.test(tab, simulate.p.value = TRUE, B = N_SIM)
        }
        tibble(n         = sum(tab),
               n_redo    = sum(redo),
               n_levels  = nrow(tab),
               fisher_p  = fisher$p.value,
               cramers_v = CramerV(tab))
}

redo_association <- expand_grid(set    = names(sets),
                                marker = names(REDO_COLUMNS),
                                factor = FACTORS) |>
        mutate(result = pmap(list(set, marker, factor),
                             \(s, m, f) test_redo(sets[[s]], REDO_COLUMNS[[m]], f))) |>
        unnest(result)

# ── Output ────────────────────────────────────────────────────────────────────
dir.create(file.path(RESULTS_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

write_csv(redo_counts,
          file.path(RESULTS_DIR, "tables", "01c_pcr_redo_counts.csv"))
write_csv(redo_association,
          file.path(RESULTS_DIR, "tables", "01c_pcr_redo_association.csv"))

cat("\nFlowers with a PCR redo, per farm x treatment\n")
redo_counts |> filter(set == "main_288") |> print(n = Inf)

cat("\nRedo vs farm and treatment (Fisher's exact test, Cramer's V)\n")
redo_association |>
        mutate(across(c(fisher_p, cramers_v), \(x) round(x, 3))) |>
        print(n = Inf)
