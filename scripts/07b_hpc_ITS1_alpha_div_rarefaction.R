# ---- Libraries ----
library(phyloseq)
library(dplyr)    # bind_rows, left_join, mutate, group_by, summarise, across, select
library(tibble)   # tibble(), rownames_to_column()
library(parallel)

# ---- Project paths (Hamilton HPC — no here()) ----
project_dir <- "/nobackup/jcnx71/cacao_flower_microbiome"
rds_dir     <- file.path(project_dir, "results", "rds")

dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Parallelism: inherit core count from SLURM allocation ----
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat(sprintf("Using %d cores\n", N_CORES))

# ---- Parameters ----
RAREFY_DEPTH     <- 5000    # primary depth: retains >99% of samples
RAREFY_DEPTH_15k <- 15000   # secondary depth: ~17 bagged + 1 unbagged drop out
N_ITER           <- 999L    # rarefaction iterations (odd number -> integer median)
BASE_SEED        <- 42

# Farm order: inside_forest -> full_sun (deforestation gradient)
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
MGMT_LEVELS <- c("inside_forest", "near_forest", "agroforest", "full_sun")

# ---- Input ----
cat(sprintf("[%s] Loading phyloseq object...\n", Sys.time()))
ps <- readRDS(file.path(rds_dir, "ps_ITS1_otu97_fungi_biosamples.rds"))
cat(sprintf("[%s] Loaded: %d samples, %d OTUs\n",
            Sys.time(), nsamples(ps), ntaxa(ps)))

################################################################################
#### 1. Helper: single rarefaction + all diversity metrics #####################
################################################################################

# Returns a tibble with one row per sample surviving rarefaction.
# Metrics:
#   Observed     = Hill q=0 (raw OTU count)
#   Shannon      = log-diversity H
#   InvSimpson   = Hill q=2 = 1/sum(p_i^2)
#   Hill_q1      = exp(Shannon) = Hill q=1
#   BergerParker = max(count)/total (dominance of most abundant OTU)
#   Hill_q0_5, Hill_q1_5, Hill_q3, Hill_q5 = intermediate Hill orders
# NOTE: Chao1 omitted — downward biased after DADA2 removes global singletons
# NOTE: FaithPD not computed — ITS1 tree unreliable due to high length variability
rarefy_alpha_once <- function(ps, depth, seed) {

    ps_rare <- suppressMessages(
        phyloseq::rarefy_even_depth(
            ps,
            sample.size = depth,
            rngseed     = seed,
            replace     = FALSE,
            trimOTUs    = TRUE,
            verbose     = FALSE
        )
    )

    # OTU matrix in samples x OTUs orientation
    otu_mat <- phyloseq::otu_table(ps_rare) |>
        as("matrix") |>
        (\(m) if (phyloseq::taxa_are_rows(ps_rare)) t(m) else m)()

    # Proportional abundances for Hill numbers at non-integer q
    prop_mat <- otu_mat / rowSums(otu_mat)

    # General Hill number D_q = (sum p_i^q)^(1/(1-q)) for q != 0, 1, inf
    hill_q_general <- function(p, q) sum(p[p > 0]^q)^(1 / (1 - q))

    bp_vals <- tibble::tibble(
        sample_id    = rownames(otu_mat),
        BergerParker = apply(otu_mat, 1, \(x) max(x) / sum(x))
    )

    hill_extra <- tibble::tibble(
        sample_id = rownames(prop_mat),
        Hill_q0_5 = apply(prop_mat, 1, \(p) hill_q_general(p, 0.5)),
        Hill_q1_5 = apply(prop_mat, 1, \(p) hill_q_general(p, 1.5)),
        Hill_q3   = apply(prop_mat, 1, \(p) hill_q_general(p, 3)),
        Hill_q5   = apply(prop_mat, 1, \(p) hill_q_general(p, 5))
    )

    phyloseq::estimate_richness(
        ps_rare,
        measures = c("Observed", "Shannon", "InvSimpson")
    ) |>
        tibble::rownames_to_column("sample_id") |>
        dplyr::mutate(Hill_q1 = exp(Shannon)) |>
        dplyr::left_join(bp_vals,    by = "sample_id") |>
        dplyr::left_join(hill_extra, by = "sample_id") |>
        dplyr::select(sample_id, Observed, Shannon, InvSimpson,
                      Hill_q1, BergerParker, Hill_q0_5, Hill_q1_5, Hill_q3, Hill_q5)
}

################################################################################
#### 1b. Test single iteration before committing to the full run ###############
################################################################################

# Initialise .Random.seed in the main process — required by rarefy_even_depth()
set.seed(BASE_SEED)

cat(sprintf("[%s] Testing single iteration...\n", Sys.time()))
test_result <- tryCatch(
    rarefy_alpha_once(ps, RAREFY_DEPTH, BASE_SEED + 1L),
    error = function(e) {
        cat(sprintf("FATAL: rarefy_alpha_once() failed: %s\n", conditionMessage(e)))
        stop(e)
    }
)
cat(sprintf("[%s] Test passed (%d samples).\n", Sys.time(), nrow(test_result)))
rm(test_result)

################################################################################
#### 1c. PSOCK cluster setup ###################################################
#### PSOCK starts fresh worker processes via sockets — no fork(), so           #
#### phyloseq objects are safely serialised to each worker.                    #
################################################################################

cat(sprintf("[%s] Creating PSOCK cluster with %d workers...\n", Sys.time(), N_CORES))
cl <- parallel::makePSOCKcluster(N_CORES)

parallel::clusterExport(cl, varlist = c("ps", "RAREFY_DEPTH", "RAREFY_DEPTH_15k",
                                        "BASE_SEED", "rarefy_alpha_once"))
parallel::clusterEvalQ(cl, {
    library(phyloseq)
    library(dplyr)
    library(tibble)
    set.seed(42)
})

cat(sprintf("[%s] Cluster ready.\n", Sys.time()))

################################################################################
#### 2a. Parallel rarefaction: 999 iterations at 5000 reads ####################
################################################################################

cat(sprintf("[%s] Starting 5k rarefaction (%d iter, %d workers)...\n",
            Sys.time(), N_ITER, N_CORES))

alpha_iters_5k <- parallel::parLapply(
    cl,
    seq_len(N_ITER),
    function(i) {
        rarefy_alpha_once(ps, RAREFY_DEPTH, BASE_SEED + i) |>
            dplyr::mutate(iteration = i)
    }
) |>
    dplyr::bind_rows() |>
    dplyr::left_join(data.frame(phyloseq::sample_data(ps)), by = "sample_id")

saveRDS(alpha_iters_5k, file.path(rds_dir, "alpha_ITS1_otu97_5k_999iters.rds"))
cat(sprintf("[%s] Saved alpha_ITS1_otu97_5k_999iters.rds\n", Sys.time()))

################################################################################
#### 2b. Aggregate 5k: median + SD per sample across 999 iterations ############
################################################################################

alpha_summary_5k <- alpha_iters_5k |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(
        dplyr::across(
            c(Observed, Shannon, InvSimpson, Hill_q1, BergerParker,
              Hill_q0_5, Hill_q1_5, Hill_q3, Hill_q5),
            list(
                median = \(x) median(x, na.rm = TRUE),
                sd     = \(x) sd(x, na.rm = TRUE)
            )
        ),
        .groups = "drop"
    ) |>
    dplyr::left_join(
        data.frame(phyloseq::sample_data(ps)),
        by = "sample_id"
    ) |>
    dplyr::mutate(
        farm_id         = factor(farm_id,         levels = FARM_LEVELS),
        management_type = factor(management_type, levels = MGMT_LEVELS),
        sample_type     = factor(sample_type)
    )

saveRDS(alpha_summary_5k, file.path(rds_dir, "alpha_ITS1_otu97_5k_999iters_summary.rds"))
cat(sprintf("[%s] Saved alpha_ITS1_otu97_5k_999iters_summary.rds\n", Sys.time()))

################################################################################
#### 3a. Parallel rarefaction: 999 iterations at 15000 reads ###################
#### NB: at 15k, ~17 bagged + 1 unbagged samples drop out.                    #
#### Results for comparison only; primary analysis uses 5k.                    #
################################################################################

cat(sprintf("[%s] Starting 15k rarefaction (%d iter, %d workers)...\n",
            Sys.time(), N_ITER, N_CORES))

alpha_iters_15k <- parallel::parLapply(
    cl,
    seq_len(N_ITER),
    function(i) {
        rarefy_alpha_once(ps, RAREFY_DEPTH_15k, BASE_SEED + i) |>
            dplyr::mutate(iteration = i)
    }
) |>
    dplyr::bind_rows() |>
    dplyr::left_join(data.frame(phyloseq::sample_data(ps)), by = "sample_id")

saveRDS(alpha_iters_15k, file.path(rds_dir, "alpha_ITS1_otu97_15k_999iters.rds"))
cat(sprintf("[%s] Saved alpha_ITS1_otu97_15k_999iters.rds\n", Sys.time()))

################################################################################
#### 3b. Aggregate 15k: median + SD per sample across 999 iterations ###########
################################################################################

alpha_summary_15k <- alpha_iters_15k |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(
        dplyr::across(
            c(Observed, Shannon, InvSimpson, Hill_q1, BergerParker,
              Hill_q0_5, Hill_q1_5, Hill_q3, Hill_q5),
            list(
                median = \(x) median(x, na.rm = TRUE),
                sd     = \(x) sd(x, na.rm = TRUE)
            )
        ),
        .groups = "drop"
    ) |>
    dplyr::left_join(
        data.frame(phyloseq::sample_data(ps)),
        by = "sample_id"
    ) |>
    dplyr::mutate(
        farm_id         = factor(farm_id,         levels = FARM_LEVELS),
        management_type = factor(management_type, levels = MGMT_LEVELS),
        sample_type     = factor(sample_type)
    )

saveRDS(alpha_summary_15k, file.path(rds_dir, "alpha_ITS1_otu97_15k_999iters_summary.rds"))
cat(sprintf("[%s] Saved alpha_ITS1_otu97_15k_999iters_summary.rds\n", Sys.time()))

parallel::stopCluster(cl)
cat(sprintf("[%s] All done.\n", Sys.time()))
