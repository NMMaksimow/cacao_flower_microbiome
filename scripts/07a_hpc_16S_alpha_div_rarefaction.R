# ---- Libraries ----
library(phyloseq)
library(dplyr)    # bind_rows, left_join, mutate, group_by, summarise, across, select
library(tibble)   # tibble(), rownames_to_column()
library(parallel)
# No extra package needed for Faith's PD — see faithPD_fast() below.

# ---- Project paths (Hamilton HPC — no here()) ----
project_dir <- "/nobackup/jcnx71/cacao_flower_microbiome"
rds_dir     <- file.path(project_dir, "results", "rds")

dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Parallelism: inherit core count from SLURM allocation ----
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat(sprintf("Using %d cores\n", N_CORES))

# ---- Parameters ----
RAREFY_DEPTH     <- 2000    # primary depth: retains >99% of samples
RAREFY_DEPTH_10k <- 10000   # secondary depth: some samples drop out
N_ITER           <- 999L    # rarefaction iterations (odd number -> integer median)
BASE_SEED        <- 42

# Farm order: inside_forest -> full_sun (deforestation gradient)
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
MGMT_LEVELS <- c("inside_forest", "near_forest", "agroforest", "full_sun")

# ---- Input ----
cat(sprintf("[%s] Loading phyloseq object...\n", Sys.time()))
ps <- readRDS(file.path(rds_dir, "ps_16S_bacteria_biosamples.rds"))
cat(sprintf("[%s] Loaded: %d samples, %d ASVs\n",
            Sys.time(), nsamples(ps), ntaxa(ps)))

################################################################################
#### 1. Faith's PD — vectorised post-order implementation (no extra packages) ##
################################################################################

# Vectorised Faith's PD for all samples in one pass.
#
# picante::pd() calls ape::drop.tip() once per sample in an R-level loop —
# O(n_samples * n_edges). For 282 samples x ~10k edges that is ~120 seconds
# per rarefaction iteration.
#
# This implementation does one post-order tree traversal over all samples
# simultaneously via matrix operations: O(n_edges * n_samples) with BLAS
# matrix-vector multiply at the end — typically < 1 second per call.
#
# otu_mat : samples x ASVs matrix; column names must match tree$tip.label.
# tree    : ape phylo (rooted); tips must be the same set as otu_mat columns.
# Returns : numeric vector of FaithPD values (one per sample row in otu_mat).
faithPD_fast <- function(otu_mat, tree) {
    tree2  <- ape::reorder.phylo(tree, "postorder")
    n_tip  <- ape::Ntip(tree2)
    n_node <- ape::Nnode(tree2)
    n_samp <- nrow(otu_mat)
    edges  <- tree2$edge  # n_edges x 2: (parent_node, child_node)

    tip_idx <- match(tree2$tip.label, colnames(otu_mat))

    # node_pres[i, n] = 1 if sample i has >= 1 present tip in the subtree of node n
    node_pres <- matrix(0L, nrow = n_samp, ncol = n_tip + n_node)
    node_pres[, seq_len(n_tip)] <- (otu_mat[, tip_idx, drop = FALSE] > 0L) * 1L

    # Post-order: propagate presence from each child node to its parent.
    # Because ape "postorder" ensures child edges precede parent edges,
    # iterating in order fully propagates before any parent is visited.
    for (k in seq_len(nrow(edges))) {
        par_node   <- edges[k, 1]
        child_node <- edges[k, 2]
        node_pres[, par_node] <- pmax(node_pres[, par_node], node_pres[, child_node])
    }

    # FaithPD[i] = sum of edge lengths where the child node of that edge
    # has at least one present descendant in sample i.
    as.vector(node_pres[, edges[, 2], drop = FALSE] %*% tree2$edge.length)
}

################################################################################
#### 2. Helper: single rarefaction + all diversity metrics #####################
################################################################################

# Returns a tibble with one row per sample surviving rarefaction.
# Metrics:
#   Observed     = Hill q=0 (raw ASV count)
#   Shannon      = log-diversity H
#   InvSimpson   = Hill q=2 = 1/sum(p_i^2)
#   Hill_q1      = exp(Shannon) = Hill q=1
#   BergerParker = max(count)/total (dominance of most abundant ASV)
#   Hill_q0_5, Hill_q1_5, Hill_q3, Hill_q5 = intermediate Hill orders
#   FaithPD      = sum of branch lengths to root (requires embedded phylogenetic tree)
# NOTE: Chao1 omitted — downward biased after DADA2 removes global singletons
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

    # OTU matrix in samples x ASVs orientation
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

    # Faith's PD: one vectorised post-order traversal over all samples at once.
    pd_vals <- tibble::tibble(
        sample_id = rownames(otu_mat),
        FaithPD   = faithPD_fast(otu_mat, phyloseq::phy_tree(ps_rare))
    )

    phyloseq::estimate_richness(
        ps_rare,
        measures = c("Observed", "Shannon", "InvSimpson")
    ) |>
        tibble::rownames_to_column("sample_id") |>
        dplyr::mutate(Hill_q1 = exp(Shannon)) |>
        dplyr::left_join(bp_vals,    by = "sample_id") |>
        dplyr::left_join(hill_extra, by = "sample_id") |>
        dplyr::left_join(pd_vals,    by = "sample_id") |>
        dplyr::select(sample_id, Observed, Shannon, InvSimpson,
                      Hill_q1, BergerParker, Hill_q0_5, Hill_q1_5, Hill_q3, Hill_q5,
                      FaithPD)
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
#### phyloseq/ape C-level pointers are safely serialised to each worker.       #
################################################################################

cat(sprintf("[%s] Creating PSOCK cluster with %d workers...\n", Sys.time(), N_CORES))
cl <- parallel::makePSOCKcluster(N_CORES)

parallel::clusterExport(cl, varlist = c("ps", "RAREFY_DEPTH", "RAREFY_DEPTH_10k",
                                        "BASE_SEED", "rarefy_alpha_once", "faithPD_fast"))
parallel::clusterEvalQ(cl, {
    library(phyloseq)
    library(dplyr)
    library(tibble)
    set.seed(42)
})

cat(sprintf("[%s] Cluster ready.\n", Sys.time()))

################################################################################
#### 2a. Parallel rarefaction: 999 iterations at 2000 reads ####################
################################################################################

cat(sprintf("[%s] Starting 2k rarefaction (%d iter, %d workers)...\n",
            Sys.time(), N_ITER, N_CORES))

alpha_iters_2k <- parallel::parLapply(
    cl,
    seq_len(N_ITER),
    function(i) {
        rarefy_alpha_once(ps, RAREFY_DEPTH, BASE_SEED + i) |>
            dplyr::mutate(iteration = i)
    }
) |>
    dplyr::bind_rows() |>
    dplyr::left_join(data.frame(phyloseq::sample_data(ps)), by = "sample_id")

saveRDS(alpha_iters_2k, file.path(rds_dir, "alpha_16S_2k_999iters.rds"))
cat(sprintf("[%s] Saved alpha_16S_2k_999iters.rds\n", Sys.time()))

################################################################################
#### 2b. Aggregate 2k: median + SD per sample across 999 iterations ############
################################################################################

alpha_summary_2k <- alpha_iters_2k |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(
        dplyr::across(
            c(Observed, Shannon, InvSimpson, Hill_q1, BergerParker,
              Hill_q0_5, Hill_q1_5, Hill_q3, Hill_q5, FaithPD),
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

saveRDS(alpha_summary_2k, file.path(rds_dir, "alpha_16S_2k_999iters_summary.rds"))
cat(sprintf("[%s] Saved alpha_16S_2k_999iters_summary.rds\n", Sys.time()))

################################################################################
#### 3a. Parallel rarefaction: 999 iterations at 10000 reads ###################
#### NB: at 10k, some samples drop out.                                        #
#### Results for comparison only; primary analysis uses 2k.                    #
################################################################################

cat(sprintf("[%s] Starting 10k rarefaction (%d iter, %d workers)...\n",
            Sys.time(), N_ITER, N_CORES))

alpha_iters_10k <- parallel::parLapply(
    cl,
    seq_len(N_ITER),
    function(i) {
        rarefy_alpha_once(ps, RAREFY_DEPTH_10k, BASE_SEED + i) |>
            dplyr::mutate(iteration = i)
    }
) |>
    dplyr::bind_rows() |>
    dplyr::left_join(data.frame(phyloseq::sample_data(ps)), by = "sample_id")

saveRDS(alpha_iters_10k, file.path(rds_dir, "alpha_16S_10k_999iters.rds"))
cat(sprintf("[%s] Saved alpha_16S_10k_999iters.rds\n", Sys.time()))

################################################################################
#### 3b. Aggregate 10k: median + SD per sample across 999 iterations ###########
################################################################################

alpha_summary_10k <- alpha_iters_10k |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(
        dplyr::across(
            c(Observed, Shannon, InvSimpson, Hill_q1, BergerParker,
              Hill_q0_5, Hill_q1_5, Hill_q3, Hill_q5, FaithPD),
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

saveRDS(alpha_summary_10k, file.path(rds_dir, "alpha_16S_10k_999iters_summary.rds"))
cat(sprintf("[%s] Saved alpha_16S_10k_999iters_summary.rds\n", Sys.time()))

parallel::stopCluster(cl)
cat(sprintf("[%s] All done.\n", Sys.time()))
