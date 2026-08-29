## ============================================================================
## 22 — Timing probe: how expensive is netCompare with permutations?
##
## The first SparCC/SPRING runs had permTest = FALSE, so netCompare cost nothing
## and the jobs finished in ~2 h. With permTest = TRUE the permutation test
## re-estimates the association matrix nPerm times per farm, which dominates
## runtime. This script measures the cost on ONE farm with a small nPerm so the
## real jobs can be sized instead of guessed.
##
## Input:
##   results/rds/20a_network_results.rds   (SparCC, 16S — already on HPC)
##   results/rds/21a_network_results.rds   (SPRING, 16S — already on HPC)
## Output:
##   logs/network/22_netcompare_timing.txt
## ============================================================================

N_PERM_PROBE <- 10L    # small probe; runtime is extrapolated linearly from this
PROBE_FARM   <- "ib"   # IB has the largest networks (123 taxa), so this is worst case
PROBE_CORES  <- c(1L, 4L)   # measure the parallel speedup as well

library(here)
library(NetCoMi)

sink(here("logs", "network", "22_netcompare_timing.txt"), split = TRUE)

cat("netCompare timing probe\n")
cat("=======================\n")
cat(sprintf("farm = %s   nPerm = %d   cores tested = %s\n\n",
            PROBE_FARM, N_PERM_PROBE, paste(PROBE_CORES, collapse = ", ")))

probe <- function(rds, label) {
        f <- here("results", "rds", rds)
        if (!file.exists(f)) {
                cat(sprintf("%-10s SKIPPED (%s not found)\n", label, rds))
                return(invisible(NULL))
        }
        r <- readRDS(f)[[PROBE_FARM]]
        if (is.null(r) || is.null(r$analyzed)) {
                cat(sprintf("%-10s SKIPPED (no analyzed network for %s)\n", label, PROBE_FARM))
                return(invisible(NULL))
        }
        n_taxa <- ncol(r$raw$adjaMat1)

        for (cs in PROBE_CORES) {
                el <- tryCatch(
                        system.time(
                                netCompare(r$analyzed, permTest = TRUE, nPerm = N_PERM_PROBE,
                                           cores = cs, verbose = FALSE, seed = 42)
                        )[["elapsed"]],
                        error = function(e) { cat("   FAILED: ", conditionMessage(e), "\n"); NA_real_ }
                )
                if (is.na(el)) next

                per_perm <- el / N_PERM_PROBE
                # 7 farms x 2 markers, both run inside one Slurm job
                total_h  <- per_perm * 1000 * 7 * 2 / 3600
                cat(sprintf("%-10s taxa=%3d  cores=%d  %6.1f s for %d perms  |  %5.2f s/perm  |  projected %5.1f h for 1000 perms x 7 farms x 2 markers\n",
                            label, n_taxa, cs, el, N_PERM_PROBE, per_perm, total_h))
        }
        cat("\n")
}

probe("20a_network_results.rds", "SparCC")
probe("21a_network_results.rds", "SPRING")

cat("Projection assumes runtime scales linearly with nPerm and that all farms\n")
cat("cost about as much as the largest one, so it is an upper bound.\n")

sink()
