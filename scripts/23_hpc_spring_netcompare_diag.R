## ============================================================================
## 23 — Diagnostic: why does netCompare fail for SPRING networks?
##
## The timing probe (script 22) hit this on farm IB for both 16S and ITS1:
##   Warning: 'nearPD()' did not converge in 100 iterations
##   Error:   task 3 failed - "NAs are not allowed in subscripted assignments"
##
## netCompare's permutation test re-estimates the association matrix on permuted
## group labels. For SPRING that goes through a rank correlation that must be
## projected to the nearest positive-definite matrix; when nearPD does not
## converge the result carries NAs and the downstream assignment fails.
##
## Questions this answers:
##   1. Does it fail on every farm, or only some?
##   2. Does it fail for both markers?
##   3. Do SparCC networks from the same farms pass? (isolates SPRING as the cause)
##
## Input:  results/rds/{20a,20b,21a,21b}_network_results.rds
## Output: logs/network/23_spring_netcompare_diag.txt
## ============================================================================

N_PERM_DIAG <- 10L   # tiny: we are testing for failure, not measuring runtime
FARMS       <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

library(here)
library(NetCoMi)

sink(here("logs", "network", "23_spring_netcompare_diag.txt"), split = TRUE)

cat("SPRING netCompare diagnostic\n")
cat("============================\n")
cat(sprintf("nPerm = %d, cores = 1, all %d farms\n\n", N_PERM_DIAG, length(FARMS)))

# Returns "OK", "SKIP", or the error message, plus whether nearPD warned
try_compare <- function(net_ana) {
        warned <- FALSE
        res <- withCallingHandlers(
                tryCatch({
                        netCompare(net_ana, permTest = TRUE, nPerm = N_PERM_DIAG,
                                   cores = 1, verbose = FALSE, seed = 42)
                        "OK"
                }, error = function(e) conditionMessage(e)),
                warning = function(w) {
                        if (grepl("nearPD", conditionMessage(w))) warned <<- TRUE
                        invokeRestart("muffleWarning")
                }
        )
        list(status = res, nearpd = warned)
}

run_set <- function(rds, label) {
        f <- here("results", "rds", rds)
        if (!file.exists(f)) {
                cat(sprintf("%-18s FILE MISSING (%s)\n\n", label, rds)); return(invisible())
        }
        results <- readRDS(f)
        cat(sprintf("── %s ──\n", label))
        for (farm in FARMS) {
                r <- results[[farm]]
                if (is.null(r) || is.null(r$analyzed)) {
                        cat(sprintf("  %-4s SKIP (no analyzed network)\n", toupper(farm))); next
                }
                out <- try_compare(r$analyzed)
                cat(sprintf("  %-4s taxa=%3d  nearPD_warning=%-5s  %s\n",
                            toupper(farm), ncol(r$raw$adjaMat1),
                            out$nearpd,
                            if (out$status == "OK") "OK" else paste("FAIL:", out$status)))
        }
        cat("\n")
}

# SPRING first (the suspect), then SparCC on the same farms as a control
run_set("21a_network_results.rds", "SPRING 16S")
run_set("21b_network_results.rds", "SPRING ITS1")
run_set("20a_network_results.rds", "SparCC 16S (control)")
run_set("20b_network_results.rds", "SparCC ITS1 (control)")

cat("If SPRING fails everywhere and SparCC passes everywhere, the permutation\n")
cat("path for rank-based association estimates is the problem, not the data.\n")

sink()
