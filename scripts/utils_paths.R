## ============================================================================
## utils_paths — data version switch (v1 / v2) and the folders it selects
##
## Sourced by every analysis script.
##
## CFM_VERSION is an environment variable, "v1" or "v2"; unset means v2.
##   RStudio:   Sys.setenv(CFM_VERSION = "v1") before sourcing a script
##   terminal:  WSLENV=CFM_VERSION CFM_VERSION=v1 Rscript.exe scripts/<script>.R
##   SLURM:     export CFM_VERSION=v1 in the sbatch script
##
## v1 = qiime2/export (QIIME 2 2024.10 processing); v2 = qiime2_v2/export
## (Snakemake v2 rerun). Results of each version go to results/<version>/.
## ============================================================================

CFM_VERSION <- Sys.getenv("CFM_VERSION", unset = "v2")
stopifnot(CFM_VERSION %in% c("v1", "v2"))

# QIIME 2 exports of this version
EXPORT_DIR <- here::here(c(v1 = "qiime2", v2 = "qiime2_v2")[[CFM_VERSION]], "export")

# Results of this version. CFM_RESULTS_DIR is for verification runs only, so a v1
# check never writes into the frozen results/v1
RESULTS_DIR <- Sys.getenv("CFM_RESULTS_DIR", unset = here::here("results", CFM_VERSION))

message("Data version ", CFM_VERSION, ": ", EXPORT_DIR, " -> ", RESULTS_DIR)
