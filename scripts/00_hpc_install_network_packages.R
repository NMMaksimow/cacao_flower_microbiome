# One-time installation of all packages needed for network analysis scripts
# (19a/b, 20a/b, 21a/b) on HPC Hamilton.
# Run from the login node — compute nodes may not have outbound internet access.
# Usage: module load r/4.6.0 && Rscript scripts/00_hpc_install_network_packages.R

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ── CRAN packages ─────────────────────────────────────────────────────────────
cran_pkgs <- c(
    "here", "dplyr", "tidyr", "readr", "tibble", "purrr", "stringr",
    "ggplot2", "ggraph", "igraph", "patchwork", "ggsci",
    "remotes", "BiocManager"
)
to_install <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(to_install)) {
    cat(sprintf("Installing CRAN packages: %s\n", paste(to_install, collapse = ", ")))
    install.packages(to_install)
} else {
    cat("CRAN packages: all present\n")
}

# ── Bioconductor packages ─────────────────────────────────────────────────────
# R 4.6 requires Bioconductor 3.23; ensure BiocManager is up to date.
if (BiocManager::version() < "3.23") {
    cat("Updating BiocManager to version 3.23 (required for R 4.6)...\n")
    BiocManager::install(version = "3.23", ask = FALSE)
}

bioc_pkgs   <- c("phyloseq", "Biostrings", "BiocParallel")
bioc_missing <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(bioc_missing)) {
    cat(sprintf("Installing Bioconductor packages: %s\n", paste(bioc_missing, collapse = ", ")))
    BiocManager::install(bioc_missing, ask = FALSE)
} else {
    cat("Bioconductor packages: all present\n")
}

# ── GitHub packages ───────────────────────────────────────────────────────────

# pulsar: STARS stability selection used by SpiecEasi.
# Not available on CRAN for conda R builds — install from GitHub.
if (!requireNamespace("pulsar", quietly = TRUE)) {
    cat("Installing pulsar from GitHub (zdk123/pulsar)...\n")
    remotes::install_github("zdk123/pulsar")
} else {
    cat("pulsar: present\n")
}

# SpiecEasi 2.0.0 from Bioconductor (NetCoMi 1.3.0 requires >= 2.0.0).
# GitHub version zdk123/SpiecEasi is 1.99.0 — too old.
needs_spieceasi <- !requireNamespace("SpiecEasi", quietly = TRUE) ||
                   packageVersion("SpiecEasi") < "2.0.0"
if (needs_spieceasi) {
    cat("Installing SpiecEasi from Bioconductor (>= 2.0.0 required)...\n")
    BiocManager::install("SpiecEasi", ask = FALSE)
} else {
    cat(sprintf("SpiecEasi: present (%s)\n", packageVersion("SpiecEasi")))
}

# SPRING: rank-based SpiecEasi variant; same Suggests caveat as above.
if (!requireNamespace("SPRING", quietly = TRUE)) {
    cat("Installing SPRING from GitHub (GraceYoon/SPRING)...\n")
    remotes::install_github("GraceYoon/SPRING")
} else {
    cat("SPRING: present\n")
}

# NetCoMi >= 1.3.0: wraps SpiecEasi/SparCC/SPRING + netCompare permutation tests.
# upgrade = "never" avoids silently bumping Bioconductor packages across a release boundary.
needs_netcomi <- !requireNamespace("NetCoMi", quietly = TRUE) ||
                 packageVersion("NetCoMi") < "1.3.0"
if (needs_netcomi) {
    cat("Installing NetCoMi from GitHub (stefpeschel/NetCoMi, ref = main)...\n")
    remotes::install_github("stefpeschel/NetCoMi",
                            ref          = "main",
                            dependencies = TRUE,
                            upgrade      = "never")
} else {
    cat(sprintf("NetCoMi: present (%s)\n", packageVersion("NetCoMi")))
}

# ── Verify ────────────────────────────────────────────────────────────────────
pkgs_check <- c("NetCoMi", "SpiecEasi", "SPRING", "phyloseq",
                "ggraph", "patchwork", "ggsci", "here")
ok  <- pkgs_check[vapply(pkgs_check, requireNamespace, logical(1), quietly = TRUE)]
bad <- setdiff(pkgs_check, ok)

cat("\n────────────────────────────────────────\n")
cat(sprintf("Installed OK (%d): %s\n", length(ok), paste(ok, collapse = ", ")))
if (length(bad))
    cat(sprintf("FAILED      (%d): %s\n", length(bad), paste(bad, collapse = ", ")))
cat("────────────────────────────────────────\n")
