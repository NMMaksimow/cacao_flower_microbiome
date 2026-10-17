#!/bin/bash
# 00_hpc_r_setup.sh — one-time R package setup on Hamilton HPC
# Run directly (NOT via sbatch) from the project root:
#   bash scripts/00_hpc_r_setup.sh

if [[ -n "$SLURM_JOB_ID" ]]; then
    echo "ERROR: run with 'bash', not 'sbatch'"
    exit 1
fi

echo "=== Hamilton HPC: R package setup ==="
echo "Date: $(date)"

module load r/4.5.1

# ---- Personal R library --------------------------------------------------
mkdir -p ~/R/library
if ! grep -q "R_LIBS_USER" ~/.Renviron 2>/dev/null; then
    echo 'R_LIBS_USER=~/R/library' >> ~/.Renviron
    echo "R_LIBS_USER added to ~/.Renviron"
else
    echo "R_LIBS_USER already set — skipping"
fi

# ---- CRAN packages -------------------------------------------------------
echo ""
echo "--- Installing CRAN packages ---"
R --no-save -e "
# scales is reinstalled explicitly: the Hamilton system R 4.5.1 library has
# scales but without its transitive deps (farver, RColorBrewer, labeling...).
# Reinstalling into ~/R/library pulls all deps and prevents cascade failures
# when zetadiv loads cowplot/forecast -> scales -> missing dep.
install.packages(
    c('scales', 'zetadiv', 'dplyr', 'tibble', 'here'),
    repos = 'https://cloud.r-project.org',
    lib   = '~/R/library'
)
"

# ---- Bioconductor packages -----------------------------------------------
echo ""
echo "--- Installing Bioconductor packages ---"
R --no-save -e "
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager',
                     repos = 'https://cloud.r-project.org',
                     lib   = '~/R/library')
# BiocGenerics installed explicitly first — it is a core Bioc dependency
# that caused 'no package called BiocGenerics' errors in the past.
BiocManager::install(c('BiocGenerics', 'phyloseq'),
                     lib    = '~/R/library',
                     update = FALSE,
                     ask    = FALSE,
                     force  = TRUE)
"

# ---- Verify --------------------------------------------------------------
echo ""
echo "--- Verifying ---"
R --no-save -e "
pkgs <- c('phyloseq', 'zetadiv', 'dplyr', 'tibble', 'here')
for (p in pkgs) {
    ok <- requireNamespace(p, quietly = TRUE)
    cat(sprintf('  %-12s  %s\n', p, ifelse(ok, 'OK', 'MISSING')))
}
"

echo ""
echo "Setup complete. Submit zeta jobs with:"
echo "  sbatch scripts/11a_hpc_16S_zeta_by_farm.sh"
echo "  sbatch scripts/11b_hpc_ITS1_zeta_by_farm.sh"

# ===========================================================================
# Part 2 — Network analysis packages (19a/b, 20a/b, 21a/b)
# Run this section after Part 1 above is complete.
# NetCoMi and SpiecEasi are GitHub packages; remotes is required.
# Installation order matters: Bioc deps → SpiecEasi → NetCoMi.
# ===========================================================================

echo ""
echo "=== Part 2: Network analysis packages ==="
echo "Date: $(date)"

# ---- CRAN packages for network analysis ------------------------------------
echo ""
echo "--- Installing CRAN network packages ---"
R --no-save -e "
install.packages(
    c('remotes',    # needed to install GitHub packages below
      'igraph',     # graph data structure and algorithms
      'ggraph',     # ggplot2-based network visualisation
      'patchwork',  # multi-panel figure assembly
      'SPRING',     # semi-parametric rank-based graphical model (21a/21b)
      'huge',       # high-dimensional undirected graph estimation (NetCoMi dep)
      'pulsar'),    # StARS model selection (SpiecEasi dep)
    repos = 'https://cloud.r-project.org',
    lib   = '~/R/library'
)
"

# ---- Bioconductor dependencies of NetCoMi ----------------------------------
echo ""
echo "--- Installing Bioconductor dependencies of NetCoMi ---"
R --no-save -e "
BiocManager::install(
    c('limma', 'edgeR', 'DESeq2', 'metagenomeSeq', 'BiocParallel'),
    lib    = '~/R/library',
    update = FALSE,
    ask    = FALSE
)
"

# ---- GitHub: SpiecEasi -----------------------------------------------------
echo ""
echo "--- Installing SpiecEasi from GitHub ---"
R --no-save -e "
remotes::install_github(
    'zdk123/SpiecEasi',
    lib          = '~/R/library',
    upgrade      = 'never',
    dependencies = TRUE,
    repos        = 'https://cloud.r-project.org'
)
"

# ---- GitHub: NetCoMi -------------------------------------------------------
echo ""
echo "--- Installing NetCoMi from GitHub ---"
R --no-save -e "
remotes::install_github(
    'stefpeschel/NetCoMi',
    lib          = '~/R/library',
    upgrade      = 'never',
    dependencies = TRUE,
    repos        = 'https://cloud.r-project.org'
)
"

# ---- Verify ----------------------------------------------------------------
echo ""
echo "--- Verifying network packages ---"
R --no-save -e "
pkgs <- c('igraph', 'ggraph', 'patchwork', 'SPRING', 'SpiecEasi', 'NetCoMi')
for (p in pkgs) {
    ok <- requireNamespace(p, quietly = TRUE)
    cat(sprintf('  %-12s  %s\n', p, ifelse(ok, 'OK', 'MISSING')))
}
"

echo ""
echo "Part 2 complete. Run network scripts with:"
echo "  Rscript scripts/19a_16S_network_by_farm.R"
echo "  Rscript scripts/21b_ITS1_network_by_farm_SPRING.R"
