#!/bin/bash
################################################################################
# setup_hamilton.sh — one-time setup on Hamilton HPC
# Run this ONCE after first login to configure R personal library
# and install required packages.
#
# Usage (from project root — run directly, NOT via sbatch):
#   bash scripts/setup_hamilton.sh
################################################################################

if [[ -n "$SLURM_JOB_ID" ]]; then
    echo "ERROR: run this script directly with 'bash', not via sbatch"
    exit 1
fi

echo "=== Hamilton HPC: one-time R setup ==="
echo "Date: $(date)"

# ---- 1. Load R module ----
module load r/4.5.2

# ---- 2. Create personal R library ----
mkdir -p ~/R/library
echo "Personal R library: ~/R/library"

# ---- 3. Set R_LIBS_USER in ~/.Renviron ----
# This tells R to look in ~/R/library for packages at every session
if ! grep -q "R_LIBS_USER" ~/.Renviron 2>/dev/null; then
    echo 'R_LIBS_USER=~/R/library' >> ~/.Renviron
    echo "R_LIBS_USER added to ~/.Renviron"
else
    echo "R_LIBS_USER already set in ~/.Renviron — skipping"
fi

# ---- 4. Install CRAN packages ----
echo ""
echo "Installing CRAN packages..."
R --no-save -e "
install.packages(
    c('zetadiv', 'dplyr', 'ggplot2', 'here'),
    repos  = 'https://cloud.r-project.org',
    lib    = '~/R/library'
)
"

# ---- 5. Install Bioconductor packages ----
echo ""
echo "Installing Bioconductor packages (phyloseq)..."
R --no-save -e "
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager',
                     repos = 'https://cloud.r-project.org',
                     lib   = '~/R/library')
BiocManager::install('phyloseq', lib = '~/R/library', update = FALSE, ask = FALSE)
"

# ---- 6. Verify ----
echo ""
echo "=== Verifying installations ==="
R --no-save -e "
pkgs <- c('phyloseq', 'zetadiv', 'dplyr', 'ggplot2', 'here')
for (p in pkgs) {
    ok <- requireNamespace(p, quietly = TRUE)
    cat(sprintf('  %-12s  %s\n', p, ifelse(ok, 'OK', 'MISSING')))
}
"

echo ""
echo "Setup complete."
