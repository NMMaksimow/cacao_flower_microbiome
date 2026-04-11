#!/bin/bash
#SBATCH --job-name=qiime_phylogeny
#SBATCH -p shared
#SBATCH --mem=32G
#SBATCH -c 16
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --output=logs/12_qiime2_phylogeny_16S.out
#SBATCH --error=logs/12_qiime2_phylogeny_16S.err
#SBATCH --mail-type=FAIL

# Script to build a phylogenetic tree for 16S data
# Step 12: de novo phylogeny via MAFFT alignment + FastTree
# Input:  qiime2/denoise/CFM_16S_dada2_repseqs.qza
# Output: qiime2/phylogeny/CFM_16S_rooted_tree.qza  (+ intermediates)
#
# Note: ITS1 data is intentionally excluded — ITS1 marker is too variable for reliable
# phylogenetics and UniFrac / Faith's PD are not meaningful for fungi.

echo "Starting 16S phylogenetic tree construction"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 32G"

module load bioinformatics
module load qiime2/2024.10amplicon

REPSEQS="qiime2/denoise/CFM_16S_dada2_repseqs.qza"
PHYLO_DIR="qiime2/phylogeny"
mkdir -p ${PHYLO_DIR}

# ------------------------------------------------------------------
# 1. Align representative sequences, mask, build unrooted + rooted tree
# ------------------------------------------------------------------
echo "Running align-to-tree-mafft-fasttree..."

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${REPSEQS} \
  --o-alignment ${PHYLO_DIR}/CFM_16S_aligned_repseqs.qza \
  --o-masked-alignment ${PHYLO_DIR}/CFM_16S_masked_aligned_repseqs.qza \
  --o-tree ${PHYLO_DIR}/CFM_16S_unrooted_tree.qza \
  --o-rooted-tree ${PHYLO_DIR}/CFM_16S_rooted_tree.qza \
  --p-n-threads ${SLURM_CPUS_PER_TASK} \
  --verbose

echo "Phylogenetic tree construction complete"
echo "Date: $(date)"

# ------------------------------------------------------------------
# 2. Export rooted tree (Newick) for use in R / Python
# ------------------------------------------------------------------
echo "Exporting rooted tree to Newick..."

EXPORT_DIR="qiime2/export/CFM_16S_rooted_tree"
mkdir -p ${EXPORT_DIR}

qiime tools export \
  --input-path ${PHYLO_DIR}/CFM_16S_rooted_tree.qza \
  --output-path ${EXPORT_DIR}

echo "Exported to ${EXPORT_DIR}/tree.nwk"
echo "All done"
echo "Date: $(date)"
