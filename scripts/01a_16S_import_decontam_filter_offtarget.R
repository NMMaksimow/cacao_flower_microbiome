# ---- Activate packages ----

# Importing data
library(readr) # (part of tidyverse) import of delimited text files (CSV, TSV)
library(biomformat) # reading and writing BIOM-format microbiome files
library(here) # constructs file paths relative to the project root

# Data handling
library(tidyverse) # ecosystem of packages for data manipulation and plotting
# library(dplyr) # (part of tidyverse) data manipulation (filtering, grouping, summarising)
# library(tidyr) # (part of tidyverse) data reshaping (pivoting, tidying)
# library(tibble) # (part of tidyverse) modern data frames with improved printing and behavior

# Core microbiome analysis
library(vegan) # ecological analyses (diversity, ordination, PERMANOVA)
library(phyloseq) # integrate analysis of microbiome data (OTU/ASV, taxonomy, metadata)
library(decontam) # identification and removal of contaminant taxa
library(ape)      # read.tree() for loading Newick phylogenetic tree

# Data version switch: CFM_VERSION, EXPORT_DIR, RESULTS_DIR
source(here("scripts", "utils_paths.R"))
# read_sample_metadata(): sample metadata TSV with column types
source(here("scripts", "utils_read_metadata.R"))
# read_qiime2_taxonomy(): taxonomy.tsv with ranks matched by prefix
source(here("scripts", "utils_read_taxonomy.R"))

# ---- Version-specific inputs ----
# Taxonomy: v1 SILVA 138 diverse-weighted; v2 SILVA 144 V4 "average"
taxonomy_dir_16S <- c(v1 = "CFM_16S_taxonomy_diverse_weighted",
                      v2 = "CFM_16S_taxonomy_average")[[CFM_VERSION]]
# Trees: v1 has one tree (all ASVs, midpoint-rooted). v2: bact (Bacteria without
# organelles, midpoint-rooted; primary) and the prok tree rooted on the archaea with
# the archaea kept (pipeline intermediate, not in the exports; copied from Hamilton)
tree_path_bacteria <- c(
        v1 = file.path(EXPORT_DIR, "CFM_16S_rooted_tree", "tree.nwk"),
        v2 = file.path(EXPORT_DIR, "CFM_16S_rooted_tree_bact", "tree.nwk")
)[[CFM_VERSION]]
tree_path_prokaryotes <- c(
        v1 = file.path(EXPORT_DIR, "CFM_16S_rooted_tree", "tree.nwk"),
        v2 = here("qiime2_v2", "phylogeny", "prok", "CFM_16S_rooted_tree_archaea.nwk")
)[[CFM_VERSION]]

# Attach a tree to a filtered object. Never at phyloseq() construction: phyloseq
# silently drops taxa missing from the tree. Here every taxon must be a tip.
attach_tree <- function(ps, tree_path) {
        tree <- ape::read.tree(tree_path)
        stopifnot(all(taxa_names(ps) %in% tree$tip.label))
        merge_phyloseq(ps, phy_tree(ape::keep.tip(tree, taxa_names(ps))))
}


################################################################################
################### Import data and phyloseq object building ###################
################################################################################
# ---- Creating phyloseq object for 16S microbiome data ----
# 1) Metadata: Import metadata from the TSV (see scripts/utils_read_metadata.R)
metadata <- read_sample_metadata() |>
        # To match row names in SampleID in phyloseq object
        column_to_rownames("sample_id") |>
        # Restore sample_id as metadata column (needed by downstream scripts)
        (\(x) { x$sample_id <- rownames(x); x })() |>
        dplyr::mutate(
                across(
                        c(
                          # Sample metadata character to factor
                          sample_type,
                          management_type,
                          farm_id,
                          tree_id,
                          # Batch technical factors
                          sublibrary_id,
                          collection_batch,
                          extraction_batch,
                          pcr_batch_its1,
                          pcr_batch_16s,
                          fwd_primer_16s,
                          rev_primer_16s,
                          fwd_primer_its1,
                          rev_primer_its1,
                          dye),
                        as.factor
                ),
                # Ordered factor (how many stigma lobes out of 5 are affected by fungal hyphae?)
                fungal_colonization_score = factor(
                        fungal_colonization_score,
                        levels = 0:5,
                        ordered = TRUE
                ),
                # Boolean variable with negative controls, is needed to run decontam() later 
                is_negative_control = sample_type %in% c(
                        "extraction_control",
                        "negative_pcr_control"
                )
        )

# Restore sample_id as a metadata column (rownames were set by column_to_rownames above)
metadata$sample_id <- rownames(metadata)

# Build phyloseq metadata table (sample_data-class) based on dataframe tibblle
sample_metadata_16S <- sample_data(metadata)

# 2) ASV feature table with counts (QIIME2 DADA2 output)
# Build otu_table-class phyloseq object from biom feature table converted to matrix
asv_table_16S <- read_biom(
        file.path(EXPORT_DIR, "CFM_16S_dada2_table", "feature-table.biom")
        ) |>
        biom_data() |> 
        as.matrix() |> 
        otu_table(taxa_are_rows = TRUE)

# 3) Taxonomy table (QIIME2 SILVA Naive-Bayes classifier output)
# Build taxonomyTable-class phyloseq object from tsv taxonomy table exported from QIIME2
taxonomy_table_16S <- read_qiime2_taxonomy(
        file.path(EXPORT_DIR, taxonomy_dir_16S, "taxonomy.tsv")
)

# 4) Create a full experiment-level phyloseq-class object
ps_16S_raw <- phyloseq(
        asv_table_16S,
        taxonomy_table_16S,
        sample_metadata_16S
        )

# 5) Representative sequences (FASTA) file path
repseqs_fp <- file.path(EXPORT_DIR, "CFM_16S_dada2_repseqs", "dna-sequences.fasta")

################################################################################
############# Decontamination and cleaning off-target sequences ################
################################################################################

# ---- Decontam to clean potential lab contamination ----
# I've run decontam() using prevalence method.
# It identifies sequences over represented in negative controls (negative PCR & extraction controls).
# I tried it on flower and mock community samples first. 
# It has no difference in terms of number of identified sequences.
# It seems methodologically dodgy to include mock communities in identification step,
# because they might be also contain lab microbiome sequences.

# I also tried two additional thresholds:
# a) 0.2 (resulted in 60 contaminant ASVs);
# b) 0.5 (resulted in 82 contaminant ASVs)
# c) 0.7 (resulted in 6147  contaminant ASVs).

# 1) Identify contamination with decontam(): 91 ASVs identified as contaminant
decontam_prev055_16S <- ps_16S_raw |> 
        subset_samples(sample_type != "mock_community") |>
        isContaminant(
                method = "prevalence",
                neg = "is_negative_control",
                # Posterior probability of contaminant identification
                threshold = 0.55
                )

# 2) Remove features identified as contamination from the phyloseq object
ps_16S_decontam <- prune_taxa(
        !decontam_prev055_16S$contaminant,
        ps_16S_raw
        )

################################################################################
################### Filter off-target sequences ################################
################################################################################

# ---- Taxonomy-based filtration: leave d__Bacteria excluding o__Chloroplast & f__Mitochondria ----
# First I removed off-target sequences before decontam().
# But then realized that it might be not the best option.
# isContamination() should perform better on dataset with all sequences including offtarget ASVs.

# I left only target domain: d__Bacteria and cleaned it from Chloroplast and Mitochondia ASVs.
# This dataset includes 25035 ASVs assigned to d__Bacteria.
# In SILVA DB taxonomy tables:
# a) Chloroplast rDNA occurs at Order rank (o__);
# b) Mitochondria rDNA occurs at Family rank (f__).

ps_16S_bacteria <- ps_16S_decontam |> 
        subset_taxa(
                Domain == "d__Bacteria" &
                # Organelle bacterial type ribosomal RNA 
                (is.na(Order)  | Order  != "o__Chloroplast") &
                (is.na(Family) | Family != "f__Mitochondria")
        ) |> 
        # prune empty samples (pibEC & pviEC)
        (\(ps) prune_samples(sample_sums(ps) > 0, ps))() |>
        # Tree after the filter: v2 bact tree has exactly these ASVs as tips
        attach_tree(tree_path_bacteria)


names(which(sample_sums(ps_16S_bacteria) == 0))


# ---- Additional taxonomically filtered subsets of potential exploratory interest ---- 
# The subset with Archaea sequences includes 83 ASVs assigned to d__Archaea:
ps_16S_archaea <- ps_16S_decontam |>
        subset_taxa(Domain == "d__Archaea")

# The prokaryote subset: Bacteria + Archaea, organelles removed as for ps_16S_bacteria.
# 515F/806R primers target both domains, so archaea are not off-target here.
# Used for functional annotation (25a), where FAPROTAX covers archaeal functions.
ps_16S_prokaryotes <- ps_16S_decontam |>
        subset_taxa(
                Domain %in% c("d__Bacteria", "d__Archaea") &
                (is.na(Order)  | Order  != "o__Chloroplast") &
                (is.na(Family) | Family != "f__Mitochondria")
        ) |>
        (\(ps) prune_samples(sample_sums(ps) > 0, ps))() |>
        # Tree after the filter, with the archaea (v2: rooted on the archaea)
        attach_tree(tree_path_prokaryotes)

# The subset with Eukaryota sequences includes 263 ASVs assigned to d__Eukaryota:
# This subset includes 2 ASVs assigned to g__Phytophthora, the major cocoa pathogen. NB!
ps_16S_eukaryota <- ps_16S_decontam |> 
        subset_taxa(Domain == "d__Eukaryota")

# The subset with Organelle sequences includes:
# a) 3032 ASVs assigned to f__Mitochondria;
# b) 455 ASVs assigned to o__Chloroplast.
ps_16S_organelle <- ps_16S_decontam |> 
        subset_taxa(Order == "o__Chloroplast" | Family == "f__Mitochondria")

# Additionally I created the subset with Unassigned sequences.
# It includes 402 ASVs unassigned at Domain level:
ps_16S_unassigned <- ps_16S_decontam |> 
        subset_taxa(Domain == "Unassigned")


################################################################################
########################### Biological samples  ################################
################################################################################

# ---- Filter only biological samples (remove all controls) ----
# Two sets per object:
# a) _with_deviations: all 294 flowers (supplementary ordination, sensitivity tests);
# b) main name: 288 flowers sampled by the standard protocol (EXCLUDE_DEVIATIONS,
#    defined in utils_read_metadata.R). Exclusion after decontam, so decontam
#    used every sample.
# Taxa with zero counts in the set (found only in controls or mocks) are dropped.
deviation_ids <- metadata$sample_id[metadata$sampling_deviation %in% EXCLUDE_DEVIATIONS]

ps_16S_bacteria_biosamples_with_deviations <- subset_samples(
        ps_16S_bacteria,
        !sample_type %in% c(
                "extraction_control",
                "negative_pcr_control",
                "mock_community"
        )
) |>
        (\(ps) prune_taxa(taxa_sums(ps) > 0, ps))()

ps_16S_bacteria_biosamples <- drop_deviations(
        ps_16S_bacteria_biosamples_with_deviations, deviation_ids
)

ps_16S_prokaryotes_biosamples_with_deviations <- subset_samples(
        ps_16S_prokaryotes,
        !sample_type %in% c(
                "extraction_control",
                "negative_pcr_control",
                "mock_community"
        )
) |>
        (\(ps) prune_taxa(taxa_sums(ps) > 0, ps))()

ps_16S_prokaryotes_biosamples <- drop_deviations(
        ps_16S_prokaryotes_biosamples_with_deviations, deviation_ids
)

################################################################################
############### Saving resulting phyloseq objects (RDS & BIOM) #################
################################################################################

dir.create(file.path(RESULTS_DIR, "rds"), showWarnings = FALSE, recursive = TRUE)

# Intermediate objects needed by 02b
saveRDS(ps_16S_raw,
        file = file.path(RESULTS_DIR, "rds", "ps_16S_raw.rds"),
        compress = "xz")

saveRDS(decontam_prev055_16S,
        file = file.path(RESULTS_DIR, "rds", "decontam_prev055_16S.rds"),
        compress = "xz")

saveRDS(ps_16S_decontam,
        file = file.path(RESULTS_DIR, "rds", "ps_16S_decontam.rds"),
        compress = "xz")

# Bacteria all samples (controls + biosamples), used for mock community QC
saveRDS(ps_16S_bacteria,
        file = file.path(RESULTS_DIR, "rds", "ps_16S_bacteria.rds"),
        compress = "xz")

# Bacteria biosamples: 288 (main analyses) and 294 (with deviations)
saveRDS(
        ps_16S_bacteria_biosamples,
        file = file.path(RESULTS_DIR, "rds", "ps_16S_bacteria_biosamples.rds"),
        compress = "xz"
)

saveRDS(
        ps_16S_bacteria_biosamples_with_deviations,
        file = file.path(RESULTS_DIR, "rds", "ps_16S_bacteria_biosamples_with_deviations.rds"),
        compress = "xz"
)

# Bacteria + Archaea biosamples, used for functional annotation (25a)
saveRDS(
        ps_16S_prokaryotes_biosamples,
        file = file.path(RESULTS_DIR, "rds", "ps_16S_prokaryotes_biosamples.rds"),
        compress = "xz"
)

saveRDS(
        ps_16S_prokaryotes_biosamples_with_deviations,
        file = file.path(RESULTS_DIR, "rds", "ps_16S_prokaryotes_biosamples_with_deviations.rds"),
        compress = "xz"
)

# ---- Export cleaned tables as BIOM ----

.taxon_string <- function(tax_matrix) {
        apply(tax_matrix, 1, function(row) {
                parts <- row[!is.na(row) & nchar(row) > 0]
                paste(parts, collapse = ";")
        })
}

export_to_biom <- function(ps, label, out_dir) {
        otu_mat <- as(otu_table(ps), "matrix")
        if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
        
        write_biom(
                make_biom(data = otu_mat),
                file.path(out_dir, paste0(label, "_feature-table.biom"))
        )
        tibble(
                `Feature ID` = rownames(as(tax_table(ps), "matrix")),
                Taxon        = .taxon_string(as(tax_table(ps), "matrix"))
        ) |>
                write_tsv(file.path(out_dir, paste0(label, "_taxonomy.tsv")))
        
        message(sprintf("%-45s  %5d taxa × %3d samples", label, ntaxa(ps), nsamples(ps)))
}

dir.create(file.path(RESULTS_DIR, "biom"), showWarnings = FALSE, recursive = TRUE)

export_to_biom(ps_16S_bacteria,            "cfm_16s_bacteria_allsamples",  file.path(RESULTS_DIR, "biom"))
export_to_biom(ps_16S_bacteria_biosamples, "cfm_16s_bacteria_biosamples",  file.path(RESULTS_DIR, "biom"))
export_to_biom(ps_16S_bacteria_biosamples_with_deviations,
               "cfm_16s_bacteria_biosamples_with_deviations",  file.path(RESULTS_DIR, "biom"))
