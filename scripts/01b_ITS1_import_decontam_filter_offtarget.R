# ---- Activate packages ----
# Importing data
library(biomformat) # reading and writing BIOM-format microbiome files
library(here)       # constructs file paths relative to the project root

# Data handling
library(tidyverse)  

# Core microbiome analysis
library(phyloseq)   # integrate analysis of microbiome data (OTU/ASV, taxonomy, metadata)
library(decontam)   # identification and removal of contaminant taxa

# Data version switch: CFM_VERSION, EXPORT_DIR, RESULTS_DIR
source(here("scripts", "utils_paths.R"))
# read_sample_metadata(): sample metadata TSV with column types
source(here("scripts", "utils_read_metadata.R"))
# read_qiime2_taxonomy(): taxonomy.tsv with ranks matched by prefix
source(here("scripts", "utils_read_taxonomy.R"))

################################################################################
################### Import data and phyloseq object building ###################
################################################################################

# ---- Creating phyloseq object for ITS1 microbiome data ----
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

# Build phyloseq metadata table (sample_data-class) based on imported metadata dataframe
sample_metadata_ITS1 <- phyloseq::sample_data(metadata)

# 2) ASV feature table with counts (QIIME2 DADA2 output)
# Build otu_table-class phyloseq object from biom feature table converted to matrix
# suppressWarnings: biomformat probes for JSON (BIOM v1) before falling back to
# HDF5 (BIOM v2); the JSON probe fails on binary HDF5 bytes → harmless warnings
asv_table_ITS1 <- suppressWarnings(
        read_biom(
                file.path(EXPORT_DIR, "CFM_ITS1_dada2_table", "feature-table.biom")
        )
) |>
        biom_data() |>
        as.matrix() |>
        otu_table(taxa_are_rows = TRUE)

# 3) Taxonomy table (QIIME2 UNITE Naive-Bayes classifier output)
# Ranks matched by prefix; v1 also has the UNITE species hypothesis (sh__ -> UNITE_SH)
taxonomy_table_ITS1 <- read_qiime2_taxonomy(
        file.path(EXPORT_DIR, "CFM_ITS1_taxonomy", "taxonomy.tsv")
)

# 4) Create a full experiment-level phyloseq-class object
ps_ITS1_raw <- phyloseq(
        asv_table_ITS1,
        taxonomy_table_ITS1,
        sample_metadata_ITS1
        )


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
# a) more strict - 0.2 (resulted in 11 contaminant features);
# b) and more relaxed - 0.7 (resulted in 2568 contaminant features).

# 1) Identify contamination with decontam(): 26 ASVs identified as contaminant (v1)
decontam_prev050_ITS1 <- ps_ITS1_raw |>
        subset_samples(sample_type != "mock_community") |>
        isContaminant(
                method = "prevalence",
                neg = "is_negative_control",
                # Posterior probability of contaminant identification
                threshold = 0.5
                )

# 2) Remove features identified as contamination from the phyloseq object
ps_ITS1_decontam <- prune_taxa(
        !decontam_prev050_ITS1$contaminant,
        ps_ITS1_raw
        )

# ---- Taxonomy-based filtration: leave k__Fungi ----
# I left only target kingdom: k__Fungi:
ps_ITS1_fungi <- ps_ITS1_decontam |> 
        subset_taxa(Kingdom == "k__Fungi") |>
        (\(ps) prune_samples(sample_sums(ps) > 0, ps))()

# ---- Additional taxonomically filtered subsets of potential exploratory interest ----
# Additionally I created the subset with Unassigned sequences:
ps_ITS1_unassigned <- ps_ITS1_decontam |> 
        subset_taxa(Kingdom == "Unassigned")


################################################################################
########################### Biological samples  ################################
################################################################################

# ---- Filter only biological samples (remove all controls) ----
ps_ITS1_fungi_biosamples <- subset_samples(
        ps_ITS1_fungi,
        !sample_type %in% c(
                "extraction_control",
                "negative_pcr_control",
                "mock_community"
        )
)

################################################################################
############### Saving resulting phyloseq objects (RDS & BIOM) #################
################################################################################
dir.create(file.path(RESULTS_DIR, "rds"), showWarnings = FALSE, recursive = TRUE)

# Intermediate objects needed by 02b
saveRDS(ps_ITS1_raw,
        file = file.path(RESULTS_DIR, "rds", "ps_ITS1_raw.rds"),
        compress = "xz")

saveRDS(decontam_prev050_ITS1,
        file = file.path(RESULTS_DIR, "rds", "decontam_prev050_ITS1.rds"),
        compress = "xz")

saveRDS(ps_ITS1_decontam,
        file = file.path(RESULTS_DIR, "rds", "ps_ITS1_decontam.rds"),
        compress = "xz")

# Fungi all samples (controls + biosamples), used for mock community QC
saveRDS(ps_ITS1_fungi,
        file = file.path(RESULTS_DIR, "rds", "ps_ITS1_fungi.rds"),
        compress = "xz")

saveRDS(
        ps_ITS1_fungi_biosamples,
        file = file.path(RESULTS_DIR, "rds", "ps_ITS1_fungi_biosamples.rds"),
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

export_to_biom(ps_ITS1_fungi,            "cfm_its1_fungi_allsamples",  file.path(RESULTS_DIR, "biom"))
export_to_biom(ps_ITS1_fungi_biosamples, "cfm_its1_fungi_biosamples",  file.path(RESULTS_DIR, "biom"))


################################################################################
################### OTU97 phyloseq object (parallel workflow) #################
################################################################################
# ITS1 reads clustered at 97% similarity with VSEARCH (open-reference).
# Same metadata, same decontam parameters, same off-target filter as the ASV
# workflow above — only the feature table and taxonomy paths differ.

# ---- Import OTU97 feature table ----
# OTU count matrix from QIIME2 vsearch cluster-features output
otu_table_ITS1_otu97 <- suppressWarnings(
        read_biom(
                file.path(EXPORT_DIR, "CFM_ITS1_otu97_table", "feature-table.biom")
        )
) |>
        biom_data() |>
        as.matrix() |>
        otu_table(taxa_are_rows = TRUE)

# ---- Import OTU97 taxonomy table ----
# Same UNITE NB-classifier output format as the ASV taxonomy
taxonomy_table_ITS1_otu97 <- read_qiime2_taxonomy(
        file.path(EXPORT_DIR, "CFM_ITS1_otu97_taxonomy", "taxonomy.tsv")
)

# ---- Build raw phyloseq object ----
# Reuse sample_metadata_ITS1 created above — same samples, same metadata
ps_ITS1_otu97_raw <- phyloseq(
        otu_table_ITS1_otu97,
        taxonomy_table_ITS1_otu97,
        sample_metadata_ITS1
        )

cat(sprintf("OTU97 raw: %d OTUs × %d samples\n",
            ntaxa(ps_ITS1_otu97_raw),
            nsamples(ps_ITS1_otu97_raw)))


################################################################################
############### Decontamination (OTU97) ########################################
################################################################################

# Same prevalence-method decontam as the ASV workflow (threshold = 0.5).
# Run on all sample types except mock communities.
decontam_prev050_ITS1_otu97 <- ps_ITS1_otu97_raw |>
        subset_samples(sample_type != "mock_community") |>
        isContaminant(
                method    = "prevalence",
                neg       = "is_negative_control",
                threshold = 0.5
                )

cat(sprintf("OTU97 decontam: %d OTUs flagged as contaminant\n",
            sum(decontam_prev050_ITS1_otu97$contaminant)))

# Remove contaminant OTUs from the full raw object
ps_ITS1_otu97_decontam <- prune_taxa(
        !decontam_prev050_ITS1_otu97$contaminant,
        ps_ITS1_otu97_raw
        )


################################################################################
############### Off-target filter (OTU97) ######################################
################################################################################

# Keep only k__Fungi (same as ASV workflow)
ps_ITS1_otu97_fungi <- ps_ITS1_otu97_decontam |>
        subset_taxa(Kingdom == "k__Fungi") |>
        (\(ps) prune_samples(sample_sums(ps) > 0, ps))()

# Unassigned subset (mirrors ASV workflow, for exploratory use)
ps_ITS1_otu97_unassigned <- ps_ITS1_otu97_decontam |>
        subset_taxa(Kingdom == "Unassigned")


################################################################################
############### Biological samples (OTU97) #####################################
################################################################################

# Remove extraction controls, negative PCR controls, and mock communities
ps_ITS1_otu97_fungi_biosamples <- subset_samples(
        ps_ITS1_otu97_fungi,
        !sample_type %in% c(
                "extraction_control",
                "negative_pcr_control",
                "mock_community"
        )
)

cat(sprintf("OTU97 biosamples (fungi only): %d OTUs × %d samples\n",
            ntaxa(ps_ITS1_otu97_fungi_biosamples),
            nsamples(ps_ITS1_otu97_fungi_biosamples)))


################################################################################
############### Save OTU97 RDS objects #########################################
################################################################################

saveRDS(ps_ITS1_otu97_raw,
        file     = file.path(RESULTS_DIR, "rds", "ps_ITS1_otu97_raw.rds"),
        compress = "xz")

saveRDS(decontam_prev050_ITS1_otu97,
        file     = file.path(RESULTS_DIR, "rds", "decontam_prev050_ITS1_otu97.rds"),
        compress = "xz")

saveRDS(ps_ITS1_otu97_decontam,
        file     = file.path(RESULTS_DIR, "rds", "ps_ITS1_otu97_decontam.rds"),
        compress = "xz")

saveRDS(ps_ITS1_otu97_fungi,
        file     = file.path(RESULTS_DIR, "rds", "ps_ITS1_otu97_fungi.rds"),
        compress = "xz")

saveRDS(ps_ITS1_otu97_fungi_biosamples,
        file     = file.path(RESULTS_DIR, "rds", "ps_ITS1_otu97_fungi_biosamples.rds"),
        compress = "xz")

# Export OTU97 cleaned tables as BIOM (mirrors ASV export above)
export_to_biom(ps_ITS1_otu97_fungi,
               "cfm_its1_otu97_fungi_allsamples",  file.path(RESULTS_DIR, "biom"))
export_to_biom(ps_ITS1_otu97_fungi_biosamples,
               "cfm_its1_otu97_fungi_biosamples",  file.path(RESULTS_DIR, "biom"))


################################################################################
############### R1-only OTU97 (v2 sensitivity branch) ##########################
################################################################################
# v2 only: forward reads denoised alone (cut to 220 nt), runs merged, chimeras
# removed, 97% OTUs, same UNITE classifier. Merging R1 with R2 loses the long
# Saccharomyces amplicon in the mock; this branch tests whether the flower data
# change (QC 09: ~1% of reads, no genus lost). Same metadata, decontam and
# off-target filter as the OTU97 workflow above.

if (CFM_VERSION == "v2") {

        # ---- Import R1-only OTU97 table and taxonomy ----
        otu_table_ITS1_r1_otu97 <- suppressWarnings(
                read_biom(
                        file.path(EXPORT_DIR, "CFM_ITS1_r1_otu97_table", "feature-table.biom")
                )
        ) |>
                biom_data() |>
                as.matrix() |>
                otu_table(taxa_are_rows = TRUE)

        taxonomy_table_ITS1_r1_otu97 <- read_qiime2_taxonomy(
                file.path(EXPORT_DIR, "CFM_ITS1_r1_otu97_taxonomy", "taxonomy.tsv")
        )

        ps_ITS1_r1_otu97_raw <- phyloseq(
                otu_table_ITS1_r1_otu97,
                taxonomy_table_ITS1_r1_otu97,
                sample_metadata_ITS1
                )

        cat(sprintf("R1 OTU97 raw: %d OTUs × %d samples\n",
                    ntaxa(ps_ITS1_r1_otu97_raw),
                    nsamples(ps_ITS1_r1_otu97_raw)))

        # ---- Decontam (threshold 0.5, mock communities excluded) ----
        decontam_prev050_ITS1_r1_otu97 <- ps_ITS1_r1_otu97_raw |>
                subset_samples(sample_type != "mock_community") |>
                isContaminant(
                        method    = "prevalence",
                        neg       = "is_negative_control",
                        threshold = 0.5
                        )

        cat(sprintf("R1 OTU97 decontam: %d OTUs flagged as contaminant\n",
                    sum(decontam_prev050_ITS1_r1_otu97$contaminant)))

        ps_ITS1_r1_otu97_decontam <- prune_taxa(
                !decontam_prev050_ITS1_r1_otu97$contaminant,
                ps_ITS1_r1_otu97_raw
                )

        # ---- Off-target filter: k__Fungi only ----
        ps_ITS1_r1_otu97_fungi <- ps_ITS1_r1_otu97_decontam |>
                subset_taxa(Kingdom == "k__Fungi") |>
                (\(ps) prune_samples(sample_sums(ps) > 0, ps))()

        # ---- Biological samples ----
        ps_ITS1_r1_otu97_fungi_biosamples <- subset_samples(
                ps_ITS1_r1_otu97_fungi,
                !sample_type %in% c(
                        "extraction_control",
                        "negative_pcr_control",
                        "mock_community"
                )
        )

        cat(sprintf("R1 OTU97 biosamples (fungi only): %d OTUs × %d samples\n",
                    ntaxa(ps_ITS1_r1_otu97_fungi_biosamples),
                    nsamples(ps_ITS1_r1_otu97_fungi_biosamples)))

        # ---- Save RDS and BIOM ----
        saveRDS(ps_ITS1_r1_otu97_raw,
                file     = file.path(RESULTS_DIR, "rds", "ps_ITS1_r1_otu97_raw.rds"),
                compress = "xz")

        saveRDS(decontam_prev050_ITS1_r1_otu97,
                file     = file.path(RESULTS_DIR, "rds", "decontam_prev050_ITS1_r1_otu97.rds"),
                compress = "xz")

        saveRDS(ps_ITS1_r1_otu97_decontam,
                file     = file.path(RESULTS_DIR, "rds", "ps_ITS1_r1_otu97_decontam.rds"),
                compress = "xz")

        saveRDS(ps_ITS1_r1_otu97_fungi,
                file     = file.path(RESULTS_DIR, "rds", "ps_ITS1_r1_otu97_fungi.rds"),
                compress = "xz")

        saveRDS(ps_ITS1_r1_otu97_fungi_biosamples,
                file     = file.path(RESULTS_DIR, "rds", "ps_ITS1_r1_otu97_fungi_biosamples.rds"),
                compress = "xz")

        export_to_biom(ps_ITS1_r1_otu97_fungi,
                       "cfm_its1_r1_otu97_fungi_allsamples",  file.path(RESULTS_DIR, "biom"))
        export_to_biom(ps_ITS1_r1_otu97_fungi_biosamples,
                       "cfm_its1_r1_otu97_fungi_biosamples",  file.path(RESULTS_DIR, "biom"))
}
