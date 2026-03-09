# ---- Activate packages ----

# Importing data
library(readr) # (part of tidyverse) import of delimited text files (CSV, TSV)
library(readxl) # import of Excel files (.xls, .xlsx)
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


################################################################################
################### Import data and phyloseq object building ###################
################################################################################

# ---- Creating phyloseq object for 16S microbiome data ----
# 1) Metadata: Import metadata from .xlsx
metadata <- read_xlsx(
        here("data", "sample_metadata_with_microscopy.xlsx"),
        sheet = "metadata",
        col_types = c("text","text","text","text","text",
                      "text", "date", "date", "date", "text",
                      "text", "numeric", "date", "date", "text",
                      "text", "numeric", "date", "logical", "text",
                      "numeric", "logical", "logical", "numeric", "logical",
                      "numeric", "logical", "numeric", "logical", "logical",
                      "numeric")
        )  |>
        # Remove qiime2:types row (QIIME2 metadata requirement)
        slice(-1) |>
        # Remove hyphen in sample-id (QIIME2 metadata requirement)
        rename(sample_id = `sample-id`) |>
        # To match row names in SampleID in phyloseq object
        column_to_rownames("sample_id") |> 
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

# Build phyloseq metadata table (sample_data-class) based on dataframe tiblle 
sample_metadata_16S <- sample_data(metadata)

# 2) ASV feature table with counts (QIIME2 DADA2 output)
# Build otu_table-class phyloseq object from biom feature table converted to matrix
asv_table_16S <- read_biom(
        here("qiime2", "export", "CFM_16S_dada2_table", "feature-table.biom")
        ) |> 
        biom_data() |> 
        as.matrix() |> 
        otu_table(taxa_are_rows = TRUE)

# 3) Taxonomy table (QIIME2 SILVA Naive-Bayes classifier output)
# Build taxonomyTable-class phyloseq object from tsv taxonomy table exported from QIIME2
taxonomy_table_16S <- read_tsv(
        here("qiime2", "export", "CFM_16S_taxonomy_diverse_weighted", "taxonomy.tsv"), 
        show_col_types = FALSE
) |> 
        transmute(
                # Fix space in the first column
                FeatureID = `Feature ID`,
                # Handle trailing spaces with regex (might be redundant)
                Taxon = str_replace_all(Taxon, "\\s+", " ")
        ) |> 
        # Split Taxon string into columns with taxonomic ranks
        separate(
                Taxon,
                # Species Hypothesis from UNITE database was preserved
                into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ";\\s*",
                # Explicit handling under-annotated taxa (by default separate rises warning and fill from the right)
                fill = "right"
        ) |> 
        # One more explicit handling blank taxonomic ranks
        mutate(across(-FeatureID, ~na_if(.x, ""))) |> 
        column_to_rownames("FeatureID") |> 
        # Phyloseq function requires matrix as input 
        as.matrix() |> 
        tax_table()

# 4) Create a full experiment-level phyloseq-class object
phyloseq_16S_raw <- phyloseq(
        asv_table_16S,
        taxonomy_table_16S,
        sample_metadata_16S
        )

# 5) Representative sequences (FASTA) file path
repseqs_fp <- here("qiime2", "export", "CFM_16S_dada2_repseqs", "dna-sequences.fasta")

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
decontam_prev055_16S <- phyloseq_16S_raw |> 
        subset_samples(sample_type != "mock_community") |>
        isContaminant(
                method = "prevalence",
                neg = "is_negative_control",
                # Posterior probability of contaminant identification
                threshold = 0.55
                )

# 2) Remove features identified as contamination from the phyloseq object
phyloseq_16S_decontam <- prune_taxa(
        !decontam_prev055_16S$contaminant,
        phyloseq_16S_raw
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

phyloseq_16S_bacteria <- phyloseq_16S_decontam |> 
        subset_taxa(
                Domain == "d__Bacteria" &
                # Organelle bacterial type ribosomal RNA 
                (is.na(Order)  | Order  != "o__Chloroplast") &
                (is.na(Family) | Family != "f__Mitochondria")
        ) |> 
        # prune empty samples (pibEC & pviEC)
        (\(ps) prune_samples(sample_sums(ps) > 0, ps))()


names(which(sample_sums(phyloseq_16S_bacteria) == 0))


# ---- Additional taxonomically filtered subsets of potential exploratory interest ---- 
# The subset with Archaea sequences includes 83 ASVs assigned to d__Archaea:
phyloseq_16S_archaea <- phyloseq_16S_decontam |> 
        subset_taxa(Domain == "d__Archaea")

# ?he subset with Eukaryota sequences includes 263 ASVs assigned to d__Eukaryota:
# This subset includes 2 ASVs assigned to g__Phytophthora, that 
phyloseq_16S_eukaryota <- phyloseq_16S_decontam |> 
        subset_taxa(Domain == "d__Eukaryota")

# The subset with Organelle sequences includes:
# a) 3032 ASVs assigned to f__Mitochondria;
# b) 455 ASVs assigned to o__Chloroplast.
phyloseq_16S_organelle <- phyloseq_16S_decontam |> 
        subset_taxa(Order == "o__Chloroplast" | Family == "f__Mitochondria")

# Additionally I created the subset with Unassigned sequences.
# It includes 402 ASVs unassigned at Domain level:
phyloseq_16S_unassigned <- phyloseq_16S_decontam |> 
        subset_taxa(Domain == "Unassigned")


################################################################################
################ Biological samples and Prevalence filtration ##################
################################################################################


# ---- Filter only biological samples (remove all controls) ----
phyloseq_16S_bacteria_biosamples <- subset_samples(
        phyloseq_16S_bacteria,
        !sample_type %in% c(
                "extraction_control",
                "negative_pcr_control",
                "mock_community"
        )
)

dir.create(here("results", "rds"), showWarnings = FALSE, recursive = TRUE)

saveRDS(
        phyloseq_16S_bacteria_biosamples,
        file = here("results", "rds", "phyloseq_16S_bacteria_biosamples.rds"),
        compress = "xz"
)

system("wmic computersystem get totalphysicalmemory", intern = TRUE)

# ---- Export cleaned tables as BIOM for gamma_diversity_metrics project ----

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
                here(out_dir, paste0(label, "_feature-table.biom"))
        )
        tibble(
                `Feature ID` = rownames(as(tax_table(ps), "matrix")),
                Taxon        = .taxon_string(as(tax_table(ps), "matrix"))
        ) |>
                write_tsv(here(out_dir, paste0(label, "_taxonomy.tsv")))
        
        message(sprintf("%-45s  %5d taxa × %3d samples", label, ntaxa(ps), nsamples(ps)))
}

dir.create(here("results", "biom"), showWarnings = FALSE, recursive = TRUE)

export_to_biom(phyloseq_16S_bacteria,            "cfm_16s_bacteria_allsamples",  here("results", "biom"))
export_to_biom(phyloseq_16S_bacteria_biosamples, "cfm_16s_bacteria_biosamples",  here("results", "biom"))
