# ---- AI agent for R ----
# 1. Install devtools (if not already installed)
install.packages("devtools")
library(devtools)

# Command to edit .Renviron file to add GitHub Classic token
usethis::edit_r_environ()

# 2. Install RgentAI package
devtools::install_github("NathanBresette/Rgent-AI", force = TRUE, upgrade = "never")

# 3. Load the library
library(rstudioai)

# 4. Launch the AI assistant
run_rgent()

# ---- Activate packages ----
# ANCOMBC2 requires R > 4.4.1 so I had to update R version and reinstall packages:
# install.packages("circlize")
# install.packages("indicspecies")
# BiocManager::install("ANCOMBC", force = TRUE)
# BiocManager::install("remotes")
# install.packages("zetadiv", repos = c("https://cran.r-project.org", "https://cloud.r-project.org"))

library(DiagrammeR) # an option to build DAG with analysis

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
library(zetadiv) # zeta-diversity analysis and zeta-decline

# differential abundance & indicator species analyses
library(ANCOMBC)  # differential abundance testing accounting for compositionality
library(indicspecies) # indicator species analysis (specificity and fidelity of association of taxa with groups)

# Linear mixed-effects models
library(lme4) # linear and generalized linear mixed-effects models
library(lmerTest) # p-values and tests for linear mixed-effects models
library(ggeffects) # marginal effects and predictions from regression models
library(glmmTMB) # flexible generalised linear mixed-effect models
library(DHARMa) # residual diagnostics for mixed and generalized models

# Data visualization
library(ggplot2) # (part of tidyverse) grammar of graphics for data visualization
library(ggh4x) # ggplot2 extension for advanced facets, scales
library(patchwork) # combining multiple ggplot objects into one figure
library(cowplot) # plot alignment and figure composition
library(ggrepel) # non-overlapping text labels in ggplot2
library(scales) # (part of tidyverse) scales transformation and formatting for plots
library(grid) # low-level graphical system used for complex layouts
library(viridis) # perceptually uniform colour paletts

# Heatmaps and complex visualization
library(pheatmap) # simple heatmap generation
library(ComplexHeatmap) # advanced and highly customizable heatmaps
library(circlize) # circular visualization utilities (used by ComplexHeatmap)
# library(metacoder) # visualization and manipulation of hierarchical taxonomic data

# renv::snreadxlrenv::snapshot()

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


# ---- Some diagnostic commands for 16S phyloseq object ----
sample_variables(phyloseq_16S_raw)
nsamples(phyloseq_16S_raw)
ntaxa(phyloseq_16S_raw)

table(tax_table(phyloseq_16S_raw)[, "Domain"],
      useNA = "ifany")

read_tsv(
        file.path("qiime2", "export", "CFM_16S_taxonomy_diverse_weighted", "taxonomy.tsv"),
        show_col_types = FALSE
) |>
        dplyr::filter(
                stringr::str_detect(Taxon, "\\s{2,}|\\s$|^\\s")
        )

otu_table(phyloseq_16S_raw) |> View()
tax_table(phyloseq_16S_raw) |> View()
ntaxa(phyloseq_16S_raw)
# I didn't save RDS, but these commands are could be used for this:
# saveRDS(phyloseq_16S_raw, "results/phyloseq_16S_raw.rds")
# phyloseq_16S_raw <- readRDS("results/phyloseq_16S_raw.rds")



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

# ---- Diagnostic commands for decontaminated phyloseq object ----
# isContaminant() produces data.frame with ASV IDs as row names.
# freq column contains ASV's frequency across whole ASV table.
# prev column contains the number of samples where ASV was identified.
# contaminant column is a Boolean vector.
class(decontam_prev_16S)
View(decontam_prev_16S)
# 82 ASVs were idenfied as contaminant
table(decontam_prev_16S$contaminant)

# Generate taxonomy table for contaminant ASVs identified by decontam()
contaminant_16S_taxa_tbl <- tax_table(phyloseq_16S_raw)[decontam_prev055_16S$contaminant, ] |>
        as.data.frame() |>
        rownames_to_column("FeatureID") |>
        # To add total count for each taxa
        left_join(
                taxa_sums(phyloseq_16S_raw) |> 
                        enframe(name = "FeatureID", value = "Total_abundance"),
                by = "FeatureID"
        ) |>
        arrange(desc(Total_abundance))


View(contaminant_16S_taxa_tbl)
# Extract Species identified as contaminants
contaminant_16S_taxa_tbl |>  pull(Species)

# Overview of taxa identified as contamination
contaminant_16S_taxa_tbl |>
        count(Domain, Phylum, Genus, Species, sort = TRUE)

# How many taxa were before and after decontamination?
ntaxa(phyloseq_16S_raw) # 29361
ntaxa(phyloseq_16S_decontam) # 29279

# Summary statistics on total read count by sample before and after decontam()?
summary(sample_sums(phyloseq_16S_raw))
summary(sample_sums(phyloseq_16S_decontam))

# How many samples have resulted in read count below 1000 and 5000?
sum(sample_sums(phyloseq_16S_raw) < 1000)
sum(sample_sums(phyloseq_16S_raw) < 5000)

sum(sample_sums(phyloseq_16S_decontam) < 1000)
sum(sample_sums(phyloseq_16S_decontam) < 5000)

# What samples have low read count?
sample_sums(phyloseq_16S_raw) |>
        (\(x) names(x[x < 1000]))()

# All 28 negative controls (negative_pcr_control, extraction_control have less than 1000 read count)
sample_sums(phyloseq_16S_decontam) |>
        (\(x) names(x[x < 1000]))()

# In addition 2 flower samples (cmt72, cmt73 have less than 5000 read count)
sample_sums(phyloseq_16S_decontam) |>
        (\(x) names(x[x < 5000]))()

# presence/absence
pa_mat <- as(otu_table(phyloseq_16S_raw), "matrix")
if (!taxa_are_rows(phyloseq_16S_raw)) pa_mat <- t(pa_mat)
pa_mat <- pa_mat > 0

neg_flag <- sample_data(phyloseq_16S_raw)$is_negative_control

prev_tbl <- tibble(
        ASV = rownames(pa_mat),
        prev_neg = rowMeans(pa_mat[, neg_flag == TRUE, drop = FALSE]),
        prev_bio = rowMeans(pa_mat[, neg_flag == FALSE, drop = FALSE])
)

ggplot(prev_tbl, aes(x = prev_neg, y = prev_bio)) +
        geom_point(alpha = 0.4, size = 1) +
        geom_abline(linetype = "dashed", color = "red") +
        scale_x_continuous("Prevalence in negative controls") +
        scale_y_continuous("Prevalence in biological samples") +
        theme_minimal()

# Less strict contaminant identification
decontam_prev050_16S <- phyloseq_16S_raw |> 
        subset_samples(sample_type != "mock_community") |>
        isContaminant(
                method = "prevalence",
                neg = "is_negative_control",
                # Posterior probability of contaminant identification
                threshold = 0.50
        )

contaminants_050 <- rownames(decontam_prev050_16S)[decontam_prev050_16S$contaminant]
contaminants_055 <- rownames(decontam_prev055_16S)[decontam_prev055_16S$contaminant]

# ???????????
contaminants_shared <- intersect(contaminants_050, contaminants_055)

# ?????? ??? 0.50
contaminants_only_050 <- setdiff(contaminants_050, contaminants_055)

# ?????? ??? 0.55
contaminants_only_055 <- setdiff(contaminants_055, contaminants_050)

# ???????
c(
        n_050 = length(contaminants_050),
        n_055 = length(contaminants_055),
        n_shared = length(contaminants_shared),
        only_050 = length(contaminants_only_050),
        only_055 = length(contaminants_only_055)
)

comparison_contaminants_tbl <-
        tibble::tibble(
                ASV_ID = union(contaminants_050, contaminants_055),
                contaminant_050 = ASV_ID %in% contaminants_050,
                contaminant_055 = ASV_ID %in% contaminants_055
        ) |>
        dplyr::left_join(
                as.data.frame(tax_table(phyloseq_16S_raw)) |>
                        tibble::rownames_to_column("ASV_ID"),
                by = "ASV_ID"
        ) |>
        dplyr::mutate(
                status = dplyr::case_when(
                        contaminant_050 & contaminant_055 ~ "both",
                        contaminant_050 & !contaminant_055 ~ "only_050",
                        !contaminant_050 & contaminant_055 ~ "only_055"
                )
        )

comparison_contaminants_tbl |> count(status)

comparison_contaminants_tbl |>
        filter(status == "only_055") |>
        dplyr::select(ASV_ID, Domain:Species) |> 
        View()

comparison_contaminants_tbl |>
        filter(status == "only_055") |>
        count(Genus, sort = TRUE)

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
        )



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

# ---- Taxonomy diagnostic commands ----
# Pipe for quick check of unique entries in taxonomy tables at any level:
tax_table(phyloseq_16S_bacteria) |> as.data.frame() |> pull(Family) |> unique()

# Check the ranks of taxonomy table where Chloroplast entry occurs
as.data.frame(tax_table(phyloseq_16S_organelle)) |>
        dplyr::summarise(
                dplyr::across(
                        everything(),
                        ~ sum(grepl("Chloroplast", ., ignore.case = TRUE), na.rm = TRUE),
                        .names = "Chloroplast_in_{.col}"
                )
        )

# Check the ranks of taxonomy table where Mitochondria entry occurs
as.data.frame(tax_table(phyloseq_16S_organelle)) |>
        dplyr::summarise(
                dplyr::across(
                        everything(),
                        ~ sum(grepl("Mitochondria", ., ignore.case = TRUE), na.rm = TRUE),
                        .names = "Mitochondria_in_{.col}"
                )
        )







################################################################################
### Quality Control: decontamination, off-target sequences & mock communities ##
################################################################################

# ---- Effect of decontamination on control samples ----
# Subset of control samples before decontamination transformed to RA.
# There're 42 samples and 265 taxa in this data subset.
phyloseq_16S_raw |> 
        subset_samples(
                sample_type %in% c("mock_community",
                                   "negative_pcr_control",
                                   "extraction_control")
        ) |> 
        # Prune taxa with zero counts
        (\(x) prune_taxa(taxa_sums(x) > 0, x))() |> 
        #transform_sample_counts(function(x) x / sum(x)) |>  
        plot_bar(fill = "Genus") +
        facet_grid(~ sample_type, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
                legend.position = "none"
        )

# Subset of control samples after decontamination transofrmed to RA.
# There're 42 samples and 174 taxa in this data subset.
phyloseq_16S_decontam |> 
        subset_samples(
                sample_type %in% c("mock_community",
                           "negative_pcr_control",
                           "extraction_control")
        ) |> 
        # Prune taxa with zero counts
        (\(x) prune_taxa(taxa_sums(x) > 0, x))() |> 
        #transform_sample_counts(function(x) x / sum(x)) |> 
        plot_bar(fill = "Genus") +
        facet_grid(~ sample_type, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
                legend.position = "none"
        )

# ---- Visualize decontamination effect on control samples ---- 
make_control_long_table <- function(physeq, label) {
        physeq |>
                subset_samples(sample_type %in% c("mock_community",
                                                  "negative_pcr_control",
                                                  "extraction_control")) |>
                (\(x) prune_taxa(taxa_sums(x) > 0, x))() |>
                psmelt() |>
                dplyr::mutate(
                        stage = label,
                        Family = ifelse(is.na(Family), "f__Unassigned", Family),
                        Genus  = ifelse(is.na(Genus),  "g__Unassigned", Genus),
                        Species = ifelse(is.na(Species), "s__Unassigned", Species),
                        Family_Genus_Species = paste(Family, Genus, Species, sep = " | ")
                )
}

control_long_tbl <-
        dplyr::bind_rows(
                make_control_long_table(phyloseq_16S_raw, "Before decontam"),
                make_control_long_table(phyloseq_16S_decontam, "After decontam")
        ) |>
        dplyr::mutate(stage = factor(stage, levels = c("Before decontam", "After decontam")))

top10_taxa <-
        control_long_tbl |>
        dplyr::group_by(Family_Genus_Species) |>
        dplyr::summarise(total_reads = sum(Abundance), .groups = "drop") |>
        dplyr::arrange(desc(total_reads)) |>
        dplyr::slice_head(n = 10) |>
        dplyr::pull(Family_Genus_Species)

control_long_tbl2 <-
        control_long_tbl |>
        dplyr::mutate(
                Family_Genus_Species = ifelse(Family_Genus_Species %in% top10_taxa,
                                              Family_Genus_Species,
                                              "Other")
        )

ggplot(control_long_tbl2,
       aes(x = Sample, y = Abundance, fill = Family_Genus_Species)) +
        geom_bar(stat = "identity") +
        facet_grid(stage ~ sample_type, scales = "free_x", space = "free_x") +
        labs(y = "Read count", x = NULL, fill = "Top taxa") +
        theme_bw() +
        theme(
                axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                legend.position = "right"
        )

# ---- Visualize expected Vs observed mock communities ----
# Expected mock community genera for D6300 ZymoBiomics DNA mixture
# NB: in mock community documentation Escherichia coli specified, but I added it as genus Escherichia-Shigella (SILVA complex)
# Salmonella haven't been properly annotated to the genus level
# create a new taxonomic rank Family | Genus to handle undeannotated Salmonella species
expected_mock <- data.frame(
        Family_Genus = c(
                "f__Listeriaceae | g__Listeria",
                "f__Bacillaceae | g__Bacillus",
                "f__Staphylococcaceae | g__Staphylococcus",
                "f__Lactobacillaceae | g__Lactobacillus",
                "f__Enterococcaceae | g__Enterococcus",
                "f__Pseudomonadaceae | g__Pseudomonas",
                "f__Enterobacteriaceae | g__Escherichia-Shigella",
                "f__Enterobacteriaceae | g__NA"  # Salmonella
        ),
        expected_RA = c(rep(12.5, 8))  # percents
)






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

# ---- Secondary filtering of biological samples ----
# Check the samples with read count per sample lower than 1000:
sample_sums(phyloseq_16S_bacteria_decontam)[
        sample_sums(phyloseq_16S_bacteria_decontam) < 1000
]

# Remove samples with read count per sample lower than 1000:
phyloseq_16S_bacteria_filt_samples <- prune_samples(
        sample_sums(phyloseq_16S_bacteria_decontam) >= 1000,
        phyloseq_16S_bacteria_decontam
)

nsamples(phyloseq_16S_bacteria_filt_samples)

hist(sample_sums(phyloseq_16S_bacteria_filt_samples))


# Prevalence filtering with 3 different parameters
filter_by_prevalence <- function(physeq, min_prev) {
        prune_taxa(
                apply(otu_table(physeq) > 0, 1, sum) > min_prev,
                physeq
        )
}

phyloseq_prev1  <- filter_by_prevalence(phyloseq_16S_bacteria_bio, 1)
phyloseq_prev2  <- filter_by_prevalence(phyloseq_16S_bacteria_bio, 2)
phyloseq_prev3  <- filter_by_prevalence(phyloseq_16S_bacteria_bio, 3)
phyloseq_prev9  <- filter_by_prevalence(phyloseq_16S_bacteria_bio, 9)


# compare total taxa present after four different prevalence filtering strategies
data.frame(
        prevalence_threshold = c("no prev filter", ">1", ">2", ">3", ">9"),
        ntaxa = c(
                ntaxa(phyloseq_16S_bacteria_bio),
                ntaxa(phyloseq_prev1),
                ntaxa(phyloseq_prev2),
                ntaxa(phyloseq_prev3),
                ntaxa(phyloseq_prev9)
        )
)

# Taxa per sample with non-zero counts after different prevalence filtering approaches
count_taxa_per_sample <- function(physeq) {
        apply(
                otu_table(physeq) > 0,
                2,
                sum
        )
}

taxa_per_sample <- data.frame(
        sample_id = sample_names(phyloseq_16S_bacteria_bio),
        no_prev = count_taxa_per_sample(phyloseq_16S_bacteria_bio),
        prev1   = count_taxa_per_sample(phyloseq_prev1),
        prev2   = count_taxa_per_sample(phyloseq_prev2),
        prev3   = count_taxa_per_sample(phyloseq_prev3),
        prev9   = count_taxa_per_sample(phyloseq_prev9)
)


summary(taxa_per_sample[, -1])

taxa_per_sample_grouped <- taxa_per_sample |>
        dplyr::left_join(
                data.frame(
                        sample_id = sample_names(phyloseq_16S_bacteria_bio),
                        sample_type = sample_data(phyloseq_16S_bacteria_bio)$sample_type
                ),
                by = "sample_id"
        )

class(taxa_per_sample_grouped)
View(taxa_per_sample_grouped)










################################################################################
##### Uneven sampling effort/sequencing depth: rarefaction Vs normalization ####
################################################################################

# There are two principle approaches to tackle uneven sampling effort/sequencing depth:
# a) Rarefaction
?vegan::rarefy()
# vegan::rrarefy()
# vegan::drarefy() calculated a
# vegan::rarecurve()
?vegan::rareslope()

# I don't understand about this function whether it performs iterative subsampling, or do it once.
?phyloseq::rarefy_even_depth()
# b) Normalization
?vegan::decostand()

# ---- Sequencing retention after decontamination and removal off-targe sequences ----
# Sequences retention during decontamination and taxonomic off-target sequence cleaning
read_retention_loss_table_16S <-
        dplyr::bind_rows(
                phyloseq_16S_raw |>
                        sample_sums() |>
                        tibble::enframe(name = "sample_id",
                                        value = "n_reads") |>
                        dplyr::mutate(stage = "Raw"),
                
                phyloseq_16S_decontam |>
                        sample_sums() |>
                        tibble::enframe(name = "sample_id",
                                        value = "n_reads") |>
                        dplyr::mutate(stage = "After decontam"),
                
                phyloseq_16S_bacteria |>
                        sample_sums() |>
                        tibble::enframe(name = "sample_id",
                                        value = "n_reads") |>
                        dplyr::mutate(stage = "Bacteria only")
        ) |>
        dplyr::mutate(sample_id = stringr::str_trim(as.character(sample_id))) |>
        tidyr::complete(
                sample_id,
                stage,
                fill = list(n_reads = 0)
        ) |> 
        dplyr::group_by(sample_id) |>
        dplyr::mutate(
                raw_reads = n_reads[stage == "Raw"],
                after_decontam_reads = n_reads[stage == "After decontam"],
                bacteria_only_reads = n_reads[stage == "Bacteria only"],
                
                fraction_removed_decontam =
                        (raw_reads - after_decontam_reads) / raw_reads,
                
                fraction_removed_off_target =
                        (after_decontam_reads - bacteria_only_reads) / raw_reads,
                
                fraction_retained =
                        bacteria_only_reads / raw_reads
        ) |> 
        dplyr::ungroup() |>
        dplyr::left_join(
                phyloseq_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        rownames_to_column(var = "sample_id") |>
                        dplyr::mutate(sample_id = stringr::str_trim(as.character(sample_id))),
                by = "sample_id"
        ) |>
        tidyr::pivot_longer(
                cols = c(
                        fraction_removed_decontam,
                        fraction_removed_off_target,
                        fraction_retained
                ),
                names_to = "component",
                values_to = "fraction_of_raw"
        ) |>
        dplyr::mutate(
                component = factor(
                        component,
                        levels = c(
                                "fraction_removed_decontam",
                                "fraction_removed_off_target",
                                "fraction_retained"
                        ),
                        labels = c(
                                "Removed by decontam",
                                "Removed as off-target",
                                "Retained (Bacteria only)"
                        )
                )
        )

# Helper function to plot decontamination and off-target filtering results
plot_read_retention_stacked_16S <-
        function(grouping_variable,
                 grouping_levels = NULL,
                 filter_function = NULL) {
                
                read_retention_loss_table_16S |>
                        (\(x) {
                                if (!is.null(filter_function)) {
                                        filter_function(x)
                                } else {
                                        x
                                }
                        })() |>
                        (\(x) {
                                if (!is.null(grouping_levels)) {
                                        dplyr::mutate(
                                                x,
                                                !!grouping_variable :=
                                                        factor(
                                                                .data[[grouping_variable]],
                                                                levels = grouping_levels
                                                        )
                                        )
                                } else {
                                        x
                                }
                        })() |>
                        dplyr::group_by(
                                .data[[grouping_variable]],
                                component
                        ) |>
                        dplyr::summarise(
                                mean_fraction_of_raw =
                                        mean(fraction_of_raw, na.rm = TRUE),
                                .groups = "drop"
                        ) |>
                        ggplot(
                                aes(
                                        x = .data[[grouping_variable]],
                                        y = mean_fraction_of_raw,
                                        fill = component
                                )
                        ) +
                        geom_bar(stat = "identity") +
                        labs(
                                x = grouping_variable,
                                y = "Mean fraction of raw reads",
                                fill = "Processing outcome"
                        ) +
                        theme_bw() +
                        theme(
                                axis.text.x = element_text(
                                        angle = 45,
                                        hjust = 1
                                ),
                                legend.position = "right"
                        )
        }

# ---- Do decontamination and off-target sequences removal affect experimentla groups and phenotypes differently? ----
# Plot decontamination and off-target filtering results by sublibrary_id
plot_read_retention_stacked_16S(grouping_variable = "sublibrary_id")

# Plot decontamination and off-target filtering results by sample type
plot_read_retention_stacked_16S(
        grouping_variable = "sample_type",
        grouping_levels = c("extraction_control",
                            "negative_pcr_control",
                            "mock_community",
                            "bagged_flower",
                            "unbagged_flower"
        )
)

# Plot decontamination and off-target filtering results by management_type (only biological samples)
plot_read_retention_stacked_16S(
        grouping_variable = "management_type",
        grouping_levels = c("inside_forest", "near_forest", "agroforest", "full_sun"),
        filter_function =
                \(x) dplyr::filter(
                        x,
                        sample_type %in% c("bagged_flower", "unbagged_flower") &
                                !is.na(sample_type)
                )
)

# Plot decontamination and off-target filtering results by farm_id (only biological samples)
plot_read_retention_stacked_16S(
        grouping_variable = "farm_id",
        grouping_levels = c("ib", "vr", "sa", "kk", "mt", "vi", "yb"),
        filter_function =
                \(x) dplyr::filter(
                        x,
                        sample_type %in% c("bagged_flower", "unbagged_flower") &
                                !is.na(sample_type)
                        )
        )

# Plot decontamination and off-target filtering results by is_pollination_clsm (only phenotyped flowers)
plot_read_retention_stacked_16S(
        grouping_variable = "is_pollination_clsm",
        grouping_levels = c(FALSE, TRUE),
        filter_function =
                \(x) dplyr::filter(
                        x,
                        !is.na(is_pollination_clsm)
                )
)

# Plot decontamination and off-target filtering results by is_germination_clsm (only phenotyped flowers)
plot_read_retention_stacked_16S(
        grouping_variable = "is_germination",
        grouping_levels = c(FALSE, TRUE),
        filter_function =
                \(x) dplyr::filter(
                        x,
                        !is.na(is_germination)
                )
)

# Does sequencing depth confound pollination intensity? No
read_retention_loss_table_16S |>
        dplyr::filter(
                !is.na(pi_clsm),
                pi_clsm > 0,
                component == "Retained (Bacteria only)"
        ) |>
        ggplot(
                aes(
                        x = pi_clsm,
                        y = fraction_of_raw
                )
        ) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "loess", se = TRUE) +
        labs(
                x = "Number of pollen grains (pi_clsm)",
                y = "Fraction of raw reads retained"
        ) +
        theme_bw()

# Does sequencing depth confound callose_plug_index? No
read_retention_loss_table_16S |>
        dplyr::filter(
                !is.na(callose_plug_index),
                callose_plug_index > 0,
                component == "Retained (Bacteria only)"
        ) |>
        ggplot(
                aes(
                        x = callose_plug_index,
                        y = fraction_of_raw
                )
        ) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "loess", se = TRUE) +
        labs(
                x = "Callose plug index",
                y = "Fraction of raw reads retained"
        ) +
        theme_bw()


# ---- Are my experimental and phenotype groups confounded by sequencing depth? ----
# sample_type
dplyr::bind_rows(
        phyloseq_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
        ) |> 
        dplyr::left_join(
                phyloseq_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        rownames_to_column(var = "sample_id"),
                by = "sample_id"
        ) |> 
        dplyr::mutate(sample_type = factor(sample_type,
                levels = c("extraction_control",
                           "negative_pcr_control",
                           "mock_community",
                           "bagged_flower",
                           "unbagged_flower"
                )
        )) |> 
        ggplot(aes(x = sample_type, y = n_reads)) +
        geom_boxplot() +
        theme_bw()

# sublibrary_id
dplyr::bind_rows(
        phyloseq_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
        ) |> 
        dplyr::left_join(
                phyloseq_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        rownames_to_column(var = "sample_id"),
                by = "sample_id"
        ) |> 
        ggplot(aes(x = sublibrary_id, y = n_reads, colour = sample_type)) +
        geom_boxplot() +
        theme_bw() +
        theme(
                axis.text.x = element_text(
                        angle = 45,
                        hjust = 1
                ),
                legend.position = "bottom"
        )

# bagged Vs unbagged
dplyr::bind_rows(
        phyloseq_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |> 
        dplyr::left_join(
                phyloseq_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        rownames_to_column(var = "sample_id"),
                by = "sample_id"
        ) |> 
        filter(sample_type %in% c("bagged_flower", "unbagged_flower")) |> 
        ggplot(aes(x = sample_type, y = n_reads)) +
        geom_boxplot() +
        theme_bw()

# farm_id coloured by management_type
dplyr::bind_rows(
        phyloseq_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
        ) |> 
        dplyr::left_join(
                phyloseq_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        rownames_to_column(var = "sample_id"),
                by = "sample_id"
        ) |> 
        filter(sample_type %in% c("bagged_flower", "unbagged_flower")) |> 
        mutate(farm_id = factor(farm_id,
                                levels = c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
                                )
               ) |> 
        ggplot(aes(x = farm_id, y = n_reads, colour = management_type)) +
        geom_boxplot() +
        theme_bw() +
        theme(
                axis.text.x = element_text(
                        angle = 45,
                        hjust = 1
                ),
                legend.position = "bottom"
        )

# is_pollination_clsm
dplyr::bind_rows(
        phyloseq_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |> 
        dplyr::left_join(
                phyloseq_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        rownames_to_column(var = "sample_id"),
                by = "sample_id"
        ) |> 
        filter(!is.na(is_pollination_clsm)) |> 
        ggplot(aes(x = is_pollination_clsm, y = n_reads)) +
        geom_boxplot() +
        theme_bw()

# is_pollination_clsm
dplyr::bind_rows(
        phyloseq_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |> 
        dplyr::left_join(
                phyloseq_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        rownames_to_column(var = "sample_id"),
                by = "sample_id"
        ) |> 
        filter(!is.na(is_germination)) |> 
        ggplot(aes(x = is_germination, y = n_reads)) +
        geom_boxplot() +
        theme_bw()
        
# ---- Is there correlation between DNA concentration and sequencing depth? ----
# DNA was measured using Qubit after locus-specific PCR, bead cleaning & pooling triplicates
# Sequencing depth here is the number of reads per sample in data set after decontamination and off-targer sequences filtering
dplyr::bind_rows(
        phyloseq_16S_bacteria |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
        ) |> 
        dplyr::left_join(
                phyloseq_16S_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        rownames_to_column(var = "sample_id"),
                by = "sample_id"
        ) |> 
        ggplot(aes(x = dna_conc_hs_16s, y = n_reads, color = sample_type)) +
        geom_point() +
        theme_bw()


# ---- Examine sequencing depth ----
# All samples, including all type of controls
# 
phyloseq_16S_bacteria_biosamples |> 
        phyloseq::sample_sums() |>
        tibble::enframe(name = "sample_id",
                        value = "n_seqs") |> 
        pull(n_seqs) |> 
        summary()

# geom_histogram()
phyloseq_16S_bacteria_biosamples |> 
        phyloseq::sample_sums() |>
        tibble::enframe(name = "sample_id",
                        value = "n_seqs") |> 
        dplyr::left_join(
                rownames_to_column(metadata, var = "sample_id"),
                by = "sample_id") |> 
        filter(sample_id == "cib73") |> 
        ggplot(aes(x = n_seqs)) +
        geom_histogram(binwidth = 500) #+
        #coord_cartesian(xlim = c(0, 5000))

# geom_jitter()
phyloseq_16S_bacteria_biosamples |> 
        phyloseq::sample_sums() |>
        tibble::enframe(name = "sample_id",
                        value = "n_seqs") |> 
        ggplot(aes(x = 1, y = n_seqs)) +
        geom_jitter() +
        geom_hline(yintercept = 3000, color = "red") +
        geom_hline(yintercept = 5000, color = "blue") +
        geom_hline(yintercept = 10000, color = "green") +
        scale_y_log10()

# geom_boxplot()
phyloseq_16S_bacteria_biosamples |> 
        phyloseq::sample_sums() |>
        tibble::enframe(name = "sample_id",
                        value = "n_seqs") |> 
        dplyr::left_join(
                rownames_to_column(metadata, var = "sample_id"),
                by = "sample_id") |> 
        mutate(management_type = factor(management_type,
                                        levels = c("inside_forest", 
                                                   "near_forest",
                                                   "agroforest",
                                                   "full_sun"))) |> 
        ggplot(aes(x = sample_type, y = n_seqs, fill = management_type)) +
        geom_boxplot()

# geom_line() on arrange data
phyloseq_16S_bacteria |> 
        phyloseq::sample_sums() |>
        tibble::enframe(name = "sample_id",
                        value = "n_seqs") |> 
        arrange(n_seqs) %>% 
        ggplot(aes(x = 1:nrow(.), y = n_seqs)) +
        geom_line()

# arrange() and print()
phyloseq_16S_bacteria |> 
        phyloseq::sample_sums() |>
        tibble::enframe(name = "sample_id",
                        value = "n_seqs") |> 
        arrange(n_seqs) %>% 
        print(n = 190)


# ---- Diagnostic rarefaction curves ----
# Create a matrix for building rarefaction curves and transpose it because vegan::rarecurve() requires samples in rows
asv_mat_16S <-  phyloseq_16S_bacteria_biosamples |> 
        otu_table() |> 
        #taxa_are_rows() |> 
        as("matrix") |> 
        # Samples should be in rows for vegan::rarefy()
        t()

# Minimal seq depth could be used for 
min_depth <- min(rowSums(asv_mat_16S)) 

# Color coding using categorical variables from metadata
color_var <- "sample_type"   
metadata_for_rarecurve <- as.data.frame(sample_data(phyloseq_16S_bacteria_biosamples))
metadata_for_rarecurve <- metadata_for_rarecurve[rownames(asv_mat_16S), , drop = FALSE]

groups <- metadata_for_rarecurve[[color_var]] |>
        as.character() |>
        dplyr::coalesce("Not phenotyped") |>
        factor()

pal <- setNames(
        RColorBrewer::brewer.pal(nlevels(groups), "Set1"),
        levels(groups)
)
pal["Not phenotyped"] <- "grey25"

curve_cols <- pal[groups]

# Rarefaction curves
vegan::rarecurve(
        asv_mat_16S,
        step   = 500,
        # sample = 10000,
        col    = curve_cols,
        label  = FALSE,
        xlab   = "Sequencing depth",
        ylab   = "Observed ASVs",
        main   = paste("Rarefaction curves coloured by", color_var)
)

legend(
        "bottomright",
        legend = levels(groups),
        col    = pal,
        lwd    = 2,
        bty    = "n"
)



################################################################################
############## Alpha Diversity: diversity within the sample ####################
######## Composition (Richness) + Structure (Dominance & Evenness) #############
################################################################################
# Richness
?vegan::specnumber()
# Other indicies: Shannon (1 - ShannonMetric), Simpson, InvSimpson
?vegan::diversity()

# Alternatives for estimating alpha, beta and gamma diversities in hierarchical settings. NB!
# Additive Diversity Partitioning and Hierarchical Null Model Testing
?vegan::adipart()
# Multiplicative Diversity Partitioning: Alpha diversity compared to the total diversity (Gamma)
?vegan::multipart()

# ---- Iterative rarefaction and alpha diversity computing ----
# ---- Helper: rarefy + alpha once ----
rarefy_and_alpha_once <- function(ps_obj, depth, seed = NULL) {
        
        if (!is.null(seed)) set.seed(seed)
        
        # Rarefy
        ps_rare <- phyloseq::rarefy_even_depth(
                ps_obj,
                sample.size = depth,
                rngseed = seed,
                replace = FALSE,
                trimOTUs = TRUE,
                verbose = FALSE
        )
        
        # Alpha diversity
        alpha <- estimate_richness(
                ps_rare,
                measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")
        ) |>
                rownames_to_column("sample_id")
        
        alpha
}

# ---- Main: iterative rarefaction ----
iterative_rarefaction_alpha <- function(ps_obj, depth, n_iter = 50, base_seed = 42) {
        
        results <- map_dfr(seq_len(n_iter), function(i) {
                
                seed_i <- base_seed + i
                
                alpha_i <- rarefy_and_alpha_once(
                        ps_obj = ps_obj,
                        depth = depth,
                        seed  = seed_i
                )
                
                alpha_i |>
                        mutate(iteration = i)
        })
        
        results
}

# ---- Iterative rarefaction and alpha diversity at multiple depths ----
iterative_rarefaction_alpha_multi_depth <- function(ps_obj,
                                                    depths,
                                                    n_iter = 50,
                                                    base_seed = 42) {
        
        purrr::map_dfr(depths, function(depth_i) {
                
                alpha_depth_i <- iterative_rarefaction_alpha(
                        ps_obj   = ps_obj,
                        depth    = depth_i,
                        n_iter   = n_iter,
                        base_seed = base_seed
                )
                
                alpha_depth_i |>
                        mutate(depth = depth_i)
        })
}

# ---- Run iterative rarefaction at multiple depths: 5k, 10k, 15k ----
set.seed(42)

alpha_iters100_16S_multi_depths <- iterative_rarefaction_alpha_multi_depth(
        ps_obj  = phyloseq_16S_bacteria_biosamples,
        depths  = c(1000, 5000, 10000, 15000),
        n_iter  = 100,
        base_seed = 42
)

# Attach metadata once
alpha_iters100_16S_multi_depths <- alpha_iters100_16S_multi_depths |>
        left_join(
                data.frame(sample_data(phyloseq_ITS1_fungi_biosamples)) |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        rownames_to_column("sample_id"),
                by = "sample_id"
        )

# Save itteratively rarefied 
saveRDS(
        alpha_iters100_16S_multi_depths,
        file = "alpha_iters100_ITS1_bundle.rds"
)

# Long format
alpha_iters100_16S_multi_long <- alpha_iters100_16S_multi_depths |>
        pivot_longer(
                cols = c(Observed, Chao1, Shannon, Simpson, InvSimpson),
                names_to = "metric",
                values_to = "value"
        )

# ---- Visualise iterative subsampling results ----
# Diagnostiс of convergence of alpha-diversity metrics in iterations 
alpha_convergence <- alpha_iters100_16S_multi_long |>
        filter(metric %in% c("Observed", "Shannon", "InvSimpson")) |>
        group_by(depth, metric, iteration) |>
        summarise(
                median_value = median(value, na.rm = TRUE),
                sd_value     = sd(value, na.rm = TRUE),
                .groups = "drop"
        )

ggplot(alpha_convergence,
       aes(x = iteration, y = median_value, colour = metric)) +
        geom_line() +
        facet_wrap(~ depth, scales = "free_y") +
        theme_bw() +
        labs(
                x = "Number of rarefaction iterations",
                y = "Median alpha diversity",
                colour = "Metric"
        )

ggplot(alpha_convergence,
       aes(x = iteration, y = sd_value, colour = metric)) +
        geom_line() +
        facet_wrap(~ depth, scales = "free_y") +
        theme_bw() +
        labs(
                x = "Number of rarefaction iterations",
                y = "SD across iterations",
                colour = "Metric"
        )

alpha_iters100_16S_multi_long |>
        dplyr::filter(metric %in% c("Observed", "Shannon", "InvSimpson")) |>
        dplyr::group_by(depth, iteration, sample_type, metric) |>
        dplyr::summarise(
                median_value = median(value, na.rm = TRUE),
                .groups = "drop"
        ) |>
        ggplot(aes(x = iteration, y = median_value, colour = sample_type)) +
        geom_line() +
        facet_grid(metric ~ depth, scales = "free_y") +
        theme_bw() +
        labs(
                x = "Number of rarefaction iterations",
                y = "Median alpha diversity",
                colour = "Sample type"
        )

alpha_iters100_16S_multi_long |>
        dplyr::filter(metric %in% c("Observed", "Shannon", "InvSimpson")) |>
        dplyr::group_by(depth, iteration, sample_type, metric) |>
        dplyr::summarise(
                sd_value = sd(value, na.rm = TRUE),
                .groups = "drop"
        ) |>
        ggplot(aes(x = iteration, y = sd_value, colour = sample_type)) +
        geom_point() +
        facet_grid(metric ~ depth, scales = "free_y") +
        theme_bw() +
        labs(
                x = "Number of rarefaction iterations",
                y = "SD of alpha diversity",
                colour = "Sample type"
        )

alpha_iters100_16S_multi_long |>
        dplyr::filter(metric %in% c("Observed", "Shannon", "InvSimpson")) |>
        dplyr::group_by(depth, iteration, sample_id, sample_type, metric) |>
        dplyr::summarise(
                median_value = median(value, na.rm = TRUE),
                .groups = "drop"
        ) |>
        ggplot(aes(x = iteration, y = median_value, colour = sample_type)) +
        geom_point(alpha = 0.4, size = 0.8) +
        facet_grid(metric ~ depth, scales = "free_y") +
        theme_bw() +
        labs(
                x = "Number of rarefaction iterations",
                y = "Median alpha diversity per sample",
                colour = "Sample type"
        )

alpha_iters100_16S_multi_long |>
        dplyr::filter(metric %in% c("Observed", "Shannon", "InvSimpson")) |>
        dplyr::filter(sample_id %in% sample(unique(sample_id), 20)) |>
        ggplot(aes(x = iteration, y = value, group = sample_id, colour = sample_type)) +
        geom_line(alpha = 0.8) +
        facet_grid(metric ~ depth, scales = "free_y") +
        theme_bw() +
        labs(
                x = "Number of rarefaction iterations",
                y = "Alpha diversity",
                colour = "Sample type"
        )


alpha_iters100_16S_multi_depths |>
        dplyr::distinct(depth, iteration, sample_id) |>
        dplyr::group_by(depth) |>
        dplyr::summarise(
                min_n_samples = min(dplyr::n()),
                max_n_samples = max(dplyr::n()),
                .groups = "drop"
        )


# ---- Rarefying for richness-based alpha metrics without iterative subsampling ----
set.seed(42)

physeq_rare_1000 <- rarefy_even_depth(
        phyloseq_prev2_ab0,
        sample.size = 1000,
        rngseed = 42,
        replace = FALSE,
        trimOTUs = TRUE,
        verbose = FALSE
)

# ---- Alpha diversity metrics ----

# 1) Richness-based metrics (require rarefaction)
alpha_richness <- estimate_richness(
        physeq_rare_1000,
        measures = c("Observed", "Chao1")
)

# 2) Evenness-sensitive metrics (can be computed on non-rarefied data)
alpha_evenness <- estimate_richness(
        phyloseq_prev2_ab0,
        measures = c("Shannon", "Simpson")
)

# 3) Combine and attach metadata
alpha_df <- alpha_richness |>
        rownames_to_column("sample_id") |>
        left_join(
                alpha_evenness |>
                        rownames_to_column("sample_id"),
                by = "sample_id"
        ) |>
        left_join(
                data.frame(sample_data(phyloseq_prev2_ab0)) |>
                        rownames_to_column("sample_id"),
                by = "sample_id"
        )

head(alpha_df)

# ---- Alpha diversity summary ----
# helper: compact summary for multiple metrics
summarise_alpha <- function(df, group_var) {
        df |>
                group_by({{ group_var }}) |>
                summarise(
                        n = n(),
                        Observed_mean = mean(Observed, na.rm = TRUE),
                        Observed_sd   = sd(Observed, na.rm = TRUE),
                        Observed_med  = median(Observed, na.rm = TRUE),
                        Shannon_mean  = mean(Shannon, na.rm = TRUE),
                        Shannon_sd    = sd(Shannon, na.rm = TRUE),
                        Shannon_med   = median(Shannon, na.rm = TRUE),
                        Simpson_mean  = mean(Simpson, na.rm = TRUE),
                        Simpson_sd    = sd(Simpson, na.rm = TRUE),
                        Simpson_med   = median(Simpson, na.rm = TRUE),
                        .groups = "drop"
                ) |>
                arrange(desc(n))
}

# A) Flower treatment (all samples in alpha_df)
alpha_by_treatment <- alpha_df |>
        mutate(flower_treatment = sample_type) |>
        summarise_alpha(flower_treatment)

View(alpha_by_treatment)

# B) CLSM pollination (exclude NA; optionally keep only is_clsm == TRUE if you want)
alpha_by_clsm <- alpha_df |>
        filter(!is.na(is_pollination_clsm)) |>
        summarise_alpha(is_pollination_clsm)

View(alpha_by_clsm)

# ---- Plot alpha diversity ----
# Prepare alpha diversity data for plotting (whole dataset)
alpha_long <- alpha_df |>
        select(
                sample_id,
                sample_type,
                Observed,
                Shannon,
                Simpson
        ) |>
        pivot_longer(
                cols = c(Observed, Shannon, Simpson),
                names_to = "metric",
                values_to = "value"
        )

alpha_long$metric <- factor(
        alpha_long$metric,
        levels = c("Observed", "Shannon", "Simpson")
        )

# Alpha diversity boxplots (whole dataset)
ggplot(alpha_long, aes(x = sample_type, y = value, fill = sample_type)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
        facet_wrap(~ metric, scales = "free_y") +
        theme_bw() +
        labs(
                title = "Alpha diversity of bacterial communities",
                subtitle = "Observed richness rarefied at 1,000 reads; Shannon and Simpson on relative abundances",
                x = "",
                y = ""
        ) +
        theme(
                legend.position = "none",
                strip.background = element_rect(fill = "grey90"),
                strip.text = element_text(face = "bold")
        )

# Prepare alpha diversity data for CLSM-verified subset
alpha_long_clsm <- alpha_df |>
        filter(!is.na(is_pollination_clsm)) |>
        select(
                sample_id,
                is_pollination_clsm,
                Observed,
                Shannon,
                Simpson
        ) |>
        pivot_longer(
                cols = c(Observed, Shannon, Simpson),
                names_to = "metric",
                values_to = "value"
        )

alpha_long_clsm$metric <- factor(
        alpha_long_clsm$metric,
        levels = c("Observed", "Shannon", "Simpson")
)

# Alpha diversity boxplots (CLSM subset)
ggplot(
        alpha_long_clsm,
        aes(x = is_pollination_clsm, y = value, fill = is_pollination_clsm)
) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
        facet_wrap(~ metric, scales = "free_y") +
        theme_bw() +
        labs(
                title = "Alpha diversity of bacterial communities (CLSM subset)",
                subtitle = "Observed richness rarefied at 1,000 reads",
                x = "Microscopy-confirmed pollination",
                y = ""
        ) +
        theme(
                legend.position = "none",
                strip.background = element_rect(fill = "grey90"),
                strip.text = element_text(face = "bold")
        )

# ---- Linear Mixed-Effects Models for Alpha Diversity ----
# I have nested design farm(7) -> tree (7) -> flower (3 per group)
# management_type -> farm_id -> tree_id -> sample_id

# Observed richness (often skewed -> log-transform)
alpha_df <- alpha_df |>
        mutate(
                log_Observed = log1p(Observed)
        )

# LMM: Observed richness
lmm_obs_treat <- lmer(
        log_Observed ~ sample_type + (1 | tree_id),
        data = alpha_df
)
summary(lmm_obs_treat)

# Shannon
lmm_shannon_treat <- lmer(
        Shannon ~ sample_type + (1 | tree_id),
        data = alpha_df
)
summary(lmm_shannon_treat)

# Simpson
lmm_simpson_treat <- lmer(
        Simpson ~ sample_type + (1 | tree_id),
        data = alpha_df
)
summary(lmm_simpson_treat)

# CLSM-verified subset
alpha_clsm_df <- alpha_df |>
        filter(!is.na(is_pollination_clsm))

# Observed
lmm_obs_clsm <- lmer(
        log_Observed ~ is_pollination_clsm + (1 | tree_id),
        data = alpha_clsm_df
)
summary(lmm_obs_clsm)

# Shannon
lmm_shannon_clsm <- lmer(
        Shannon ~ is_pollination_clsm + (1 | tree_id),
        data = alpha_clsm_df
)
summary(lmm_shannon_clsm)

# Simpson
lmm_simpson_clsm <- lmer(
        Simpson ~ is_pollination_clsm + (1 | tree_id),
        data = alpha_clsm_df
)
summary(lmm_simpson_clsm)

# ---- Strict subset !! ----
alpha_df <- alpha_df |>
        mutate(
                visitation_pollination_group = case_when(
                        sample_type == "bagged_flower" &
                                is_pollination_clsm == FALSE ~ "bagged_unpollinated",
                        
                        sample_type == "unbagged_flower" &
                                is_pollination_clsm == TRUE ~ "open_pollinated",
                        
                        TRUE ~ NA_character_
                )
        )

alpha_strict <- alpha_df |>
        filter(!is.na(visitation_pollination_group))

table(alpha_strict$visitation_pollination_group)

alpha_strict |>
        count(visitation_pollination_group) |>
        mutate(prop = n / sum(n))

View(metadata)
unique(metadata$extraction_batch)
unique(metadata$`16s_pcr_batch`)
View(phyloseq_prev2_ab0)

# ---- Batch vs group association checks ----
meta <- data.frame(sample_data(phyloseq_prev2_ab0)) |>
        rownames_to_column("sample_id") |>
        mutate(
                extraction_batch = as.factor(extraction_batch),
                `X16s_pcr_batch`  = as.factor(`X16s_pcr_batch`),
                farm_id = as.factor(farm_id),
                sample_type = as.factor(sample_type)
        )

# Helper: Cramer's V
cramers_v <- function(tab) {
        chisq <- suppressWarnings(chisq.test(tab))
        n <- sum(tab)
        k <- min(nrow(tab), ncol(tab))
        sqrt(as.numeric(chisq$statistic) / (n * (k - 1)))
}

# 1) Extraction batch vs farm_id / sample_type
tab_ext_farm  <- table(meta$extraction_batch, meta$farm_id)
tab_ext_type  <- table(meta$extraction_batch, meta$sample_type)

chisq_ext_farm <- suppressWarnings(chisq.test(tab_ext_farm))
chisq_ext_type <- suppressWarnings(chisq.test(tab_ext_type))

v_ext_farm <- cramers_v(tab_ext_farm)
v_ext_type <- cramers_v(tab_ext_type)

list(
        chisq_ext_farm = chisq_ext_farm,
        cramersV_ext_farm = v_ext_farm,
        chisq_ext_type = chisq_ext_type,
        cramersV_ext_type = v_ext_type
)

# 2) PCR batch vs farm_id / sample_type
tab_pcr_farm <- table(meta$`X16s_pcr_batch`, meta$farm_id)
tab_pcr_type <- table(meta$`X16s_pcr_batch`, meta$sample_type)

chisq_pcr_farm <- suppressWarnings(chisq.test(tab_pcr_farm))
chisq_pcr_type <- suppressWarnings(chisq.test(tab_pcr_type))

v_pcr_farm <- cramers_v(tab_pcr_farm)
v_pcr_type <- cramers_v(tab_pcr_type)

list(
        chisq_pcr_farm = chisq_pcr_farm,
        cramersV_pcr_farm = v_pcr_farm,
        chisq_pcr_type = chisq_pcr_type,
        cramersV_pcr_type = v_pcr_type
)

# 3) Readable proportion tables (good for reporting)
prop_ext_type <- prop.table(tab_ext_type, margin = 2)  # within each sample_type
prop_pcr_type <- prop.table(tab_pcr_type, margin = 2)

prop_ext_type
prop_pcr_type





################################################################################
################# Beta Diversity: diversity between samples ####################
################################################################################

# There are 2 functions to compute distances in vegan package:
# a) vegan::vegdist() does not rarefy data and allows to specify variety of dissimilarity metrics/distances:
?vegan::vegdist()
distanceMethodList
# b) vegan::avgdist() performs rarefaction to control uneven sampling effort: 
# b.1) sample parameter specify rarefaction depth; 
# b.2) iterations of rarefaction (def. 100) to calculate average or median dissimilarity metric.
?vegan::avgdist()

# Function distance() from phyloseq package depends on vegdist() function from vegan.
# It also utilise dist() function from stats package. 
# It operates on phyloseq-class object.
?phyloseq::distance()
?stats::dist()

# Bray-Curtis performed after a ‘double standardization’ procedure, 
# which gives equal weight to all species and samples (Goodall 1978).
# According to Bray & Curtis (1957), values for species i are first 
# standardized by scoring abundances as a percentage of the maximum value
# attained by species i over all samples. For the second standardization,
# the species scores are then rescored as a percentage of the sample totals. 

# in vegan documentation they mention alternative names for method = bray:
# Czekanowski or Steinhaus (standartization only by sample sums)
distanceMethodList

# How to account for uneven sampling depth?
# Rarefaction or Normalization

# ---- Beta-diversity: BC dissimilarity ----
# bounded between 0 (the same communities) and 1 (no matches)
# relative abundance based => use Jaccard for presence/abscence analysis
# non-phylogenetic => use UniFrac after MSA & phylogeny building
# BC dissimilarity is directly related to the quantitative Sorensen similarity index 

compute_bray_curtis <- function(physeq) {
        physeq_rel <- transform_sample_counts(
                physeq,
                function(x) x / sum(x)
        )

        phyloseq::distance(physeq_rel, method = "bray")
}

bc_dissim_no_prev <- compute_bray_curtis(phyloseq_16S_bacteria_bio)
bc_dissim_prev1   <- compute_bray_curtis(phyloseq_prev1)
bc_dissim_prev2   <- compute_bray_curtis(phyloseq_prev2)
bc_dissim_prev3   <- compute_bray_curtis(phyloseq_prev3)
bc_dissim_prev9   <- compute_bray_curtis(phyloseq_prev9)

# Ordination: Principle Coordinate Analysis
ord_no_prev <- ordinate(phyloseq_16S_bacteria_bio, method = "PCoA", distance = bc_dissim_no_prev)
ord_prev1   <- ordinate(phyloseq_prev1, method = "PCoA", distance = bc_dissim_prev1)
ord_prev2   <- ordinate(phyloseq_prev2, method = "PCoA", distance = bc_dissim_prev2)
ord_prev3   <- ordinate(phyloseq_prev3, method = "PCoA", distance = bc_dissim_prev3)
ord_prev9   <- ordinate(phyloseq_prev9, method = "PCoA", distance = bc_dissim_prev9)

# Visualization PCoA
plot_ordination(
        phyloseq_16S_bacteria_bio,
        ord_no_prev,
        color = "sample_type"
) + ggtitle("Bray?Curtis: no prevalence filter")

plot_ordination(
        phyloseq_prev1,
        ord_prev1,
        color = "sample_type"
) + ggtitle("Bray?Curtis: prevalence > 1")

plot_ordination(
        phyloseq_prev2,
        ord_prev2,
        color = "sample_type"
) + ggtitle("Bray?Curtis: prevalence > 2")

plot_ordination(
        phyloseq_prev3,
        ord_prev3,
        color = "sample_type"
) + ggtitle("Bray?Curtis: prevalence > 3")

plot_ordination(
        phyloseq_prev9,
        ord_prev9,
        color = "sample_type"
) + ggtitle("Bray?Curtis: prevalence > 9")

# Mantel test: BC dissimilarity robustness after different prevalence filtration
mantel_prev1 <- vegan::mantel(bc_dissim_no_prev, bc_dissim_prev1, method = "spearman")
mantel_prev2 <- vegan::mantel(bc_dissim_no_prev, bc_dissim_prev2, method = "spearman")
mantel_prev3 <- vegan::mantel(bc_dissim_no_prev, bc_dissim_prev3, method = "spearman")
mantel_prev9 <- vegan::mantel(bc_dissim_no_prev, bc_dissim_prev9, method = "spearman")

data.frame(
        comparison = c("no vs >1", "no vs >2", "no vs >3", "no vs >9"),
        rho = c(
                mantel_prev1$statistic,
                mantel_prev2$statistic,
                mantel_prev3$statistic,
                mantel_prev9$statistic
        ),
        p_value = c(
                mantel_prev1$signif,
                mantel_prev2$signif,
                mantel_prev3$signif,
                mantel_prev9$signif
        )
)

# ---- Read count filtering ----
# Filter taxa by total read count across all samples
filter_by_abundance <- function(physeq, min_reads) {
        prune_taxa(
                taxa_sums(physeq) >= min_reads,
                physeq
        )
}

phyloseq_prev2_ab0 <- phyloseq_prev2
phyloseq_prev2_ab10 <- filter_by_abundance(phyloseq_prev2, 10)
phyloseq_prev2_ab40 <- filter_by_abundance(phyloseq_prev2, 40)
phyloseq_prev2_ab150 <- filter_by_abundance(phyloseq_prev2, 150)

# Compare total number of taxa
data.frame(
        filtering = c("prev>2, no count", "prev>2, >=10", "prev>2, >=40", "prev>2, >=150"),
        ntaxa = c(
                ntaxa(phyloseq_prev2_ab0),
                ntaxa(phyloseq_prev2_ab10),
                ntaxa(phyloseq_prev2_ab40),
                ntaxa(phyloseq_prev2_ab150)
        )
)

# ---- Taxa per sample after different abundance filters ----
# ---- Prevalence filter (fixed) ----
count_taxa_per_sample <- function(physeq) {
        apply(
                otu_table(physeq) > 0,
                2,
                sum
        )
}

taxa_per_sample_abundance <- data.frame(
        sample_id = sample_names(phyloseq_prev2_ab0),
        ab0   = count_taxa_per_sample(phyloseq_prev2_ab0),
        ab10  = count_taxa_per_sample(phyloseq_prev2_ab10),
        ab40  = count_taxa_per_sample(phyloseq_prev2_ab40),
        ab150 = count_taxa_per_sample(phyloseq_prev2_ab150)
)

summary(taxa_per_sample_abundance[, -1])

# ---- Bray?Curtis dissimilarity ----
compute_bray_curtis <- function(physeq) {
        physeq_rel <- transform_sample_counts(
                physeq,
                function(x) x / sum(x)
        )
        phyloseq::distance(physeq_rel, method = "bray")
}

bc_prev2_ab0   <- compute_bray_curtis(phyloseq_prev2_ab0)
bc_prev2_ab10  <- compute_bray_curtis(phyloseq_prev2_ab10)
bc_prev2_ab40  <- compute_bray_curtis(phyloseq_prev2_ab40)
bc_prev2_ab150 <- compute_bray_curtis(phyloseq_prev2_ab150)

# ---- Mantel test: robustness to abundance filtering ----
mantel_ab10  <- vegan::mantel(bc_prev2_ab0, bc_prev2_ab10,  method = "spearman")
mantel_ab40  <- vegan::mantel(bc_prev2_ab0, bc_prev2_ab40,  method = "spearman")
mantel_ab150 <- vegan::mantel(bc_prev2_ab0, bc_prev2_ab150, method = "spearman")

data.frame(
        comparison = c("ab0 vs ab10", "ab0 vs ab40", "ab0 vs ab150"),
        rho = c(
                mantel_ab10$statistic,
                mantel_ab40$statistic,
                mantel_ab150$statistic
        ),
        p_value = c(
                mantel_ab10$signif,
                mantel_ab40$signif,
                mantel_ab150$signif
        )
)

# ---- PCoA ordinations ----
ord_prev2_ab0   <- ordinate(phyloseq_prev2_ab0,   "PCoA", bc_prev2_ab0)
ord_prev2_ab10  <- ordinate(phyloseq_prev2_ab10,  "PCoA", bc_prev2_ab10)
ord_prev2_ab40  <- ordinate(phyloseq_prev2_ab40,  "PCoA", bc_prev2_ab40)
ord_prev2_ab150 <- ordinate(phyloseq_prev2_ab150, "PCoA", bc_prev2_ab150)

# ---- PCoA visualization ----
plot_ordination(
        phyloseq_prev2_ab0,
        ord_prev2_ab0,
        color = "sample_type"
) + ggtitle("BC PCoA: prev>2, no count filter")

plot_ordination(
        phyloseq_prev2_ab10,
        ord_prev2_ab10,
        color = "sample_type"
) + ggtitle("BC PCoA: prev>2, count >=10")

plot_ordination(
        phyloseq_prev2_ab40,
        ord_prev2_ab40,
        color = "sample_type"
) + ggtitle("BC PCoA: prev>2, count >=40")

plot_ordination(
        phyloseq_prev2_ab150,
        ord_prev2_ab150,
        color = "sample_type"
) + ggtitle("BC PCoA: prev>2, count >=150")



# ---- Bray?Curtis (relative abundance based) and Jaccard (presence/absence) dissimilarity ----
compute_beta_distance <- function(physeq, method = c("bray", "jaccard")) {
        method <- match.arg(method)
        
        if (method == "bray") {
                physeq_rel <- transform_sample_counts(
                        physeq,
                        function(x) x / sum(x)
                )
                return(phyloseq::distance(physeq_rel, method = "bray"))
        }
        
        if (method == "jaccard") {
                physeq_pa <- transform_sample_counts(
                        physeq,
                        function(x) {
                                y <- as.integer(x > 0)
                                names(y) <- names(x)   # ????? ? ???? ????? ?????? ?????
                                y
                        }
                )
                return(phyloseq::distance(physeq_pa, method = "jaccard"))
        }
}

bc_prev2_ab0  <- compute_beta_distance(phyloseq_prev2_ab0, "bray")
jac_prev2_ab0 <- compute_beta_distance(phyloseq_prev2_ab0, "jaccard")

ord_bc  <- ordinate(phyloseq_prev2_ab0, "PCoA", bc_prev2_ab0)
ord_jac <- ordinate(phyloseq_prev2_ab0, "PCoA", jac_prev2_ab0)

plot_ordination(phyloseq_prev2_ab0, ord_bc,  color = "sample_type") +
        ggtitle("PCoA: Bray?Curtis (RA)")

plot_ordination(phyloseq_prev2_ab0, ord_jac, color = "sample_type") +
        ggtitle("PCoA: Jaccard (presence/absence)")

# PCoA by major variables
plot_pcoa_by <- function(physeq, ord, color_var, title_suffix) {
        plot_ordination(
                physeq,
                ord,
                color = color_var
        ) +
                theme_bw() +
                ggtitle(paste("PCoA (Bray?Curtis):", title_suffix)) +
                theme(
                        legend.title = element_text(size = 10),
                        legend.text  = element_text(size = 9)
                )
}

# ---- PCoA: categorical variables ----
plot_pcoa_by(
        phyloseq_prev2_ab0,
        ord_bc,
        "sample_type",
        "flower visitation (bagged vs open)"
)

plot_pcoa_by(
        phyloseq_prev2_ab0,
        ord_bc,
        "is_pollination_clsm",
        "CLSM-confirmed pollination"
)

plot_pcoa_by(
        phyloseq_prev2_ab0,
        ord_bc,
        "management_type",
        "management type"
)

plot_pcoa_by(
        phyloseq_prev2_ab0,
        ord_bc,
        "farm_id",
        "farm ID"
)
# Nested design
plot_ordination(
        phyloseq_prev2_ab0,
        ord_bc,
        color = "tree_id"
) +
        theme_bw() +
        ggtitle("PCoA (Bray?Curtis): tree-level clustering") +
        guides(color = guide_legend(ncol = 2))


# ---- Pollination intensity: Ordination dataframe + metadata ----
ord_df <- plot_ordination(
        phyloseq_prev2_ab0,
        ord_bc,
        justDF = TRUE
)

if (!"SampleID" %in% names(ord_df)) {
        ord_df <- ord_df |>
                rownames_to_column("SampleID")
}

meta_df <- data.frame(sample_data(phyloseq_prev2_ab0)) |>
        rownames_to_column("SampleID")

ord_df <- ord_df |>
        left_join(meta_df, by = "SampleID")

# ---- Gradient plot ----
# Pollination intensity
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = log1p(pi_clsm.x))) +
        geom_point(size = 2, alpha = 0.8) +
        scale_color_viridis_c(na.value = "grey80") +
        theme_bw() +
        labs(
                title = "PCoA (Bray?Curtis)",
                subtitle = "Gradient coloured by pollen intensity (pi_clsm)",
                x = "PCoA1",
                y = "PCoA2",
                color = "log1(pi_clsm)"
        )

# Callose plug index
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = germination_rate.x)) +
        geom_point(size = 2, alpha = 0.8) +
        scale_color_viridis_c(na.value = "grey80") +
        theme_bw() +
        labs(
                title = "PCoA (Bray?Curtis)",
                subtitle = "Gradient coloured by callose plug index ",
                x = "PCoA1",
                y = "PCoA2",
                color = "callose_plug_index"
        )

# ---- PERMANOVA Bray-Curtis & Jaccard (flower_treatment)----
# Metadata
meta_prev2 <- data.frame(sample_data(phyloseq_prev2_ab0))

# Beta-diversity distances
bc_prev2  <- compute_beta_distance(phyloseq_prev2_ab0, method = "bray")
jac_prev2 <- compute_beta_distance(phyloseq_prev2_ab0, method = "jaccard")

# PERMANOVA: Bray?Curtis
adonis_bc_sample_type <- adonis2(
        bc_prev2 ~ sample_type,
        data = meta_prev2,
        permutations = 999,
        by = "margin"
)

# PERMANOVA: Jaccard
adonis_jac_sample_type <- adonis2(
        jac_prev2 ~ sample_type,
        data = meta_prev2,
        permutations = 999,
        by = "margin"
)

adonis_bc_sample_type
adonis_jac_sample_type

# Homogeneity of dispersions (Bray?Curtis)
bd_bc <- betadisper(bc_prev2, meta_prev2$sample_type)

# Permutation test
permutest(bd_bc)

# Optional: ANOVA table
anova(bd_bc)

# Homogeneity of dispersions (Jaccard)
bd_jac <- betadisper(jac_prev2, meta_prev2$sample_type)

# Permutation test
permutest(bd_jac)

# Optional: ANOVA table
anova(bd_jac)

# Bray?Curtis dispersion
plot(bd_bc)
boxplot(bd_bc, main = "Bray?Curtis dispersion by sample_type")

# Jaccard dispersion
plot(bd_jac)
boxplot(bd_jac, main = "Jaccard dispersion by sample_type")

# ---- Microscopy-guided subset (CLSM) ----
# Start from filtered dataset: prevalence > 2, no abundance filter
phyloseq_prev2_ab0

# Keep only samples with CLSM metadata and non-missing pollination label
phyloseq_prev2_ab0_clsm <- phyloseq_prev2_ab0 |>
        subset_samples(is_clsm == TRUE) |> # this filtration layer could be omitted
        subset_samples(!is.na(is_pollination_clsm))

# Drop taxa with zero counts after subsetting
phyloseq_prev2_ab0_clsm <- prune_taxa(
        taxa_sums(phyloseq_prev2_ab0_clsm) > 0,
        phyloseq_prev2_ab0_clsm
)

# Quick checks
nsamples(phyloseq_prev2_ab0_clsm)
table(sample_data(phyloseq_prev2_ab0_clsm)$is_pollination_clsm, useNA = "ifany")
table(sample_data(phyloseq_prev2_ab0_clsm)$sample_type, useNA = "ifany")


# ---- Beta distances (Bray RA; Jaccard presence/absence) ----
compute_beta_distance <- function(physeq, method = c("bray", "jaccard")) {
        method <- match.arg(method)
        
        if (method == "bray") {
                physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))
                return(phyloseq::distance(physeq_rel, method = "bray"))
        }
        
        if (method == "jaccard") {
                # vegan::vegdist with binary=TRUE treats any positive count as presence
                return(phyloseq::distance(physeq, method = "jaccard", binary = TRUE))
        }
}

bc_clsm  <- compute_beta_distance(phyloseq_prev2_ab0_clsm, method = "bray")
jac_clsm <- compute_beta_distance(phyloseq_prev2_ab0_clsm, method = "jaccard")


# ---- Ordination (PCoA) ----
ord_bc_clsm  <- ordinate(phyloseq_prev2_ab0_clsm, method = "PCoA", distance = bc_clsm)
ord_jac_clsm <- ordinate(phyloseq_prev2_ab0_clsm, method = "PCoA", distance = jac_clsm)

# PCoA plots
plot_ordination(
        phyloseq_prev2_ab0_clsm,
        ord_bc_clsm,
        color = "is_pollination_clsm",
        shape = "sample_type"
) + ggtitle("PCoA (Bray?Curtis, RA): CLSM subset")

plot_ordination(
        phyloseq_prev2_ab0_clsm,
        ord_jac_clsm,
        color = "is_pollination_clsm",
        shape = "sample_type"
) + ggtitle("PCoA (Jaccard, presence/absence): CLSM subset")


# ---- PERMANOVA + betadisper ----
meta_clsm <- data.frame(sample_data(phyloseq_prev2_ab0_clsm))

# 1) Simple model: pollination status only
adonis_bc_poll <- vegan::adonis2(
        bc_clsm ~ is_pollination_clsm,
        data = meta_clsm,
        permutations = 999,
        by = "margin"
)

adonis_jac_poll <- vegan::adonis2(
        jac_clsm ~ is_pollination_clsm,
        data = meta_clsm,
        permutations = 999,
        by = "margin"
)

adonis_bc_poll
adonis_jac_poll

# Dispersion check (assumption check for PERMANOVA interpretation)
bd_bc_poll <- vegan::betadisper(bc_clsm, meta_clsm$is_pollination_clsm)
bd_jac_poll <- vegan::betadisper(jac_clsm, meta_clsm$is_pollination_clsm)

vegan::permutest(bd_bc_poll)
vegan::permutest(bd_jac_poll)

boxplot(bd_bc_poll, main = "Bray?Curtis dispersion by is_pollination_clsm")
boxplot(bd_jac_poll, main = "Jaccard dispersion by is_pollination_clsm")


# ---- Optional: adjust for sample_type (bagged vs open) ----
# This tests whether microscopy-confirmed pollination explains differences
# beyond the bagged/open label.
adonis_bc_poll_adj <- vegan::adonis2(
        bc_clsm ~ sample_type + is_pollination_clsm,
        data = meta_clsm,
        permutations = 999,
        by = "margin"
)

adonis_jac_poll_adj <- vegan::adonis2(
        jac_clsm ~ sample_type + is_pollination_clsm,
        data = meta_clsm,
        permutations = 999,
        by = "margin"
)

adonis_bc_poll_adj
adonis_jac_poll_adj


# ---- Optional: restrict permutations within tree_id (repeated flowers per tree) ----
# Use this if tree_id exists and you have multiple flowers per tree.
if ("tree_id" %in% names(meta_clsm)) {
        adonis_bc_poll_tree <- vegan::adonis2(
                bc_clsm ~ is_pollination_clsm,
                data = meta_clsm,
                permutations = 999,
                by = "margin",
                strata = meta_clsm$tree_id
        )
        
        adonis_jac_poll_tree <- vegan::adonis2(
                jac_clsm ~ is_pollination_clsm,
                data = meta_clsm,
                permutations = 999,
                by = "margin",
                strata = meta_clsm$tree_id
        )
        
        adonis_bc_poll_tree
        adonis_jac_poll_tree
}

# ---- Optional: add pi_clsm as a quantitative predictor ----
# This is useful if pi_clsm is numeric and not too sparse.
if ("pi_clsm" %in% names(meta_clsm)) {
        meta_clsm$pi_clsm <- as.numeric(meta_clsm$pi_clsm)
        
        adonis_bc_pi <- vegan::adonis2(
                bc_clsm ~ pi_clsm,
                data = meta_clsm,
                permutations = 999,
                by = "margin"
        )
        
        adonis_jac_pi <- vegan::adonis2(
                jac_clsm ~ pi_clsm,
                data = meta_clsm,
                permutations = 999,
                by = "margin"
        )
        
        adonis_bc_pi
        adonis_jac_pi
}

# ---- Alpha Diversity: diversity within the sample 
# Composition (Richness/Observed ASV) + structure (Dominance/Evenness)




################################################################################
################# Zeta Diversity: diversity between samples ####################
################################################################################
# I used zeta-diversity decline to answer the question:
# Do flower visitors homogenize or heterogenize flower bacterial microbiome?
# I expect zeta-decline curve for unbagged flowers fit:
# a) exp law in case of heterogenization effect of visitors;
# b) power law in case of homogenization effect of visitors.
# c) If curves are parallel the effect of visitors is neglegible and the same processes drive community assembly.

# There're two functions in zetadiv package for zeta-diversity decline:
?Zeta.decline.ex() # computes the expectation of zeta diversity fo
?Zeta.decline.mc()
?Zeta.decline()

# I transpose feature table and convert it into incidence matrix
# Convert ASV table into presence/absence incidence matrix
incidence_df_16S <-
        phyloseq_16S_bacteria_biosamples |>
        otu_table() |>
        as("matrix") |>
        t() |>
        (\(x) {
                x[x > 0] <- 1
                as.data.frame(x)
        })()

# Prepare metadata as dataframe to analyse zeta-decline by experimental and phenotypical groups
metadata_16S <-
        sample_data(phyloseq_16S_bacteria_biosamples) |>
        as("data.frame")

all(rownames(incidence_df_16S) ==
            rownames(metadata_16S))

# ---- Compute zeta-diversity decline (z1 - z10 components) on ASVs ----
# using Monte Carlo sampling method (1000 samples)
zeta_bagged_16S <-
        incidence_df_16S[
                metadata_16S$sample_type == "bagged_flower",
        ] |>
        Zeta.decline.mc(
                orders = 1:10,
                sam = 1000,
                plot = FALSE
        )

zeta_unbagged_16S <-
        incidence_df_16S[
                metadata_16S$sample_type == "unbagged_flower",
        ] |>
        Zeta.decline.mc(
                orders = 1:10,
                sam = 1000,
                plot = FALSE
        )
# using Monte Carlo sampling method (10000 samples)
zeta_bagged_16S_10k <-
        incidence_df_16S[
                metadata_16S$sample_type == "bagged_flower",
        ] |>
        Zeta.decline.mc(
                orders = 1:10,
                sam = 10000,
                plot = FALSE
        )

zeta_unbagged_16S_10k <-
        incidence_df_16S[
                metadata_16S$sample_type == "unbagged_flower",
        ] |>
        Zeta.decline.mc(
                orders = 1:10,
                sam = 10000,
                plot = FALSE
        )

# ---- Agregate zeta-diversity values for bagged and unbagged flowers from ASVs ----
zeta_plot_tbl_16S_10k <-
        dplyr::bind_rows(
                data.frame(
                        order = zeta_bagged_16S_10k$zeta.order,
                        zeta = zeta_bagged_16S_10k$zeta.val,
                        zeta_sd = zeta_bagged_16S_10k$zeta.val.sd,
                        group = "bagged_flower"
                ),
                data.frame(
                        order = zeta_unbagged_16S_10k$zeta.order,
                        zeta = zeta_unbagged_16S_10k$zeta.val,
                        zeta_sd = zeta_unbagged_16S_10k$zeta.val.sd,
                        group = "unbagged_flower"
                )
        )

#  ---- Visualize zeta-diversity decline ASVs ----
ggplot(zeta_plot_tbl_16S_10k,
       aes(x = order, y = zeta, colour = group)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 2) +
        geom_errorbar(
                aes(ymin = pmax(zeta - zeta_sd, .Machine$double.eps),
                    ymax = zeta + zeta_sd),
                width = 0.15,
                alpha = 0.9
        ) +
        theme_bw() +
        scale_y_log10() +
        scale_x_continuous(breaks = 1:10) +
        labs(
                x = "Zeta order (number of flowers)",
                y = "Zeta diversity (average number of shared ASVs)",
                title = "Zeta diversity decline curves for 16S ASVs",
                colour = "Flower treatment"
        )

# ---- Diagnostic commands for Zetadiv.decline.mc() output lists ----
# Compare coefficients of determination for power and exponential law fit:
# how well the fitted model explains the variance in the observed zeta decline
summary(zeta_bagged_16S$zeta.pl)$r.squared
summary(zeta_bagged_16S$zeta.exp)$r.squared

summary(zeta_unbagged_16S$zeta.pl)$r.squared
summary(zeta_unbagged_16S$zeta.exp)$r.squared

extract_zeta_fits <- function(zeta_result) {
        
        power_r2 <- summary(zeta_result$zeta.pl)$r.squared
        exp_r2   <- summary(zeta_result$zeta.exp)$r.squared
        
        power_coef <- coef(zeta_result$zeta.pl)
        exp_coef   <- coef(zeta_result$zeta.exp)
        
        list(
                power_law = list(
                        r2 = power_r2,
                        a_intercept  = 10^(power_coef[["(Intercept)"]]),
                        b_coef  = power_coef[["log10(c(orders))"]]
                ),
                exponential = list(
                        r2 = exp_r2,
                        a_intercept  = 10^(exp_coef[["(Intercept)"]]),
                        b_coef = exp_coef[["c(orders)"]]
                ),
                aic = zeta_result$aic
        )
}

extract_zeta_fits(zeta_bagged)
extract_zeta_fits(zeta_unbagged)


# ---- Agregate zeta-diversity values for bagged and unbagged flowers from ASVs ----
zeta_plot_tbl_16S <-
        dplyr::bind_rows(
                data.frame(
                        order = zeta_bagged_16S$zeta.order,
                        zeta = zeta_bagged_16S$zeta.val,
                        zeta_sd = zeta_bagged_16S$zeta.val.sd,
                        group = "bagged_flower"
                ),
                data.frame(
                        order = zeta_unbagged_16S$zeta.order,
                        zeta = zeta_unbagged_16S$zeta.val,
                        zeta_sd = zeta_unbagged_16S$zeta.val.sd,
                        group = "unbagged_flower"
                )
        )

#  ---- Visualize zeta-diversity decline ASVs ----
ggplot(zeta_plot_tbl_16S,
       aes(x = order, y = zeta, colour = group)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 2) +
        geom_errorbar(
                aes(ymin = pmax(zeta - zeta_sd, .Machine$double.eps),
                    ymax = zeta + zeta_sd),
                width = 0.15,
                alpha = 0.9
        ) +
        theme_bw() +
        scale_y_log10() +
        scale_x_continuous(breaks = 1:10) +
        labs(
                x = "Zeta order (number of flowers)",
                y = "Zeta diversity (average number of shared ASVs)",
                title = "Zeta diversity decline curves for 16S ASVs",
                colour = "Flower treatment"
        )

# ---- Zeta diversity decline on Genus level (sensitivity test & ASV noise reduction) ----
# tax_glom() with handling NAs to preserve matrix topology 
phyloseq_16S_bacteria_genus <-
        phyloseq_16S_bacteria_biosamples |>
        (\(physeq_16S) {
                
                # Extract taxonomy as matrix for direct manipulation
                taxonomy_matrix_16S <-
                        tax_table(physeq_16S) |>
                        as.matrix()
                
                # Identify ASVs with missing genus annotation
                genus_is_na <-
                        is.na(taxonomy_matrix_16S[, "Genus"])
                
                # Assign unique artificial genus labels to unclassified ASVs
                taxonomy_matrix_16S[genus_is_na, "Genus"] <-
                        paste0(
                                "Unclassified_g__ASV_",
                                rownames(taxonomy_matrix_16S)[genus_is_na]
                        )
                
                # Write modified taxonomy back into phyloseq object
                tax_table(physeq_16S) <-
                        tax_table(taxonomy_matrix_16S)
                
                # Agglomerate at genus level (counts are summed internally)
                phyloseq::tax_glom(
                        physeq_16S,
                        taxrank = "Genus",
                        NArm = FALSE
                )
                
        })()

# Samples
nsamples(phyloseq_16S_bacteria_genus)

tax_table(phyloseq_16S_bacteria_genus) |> View()

# How many ASVs were aglomerated into genera?
ntaxa(phyloseq_16S_bacteria_biosamples)  # ASV count
ntaxa(phyloseq_16S_bacteria_genus)                # genus count

# How efficient annotation by NB classifier QIIME2 (SILVA DB)?
sum(grepl("^Unclassified_g__ASV_", 
          tax_table(phyloseq_16S_bacteria_genus)[, "Genus"]))

# 
table(tax_table(phyloseq_16S_bacteria_biosamples)[, "Genus"])

# Convert ASV table into presence/absence incidence matrix
incidence_df_16S_genus <-
        phyloseq_16S_bacteria_genus |>
        otu_table() |>
        as("matrix") |>
        t() |>
        (\(x) {
                x[x > 0] <- 1
                as.data.frame(x)
        })()

# ---- Compute zeta-diversity decline (z1 - z10 components) on Genus level ----
# using Monte Carlo sampling method (1000 samples)
zeta_bagged_16S_genus <-
        incidence_df_16S_genus[
                metadata_16S$sample_type == "bagged_flower",
        ] |>
        Zeta.decline.mc(
                orders = 1:10,
                sam = 1000,
                plot = FALSE
        )

zeta_unbagged_16S_genus <-
        incidence_df_16S_genus[
                metadata_16S$sample_type == "unbagged_flower",
        ] |>
        Zeta.decline.mc(
                orders = 1:10,
                sam = 1000,
                plot = FALSE
        )

str(zeta_unbagged_16S_genus)
str(zeta_bagged_16S_genus)

# ---- Agregate zeta-diversity values for bagged and unbagged flowers from Genus ----
zeta_plot_tbl_16S_genus <-
        dplyr::bind_rows(
                data.frame(
                        order = zeta_bagged_16S_genus$zeta.order,
                        zeta = zeta_bagged_16S_genus$zeta.val,
                        zeta_sd = zeta_bagged_16S_genus$zeta.val.sd,
                        group = "bagged_flower"
                ),
                data.frame(
                        order = zeta_unbagged_16S_genus$zeta.order,
                        zeta = zeta_unbagged_16S_genus$zeta.val,
                        zeta_sd = zeta_unbagged_16S_genus$zeta.val.sd,
                        group = "unbagged_flower"
                )
        )

# ---- Visualize zeta-diversity decline Genus ----
ggplot(zeta_plot_tbl_16S_genus,
       aes(x = order, y = zeta, colour = group)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 2) +
        geom_errorbar(
                aes(ymin = pmax(zeta - zeta_sd, .Machine$double.eps),
                    ymax = zeta + zeta_sd),
                width = 0.15,
                alpha = 0.9
        ) +
        theme_bw() +
        scale_y_log10() +
        scale_x_continuous(breaks = 1:10) +
        labs(
                x = "Zeta order (number of flowers)",
                y = "Zeta diversity (average number of shared genera)",
                title = "Zeta diversity decline curves for 16S Genus",
                colour = "Flower treatment"
        )

################################################################################
############## A part of dirty wet exploratort code with some heatmaps #########
################################################################################
# ---- Heatmap  ----
# ---- Extract OTU matrix ----
otu_mat <- as(otu_table(phyloseq_prev2_ab0), "matrix")

if (taxa_are_rows(phyloseq_prev2_ab0)) {
        otu_mat <- t(otu_mat)
}

otu_df <- as.data.frame(otu_mat) |>
        rownames_to_column("sample_id")

# ---- Metadata ----
meta_df <- data.frame(sample_data(phyloseq_prev2_ab0)) |>
        rownames_to_column("sample_id") |>
        select(sample_id, sample_type)

# ---- Taxonomy (use Genus; switch to Family if needed) ----
tax_df_use <- tax_df_clean |>
        rownames_to_column("taxon_id") |>
        select(taxon_id, Genus)


# ---- Long format ----
long_df <- otu_df |>
        left_join(meta_df, by = "sample_id") |>
        pivot_longer(
                cols = -c(sample_id, sample_type),
                names_to = "taxon_id",
                values_to = "count"
        ) |>
        left_join(tax_df_use, by = "taxon_id")

# ---- Relative abundance per sample ----
long_df <- long_df |>
        group_by(sample_id) |>
        mutate(rel_abund = count / sum(count)) |>
        ungroup()

# ---- Keep most abundant genera for readability ----
top_genera <- mean_ra_genus |>
        group_by(Genus) |>
        summarise(total_ra = sum(mean_ra)) |>
        slice_max(total_ra, n = 50) |>
        pull(Genus)

mean_ra_genus <- mean_ra_genus |>
        filter(Genus %in% top_genera)

# ---- Heatmap ----
ggplot(mean_ra_genus,
       aes(x = sample_type,
           y = fct_reorder(Genus, mean_ra, .fun = max),
           fill = mean_ra)) +
        geom_tile(color = "white") +
        scale_fill_viridis_c(
                trans = "log10",
                labels = scientific,
                name = "Mean relative abundance"
        ) +
        theme_bw() +
        labs(
                title = "Taxonomic structure of bacterial communities",
                subtitle = "Top genera by mean relative abundance",
                x = "",
                y = "Genus"
        ) +
        theme(
                panel.grid = element_blank(),
                axis.text.y = element_text(size = 9)
        )

# ---- Taxonomic heat trees with metacoder ----
View(tax_table(phyloseq_prev2_ab0))


################## I tried metacoder for taxonomic tree ########################
# ---- Heat tree: bagged vs open flowers (metacoder) ----
library(metacoder)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 1) Prepare OTU table (samples x ASVs)
# ------------------------------------------------------------
otu_mat <- as(otu_table(phyloseq_prev2_ab0), "matrix")
if (taxa_are_rows(phyloseq_prev2_ab0)) {
        otu_mat <- t(otu_mat)
}

otu_df <- as.data.frame(otu_mat) |>
        rownames_to_column("sample_id")

# ------------------------------------------------------------
# 2) Minimal metadata
# ------------------------------------------------------------
meta_df <- data.frame(sample_data(phyloseq_prev2_ab0)) |>
        rownames_to_column("sample_id") |>
        select(sample_id, sample_type)

# ------------------------------------------------------------
# 3) Taxonomy table (Genus level only)
# ------------------------------------------------------------
tax_df_genus <- tax_df_clean |>
        rownames_to_column("taxon_id") |>
        select(taxon_id, Domain, Phylum, Class, Order, Family, Genus) |>
        mutate(
                Genus = ifelse(is.na(Genus) | Genus == "", "unclassified", Genus)
        )

# ------------------------------------------------------------
# 4) Mean relative abundance per Genus per group
# ------------------------------------------------------------
genus_ra <- otu_df |>
        left_join(meta_df, by = "sample_id") |>
        pivot_longer(
                cols = -c(sample_id, sample_type),
                names_to = "taxon_id",
                values_to = "count"
        ) |>
        left_join(tax_df_genus, by = "taxon_id") |>
        group_by(sample_type, Genus) |>
        summarise(mean_count = mean(count), .groups = "drop") |>
        group_by(sample_type) |>
        mutate(mean_ra = mean_count / sum(mean_count)) |>
        ungroup() |>
        select(sample_type, Genus, mean_ra) |>
        pivot_wider(
                names_from = sample_type,
                values_from = mean_ra,
                values_fill = 0
        )

# ------------------------------------------------------------
# 5) Build taxmap object (Genus hierarchy)
# ------------------------------------------------------------
taxmap_obj <- parse_tax_data(
        tax_df_genus,
        class_cols = c(
                Domain = "Domain",
                Phylum = "Phylum",
                Class  = "Class",
                Order  = "Order",
                Family = "Family",
                Genus  = "Genus"
        )
)

# ------------------------------------------------------------
# 6) Attach genus-level abundance data
# ------------------------------------------------------------
taxmap_obj$data$genus_ra <- genus_ra

# ------------------------------------------------------------
# 7) Heat tree: bagged flowers
# ------------------------------------------------------------
heat_tree(
        taxmap_obj,
        data = "genus_ra",
        node_label = taxon_names,
        node_size = unbagged_flower + bagged_flower,
        node_color = bagged_flower,
        node_color_range = c("grey90", "blue"),
        node_color_trans = "log10",
        node_size_range = c(0.01, 0.06),
        layout = "davidson-harel",
        title = "Genus-level composition: bagged flowers"
)

# ------------------------------------------------------------
# 8) Heat tree: openly pollinated flowers
# ------------------------------------------------------------
heat_tree(
        taxmap_obj,
        data = "genus_ra",
        node_label = taxon_names,
        node_size = unbagged_flower + bagged_flower,
        node_color = unbagged_flower,
        node_color_range = c("grey90", "red"),
        node_color_trans = "log10",
        node_size_range = c(0.01, 0.06),
        layout = "davidson-harel",
        title = "Genus-level composition: openly pollinated flowers"
)

install.packages("pheatmap")
library(pheatmap)

# Top genera by mean abundance
genus_mat <- genus_ra |>
        column_to_rownames("Genus") |>
        select(bagged_flower, unbagged_flower)

# Log-transform for visualization
genus_mat_log <- log10(genus_mat + 1e-6)

# Select top 40 genera by total abundance
top_genera <- rowSums(genus_mat) |>
        sort(decreasing = TRUE) |>
        head(50) |>
        names()

pheatmap(
        genus_mat_log[top_genera, ],
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        main = "Genus-level composition: bagged vs open flowers",
        fontsize_row = 8
)

# ---- Visitation footrpint ----
genus_fc <- genus_ra |>
        mutate(
                log2_fc_open_vs_bagged =
                        log2((unbagged_flower + 1e-6) /
                                     (bagged_flower + 1e-6))
        ) |>
        arrange(desc(log2_fc_open_vs_bagged))
View(genus_fc)

#

heat_mat_fc <- fc_df |>
        select(Genus, log2FC) |>
        column_to_rownames("Genus") |>
        as.matrix()

pheatmap(
        heat_mat_fc,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-max(abs(heat_mat_fc)),
                     max(abs(heat_mat_fc)),
                     length.out = 101),
        main = "Differentially enriched genera (log2FC)",
        fontsize_row = 7,
        border_color = NA
)

#
top_fc <- fc_df |>
        slice_max(order_by = abs(log2FC), n = 80)

heat_mat <- top_fc |>
        select(Genus, bagged_flower, unbagged_flower) |>
        column_to_rownames("Genus") |>
        mutate(across(everything(), ~ log10(.x + 1e-6))) |>
        as.matrix()

pheatmap(
        heat_mat,
        scale = "none",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("grey90", "orange", "red"))(100),
        main = "Top 80 genera by |log2FC|",
        fontsize_row = 8
)


# ============================================================
# Heatmap across samples: top50 enriched in OPEN + top50 in BAGGED
# (Genus level; columns = samples)
# ============================================================

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)

# Heatmap packages
library(pheatmap)

# Optional (for bold row labels): ComplexHeatmap
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(grid)

# ----------------------------
# 0) Choose phyloseq + groups
# ----------------------------
ps <- phyloseq_prev2_ab0

# Make sure sample_type exists and is exactly the two groups you want
# (adjust names if needed)
meta <- data.frame(sample_data(ps)) |>
        rownames_to_column("sample_id") |>
        mutate(sample_type = as.character(sample_type))

stopifnot(all(c("bagged_flower", "unbagged_flower") %in% unique(meta$sample_type)))

# -----------------------------------------
# 1) Agglomerate to Genus + Relative Abundance
# -----------------------------------------
# If your taxonomy rank is named differently, adjust "Genus"
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)

ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

# Extract genus-by-sample matrix (rows = taxa, cols = samples)
mat <- as(otu_table(ps_genus_rel), "matrix")
if (!taxa_are_rows(ps_genus_rel)) mat <- t(mat)

# Get genus names for rows
tax_genus <- as.data.frame(tax_table(ps_genus_rel)) |>
        rownames_to_column("taxa_id")

# Replace missing genus with a readable label
# (keeps unique taxa_id to avoid duplicated rownames)
tax_genus <- tax_genus |>
        mutate(
                Genus = as.character(Genus),
                Genus = ifelse(is.na(Genus) | Genus == "" | Genus == "unclassified", "unclassified", Genus),
                Genus_label = paste0(Genus, " [", taxa_id, "]")
        )

# Align labels to matrix rows
tax_genus <- tax_genus |> filter(taxa_id %in% rownames(mat))
mat <- mat[tax_genus$taxa_id, , drop = FALSE]
rownames(mat) <- tax_genus$Genus_label

View(tax_genus)

# -----------------------------------------
# 2) Compute mean RA per genus in each group + log2FC
# -----------------------------------------
sample_type_vec <- meta$sample_type
names(sample_type_vec) <- meta$sample_id
sample_type_vec <- sample_type_vec[colnames(mat)]

open_samps   <- names(sample_type_vec)[sample_type_vec == "unbagged_flower"]
bagged_samps <- names(sample_type_vec)[sample_type_vec == "bagged_flower"]

eps <- 1e-6

mean_open   <- rowMeans(mat[, open_samps, drop = FALSE], na.rm = TRUE)
mean_bagged <- rowMeans(mat[, bagged_samps, drop = FALSE], na.rm = TRUE)

log2fc_open_vs_bagged <- log2((mean_open + eps) / (mean_bagged + eps))

fc_df <- tibble(
        Genus_label = names(log2fc_open_vs_bagged),
        mean_open = mean_open,
        mean_bagged = mean_bagged,
        log2FC = log2fc_open_vs_bagged
)

# -----------------------------------------
# 3) Pick top50 enriched in OPEN and top50 enriched in BAGGED
# -----------------------------------------
top_open <- fc_df |> arrange(desc(log2FC)) |> slice_head(n = 50)
top_bag  <- fc_df |> arrange(log2FC)       |> slice_head(n = 50)

selected <- bind_rows(top_open, top_bag) |> distinct(Genus_label) |> pull(Genus_label)

mat_sel <- mat[selected, , drop = FALSE]

# Optional: filter out extremely rare rows to improve readability
# e.g. keep only genera that reach at least 0.1% in at least one sample
# (tune threshold if you want)
keep <- apply(mat_sel, 1, function(x) max(x, na.rm = TRUE) >= 0.001)
mat_sel <- mat_sel[keep, , drop = FALSE]

# -----------------------------------------
# 4) Transform for plotting (recommended)
# -----------------------------------------
# With microbiome RA, log transform usually helps a lot:
mat_plot <- log10(mat_sel + eps)

# -----------------------------------------
# 5) Column annotation (sample_type + optional farm/tree/etc.)
# -----------------------------------------
ann_col <- meta |>
        filter(sample_id %in% colnames(mat_plot)) |>
        distinct(sample_id, .keep_all = TRUE) |>
        select(sample_id, sample_type, farm_id, tree_id, management_type, is_pollination_clsm, pi_clsm) |>
        column_to_rownames("sample_id")

# Order columns: bagged first then open (optional)
ord_cols <- c(bagged_samps, open_samps)
ord_cols <- ord_cols[ord_cols %in% colnames(mat_plot)]
mat_plot <- mat_plot[, ord_cols, drop = FALSE]
ann_col  <- ann_col[ord_cols, , drop = FALSE]

# -----------------------------------------
# 6A) Heatmap with pheatmap (fast, simple)
# -----------------------------------------
# IMPORTANT: Do NOT use scale="row" if you want absolute differences; 
# use scale="row" if you want pattern-based comparison across samples.
pheatmap(
        mat_plot,
        annotation_col = ann_col["sample_type", drop = FALSE],
        show_colnames = FALSE,
        fontsize_row = 6,
        border_color = NA,
        clustering_method = "ward.D2",
        main = "Top genera enriched in OPEN (top50) and BAGGED (top50) across samples\nValues: log10(RA + 1e-6)"
)

# -----------------------------------------
# 6B) Heatmap with ComplexHeatmap + bold insect-associated genera (optional)
# -----------------------------------------
# Put here the genera you consider insect-associated (you can expand later)
# NOTE: rownames are "Genus [taxa_id]"; so we match by partial string.
insect_genera <- c(
        "Wolbachia", "Arsenophonus", "Rickettsia", "Rickettsiella",
        "Candidatus_Tremblaya", "Candidatus_Carsonella", "Gilliamella"
)

is_insect <- sapply(rownames(mat_plot), function(x) any(grepl(paste0("\\b", insect_genera, "\\b"), x)))

row_lab_gp <- gpar(fontsize = 7, fontface = ifelse(is_insect, "bold", "plain"))

ha <- HeatmapAnnotation(
        sample_type = ann_col$sample_type,
        annotation_name_side = "left"
)

col_fun <- colorRamp2(
        c(quantile(as.vector(mat_plot), 0.02, na.rm = TRUE),
          median(as.vector(mat_plot), na.rm = TRUE),
          quantile(as.vector(mat_plot), 0.98, na.rm = TRUE)),
        c("blue", "white", "red")
)

Heatmap(
        mat_plot,
        name = "log10(RA)",
        top_annotation = ha,
        col = col_fun,
        show_column_names = FALSE,
        row_names_gp = row_lab_gp,
        column_title = "Top50 OPEN-enriched + Top50 BAGGED-enriched genera (across samples)",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2"
)


# ---- Enodsymbioints check ----
# Find ASVs assigned to genus "endosymbionts"
endo_asvs <- taxa_names(phyloseq_prev2_ab0)[
        tax_table(phyloseq_prev2_ab0)[, "Genus"] %in% c("endosymbionts", "g__endosymbionts")
]
length(endo_asvs)


# Relative abundance phyloseq
physeq_rel <- transform_sample_counts(
        phyloseq_prev2_ab0,
        function(x) x / sum(x)
)

# Sum RA of endosymbionts per sample
endo_ra <- otu_table(physeq_rel)[endo_asvs, , drop = FALSE] |>
        as.matrix() |>
        colSums()

endo_df <- data.frame(
        sample_id = names(endo_ra),
        endosymbionts_ra = as.numeric(endo_ra)
)

endo_df <- endo_df |>
        left_join(
                data.frame(sample_data(phyloseq_prev2_ab0)) |>
                        rownames_to_column("sample_id"),
                by = "sample_id"
        ) |>
        arrange(desc(endosymbionts_ra))
View(endo_df)
