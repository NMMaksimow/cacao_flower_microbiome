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
library(phyloseq) # integrate analysis of microbiome data (OTU/ASV, taxonomy, metadata)
library(decontam) # identification and removal of contaminant taxa
library(vegan) # ecological analyses (diversity, ordination, PERMANOVA)

# differential abundance & indicator species analyses
library(ANCOMBC)  # differential abundance testing accounting for compositionality
library(indicspecies) # indicator species analysis (specificity and fidelity of association of taxa with groups)
install.packages("gsl", type = "source")

# Linear mixed-effects models
library(lme4) # linear and generalized linear mixed-effects models
library(lmerTest) # p-values and tests for linear mixed-effects models
library(ggeffects) # marginal effects and predictions from regression models
library(glmmTMB) # flexible generalised linear mixed-effect models
library(DHARMa) # residual diagnostics for mixed and generalized models

install.packages("gsl", type = "source")
sessionInfo()

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
library(metacoder) # visualization and manipulation of hierarchical taxonomic data

# renv::snreadxlrenv::snapshot()

################################################################################
################### Import data and phyloseq object building ###################
################################################################################

# ---- Creating phyloseq object for ITS1 microbiome data ----
# 1) Metadata: Import metadata from .xlsx
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

# Build phyloseq metadata table (sample_data-class) based on imported metadata dataframe
sample_metadata_ITS1 <- phyloseq::sample_data(metadata)

# 2) ASV feature table with counts (QIIME2 DADA2 output)
# Build otu_table-class phyloseq object from biom feature table converted to matrix
asv_table_ITS1 <- read_biom(
        file.path("qiime2", "export", "CFM_ITS1_dada2_table", "feature-table.biom")
        ) |> 
        biom_data() |> 
        as.matrix() |> 
        otu_table(taxa_are_rows = TRUE)

# 3) Taxonomy table (QIIME2 UNITE Naive-Bayes classifier output)
# Build taxonomyTable-class phyloseq object from tsv taxonomy table exported from QIIME2
taxonomy_table_ITS1 <- read_tsv(
        file.path("qiime2", "export", "CFM_ITS1_taxonomy", "taxonomy.tsv"), 
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
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "UNITE_SH"),
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
phyloseq_ITS1_raw <- phyloseq(
        asv_table_ITS1,
        taxonomy_table_ITS1,
        sample_metadata_ITS1
        )

# 5) Representative sequences (FASTA) file path
repseqs_fp <- file.path("qiime2", "export", "CFM_ITS1_dada2_repseqs", "dna-sequences.fasta")

# ---- Some diagnostic commands for ITS1 phyloseq object ----
# Overview of phyloseq object
phyloseq_ITS1_raw

# Exploring metadata
View(sample_data(phyloseq_ITS1_raw))
colnames(sample_data((phyloseq_ITS1_raw)))
rownames(sample_data(phyloseq_ITS1_raw))
dim.data.frame(metadata)
sample_variables(phyloseq_ITS1_raw)
nsamples(phyloseq_ITS1_raw)
table(sample_data(phyloseq_ITS1_raw)$sample_type)
summary(sample_data(phyloseq_ITS1_raw)$pi_clsm)

# Exploring feature table with counts
# NB! Heavy sparce matrix with features in raws and sample_id in columns
View(otu_table(phyloseq_ITS1_raw))

# sanity check wether sample names identical in feature table with counts and metadata
identical(colnames(otu_table(phyloseq_ITS1_raw)), rownames(sample_data(phyloseq_ITS1_raw)))

# read count summary statistics
sample_sums(phyloseq_ITS1_raw) |> summary() 
sample_sums(phyloseq_ITS1_raw) |> hist() 

# Exploring taxonomy annotation
View(tax_table(phyloseq_ITS1_raw))
ntaxa(phyloseq_ITS1_raw)
table(tax_table(phyloseq_ITS1_raw)[, "Kingdom"], useNA = "ifany")
# summary statistics on counts by taxa
taxa_sums(phyloseq_ITS1_raw) |> summary() 
taxa_sums(phyloseq_ITS1_raw) |> hist() # highly skewed (most taxa are singletons and low abundance)

# Check Taxon column lengths in row TSV annotation file
taxonomy_raw <- read_tsv(
        file.path("qiime2", "export", "CFM_ITS1_taxonomy", "taxonomy.tsv"),
        show_col_types = FALSE
)
tax_split_lengths <- lengths(strsplit(taxonomy_raw$Taxon, ";"))
table(tax_split_lengths)
taxonomy_raw$Taxon[tax_split_lengths > 7] |> head(10)




# ---- Save intermidiate results ---
# saveRDS(phyloseq_ITS1_raw, "results/phyloseq_ITS1_raw.rds")
# phyloseq_ITS1_raw <- readRDS("results/phyloseq_ITS1_raw.rds")



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

# 1) Identify contamination with decontam(): 27 ASVs identified as contaminant
decontam_prev05_ITS1 <- phyloseq_ITS1_raw |>
        subset_samples(sample_type != "mock_community") |>
        isContaminant(
                method = "prevalence",
                neg = "is_negative_control",
                # Posterior probability of contaminant identification
                threshold = 0.5
                )

# 2) Remove features identified as contamination from the phyloseq object
phyloseq_ITS1_decontam <- prune_taxa(
        !decontam_prev05_ITS1$contaminant,
        phyloseq_ITS1_raw
        )

# ---- Diagnostic commands for decontaminated phyloseq object ----
# Diagnostic of identified contamination
# 27 ASVs were identified as contamination
table(decontam_prev_ITS1$contaminant) 

# Generate taxonomy table for contaminant ASVs identified by decontam()
contaminant_ITS1_taxa_tbl <- tax_table(phyloseq_ITS1_raw)[decontam_prev_ITS1$contaminant, ] |>
        as.data.frame() |>
        rownames_to_column("FeatureID") |>
        # To add total count for each taxa
        left_join(
                taxa_sums(phyloseq_ITS1_raw) |> 
                        enframe(name = "FeatureID", value = "Total_abundance"),
                by = "FeatureID"
        ) |>
        arrange(desc(Total_abundance))

View(contaminant_ITS1_taxa_tbl)
# s__Cutaneotrichosporon_curvatum => human mycobiome
# s__Apiotrichum_domesticum => human mycobiome
# s__Hannaella_taiwanensis ?

# Overview of taxa identified as contamination
contaminant_ITS1_taxa_tbl |>
        count(Kingdom, Phylum, Genus, Species, sort = TRUE)

# How many taxa was removed by decontam()?
ntaxa(phyloseq_ITS1_raw) # 12130
ntaxa(phyloseq_ITS1_decontam) # 12103

# Summary statistics on total read count by sample before and after decontam()?
summary(sample_sums(phyloseq_ITS1_raw))
summary(sample_sums(phyloseq_ITS1_decontam))

# How many samples have resulted in read count below 1000 and 5000?
sum(sample_sums(phyloseq_ITS1_raw) < 1000)
sum(sample_sums(phyloseq_ITS1_raw) < 5000)

sum(sample_sums(phyloseq_ITS1_decontam) < 1000)
sum(sample_sums(phyloseq_ITS1_decontam) < 5000)

# What samples have low read count?
sample_sums(phyloseq_ITS1_raw) |>
        (\(x) names(x[x < 1000]))()

# All 28 negative controls (negative_pcr_control, extraction_control have less than 1000 read count)
sample_sums(phyloseq_ITS1_decontam) |>
        (\(x) names(x[x < 1000]))()

# In addition 2 flower samples (cmt72, cmt73 have less than 5000 read count)
sample_sums(phyloseq_ITS1_decontam) |>
        (\(x) names(x[x < 5000]))()

# ---- Taxonomy-based filtration: leave k__Fungi ----
# I left only target kingdom: k__Fungi:
phyloseq_ITS1_fungi <- phyloseq_ITS1_decontam |> 
        subset_taxa(Kingdom == "k__Fungi")

# ---- Additional taxonomically filtered subsets of potential exploratory interest ----
# Additionally I created the subset with Unassigned sequences:
phyloseq_ITS1_unassigned <- phyloseq_ITS1_decontam |> 
        subset_taxa(Kingdom == "Unassigned")

# ---- Diagnostic commands after cleaning off-target sequences ----

# Before cleanin ITS1 sequences represent 6 eukaryotic kingdoms and Eukaryotes incertae sedis as well as 122 Unassigned sequences.
table(tax_table(phyloseq_ITS1_decontam)[, "Kingdom"], useNA = "ifany")

# After cleaning the only kingdom is Fungi
table(tax_table(phyloseq_ITS1_fungi)[, "Kingdom"], useNA = "ifany")

table(tax_table(phyloseq_ITS1_unassigned)[, "Kingdom"], useNA = "ifany")

ntaxa(phyloseq_ITS1_raw)
ntaxa(phyloseq_ITS1_decontam)
ntaxa(phyloseq_ITS1_fungi)
ntaxa(phyloseq_ITS1_unassigned)



################################################################################
### Quality Control: decontamination, off-target sequences & mock communities ##
################################################################################


# ---- Visualizing mock community using taxa barplots ----
phyloseq_ITS1_fungi |>
        subset_samples(sample_type == "mock_community") |>
        (\(ps) prune_taxa(taxa_sums(ps) > 0, ps))() |>
        tax_glom(taxrank = "Species", NArm = FALSE) |>
        transform_sample_counts(\(x) x / sum(x)) |>
        plot_bar(fill = "Species") +
        theme_bw() +
        theme(
                axis.text.x = element_text(angle = 90, hjust = 1),
                legend.position = "right"
        )

# Taxonomic barplot visualization using phyloseq::plot_bar()
# NB! Normalization before taxa agglomeration => which result in 100% barplots
phyloseq_ITS1_fungi |> 
        subset_samples(pi_clsm >= 50 ) |>
        (\(ps) {
                # Transform counts to relative abundances (RA)
                ps_rel <- transform_sample_counts(ps, \(x) x / sum(x))
                
                mean_ra <- taxa_sums(ps_rel) / nsamples(ps_rel)
                
                # Plot only taxa with RA >= mean for this taxa across the whole sample subset
                keep <- mean_ra >= 0.000001
                prune_taxa(keep, ps_rel)
        })() |>
        # Agglomeration of ASVs and counts by taxa at taxonomic rank == taxrank
        tax_glom(taxrank = "Phylum", NArm = FALSE) |>
        plot_bar(fill = "Phylum") +
        theme_bw() +
        theme(
                axis.text.x = element_text(angle = 90, hjust = 1)
        ) +
        labs(title = "XXX", subtitle = "YYY", x = NULL, y = "Relative abundance")

# Taxonomic barplot visualization using ggplot()
phyloseq_ITS1_raw |>
        # Subset samples using metadata variables
        subset_samples(sample_type %in% c("mock_community", "extraction_control", "negative_pcr_control")) |>
        # Agglomerate taxa at the defined taxonomic rank
        tax_glom(taxrank = "Species", NArm = FALSE) |>
        transform_sample_counts(\(x) x / sum(x)) |>
        # Transform phyloseq object to long data.frame for gggplot
        psmelt() |>
        # All taxa with RA < RA_threshold per sample are aggregated to Other
        mutate(
                Species = ifelse(Abundance < 0.01, "Other", Species)
        ) |>
        group_by(Sample, Species) |>
        summarise(Abundance = sum(Abundance), .groups = "drop") |>
        ggplot(aes(x = Sample, y = Abundance, fill = Species)) +
        geom_col() +
        scale_y_continuous(labels = scales::percent_format()) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = "XXX", subtitle = "YYY", x = NULL, y = "Relative abundance")

# ---- Visualize expected Vs observed mock communities ----
# ZymoBiomics D6300 DNA mock community include two fungi:
# a) Saccharomyces cerevisiae - 2%, and 
# b) Cryptococcus neoformans - 2%
expected_mock_ITS1 <- data.frame(
        Family = c(
                "g__Cryptococcus",
                "g__Saccharomyces"
        ),
        expected_RA = c(rep(50, 2))
)


# Expected mock composition was rescaled to fungal taxa only (50% each), as fungil DNA is not expected to amplify with ITS1 primers.
expected_mock_zymoD6300_ITS1 <- tibble::tibble(
        Genus = c("g__Cryptococcus", "g__Saccharomyces"),
        Abundance = c(0.5, 0.5),
        Type = "Expected",
        Sample = "Expected mock"
        )

observed_mock <- phyloseq_ITS1_fungi |>
        subset_samples(sample_type == "mock_community") |>
        tax_glom(taxrank = "Genus", NArm = FALSE) |>
        transform_sample_counts(function(x) x / sum(x)) |>
        psmelt() |>
        filter(Abundance > 0) |> 
        mutate(
                Genus = ifelse(is.na(Genus), "Unassigned", Genus),
                Type = "Observed",
                Sample = "Mock samples"
        )

plot_data <- bind_rows(
        observed_mock,
        expected_mock_zymoD6300_ITS1
        )

ggplot(plot_data,
       aes(x = Sample, y = Abundance, fill = Genus)) +
        geom_col(width = 0.7, color = "white") +
        facet_wrap(~Type, scales = "free_x") +
        scale_y_continuous(labels = scales::percent_format()) +
        labs(
                title = "Mock community composition (ITS1)",
                y = "Relative abundance",
                x = ""
        ) +
        theme(
                axis.text.x = element_text(angle = 90, hjust = 1),
                legend.position = "right"
        )




################################################################################
################ Biological samples and Prevalence filtration ##################
################################################################################

# ---- Filter only biological samples (remove all controls) ----
phyloseq_ITS1_fungi_biosamples <- subset_samples(
        phyloseq_ITS1_fungi,
        !sample_type %in% c(
                "extraction_control",
                "negative_pcr_control",
                "mock_community"
        )
)

# ---- ASV prevalence diagnostics ----
# Empirical Cumulative Distribution Function (ECDF)
phyloseq_ITS1_fungi_biosamples |>
        (\(physeq) {
                
                # extract OTU table and ensure ASVs are rows
                otu_table(physeq) |>
                        (\(otu) {
                                if (!taxa_are_rows(physeq)) t(otu) else otu
                        })() |>
                        # convert to presence/absence
                        as("matrix") |>
                        (\(x) x > 0)() |>
                        # calculate prevalence per ASV
                        rowMeans() |>
                        tibble::enframe(
                                name  = "ASV_ID",
                                value = "prevalence"
                        )
                
        })() |>
        ggplot(aes(x = prevalence)) +
        stat_ecdf(size = 1) +
        scale_x_continuous(
                labels = scales::percent_format(accuracy = 1)
        ) +
        geom_vline(
                xintercept = c(0.01, 0.05, 0.10),
                linetype = "dashed",
                colour = "grey40"
        ) +
        labs(
                x = "ASV prevalence across biological samples",
                y = "Cumulative proportion of ASVs",
                title = "Empirical cumulative distribution of ASV prevalence"
        ) +
        theme_minimal()

c(0, 0.01, 0.05, 0.10, 0.5) |>
        purrr::map_dfr(\(prevalence_threshold) {
                
                phyloseq_ITS1_fungi_biosamples |>
                        (\(physeq) {
                                
                                otu_table(physeq) |>
                                        (\(otu) {
                                                if (!taxa_are_rows(physeq)) t(otu) else otu
                                        })() |>
                                        as("matrix") |>
                                        (\(x) x > 0)() |>
                                        rowMeans()
                                
                        })() |>
                        (\(prevalence_vector) {
                                tibble::tibble(
                                        prevalence_threshold = prevalence_threshold,
                                        proportion_ASVs_below_threshold =
                                                mean(prevalence_vector < prevalence_threshold),
                                        n_ASVs_below_threshold =
                                                sum(prevalence_vector < prevalence_threshold),
                                        total_ASVs =
                                                length(prevalence_vector)
                                )
                        })()
                
        })

phyloseq_ITS1_fungi_biosamples |>
        (\(physeq) {
                otu_table(physeq) |>
                        (\(otu) {
                                if (!taxa_are_rows(physeq)) t(otu) else otu
                        })() |>
                        as("matrix") |>
                        (\(x) x > 0)() |>
                        rowMeans()
        })() |>
        sort() |>
        (\(prevalence_vector) {
                tibble::tibble(
                        prevalence = prevalence_vector,
                        ecdf = seq_along(prevalence_vector) / length(prevalence_vector)
                )
        })() |>
        (\(df) {
                df |>
                        dplyr::mutate(
                                d1 = c(NA, diff(ecdf)),
                                d2 = c(NA, diff(d1))
                        )
        })() |>
        dplyr::slice_max(abs(d2), n = 1) |> View()

# ---- Prevalence filtration ----
phyloseq_ITS1_fungi_biosample_prev_ge_1percent <-
        phyloseq_ITS1_fungi_biosamples |>
        (\(physeq) {
                
                # number of biological samples
                n_samples <- nsamples(physeq)
                
                # minimum number of samples corresponding to 1% prevalence
                min_samples_required <- ceiling(0.01 * n_samples)
                
                prune_taxa(
                        apply(
                                as(
                                        if (taxa_are_rows(physeq))
                                                otu_table(physeq)
                                        else
                                                t(otu_table(physeq)),
                                        "matrix"
                                ) > 0,
                                1,
                                sum
                        ) >= min_samples_required,
                        physeq
                )
                
        })()

ntaxa(phyloseq_ITS1_fungi_biosample_prev_ge_1percent)
ntaxa(phyloseq_ITS1_fungi_biosamples_prev3)




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


# ---- Effect of decontamination on control samples ----
# Subset of control samples before decontamination transformed to RA.
# There're 42 samples and 265 taxa in this data subset.
phyloseq_ITS1_raw |> 
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
phyloseq_ITS1_decontam |> 
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
                make_control_long_table(phyloseq_ITS1_raw, "Before decontam"),
                make_control_long_table(phyloseq_ITS1_decontam, "After decontam")
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


# ---- Sequencing retention after decontamination and removal off-targe sequences ----
# Sequences retention during decontamination and taxonomic off-target sequence cleaning
read_retention_loss_table_ITS1 <-
        dplyr::bind_rows(
                phyloseq_ITS1_raw |>
                        sample_sums() |>
                        tibble::enframe(name = "sample_id",
                                        value = "n_reads") |>
                        dplyr::mutate(stage = "Raw"),
                
                phyloseq_ITS1_decontam |>
                        sample_sums() |>
                        tibble::enframe(name = "sample_id",
                                        value = "n_reads") |>
                        dplyr::mutate(stage = "After decontam"),
                
                phyloseq_ITS1_fungi |>
                        sample_sums() |>
                        tibble::enframe(name = "sample_id",
                                        value = "n_reads") |>
                        dplyr::mutate(stage = "Fungi only")
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
                fungi_only_reads = n_reads[stage == "Fungi only"],
                
                fraction_removed_decontam =
                        (raw_reads - after_decontam_reads) / raw_reads,
                
                fraction_removed_off_target =
                        (after_decontam_reads - fungi_only_reads) / raw_reads,
                
                fraction_retained =
                        fungi_only_reads / raw_reads
        ) |> 
        dplyr::ungroup() |>
        dplyr::left_join(
                phyloseq_ITS1_raw |>
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
                                "Retained (Fungi only)"
                        )
                )
        )

# Helper function to plot decontamination and off-target filtering results
plot_read_retention_stacked_ITS1 <-
        function(grouping_variable,
                 grouping_levels = NULL,
                 filter_function = NULL) {
                
                read_retention_loss_table_ITS1 |>
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
plot_read_retention_stacked_ITS1(grouping_variable = "sublibrary_id")

# Plot decontamination and off-target filtering results by sample type
plot_read_retention_stacked_ITS1(
        grouping_variable = "sample_type",
        grouping_levels = c("extraction_control",
                            "negative_pcr_control",
                            "mock_community",
                            "bagged_flower",
                            "unbagged_flower"
        )
)

# Plot decontamination and off-target filtering results by management_type (only biological samples)
plot_read_retention_stacked_ITS1(
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
plot_read_retention_stacked_ITS1(
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
plot_read_retention_stacked_ITS1(
        grouping_variable = "is_pollination_clsm",
        grouping_levels = c(FALSE, TRUE),
        filter_function =
                \(x) dplyr::filter(
                        x,
                        !is.na(is_pollination_clsm)
                )
)

# Plot decontamination and off-target filtering results by is_germination_clsm (only phenotyped flowers)
plot_read_retention_stacked_ITS1(
        grouping_variable = "is_germination",
        grouping_levels = c(FALSE, TRUE),
        filter_function =
                \(x) dplyr::filter(
                        x,
                        !is.na(is_germination)
                )
)

# Does sequencing depth confound pollination intensity? No
read_retention_loss_table_ITS1 |>
        dplyr::filter(
                !is.na(pi_clsm),
                pi_clsm > 0,
                component == "Retained (Fungi only)"
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
read_retention_loss_table_ITS1 |>
        dplyr::filter(
                !is.na(callose_plug_index),
                callose_plug_index > 0,
                component == "Retained (Fungi only)"
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
        phyloseq_ITS1_fungi |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
        ) |> 
        dplyr::left_join(
                phyloseq_ITS1_raw |>
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
                                                      "unbagged_flower")
                                           )
                      ) |> 
        ggplot(aes(x = sample_type, y = n_reads)) +
        geom_boxplot() +
        theme_bw()

# sublibrary_id
dplyr::bind_rows(
        phyloseq_ITS1_fungi |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
        ) |> 
        dplyr::left_join(
                phyloseq_ITS1_raw |>
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
        phyloseq_ITS1_fungi |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
        ) |> 
        dplyr::left_join(
                phyloseq_ITS1_raw |>
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
        phyloseq_ITS1_fungi |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
        ) |> 
        dplyr::left_join(
                phyloseq_ITS1_raw |>
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
        phyloseq_ITS1_fungi |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |> 
        dplyr::left_join(
                phyloseq_ITS1_raw |>
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
        phyloseq_ITS1_fungi |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |> 
        dplyr::left_join(
                phyloseq_ITS1_raw |>
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
        phyloseq_ITS1_fungi |>
                sample_sums() |>
                tibble::enframe(name = "sample_id",
                                value = "n_reads")
) |> 
        dplyr::left_join(
                phyloseq_ITS1_raw |>
                        sample_data() |>
                        as.data.frame() |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        rownames_to_column(var = "sample_id"),
                by = "sample_id"
        ) |> 
        ggplot(aes(x = dna_conc_hs_its1, y = n_reads, color = sample_type)) +
        geom_point() +
        theme_bw()



# ---- Sequencing depth summary ----
# All samples, including all type of controls
depths_fungi_all_samples <- sample_sums(phyloseq_ITS1_fungi)
summary(depths_fungi_all_samples)
quantile(depths_fungi_all_samples)
quantile(depths_fungi_all_samples, probs = c(0.05, 0.1))

# Only biological samples here:
depths_fungi_bio <- sample_sums(phyloseq_ITS1_fungi_biosamples)
summary(depths_fungi_bio)
quantile(depths_fungi_bio)
quantile(depths_fungi_bio, probs = c(0.05, 0.1))
nsamples(phyloseq_ITS1_fungi_biosamples)*0.1

# Depth distribution
hist(depths_fungi_all_samples,
     breaks = 50,
     main = "Sequencing depth per sample after deontamination, filtation off-target sequences\n(all samples)",
     xlab = "Read count"
)

hist(depths_fungi_bio,
     breaks = 50,
     main = "Sequencing depth per sample after deontamination, filtation off-target sequences\n(only biological samples)",
     xlab = "Read count"
)

# ---- Diagnostic rarefaction curves ----
# Create a matrix for building rarefaction curves and transpose it because vegan::rarecurve() requires samples in rows
asv_mat_ITS1 <-  phyloseq_ITS1_fungi_biosamples |> 
        otu_table() |> 
        #taxa_are_rows() |> 
        as("matrix") |> 
        t()

# Minimal seq depth could be used for 
min_depth <- min(rowSums(asv_mat_ITS1)) 

# Color coding using categorical variables from metadata
color_var <- "sample_type"   
metadata_for_rarecurve <- as.data.frame(sample_data(phyloseq_ITS1_fungi_biosamples))
metadata_for_rarecurve <- metadata_for_rarecurve[rownames(asv_mat_ITS1), , drop = FALSE]

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
        asv_mat_ITS1,
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

# ---- Diagnostic QC comparison of alpha metrics between flower samples and mock community ----
# Here I tried to rarefy and calculate alpha diversity for all samples, incl. all controls. 
# All negative_pcr_control and extraction_control samples didn't pass rarefaction at 10 000 reads. 
# All mock_community, with ASVs primarily from Cryptococcus, samples as expected had low alpha diversity metrics compared to flower communities.
phyloseq_ITS1_fungi |> 
        rarefy_even_depth(
                sample.size = 10000,
                rngseed = 42,
                replace = FALSE,
                trimOTUs = TRUE,
                verbose = TRUE) |> 
        (\(ps) {
        estimate_richness(
                ps,
                measures = c("Observed", "Shannon", "InvSimpson")
        ) |> 
        rownames_to_column("sample_id") |> 
        left_join(
                data.frame(
                        sample_data(phyloseq_ITS1_fungi)) |>
                        rownames_to_column("sample_id"),
                by = "sample_id"
        )
                }
        )() |> 
        dplyr::select(
                sample_id,
                Observed,
                Shannon,
                InvSimpson,
                sample_type
                ) |>
        mutate(
                sample_type_groupped = case_when(
                        sample_type == "mock_community" ~ "mock_community",
                        sample_type %in% c("bagged_flower", "unbagged_flower") ~ "flower_sample",
                        TRUE ~ NA_character_
                ),
                sample_type_groupped = factor(sample_type_groupped, levels = c("mock_community", "flower_sample"))
        ) |> 
        pivot_longer(
                cols = c(Observed, Shannon, InvSimpson),
                names_to = "metric",
                values_to = "value"
        ) |>
        mutate(
                metric = factor(metric,
                                levels = c("Observed",
                                           "Shannon",
                                           "InvSimpson")
                )
        ) |> 
        ggplot(aes(x = sample_type_groupped, y = value)) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
        facet_wrap(~ metric, scales = "free_y") +
        theme_bw() +
        labs(
                title = "Alpha diversity of flower fungal communities and positive control mock communities",
                subtitle = "Alpha diversity metrics were calculated on an ASV table rarefied to 10,000 reads per sample",
                x = "",
                y = ""
        ) +
        theme(
                legend.position = "none",
                strip.background = element_rect(fill = "grey90"),
                strip.text = element_text(face = "bold")
        )

otu_table(phyloseq_16S_bacteria) |> View()

        
# ---- Rarefaction of flower samples to calculate alpha diversity metrics ----
set.seed(42)
# Rarefaction of biological samples only:
phyloseq_ITS1_rare_10k <- phyloseq::rarefy_even_depth(
        phyloseq_ITS1_fungi_biosamples,
        sample.size = 10000,
        rngseed = 42,
        replace = FALSE,
        trimOTUs = TRUE,
        verbose = TRUE
)

# Five samples were removed at rarefaction to 10 000 reads (cmt12, cmt71, cmt72, cmt73, cyb32).
# All five samples represent bagged_flower sample class, causing slight imbalance.
# At tree level, one tree was lost due to rarefaction (mt7).
# 10k rarefied data set contain 289 samples.

nsamples(phyloseq_ITS1_rare_10k)

# ---- Calculate alpha diversity metrics for flower samples ----
alpha_diversity_ITS1_rare_10k <- estimate_richness(
        phyloseq_ITS1_rare_10k,
        measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")
        ) |> 
        rownames_to_column("sample_id") |> 
        left_join(
                data.frame(
                        sample_data(phyloseq_ITS1_rare_10k)) |>
                        rownames_to_column("sample_id"),
                by = "sample_id"
        ) |> 
        dplyr::select(
                sample_id,
                Observed,
                Chao1,
                se.chao1,
                Shannon,
                Simpson,
                InvSimpson,
                sublibrary_id,
                sample_type,
                management_type,
                farm_id,
                tree_id,
                collection_batch,
                extraction_batch,
                pcr_batch_its1,
                is_pollination_clsm,
                pi_clsm,
                is_germination,
                callose_plug_index,
                is_hyphae_present,
                is_spore_present,
                fungal_colonization_score
        )

# Transform dataframe with alpha diversity and metadata to long format
alpha_diversity_ITS1_rare_10k_long <- alpha_diversity_ITS1_rare_10k |>
        # Transform to long format for plotting
        pivot_longer(
                cols = c(Observed, Chao1, `se.chao1`, Shannon, Simpson, InvSimpson),
                names_to = "metric",
                values_to = "value"
        ) |>
        mutate(
                metric = factor(metric,
                                levels = c("Observed",
                                           "Chao1",
                                           "se.chao1",
                                           "Shannon",
                                           "Simpson",
                                           "InvSimpson")
                                )
        )

# ---- Visualize alpha diversity metrics by metadata groups ----
# Plot alpha diversity metrics by flower treatment (sample_type: bagged Vs openly pollinated)
alpha_diversity_ITS1_rare_10k_long |> 
        filter(metric %in% c("Observed", "Shannon", "InvSimpson")) |>
        ggplot(aes(x = sample_type, y = value, fill = sample_type)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
        facet_wrap(~ metric, scales = "free_y") +
        theme_bw() +
        labs(
                title = "Alpha diversity of fungal communities (ITS1) by flower treatment groups",
                subtitle = "Alpha diversity metrics were calculated on an ASV table rarefied to 10,000 reads per sample",
                x = "",
                y = ""
        ) +
        theme(
                legend.position = "none",
                strip.background = element_rect(fill = "grey90"),
                strip.text = element_text(face = "bold")
        )

# Plot alpha diversity metrics by is_pollination_clsm verified by CLSM
alpha_diversity_ITS1_rare_10k_long |> 
        filter(metric %in% c("Observed", "Shannon", "InvSimpson")) |>
        filter(!is.na(is_pollination_clsm)) |> 
        ggplot(aes(x = is_pollination_clsm, y = value, fill = is_pollination_clsm)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
        facet_wrap(~ metric, scales = "free_y") +
        theme_bw() +
        labs(
                title = "Alpha diversity of fungal communities (ITS1) for unpollinated and pollinated flowers\n(verified by confocal microscopy)",
                subtitle = "Alpha iversity metrix on samples rarefied at 10,000 reads",
                x = "",
                y = ""
        ) +
        theme(
                legend.position = "bottom",
                strip.background = element_rect(fill = "grey90"),
                strip.text = element_text(face = "bold")
        )

# Plot alpha diversity metrics by is_pollination_clsm verified by CLSM
alpha_diversity_ITS1_rare_10k_long |> 
        filter(metric %in% c("Observed", "Shannon", "InvSimpson")) |>
        filter(!is.na(is_pollination_clsm)) |> 
        ggplot(aes(x = sample_type, y = value, fill = is_pollination_clsm)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(
                position = position_jitterdodge(
                        jitter.width = 0.15,
                        dodge.width = 0.75),
                alpha = 0.3, size = 0.8) +
        facet_wrap(~ metric, scales = "free_y") +
        theme_bw() +
        labs(
                title = "Alpha diversity of fungal communities (ITS1) for unpollinated and pollinated flowers\n(verified by confocal microscopy)",
                subtitle = "Alpha iversity metrix on samples rarefied at 10,000 reads",
                x = "",
                y = ""
        ) +
        theme(
                legend.position = "bottom",
                strip.background = element_rect(fill = "grey90"),
                strip.text = element_text(face = "bold")
        )

# ---- Visualize alpha diversity by farm_id ----
# Plot alpha diversity metrics by sample_type x farm_id (ordered by management_type)
alpha_diversity_ITS1_rare_10k_long |> 
        filter(metric %in% c("Observed", "Shannon", "InvSimpson")) |>
        mutate(
                farm_id = factor(
                        farm_id,
                        levels = c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
                )
        ) |> 
        ggplot(aes(x = farm_id, y = value, fill = sample_type)) +
        geom_boxplot(
                outlier.shape = NA,
                alpha = 0.7,
                position = position_dodge(width = 0.75)
                ) +
        geom_jitter(aes(colour = is_pollination_clsm),
                position = position_jitterdodge(
                        jitter.width = 0.05,
                        dodge.width = 0.05),
                alpha = 0.3,
                size = 0.8) +
        # Add sample size above the geom_jitter()
        #geom_text(
        #        aes(label = n),
        #        stat = "summary",
        #        fun = max,
        #        position = position_dodge(width = 0.75),
        #        vjust = -0.5,
        #        size = 3
        #) +
        facet_wrap(~ metric, scales = "free_y") +
        theme_bw() +
        labs(
                title = "Alpha diversity of fungal communities (ITS1) across farms",
                subtitle = "Alpha iversity metrix on samples rarefied at 10,000 reads",
                x = "",
                y = ""
        ) +
        theme(
                legend.position = "bottom",
                strip.background = element_rect(fill = "grey90"),
                strip.text = element_text(face = "bold")
        )

# ---- Alpha diversity boxplots by farm_id and management_type ----
plot_alpha_metric <- function(metric_name, y_label, show_legend = FALSE) {
        
alpha_diversity_ITS1_rare_10k_long |>
                filter(metric == metric_name) |>
                mutate(
                        farm_id = factor(
                                farm_id,
                                levels = c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
                        ),
                        management_type = factor(
                                management_type,
                                levels = c("inside_forest", "near_forest", "agroforest", "full_sun")
                        )
                ) |>
                ggplot(aes(x = farm_id, y = value, fill = sample_type)) +
                geom_boxplot(
                        outlier.shape = NA,
                        alpha = 0.7,
                        position = position_dodge(width = 0.75)
                ) +
                geom_jitter(
                        position = position_jitterdodge(
                                jitter.width = 0.05,
                                dodge.width = 0.15
                        ),
                        alpha = 0.3,
                        size = 0.8,
                        stroke = 0.2
                ) +
                ggh4x::facet_nested(
                        ~ management_type,
                        scales = "free_x",
                        space = "free_x",
                        switch = "x",
                        labeller = as_labeller(c(
                                inside_forest = "inside\nforest",
                                near_forest   = "near\nforest",
                                agroforest    = "agroforest",
                                full_sun      = "full\nsun")
                                )
                ) +
                theme_bw() +
                labs(
                        x = "",
                        y = y_label
                ) +
                theme(
                        strip.placement = "outside",
                        panel.border = element_blank(),
                        panel.spacing.x = unit(0.2, "lines"),
                        panel.spacing.y = unit(0.25, "lines"),
                        strip.padding = unit(0.25, "lines"),
                        legend.position = "none",
                        strip.background = element_blank(),
                        strip.text.x = element_text(size = 9, 
                                                    margin = margin(t = 5, b = 5))
                )
}

legend_source <-
        alpha_diversity_ITS1_rare_10k_long |>
        filter(metric == "Observed") |>
        ggplot(aes(x = farm_id, y = value, fill = sample_type)) +
        geom_boxplot(alpha = 0.7) +
        theme_bw() +
        theme(
                legend.position = "bottom",
                legend.title = element_blank()
        )

alpha_legend <- cowplot::get_legend(legend_source)

p1_observed    <- plot_alpha_metric("Observed", "Observed ASVs")
p2_shannon     <- plot_alpha_metric("Shannon", "Shannon")
p3_inv_simpson <- plot_alpha_metric("InvSimpson", "InvSimpson")

final_alpha_slide <-
        ( p1_observed | p2_shannon | p3_inv_simpson ) /
        alpha_legend +
        plot_layout(heights = c(1, 0.001))
        
# ---- Alpha diversity metrics and quantitative metadata visualizations ----
# I've tried to visulize alpha diversity against quantitative microscopy variables.
# There are few samples with high pollination intensities (> 50) so these data point have large leverage and make regression dodgy.
alpha_diversity_ITS1_rare_10k_long |>
        dplyr::filter(
                metric %in% c("Observed", "Shannon", "InvSimpson"),
                !is.na(pi_clsm),
                sample_type %in% c("bagged_flower", "unbagged_flower")
        ) |>
        ggplot(
                aes(
                        x = pi_clsm,
                        y = value,
                        colour = sample_type
                )
        ) +
        geom_point(
                alpha = 0.9,
                size = 2
        ) +
        geom_smooth(
                aes(fill = sample_type),
                method = "lm",
                se = TRUE
        ) +
        facet_wrap(. ~ metric, scales = "free_y") +
        theme_bw() +
        labs(
                x = "Pollination intensity",
                y = "",
                title = "Alpha diversity vs pollination intensity"
        )

# Alpha diversity metrics vs callose plug index
alpha_diversity_ITS1_rare_10k_long |>
        dplyr::filter(
                metric %in% c("Observed", "Shannon", "InvSimpson"),
                !is.na(is_germination),
                sample_type %in% c("bagged_flower", "unbagged_flower")
        ) |>
        ggplot(
                aes(
                        x = callose_plug_index,
                        y = value,
                        colour = sample_type
                )
        ) +
        geom_point(
                alpha = 0.9,
                size = 2
        ) +
        geom_smooth(
                aes(fill = sample_type),
                method = "lm",
                se = TRUE
        ) +
        facet_wrap(. ~ metric, scales = "free_y") +
        theme_bw() +
        labs(
                x = "Callose plug index",
                y = "",
                title = "Alpha diversity vs callose plug index"
        )



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

alpha_iters100_ITS1_multi <- iterative_rarefaction_alpha_multi_depth(
        ps_obj  = phyloseq_ITS1_fungi_biosamples,
        depths  = c(5000, 10000, 15000),
        n_iter  = 100,
        base_seed = 42,
        timestamp    = Sys.time(),
        session_info = sessionInfo()
)

# attach metadata once
alpha_iters100_ITS1_multi <- alpha_iters100_ITS1_multi |>
        left_join(
                data.frame(sample_data(phyloseq_ITS1_fungi_biosamples)) |>
                        (\(x) { class(x) <- "data.frame"; x })() |>
                        rownames_to_column("sample_id"),
                by = "sample_id"
        )

saveRDS(
        alpha_iters100_ITS1_multi,
        file = "alpha_iters100_ITS1_bundle.rds"
)

# long format
alpha_iters100_ITS1_multi_long <- alpha_iters100_ITS1_multi |>
        pivot_longer(
                cols = c(Observed, Chao1, Shannon, Simpson, InvSimpson),
                names_to = "metric",
                values_to = "value"
        )

# ---- Visualise iterative subsampling results ----
# Diagnosti� of convergence of alpha-diversity metrics in iterations 
alpha_convergence <- alpha_iters100_ITS1_multi_long |>
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

alpha_iters100_ITS1_multi_long |>
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

alpha_iters100_ITS1_multi_long |>
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

alpha_iters100_ITS1_multi_long |>
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

alpha_iters100_ITS1_multi_long |>
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


alpha_iters100_ITS1_multi |>
        dplyr::distinct(depth, iteration, sample_id) |>
        dplyr::group_by(depth) |>
        dplyr::summarise(
                min_n_samples = min(dplyr::n()),
                max_n_samples = max(dplyr::n()),
                .groups = "drop"
        )



# ---- Aggregate across iterations ----
alpha_ITS1_10k_summary <- alpha_iters_ITS1_10k |>
        pivot_longer(
                cols = c(Observed, Chao1, Shannon, Simpson, InvSimpson),
                names_to = "metric",
                values_to = "value"
        ) |>
        group_by(sample_id, metric) |>
        summarise(
                median_value = median(value, na.rm = TRUE),
                mean_value   = mean(value, na.rm = TRUE),
                sd_value     = sd(value, na.rm = TRUE),
                .groups = "drop"
        ) |>
        left_join(
                alpha_iters_ITS1_10k |>
                        distinct(sample_id, .keep_all = TRUE) |>
                        select(
                                sample_id,
                                sublibrary_id,
                                sample_type,
                                management_type,
                                farm_id,
                                tree_id,
                                collection_batch,
                                extraction_batch,
                                pcr_batch_its1,
                                is_pollination_clsm,
                                pi_clsm,
                                is_germination,
                                callose_plug_index,
                                is_hyphae_present,
                                is_spore_present,
                                fungal_colonization_score
                        ),
                by = "sample_id"
        )





# ---- Rarefy for richness metrics and no rarefaction for evenness metrics approach ----
# 1) Richness-based metrics (require rarefaction)
alpha_richness <- estimate_richness(
        phyloseq_ITS1_rare_10k,
        measures = c("Observed", "Chao1")
)

View(alpha_richness)

# 2) Evenness-sensitive metrics (can be computed on non-rarefied data)
alpha_evenness <- estimate_richness(
        phyloseq_ITS1_fungi_biosamples,
        measures = c("Shannon", "Simpson")
)

View(alpha_evenness)

# 3) Combine and attach metadata
alpha_df <- alpha_richness |>
        rownames_to_column("sample_id") |>
        left_join(
                alpha_evenness |>
                        rownames_to_column("sample_id"),
                by = "sample_id"
        ) |>
        left_join(
                data.frame(sample_data(phyloseq_ITS1_fungi_biosamples)) |>
                        rownames_to_column("sample_id"),
                by = "sample_id"
        )

View(alpha_df)

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

# ---- Linear Mixed-Effects Models for Alpha Diversity ----
# I have nested design farm(7) -> tree (7) -> flower (3 per group)
# management_type -> farm_id -> tree_id -> sample_id
fit_alpha_lmm <- function(response, data, fixed, random) {
        formula <- as.formula(
                paste(response, "~", fixed, "+", random)
        )
        lmer(formula, data = data)
        }

# ---- LMM for Observed ASVs (Community Richness) ----
# Does flower visitor exclusion (bagging) affect fungal ASV richness per flower?

# Variance for farm_id is 0, so this model is singular:
lmm_observed_ITS1_fixed_treatment_rand_farm_tree <- fit_alpha_lmm(
        response = "Observed",
        fixed = "sample_type",
        random = "(1 | farm_id / tree_id)",
        data = alpha_diversity_ITS1_rare_10k
        )

summary(lmm_observed_ITS1_fixed_treatment_rand_farm_tree)

# I removed farm_id from the model:
lmm_observed_ITS1_fixed_treatment_rand_tree <- fit_alpha_lmm(
        response = "Observed",
        fixed = "sample_type",
        random = "(1 | tree_id)",
        data = alpha_diversity_ITS1_rare_10k
        )

summary(lmm_observed_ITS1_fixed_treatment_rand_tree)

# I added management_type as fixed effect:
lmm_observed_ITS1_fixed_treatment_management_rand_farm_tree <- fit_alpha_lmm(
        response = "Observed",
        fixed = "sample_type + management_type",
        random = "(1 | farm_id / tree_id)",
        data = alpha_diversity_ITS1_rare_10k
        )

summary(lmm_observed_ITS1_fixed_treatment_management_rand_farm_tree)

# Flower visitor access does not increase the number of detectable fungal ASVs per flower after rarefaction.

# Assumptions check:
performance::check_normality(lmm_observed_ITS1_fixed_treatment_rand_tree)
performance::check_heteroscedasticity(lmm_observed_ITS1_fixed_treatment_rand_tree)
performance::check_outliers(lmm_observed_ITS1_fixed_treatment_rand_tree)

# LMM didn't pass assumptions check, Observed ASVs is a count data

# Generalised linear mixed effect model for count data with Negative Binomial distribution
glmm_observed_ITS1_fixed_treatment_rand_tree <- glmmTMB::glmmTMB(
        Observed ~ sample_type + (1|tree_id),
        family = nbinom2,
        data = alpha_diversity_ITS1_rare_10k
)

summary(glmm_observed_ITS1_fixed_treatment_rand_tree)

# DHARMa simmulation for GLMM
sim_res <- simulateResiduals(
        fittedModel = glmm_observed_ITS1_fixed_treatment_rand_tree,
        n = 1000
)

plot(sim_res)

testUniformity(sim_res)
testDispersion(sim_res)
testOutliers(sim_res)

# ---- LMM for Shannon (Richness & Evenness) ----
# Does access of flower visitors alter fungal community diversity (richness and evenness) on flowers?

# Variance for both random effects is non zero:
lmm_shannon_ITS1_fixed_treatment_rand_farm_tree <- fit_alpha_lmm(
        response = "Shannon",
        fixed = "sample_type",
        random = "(1 | farm_id / tree_id)",
        data = alpha_diversity_ITS1_rare_10k
        )

summary(lmm_shannon_ITS1_fixed_treatment_rand_farm_tree)

# I removed farm_id from the model:
lmm_shannon_ITS1_fixed_treatment_rand_tree <- fit_alpha_lmm(
        response = "Shannon",
        fixed = "sample_type",
        random = "(1 | tree_id)",
        data = alpha_diversity_ITS1_rare_10k
        )

summary(lmm_shannon_ITS1_fixed_treatment_rand_tree)

# I added management_type as fixed effect:
lmm_shannon_ITS1_fixed_treatment_management_rand_farm_tree <- fit_alpha_lmm(
        response = "Shannon",
        fixed = "sample_type + management_type",
        random = "(1 | farm_id / tree_id)",
        data = alpha_diversity_ITS1_rare_10k
        )

summary(lmm_shannon_ITS1_fixed_treatment_management_rand_farm_tree)

# Assumptions check:
performance::check_normality(lmm_shannon_ITS1_fixed_treatment_rand_tree)
performance::check_heteroscedasticity(lmm_shannon_ITS1_fixed_treatment_rand_tree)
performance::check_outliers(lmm_shannon_ITS1_fixed_treatment_rand_tree)

# ---- LMM for InvSimpson (Effective number of dominant ASVs) ----
# Does access of flower visitors increase dominance of specific fungal ASV in floral communities?

# Variance for both random effects is non zero:
lmm_InvSimpson_ITS1_fixed_treatment_rand_farm_tree <- fit_alpha_lmm(
        response = "InvSimpson",
        fixed = "sample_type",
        random = "(1 | farm_id / tree_id)",
        data = alpha_diversity_ITS1_rare_10k
        )

summary(lmm_InvSimpson_ITS1_fixed_treatment_rand_farm_tree)

# I removed farm_id from the model:
lmm_InvSimpson_ITS1_fixed_treatment_rand_tree <- fit_alpha_lmm(
        response = "InvSimpson",
        fixed = "sample_type",
        random = "(1 | tree_id)",
        data = alpha_diversity_ITS1_rare_10k
        )

summary(lmm_InvSimpson_ITS1_fixed_treatment_rand_tree)

# I added management_type as fixed effect:
lmm_InvSimpson_ITS1_fixed_treatment_management_rand_farm_tree <- fit_alpha_lmm(
        response = "InvSimpson",
        fixed = "sample_type + management_type",
        random = "(1 | farm_id / tree_id)",
        data = alpha_diversity_ITS1_rare_10k
        )

summary(lmm_InvSimpson_ITS1_fixed_treatment_management_rand_farm_tree)

# Assumptions check:
performance::check_normality(lmm_InvSimpson_ITS1_fixed_treatment_rand_tree)
performance::check_heteroscedasticity(lmm_InvSimpson_ITS1_fixed_treatment_rand_tree)
performance::check_outliers(lmm_InvSimpson_ITS1_fixed_treatment_rand_tree)

# ---- Models for alpha diversity and pic_clsm on CLSM-verified subset ----

# NB!
# These models better to try on unbagged flowers only (natural pollination):
# subset = sample_type == "unbagged_flower"

# Is fungal ASV richness associated with pollination intensity measured by CLSM?
lmm_observed_ITS1_fixed_pi_management_rand_farm_tree <- lmer(
        Observed ~ pi_clsm + management_type + (1 | farm_id / tree_id),
        data = alpha_diversity_ITS1_rare_10k,
        subset = !is.na(pi_clsm)
        )

summary(lmm_observed_ITS1_fixed_pi_management_rand_farm_tree)

# Is fungal Shannon diversity associated with pollination intensity?
lmm_shannon_ITS1_fixed_pi_management_rand_farm_tree <- lmer(
        Shannon ~ pi_clsm + management_type + (1 | farm_id / tree_id),
        data = alpha_diversity_ITS1_rare_10k,
        subset = !is.na(pi_clsm)
        )

summary(lmm_shannon_ITS1_fixed_pi_management_rand_farm_tree)

# Simpson
lmm_simpson_ITS1_fixed_pi_management_rand_farm_tree <- lmer(
        Shannon ~ pi_clsm + management_type + (1 | farm_id / tree_id),
        data = alpha_diversity_ITS1_rare_10k,
        subset = !is.na(pi_clsm)
        )

summary(lmm_simpson_ITS1_fixed_pi_management_rand_farm_tree)

# ---- LMM with interaction of bagging and pollination intensity ----
# Does pollination intensity affect fungal ASV richness differently in bagged and unbagged flowers?
# LMM: Observed ASVs ~ bagging Г— pollination intensity
lmm_observed_interaction <- lmer(
        Observed ~ sample_type * pi_clsm + management_type + (1 | tree_id),
        data = alpha_diversity_ITS1_rare_10k,
        subset = !is.na(pi_clsm)
)

summary(lmm_observed_interaction)

# Does pollination intensity modify the effect of flower visitor access on fungal community evenness?
# LMM: Shannon ~ bagging Г— pollination intensity
lmm_shannon_interaction <- lmer(
        Shannon ~ sample_type * pi_clsm + management_type + (1 | tree_id),
        data = alpha_diversity_ITS1_rare_10k,
        subset = !is.na(pi_clsm)
)

summary(lmm_shannon_interaction)

# Does pollination intensity influence dominance structure of fungal communities differently in bagged and unbagged flowers?
# LMM: InvSimpson ~ bagging Г— pollination intensity
lmm_invsimpson_interaction <- lmer(
        InvSimpson ~ sample_type * pi_clsm + management_type + (1 | tree_id),
        data = alpha_diversity_ITS1_rare_10k,
        subset = !is.na(pi_clsm)
)

summary(lmm_invsimpson_interaction)

 ---- Visulize predicted marginal means from the LMM ----
# What is the average expected alpha diversity across farms and trees, given a treatment?
plot_marginal_means <- function(model, response_label, title) {
        
        pred <- ggpredict(model, terms = "sample_type")
        
        ggplot(pred, aes(x = x, y = predicted)) +
                geom_point(size = 3) +
                geom_errorbar(
                        aes(ymin = conf.low, ymax = conf.high),
                        width = 0.15
                ) +
                theme_bw() +
                labs(
                        x = "",
                        y = response_label,
                        title = title
                )
}

plot_marginal_means <- function(model, response_label, title) {
  
  pred <- ggpredict(model, terms = "sample_type")
  
  ggplot(pred, aes(x = x, y = predicted)) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high),
      width = 0.15
    ) +
    theme_bw() +
    labs(
      x = "",
      y = response_label,
      title = title
    )
}

p_obs  <- plot_marginal_means(
        lmm_observed_ITS1_fixed_treatment_rand_tree,
        "Observed ASVs",
        ""
)

p_sha  <- plot_marginal_means(
        lmm_shannon_ITS1_fixed_treatment_rand_tree,
        "Shannon diversity",
        ""
)

p_inv  <- plot_marginal_means(
        lmm_InvSimpson_ITS1_fixed_treatment_rand_tree,
        "Inverse Simpson",
        ""
)

(p_obs | p_sha | p_inv) +
        plot_layout(widths = c(1, 1, 1))


# ---- Visulize LMM for alpha diversity with interaction of sample_type x pi_clsm ----
# InvSimpson
pred_invsimpson_treat <- ggpredict(
       pred_shannon_pi <- ggpredict(
        lmm_shannon_interaction,
        terms = c("pi_clsm", "sample_type")
        )

ggplot(pred_shannon_pi,
       aes(x = x, y = predicted,
           colour = group, fill = group)) +
        geom_line(size = 1) +
        geom_ribbon(
                aes(ymin = conf.low, ymax = conf.high),
                alpha = 0.2,
                colour = NA
        ) +
        theme_bw() +
        labs(
                x = "Pollination intensity (CLSM)",
                y = "Shannon diversity",
                colour = "Flower treatment",
                fill = "Flower treatment",
                title = "Model-predicted relationship between pollination intensity and fungal diversity"
        )
,
        terms = "sample_type"
)

ggplot(pred_invsimpson_treat,
       aes(x = x, y = predicted)) +
        geom_point(size = 3) +
        geom_errorbar(
                aes(ymin = conf.low, ymax = conf.high),
                width = 0.15
        ) +
        theme_bw() +
        labs(
                x = "",
                y = "Inverse Simpson",
                title = "Effect of flower visitor access on fungal dominance"
        )

# Shannon
pred_shannon_pi <- ggpredict(
        lmm_shannon_interaction,
        terms = c("pi_clsm", "sample_type")
)

ggplot(pred_shannon_pi,
       aes(x = x, y = predicted,
           colour = group, fill = group)) +
        geom_line(size = 1) +
        geom_ribbon(
                aes(ymin = conf.low, ymax = conf.high),
                alpha = 0.2,
                colour = NA
        ) +
        theme_bw() +
        labs(
                x = "Pollination intensity (CLSM)",
                y = "Shannon diversity",
                colour = "Flower treatment",
                fill = "Flower treatment",
                title = "Model-predicted relationship between pollination intensity and fungal diversity"
        )

# InvSimpson
pred_invsimpson_pi <- ggpredict(
        lmm_invsimpson_interaction,
        terms = c("pi_clsm", "sample_type")
)

ggplot(pred_invsimpson_pi,
       aes(x = x, y = predicted,
           colour = group, fill = group)) +
        geom_line(size = 1) +
        geom_ribbon(
                aes(ymin = conf.low, ymax = conf.high),
                alpha = 0.2,
                colour = NA
        ) +
        theme_bw() +
        labs(
                x = "Pollination intensity (CLSM)",
                y = "Inverse Simpson",
                colour = "Flower treatment",
                fill = "Flower treatment",
                title = "Pollination intensity and dominance structure of fungal communities"
        )





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
unique(metadata$`ITS1_pcr_batch`)

# ---- Batch vs group association checks ----
meta <- data.frame(sample_data(phyloseq_prev2_ab0)) |>
        rownames_to_column("sample_id") |>
        mutate(
                extraction_batch = as.factor(extraction_batch),
                `XITS1_pcr_batch`  = as.factor(`XITS1_pcr_batch`),
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
tab_pcr_farm <- table(meta$`XITS1_pcr_batch`, meta$farm_id)
tab_pcr_type <- table(meta$`XITS1_pcr_batch`, meta$sample_type)

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

# ---- Filter samples by sample total read count (sequencing depth) ----
# Check the samples with total read count per sample lower than 5000:
sample_sums(phyloseq_ITS1_fungi)[
        sample_sums(phyloseq_ITS1_fungi) < 10000
        ]

# Remove samples with read count per sample lower than 5000:
phyloseq_ITS1_fungi_high_qual_samples <- prune_samples(
        sample_sums(phyloseq_ITS1_fungi_biosamples) >= 5000,
        phyloseq_ITS1_fungi_biosamples
        )

nsamples(phyloseq_ITS1_raw)
nsamples(phyloseq_ITS1_decontam)
nsamples(phyloseq_ITS1_fungi)

nsamples(phyloseq_ITS1_fungi_biosamples)
# These removes only 3 flowers issued from the same tree (namely cmt71, cmt72, cmt73)
nsamples(phyloseq_ITS1_fungi_high_qual_samples)


# ---- Filter taxa by prevalence across biological samples ----
# Different parameters for prevalence filtering were applied 
# apply(abscence/presence_matrix, MARGIN == 1 rows, sum) sum by rows result in number of samples containing ASV
filter_by_prevalence <- function(phyloseq, min_prev) {
        prune_taxa(
                apply(otu_table(phyloseq) > 0, 1, sum) > min_prev,
                phyloseq
        )
}

taxa_are_rows(otu_table(phyloseq_ITS1_fungi_biosamples))

# How many samples have this taxa?
phyloseq_ITS1_fungi_biosamples_prev1  <- filter_by_prevalence(phyloseq_ITS1_fungi_biosamples, 1)
phyloseq_ITS1_fungi_biosamples_prev2  <- filter_by_prevalence(phyloseq_ITS1_fungi_biosamples, 2)
phyloseq_ITS1_fungi_biosamples_prev3  <- filter_by_prevalence(phyloseq_ITS1_fungi_biosamples, 3)
phyloseq_ITS1_fungi_biosamples_prev9  <- filter_by_prevalence(phyloseq_ITS1_fungi_biosamples, 9)

# compare total taxa present after four different prevalence filtering strategies
data.frame(
        prevalence_threshold = c("no prev filter", ">1", ">2", ">3", ">9"),
        ntaxa = c(
                ntaxa(phyloseq_ITS1_fungi_biosamples),
                ntaxa(phyloseq_ITS1_fungi_biosamples_prev1),
                ntaxa(phyloseq_ITS1_fungi_biosamples_prev2),
                ntaxa(phyloseq_ITS1_fungi_biosamples_prev3),
                ntaxa(phyloseq_ITS1_fungi_biosamples_prev9)
        )
)

# Taxa per sample with non-zero counts after different prevalence filtering approaches
count_taxa_per_sample <- function(phyloseq) {
        apply(
                otu_table(phyloseq) > 0,
                2,
                sum
        )
        }

taxa_per_sample <- data.frame(
        sample_id = sample_names(phyloseq_ITS1_fungi_biosamples),
        no_prev = count_taxa_per_sample(phyloseq_ITS1_fungi_biosamples),
        prev1   = count_taxa_per_sample(phyloseq_ITS1_fungi_biosamples_prev1),
        prev2   = count_taxa_per_sample(phyloseq_ITS1_fungi_biosamples_prev2),
        prev3   = count_taxa_per_sample(phyloseq_ITS1_fungi_biosamples_prev3),
        prev9   = count_taxa_per_sample(phyloseq_ITS1_fungi_biosamples_prev9)
        )


summary(taxa_per_sample[, -1])

taxa_per_sample_grouped <- taxa_per_sample |>
        dplyr::left_join(
                data.frame(
                        sample_id = sample_names(phyloseq_ITS1_fungi_biosamples),
                        sample_type = sample_data(phyloseq_ITS1_fungi_biosamples)$sample_type
                ),
                by = "sample_id"
        )

class(taxa_per_sample_grouped)
View(taxa_per_sample_grouped)

# ---- Robustness of BC dissimilarity after different prevalence filtering ----
# bounded between 0 (the same communities) and 1 (no matches)
# relative abundance based => use Jaccard for presence/abscence analysis
# non-phylogenetic => use UniFrac after MSA & phylogeny building
# BC dissimilarity is directly related to the quantitative Sorensen similarity index 

compute_bray_curtis <- function(phyloseq) {
        phyloseq_rel <- transform_sample_counts(
                phyloseq,
                function(x) x / sum(x)
        )

        phyloseq::distance(phyloseq_rel, method = "bray")
}

bc_dissim_no_prev <- compute_bray_curtis(phyloseq_ITS1_fungi_bio)
bc_dissim_prev1   <- compute_bray_curtis(phyloseq_prev1)
bc_dissim_prev2   <- compute_bray_curtis(phyloseq_prev2)
bc_dissim_prev3   <- compute_bray_curtis(phyloseq_prev3)
bc_dissim_prev9   <- compute_bray_curtis(phyloseq_prev9)

# Ordination: Principle Coordinate Analysis
ord_no_prev <- ordinate(phyloseq_ITS1_fungi_bio, method = "PCoA", distance = bc_dissim_no_prev)
ord_prev1   <- ordinate(phyloseq_prev1, method = "PCoA", distance = bc_dissim_prev1)
ord_prev2   <- ordinate(phyloseq_prev2, method = "PCoA", distance = bc_dissim_prev2)
ord_prev3   <- ordinate(phyloseq_prev3, method = "PCoA", distance = bc_dissim_prev3)
ord_prev9   <- ordinate(phyloseq_prev9, method = "PCoA", distance = bc_dissim_prev9)

# Visualization PCoA
plot_ordination(
        phyloseq_ITS1_fungi_bio,
        ord_no_prev,
        color = "sample_type"
) + ggtitle("BrayвЂ“Curtis: no prevalence filter")

plot_ordination(
        phyloseq_prev1,
        ord_prev1,
        color = "sample_type"
) + ggtitle("BrayвЂ“Curtis: prevalence > 1")

plot_ordination(
        phyloseq_prev2,
        ord_prev2,
        color = "sample_type"
) + ggtitle("BrayвЂ“Curtis: prevalence > 2")

plot_ordination(
        phyloseq_prev3,
        ord_prev3,
        color = "sample_type"
) + ggtitle("BrayвЂ“Curtis: prevalence > 3")

plot_ordination(
        phyloseq_prev9,
        ord_prev9,
        color = "sample_type"
) + ggtitle("BrayвЂ“Curtis: prevalence > 9")

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

# ---- Filter taxa by read count per samples (abundance) ----
# Filter taxa by total read count across all samples
filter_by_abundance <- function(phyloseq, min_reads) {
        prune_taxa(
                taxa_sums(phyloseq) >= min_reads,
                phyloseq
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
count_taxa_per_sample <- function(phyloseq) {
        apply(
                otu_table(phyloseq) > 0,
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

# ---- Bray-Curtis dissimilarity after different prevalence filtering approaches ----
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





# ---- Compute Bray-Curtis (relative abundance based) and Jaccard (presence/absence) dissimilarity matricies ----
compute_beta_dissimilarity <- function(phyloseq, method = c("bray", "jaccard")) {
        method <- match.arg(method)
        
        if (method == "bray") {
                phyloseq_rel <- transform_sample_counts(
                        phyloseq,
                        function(x) x / sum(x)
                )
                return(phyloseq::distance(phyloseq_rel, method = "bray"))
        }
        
        if (method == "jaccard") {
                phyloseq_pa <- transform_sample_counts(
                        phyloseq,
                        function(x) {
                                y <- as.integer(x > 0)
                                names(y) <- names(x)   # РёРЅР°С‡Рµ Сѓ РјРµРЅСЏ Р·РґРµСЃСЊ РїР°РґР°СЋС‚ РёРјРµРЅР°
                                y
                        }
                )
                return(phyloseq::distance(phyloseq_pa, method = "jaccard"))
        }
}

# 1) Calculate dissimilarity matrices for the whole dataset (unfiltered by prevalence/abundance, not rarefied).
# Comparison of visitor effect: bagged Vs unbagged
bc_ITS1_fungi_biosamples <- compute_beta_dissimilarity(phyloseq_ITS1_fungi_biosamples, "bray")
jac_ITS1_fungi_biosamples <- compute_beta_dissimilarity(phyloseq_ITS1_fungi_biosamples, "jaccard")

# 2) Calculate dissimilarity matrices for the unbagged_flower with CLSM phenotyping dataset (unfiltered by prevalence/abundance, not rarefied).
# Comparison of pollinator effect in natural condition: 

# Filter phyloseq object unbagged with available CLSM data on pollination
phyloseq_ITS1_unbagged_clsm <- phyloseq_ITS1_fungi_biosamples |> 
        subset_samples(
                sample_type == "unbagged_flower" &
                        !is.na(is_pollination_clsm)
        ) |> 
        # Prune taxa with zero counts in this sample subset
        (\(x) prune_taxa(taxa_sums(x) > 0, x))()

bc_ITS1_fungi_unbagged_clsm  <- compute_beta_dissimilarity(phyloseq_ITS1_unbagged_clsm, method = "bray")
jac_ITS1_fungi_unbagged_clsm <- compute_beta_dissimilarity(phyloseq_ITS1_unbagged_clsm, method = "jaccard")

# 3) Calculate dissimilarity matrices for samples with CLSM phenotyping of stigma fungal infection (unfiltered by prevalence/abundance, not rarefied).
# Comparison of fungal flower communities with different intensity of microscopy validated fungal stigma infection.
phyloseq_ITS1_clsm_hyphae <- phyloseq_ITS1_fungi_biosamples |>
        subset_samples(!is.na(is_hyphae_present)) |>
        (\(x) prune_taxa(taxa_sums(x) > 0, x))()

bc_ITS1_fungi_clsm_hyphae <- compute_beta_dissimilarity(phyloseq_ITS1_clsm_hyphae, method = "bray")
jac_ITS1_fungi_clsm_hyphae <- compute_beta_dissimilarity(phyloseq_ITS1_clsm_hyphae, method = "jaccard")


# ---- Some diagnostic commands for dist class objects -- matrices of pairwaise sample dissimilarity ----
class(bc_ITS1_fungi_biosamples)
typeof(jac_ITS1_fungi_unbagged_clsm)
dim(jac_ITS1_fungi_unbagged_clsm)
as.matrix(bc_ITS1_fungi_biosamples)[1:5, 1:5]
as.matrix(jac_ITS1_fungi_unbagged_clsm)[1:5, 1:5]
bc_ITS1_fungi_biosamples["cib11", "cib12"]
jac_ITS1_fungi_unbagged_clsm["pmt11", "pmt12"]


# ---- Ordination (dimensionality reduction): PCoA (MDS) & NMDS ----
# PCoA (MDS)
ord_pcoa_bc_ITS1_fungi_biosamples  <- ordinate(phyloseq_ITS1_fungi_biosamples, "PCoA", bc_ITS1_fungi_biosamples)
ord_pcoa_jac_ITS1_fungi_biosamples <- ordinate(phyloseq_ITS1_fungi_biosamples, "PCoA", jac_ITS1_fungi_biosamples)

ord_pcoa_bc_ITS1_unbagged_clsm <- ordinate(phyloseq_ITS1_unbagged_clsm, "PCoA", bc_ITS1_fungi_unbagged_clsm)
ord_pcoa_jac_ITS1_unbagged_clsm <- ordinate(phyloseq_ITS1_unbagged_clsm, "PCoA", jac_ITS1_fungi_unbagged_clsm)

ord_pcoa_bc_ITS1_fungi_clsm_hyphae <- ordinate(phyloseq_ITS1_clsm_hyphae, "PCoA", bc_ITS1_fungi_clsm_hyphae)


# How much of variation is explained by Ax1, Ax2 and Ax3? 
ord_pcoa_bc_ITS1_fungi_biosamples$values$Relative_eig[1:3]
ord_pcoa_jac_ITS1_fungi_biosamples$values$Relative_eig[1:3]

ord_pcoa_bc_ITS1_unbagged_clsm$values$Relative_eig[1:3]
ord_pcoa_jac_ITS1_unbagged_clsm$values$Relative_eig[1:3]

ord_pcoa_bc_ITS1_fungi_clsm_hyphae$values$Relative_eig[1:3]


# ---- Optional: NMDS ordination ----
ord_nmds_bc_ITS1_fungi_biosamples  <- ordinate(phyloseq_ITS1_fungi_biosamples, "NMDS", bc_ITS1_fungi_biosamples)
ord_nmds_jac_ITS1_fungi_biosamples <- ordinate(phyloseq_ITS1_fungi_biosamples, "NMDS", jac_ITS1_fungi_biosamples)

ord_nmds_bc_ITS1_unbagged_clsm <- ordinate(phyloseq_ITS1_unbagged_clsm, "NMDS", bc_ITS1_fungi_unbagged_clsm)
ord_nmds_jac_ITS1_unbagged_clsm <- ordinate(phyloseq_ITS1_unbagged_clsm, "NMDS", jac_ITS1_fungi_unbagged_clsm)

# ---- Visualize PCoA by major categorical variables ----
plot_pcoa_by <- function(phyloseq, ord, color_var, title_suffix) {
        plot_ordination(
                phyloseq,
                ord,
                color = color_var
        ) +
                theme_bw() +
                ggtitle(paste("PCoA (BrayвЂ“Curtis):", title_suffix)) +
                theme(
                        legend.title = element_text(size = 10),
                        legend.text  = element_text(size = 9)
                )
}

# ---- Bray Curtis dissimilarity PCoA: categorical variables ----
# PCoA visualizations by key experimental variables in the whole dataset
# Bagged Vs Unbagged 
plot_pcoa_by(
        phyloseq_ITS1_fungi_biosamples,
        ord_pcoa_bc_ITS1_fungi_biosamples,
        "sample_type",
        "flower visitation (bagged vs open)"
)

# Management types
plot_pcoa_by(
        phyloseq_ITS1_fungi_biosamples,
        ord_pcoa_bc_ITS1_fungi_biosamples,
        "management_type",
        "management type"
)

# Farm ID
plot_pcoa_by(
        phyloseq_ITS1_fungi_biosamples,
        ord_pcoa_bc_ITS1_fungi_biosamples,
        "farm_id",
        "farm ID"
)

# Visualize fungal community on stigma with available CLSM data on stigma fungal infection
plot_pcoa_by(
        phyloseq_ITS1_clsm_hyphae,
        ord_pcoa_bc_ITS1_fungi_clsm_hyphae,
        "fungal_colonization_score",
        "Stigma fungal infection score"
        ) +
        # Add red outline for samples with microscopically visible fungal spores (is_spore_present flag in metadata)
        geom_point(
                data =
                plot_ordination(
                        phyloseq_ITS1_clsm_hyphae,
                        ord_pcoa_bc_ITS1_fungi_clsm_hyphae,
                        justDF = TRUE
                ) |>
                dplyr::filter(is_spore_present),
                aes(x = Axis.1, y = Axis.2),
                shape = 21,        # allows outline
                colour = "red",    # outline colour
                fill = NA,         # keep original fill from below
                size = 3,
                stroke = 1
                )

# ---- Diagnostic PCoA visualizations (technical variables) ----
plot_pcoa_by(
        phyloseq_ITS1_fungi_biosamples,
        ord_pcoa_bc_ITS1_fungi_biosamples,
        "is_pollination_clsm",
        "CLSM-confirmed pollination"
)

plot_pcoa_by(
        phyloseq_ITS1_fungi_biosamples,
        ord_pcoa_bc_ITS1_fungi_biosamples,
        "collection_batch",
        "Collection batch"
)

plot_pcoa_by(
        phyloseq_ITS1_fungi_biosamples,
        ord_pcoa_bc_ITS1_fungi_biosamples,
        "extraction_batch",
        "Extraction batch"
)

plot_pcoa_by(
        phyloseq_ITS1_fungi_biosamples,
        ord_pcoa_bc_ITS1_fungi_biosamples,
        "its1_pcr_batch",
        "ITS1 PCR batch"
)

plot_pcoa_by(
        phyloseq_ITS1_fungi_biosamples,
        ord_pcoa_bc_ITS1_fungi_biosamples,
        "sublibrary_id",
        "Sublibrary ID (the same external barcodes)"
)

# Nested design
plot_ordination(
        phyloseq_ITS1_fungi_biosamples,
        ord_pcoa_bc_ITS1_fungi_biosamples,
        color = "tree_id"
) +
        theme_bw() +
        ggtitle("PCoA (BrayвЂ“Curtis): tree-level clustering") +
        guides(color = guide_legend(ncol = 2))

# ---- Visualize PCoA with phenotyping data (unbagged flowers, CLSM only) ----
# Add transformed variable to sample_data
sample_data(phyloseq_ITS1_unbagged_clsm)$log1p_pi_clsm <-
        log1p(sample_data(phyloseq_ITS1_unbagged_clsm)$pi_clsm)

# Visualize unbagged flower fungal communties with different pollination intensities
plot_ordination(
        phyloseq_ITS1_unbagged_clsm,
        ord_pcoa_bc_ITS1_unbagged_clsm,
        color = "log1p_pi_clsm"
) +
        scale_color_viridis_c(
                na.value = "grey80",
                name = expression(log(1 + PI[CLSM]))
        ) +
        theme_bw() +
        labs(
                title = "PCoA (BrayвЂ“Curtis dissimilarity)",
                subtitle = "Unbagged flowers with CLSM phenotyping\nColour indicates pollination intensity"
        ) +
        theme(
                legend.title = element_text(size = 10),
                legend.text  = element_text(size = 9)
        )

# Callose plug index
plot_ordination(
        phyloseq_ITS1_unbagged_clsm,
        ord_pcoa_bc_ITS1_unbagged_clsm,
        color = "callose_plug_index"
) +
        scale_color_viridis_c(
                na.value = "grey80",
                name = "Callose plug index"
        ) +
        theme_bw() +
        labs(
                title = "PCoA (BrayвЂ“Curtis dissimilarity)",
                subtitle = "Unbagged flowers with CLSM phenotyping\nColour indicates callose plug abundance"
        ) +
        theme(
                legend.title = element_text(size = 10),
                legend.text  = element_text(size = 9)
        )


# ---- Betadisper and PERMANOVA (flower treatment: visitor effect) ----
# ---- 1) Visitor effect ----
# H0: Composition of the fungal microbiome does not differ between bagged and unbagged cocoa flowers.
# H1: Exposure of flowers to insect visitors (unbagged treatment) is associated with a shift in fungal community composition relative to bagged flowers.

# Metadata dataframe
metadata_biosamples <- data.frame(sample_data(phyloseq_ITS1_fungi_biosamples))

# ---- Variance homogeneity check (BC) ----
betadisper_bc_sample_type <- betadisper(
        bc_ITS1_fungi_biosamples,
        group = metadata_biosamples$sample_type
        )

permutest(betadisper_bc_sample_type, permutations = 999)

boxplot(
        betadisper_bc_sample_type,
        main = "BrayвЂ“Curtis dispersion by sample type",
        ylab = "Distance to centroid"
        )

# ---- PERMANOVA: visitor effect (BC) ----
adonis_bc_visitors <- adonis2(
        bc_ITS1_fungi_biosamples ~ sample_type,
        data = metadata_biosamples,
        permutations = 999,
        strata = metadata_biosamples$tree_id
)

adonis_bc_visitors

# ---- Variance homogeneity check (Jaccard) ----
betadisper_jaccard_sample_type <- betadisper(
        jac_ITS1_fungi_biosamples,
        group = metadata_biosamples$sample_type
)

permutest(
        betadisper_jaccard_sample_type,
        permutations = 999
)

boxplot(
        betadisper_jaccard_sample_type,
        main = "Jaccard dispersion by sample type",
        ylab = "Distance to centroid"
)

# ---- PERMANOVA: visitor effect (Jaccard) ----
adonis_jaccard_visitors <- adonis2(
        jac_ITS1_fungi_biosamples ~ sample_type,
        data = metadata_biosamples,
        permutations = 999,
        strata = metadata_biosamples$tree_id
)

adonis_jaccard_visitors

# ---- Envfit: linking flower fungal community structure to pollination intensity and pollen performance ----

# I use envfit to test whether continuous CLSM-derived phenotypes show a significant directional association with flower fungal community structure in a reduced multidimensional (PCoA) space.

# Extract sample metadata for unbagged flowers with CLSM phenotyping
metadata_unbagged_clsm <- data.frame(
        sample_data(phyloseq_ITS1_unbagged_clsm)
)

# Extract the first three PCoA axes to capture more community variation
ordination_scores_pcoa_bc_unbagged_clsm <-
        ord_pcoa_bc_ITS1_unbagged_clsm$vectors[, 1:3]


# ---- Envfit: pollination intensity (pi_clsm) ----
# Tests whether variation in flower fungal community structure is aligned
# with a gradient of pollination intensity quantified by CLSM.

envfit_pollination_intensity <- envfit(
        ordination_scores_pcoa_bc_unbagged_clsm,
        metadata_unbagged_clsm$pi_clsm,
        permutations = 999,
        na.rm = TRUE
)

envfit_pollination_intensity


# ---- Envfit: pollen tube performance (callose_plug_index) ----
# Tests whether flower fungal community structure covaries with
# pollen tube growth performance as reflected by callose plug index (blob detected by LoG algorithm)

envfit_callose_plug_index <- envfit(
        ordination_scores_pcoa_bc_unbagged_clsm,
        metadata_unbagged_clsm$callose_plug_index,
        permutations = 999,
        na.rm = TRUE
)

envfit_callose_plug_index

# ---- Betadisper and PERMANOVA (stigma fungal infection: CLSM hyphae) ----
# ---- 3) Fungal infection effect ----
# H0: Composition of the fungal microbiome does not differ among levels of stigma fungal colonization.
# H1: Increasing fungal colonization severity (0вЂ“5 lobes affected) is associated with shifts in fungal community composition.

# Metadata dataframe (CLSM hyphae subset)
metadata_clsm_hyphae <- data.frame(sample_data(phyloseq_ITS1_clsm_hyphae))

# Ensure fungal colonization score is an ordered factor (safety check)
metadata_clsm_hyphae$fungal_colonization_score <- factor(
        metadata_clsm_hyphae$fungal_colonization_score,
        levels = 0:5,
        ordered = TRUE
)

# ---- Variance homogeneity check (BC) ----
betadisper_bc_fungal_colonization_score <- betadisper(
        bc_ITS1_fungi_clsm_hyphae,
        group = metadata_clsm_hyphae$fungal_colonization_score
)

permutest(
        betadisper_bc_fungal_colonization_score,
        permutations = 999
)

boxplot(
        betadisper_bc_fungal_colonization_score,
        main = "BrayвЂ“Curtis dispersion by fungal colonization score",
        ylab = "Distance to centroid"
)

# ---- PERMANOVA: fungal infection effect (BC) ----
adonis_bc_fungal_colonization_score <- adonis2(
        bc_ITS1_fungi_clsm_hyphae ~ fungal_colonization_score,
        data = metadata_clsm_hyphae,
        permutations = 999,
        strata = metadata_clsm_hyphae$tree_id
)

adonis_bc_fungal_colonization_score

# ---- Variance homogeneity check (Jaccard) ----
betadisper_jaccard_fungal_colonization_score <- betadisper(
        jac_ITS1_fungi_clsm_hyphae,
        group = metadata_clsm_hyphae$fungal_colonization_score
)

permutest(
        betadisper_jaccard_fungal_colonization_score,
        permutations = 999
)

boxplot(
        betadisper_jaccard_fungal_colonization_score,
        main = "Jaccard dispersion by fungal colonization score",
        ylab = "Distance to centroid"
)

# ---- PERMANOVA: fungal infection effect (Jaccard) ----
adonis_jaccard_fungal_colonization_score <- adonis2(
        jac_ITS1_fungi_clsm_hyphae ~ fungal_colonization_score,
        data = metadata_clsm_hyphae,
        permutations = 999,
        strata = metadata_clsm_hyphae$tree_id
)

adonis_jaccard_fungal_colonization_score


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
incidence_df_ITS1 <-
        phyloseq_ITS1_fungi_biosamples |>
        otu_table() |>
        as("matrix") |>
        t() |>
        (\(x) {
                x[x > 0] <- 1
                as.data.frame(x)
        })()

# Prepare metadata as dataframe to analyse zeta-decline by experimental and phenotypical groups
metadata_ITS1 <-
        sample_data(phyloseq_ITS1_fungi_biosamples) |>
        as("data.frame")

all(rownames(incidence_df_ITS1) ==
            rownames(metadata_ITS1))

# Compute zeta-diversity decline (z1 - z10 components)
# using Monte Carlo sampling method (1000 samples)
zeta_bagged_ITS1 <-
        incidence_df_ITS1[
                metadata_ITS1$sample_type == "bagged_flower",
        ] |>
        Zeta.decline.mc(
                orders = 1:10,
                sam = 1000,
                plot = FALSE
        )

zeta_unbagged_ITS1 <-
        incidence_df_ITS1[
                metadata_ITS1$sample_type == "unbagged_flower",
        ] |>
        Zeta.decline.mc(
                orders = 1:10,
                sam = 1000,
                plot = FALSE
        )

# ---- Diagnostic commands for Zetadiv.decline.mc() output lists ----
# Compare coefficients of determination for power and exponential law fit:
# how well the fitted model explains the variance in the observed zeta decline
summary(zeta_bagged_ITS1$zeta.pl)$r.squared
summary(zeta_bagged_ITS1$zeta.exp)$r.squared

summary(zeta_unbagged_ITS1$zeta.pl)$r.squared
summary(zeta_unbagged_ITS1$zeta.exp)$r.squared

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

extract_zeta_fits(zeta_bagged_ITS1)
extract_zeta_fits(zeta_unbagged_ITS1)


# Agregate zeta-diversity values for bagged and unbagged flowers
zeta_plot_tbl_ITS1 <-
        dplyr::bind_rows(
                data.frame(
                        order = zeta_bagged_ITS1$zeta.order,
                        zeta = zeta_bagged_ITS1$zeta.val,
                        zeta_sd = zeta_bagged_ITS1$zeta.val.sd,
                        group = "bagged_flower"
                ),
                data.frame(
                        order = zeta_unbagged_ITS1$zeta.order,
                        zeta = zeta_unbagged_ITS1$zeta.val,
                        zeta_sd = zeta_unbagged_ITS1$zeta.val.sd,
                        group = "unbagged_flower"
                )
        )

# Visualize zeta-diversity decline
ggplot(zeta_plot_tbl_ITS1,
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
                y = "Zeta diversity component (average number of shared ASVs)",
                title = "Zeta diversity decline curves for ITS1 ASVs",
                colour = "Flower treatment"
        )


################################################################################
## Indicator species analysis (ISA), ANCOM-BC & Differential Abundance #########
################################################################################
# ---- ISA: bagged vs unbagged (fungi associated with insect visitors) ----
# ---- Sanity checks prior ISA ----
table(sample_data(phyloseq_ITS1_fungi_biosamples)$sample_type)

boxplot(sample_sums(phyloseq_ITS1_fungi_biosamples) ~ 
                sample_data(phyloseq_ITS1_fungi_biosamples)$sample_type)

# ---- Presence/abscence ISA without prevalence or abundance filtering on unrarefied data ----
set.seed(42)
indval_bagged_vs_unbagged_result <-
        phyloseq_ITS1_fungi_biosamples |>
        (\(x)
         multipatt(
                 # Get ASV count matrix from phyloseq object
                 # Transpose count matrix to get: samples x ASVs
                 # Convert count matrix to presence/absence matrix
                 t(as(otu_table(x), "matrix")) > 0,
                 # Extract grouping variable values for the samples from the same phyloseq object
                 sample_data(x)$sample_type,
                 func = "IndVal.g",
                 control = how(nperm = 999)
         )
        )()

# ---- Robustness of ISA on pre-filtered phyloseq dataset ----
# Remove:
# 1) global abandunce singleton ASV, 2) ASVs with low prevalence across samples
# This steps reduces the number of permutation tests?

set.seed(42)
indval_bagged_vs_unbagged_result_no_singletons_prev1 <-
        phyloseq_ITS1_fungi_biosamples |>
        # Pruned global singletons (read count == 1 across all dataset; 77 taxa pruned)
        # As expected it didn't affect ISA (singleton taxa couldn't be indicator taxa)
        (\(x) prune_taxa(taxa_sums(x) > 1, x) )() |> 
        # Filter ASVs with prevalence == 1 (found only in one sample and also couldn't be indicator)
        (\(x) prune_taxa(apply(otu_table(x) > 0, 1, sum) > 1, x))()|> 
        (\(x)
         multipatt(
                 # Get ASV count matrix from phyloseq object
                 # Transpose count matrix to get: samples x ASVs
                 # Convert count matrix to presence/absence matrix
                 t(as(otu_table(x), "matrix")) > 0,
                 # Extract grouping variable values for the samples from the same phyloseq object as vector
                 sample_data(x)$sample_type,
                 func = "IndVal.g",
                 control = how(nperm = 999)
         )
        )()


summary(indval_bagged_vs_unbagged_result)
summary.multipatt(indval_bagged_vs_unbagged_result)

# ---- Processin ISA results bagged Vs unbagged ----
# Extract IndVal results as dataframe
indval_bagged_vs_unbagged_df <- indval_bagged_vs_unbagged_result$sign

# Filter ISA results: ASVs significantly associated with unbagged flower group
indic_asv_bagged_vs_unbagged <- rownames(
        indval_bagged_vs_unbagged_df[
                #indval_bagged_vs_unbagged_df$index == 2 &
                        !is.na(indval_bagged_vs_unbagged_df$p.value) &
                        indval_bagged_vs_unbagged_df$p.value <= 0.01 &
                                indval_bagged_vs_unbagged_df$stat >= 0.3,
        ]
)

# Map taxonomy of indicator species 
taxonomy_indic_bagged_vs_unbagged <- as.data.frame(taxonomy_table_ITS1)[indic_asv_bagged_vs_unbagged,]

# Build master data frame for visualization of ISA results
filtered_highIndVal_bagged_vs_unbagged <- indval_bagged_vs_unbagged_result$sign |>
        as.data.frame() |>
        dplyr::filter(!is.na(p.value), p.value <= 0.01, stat >= 0.3) |>
        (\(df) {
                # Extract otu table for selected ASVs
                pa <- otu_table(phyloseq_ITS1_fungi_biosamples)[rownames(df), ]
                # Convert OTU table to presence/absence matrix to calculate prevalence of ASV by group
                pa <- as(pa, "matrix") > 0
                
                # Sample grouping variable (aligned with otu_table column names)
                grp <- sample_data(phyloseq_ITS1_fungi_biosamples)$sample_type
                
                # Ensure ASVs are rows
                if (!taxa_are_rows(phyloseq_ITS1_fungi_biosamples)) {
                        pa <- t(pa)
                }
                
                tibble::tibble(
                        ASV = rownames(pa),
                        # Group for which IndVal was calculated
                        indicator_group = dplyr::recode(
                                as.character(df$index),
                                `1` = "bagged_flower",
                                `2` = "unbagged_flower"
                        ),
                        # IndVal strength of association of indicator ASV
                        stat = df$stat,
                        # calculate prevalence of ASVs across groups
                        prevalence_bagged_flower =
                                rowMeans(pa[, grp == "bagged_flower", drop = FALSE]),
                        prevalence_unbagged_flower =
                                rowMeans(pa[, grp == "unbagged_flower", drop = FALSE])
                )
        })() |> 
        dplyr::left_join(
                as.data.frame(tax_table(phyloseq_ITS1_fungi_biosamples)) |>
                        tibble::rownames_to_column("ASV"),
                by = "ASV"
        )

# ---- Diagnostic commands for ISA results ----
# Collect summart on multipatt() output:
summary(indval_bagged_vs_unbagged_result)
summary.multipatt(indval_bagged_vs_unbagged_result)

# Visualize distribution of stat (IndVal.g)
indval_bagged_vs_unbagged_result$sign |>
        as.data.frame() |>
        dplyr::filter(!is.na(stat)) |>
        ggplot(aes(x = stat)) +
        geom_histogram(bins = 40) +
        labs(
                x = "IndVal statistic",
                y = "Number of ASVs",
                title = "Distribution of IndVal statistics"
        ) +
        theme_minimal()

# Explore dataframe with filtered IndVal results
dim(filtered_highIndVal_bagged_vs_unbagged)
class(filtered_highIndVal_bagged_vs_unbagged)

sample_data(phyloseq_ITS1_fungi_biosamples)$sample_type 

View(taxonomy_table_ITS1)
View(taxonomy_indic_bagged_vs_unbagged)

taxonomy_indic_bagged_vs_unbagged$Species
taxonomy_indic_bagged_vs_unbagged$Genus

# Concordance of indicator ASVs mapping to the same taxa (Genus)
indval_bagged_vs_unbagged_result$sign |>
        as.data.frame() |>
        dplyr::filter(!is.na(p.value), p.value <= 0.01 & stat >= 0.3) |>
        (\(df)
         dplyr::bind_cols(
                 df,
                 as.data.frame(taxonomy_table_ITS1)[rownames(df), c("Species", "Genus", "Family")]
         )
        )() |>
        dplyr::filter(!is.na(Genus)) |>
        dplyr::mutate(
                group = dplyr::recode(
                        as.character(index),
                        `1` = "bagged_flower",
                        `2` = "unbagged_flower"
                )
        ) |>
        dplyr::group_by(group, Genus) |>
        dplyr::summarise(
                n_asv = dplyr::n(),
                max_stat = max(stat, na.rm = TRUE),
                mean_stat = mean(stat),
                .groups = "drop"
        ) |>
        dplyr::arrange(desc(n_asv), desc(max_stat)) |>
        View()


# ---- Heat map visualization for unbagged Vs bagged indicators (the best solution so far) ----
make_indicator_heatmap <- function(df, indicator) {
        
        df |>
                dplyr::filter(indicator_group == indicator) |>
                dplyr::filter(!is.na(Species)) |>
                dplyr::group_by(Species) |>
                dplyr::summarise(
                        n_ASV = dplyr::n(),
                        mean_prev_bagged   = mean(prevalence_bagged_flower, na.rm = TRUE),
                        mean_prev_unbagged = mean(prevalence_unbagged_flower, na.rm = TRUE),
                        ordering =
                                mean(prevalence_unbagged_flower - prevalence_bagged_flower, na.rm = TRUE),
                        .groups = "drop"
                ) |>
                # Cut-off for visualizing number of aggregated ASVs into the Taxon
                dplyr::filter(n_ASV >= 1) |> 
                dplyr::arrange(ordering) |>
                dplyr::mutate(
                        Species_label = paste0(Species, " (", n_ASV, " ASVs)"),
                        Species_label = factor(Species_label, levels = unique(Species_label))
                ) |>
                tidyr::pivot_longer(
                        cols = c(mean_prev_bagged, mean_prev_unbagged),
                        names_to = "group",
                        values_to = "prevalence"
                ) |>
                dplyr::mutate(
                        group = dplyr::recode(
                                group,
                                mean_prev_bagged   = "bagged_flower",
                                mean_prev_unbagged = "unbagged_flower"
                        ),
                        group = factor(group, levels = c("unbagged_flower", "bagged_flower"))
                ) |>
                ggplot(aes(x = group, y = Species_label, fill = prevalence)) +
                geom_tile(color = "grey90") +
                scale_fill_viridis_c(
                        option = "magma",
                        limits = c(0, 1)
                ) +
                guides(
                        fill = guide_colorbar(
                                title.position = "top",
                                title.hjust = 0.5
                        )
                ) +
                labs(
                        x = NULL,
                        y = NULL,
                        fill = "Proportion of flowers\nwith taxon detected",
                        title = paste("Indicator Species for", indicator, 
                                      "(IndVal.g ≥ 0.3, permutation test p < 0.01)")
                ) +
                theme_minimal() +
                theme(panel.grid = element_blank(),
                      plot.title = element_text(hjust = 0)
                      )
}

# Creating patches
p_unbagged <- make_indicator_heatmap(
        filtered_highIndVal_bagged_vs_unbagged,
        "unbagged_flower"
        )

p_bagged <- make_indicator_heatmap(
        filtered_highIndVal_bagged_vs_unbagged,
        "bagged_flower"
)

# Assemble patches in one figure adjusting their hights to number of rows
n_unbagged <- n_distinct(
        filtered_highIndVal_bagged_vs_unbagged |>
                dplyr::filter(indicator_group == "unbagged_flower", !is.na(Species)) |>
                pull(Family)
)

n_bagged <- n_distinct(
        filtered_highIndVal_bagged_vs_unbagged |>
                dplyr::filter(indicator_group == "bagged_flower", !is.na(Species)) |>
                pull(Family)
)

# Collect patchwork
(p_unbagged / p_bagged) +
        plot_layout(
                heights = c(n_unbagged, n_bagged),
                guides = "collect"
        ) &
        theme(legend.position = "right")

# ---- Make_indicator_heatmap2 Function ----
make_indicator_heatmap2 <- function(df,
                                    indicator,
                                    rank,
                                    min_asv = 1) {
        # rank: taxonomic rank used for aggregation (e.g. "Family", "Genus", "Species")
        # indicator: indicator_group to plot (e.g. "unbagged_flower")
        # min_asv: minimum number of ASVs aggregated into a taxon to keep it
        
        rank_sym <- rlang::sym(rank)
        
        df |>
                # keep indicators of the selected group only
                dplyr::filter(indicator_group == indicator) |>
                
                # aggregate strictly at the chosen taxonomic rank
                dplyr::group_by(!!rank_sym) |>
                dplyr::summarise(
                        n_ASV = dplyr::n(),
                        
                        mean_prev_bagged =
                                mean(prevalence_bagged_flower, na.rm = TRUE),
                        
                        mean_prev_unbagged =
                                mean(prevalence_unbagged_flower, na.rm = TRUE),
                        
                        # ordering metric: contrast between groups
                        ordering =
                                mean(prevalence_unbagged_flower -
                                             prevalence_bagged_flower,
                                     na.rm = TRUE),
                        
                        .groups = "drop"
                ) |>
                
                # optional cut-off on the number of aggregated ASVs
                dplyr::filter(n_ASV >= min_asv) |>
                
                # order taxa for plotting
                dplyr::arrange(ordering) |>
                
                # construct display label with fallback for Incertae sedis
                dplyr::mutate(
                        Taxon_label = dplyr::case_when(
                                !is.na(!!rank_sym) &
                                        !grepl("Incertae_sedis", !!rank_sym) ~
                                        as.character(!!rank_sym),
                                
                                !is.na(Order)  &
                                        !grepl("Incertae_sedis", Order)  ~
                                        as.character(Order),
                                
                                !is.na(Class)  &
                                        !grepl("Incertae_sedis", Class)  ~
                                        as.character(Class),
                                
                                !is.na(Phylum) &
                                        !grepl("Incertae_sedis", Phylum) ~
                                        as.character(Phylum),
                                
                                TRUE ~ as.character(!!rank_sym)
                        ),
                        
                        Taxon_label =
                                paste0(Taxon_label, " (", n_ASV, " ASVs)"),
                        
                        Taxon_label =
                                factor(Taxon_label,
                                       levels = unique(Taxon_label))
                ) |>
                
                # reshape for heatmap
                tidyr::pivot_longer(
                        cols = c(mean_prev_bagged, mean_prev_unbagged),
                        names_to = "group",
                        values_to = "prevalence"
                ) |>
                
                dplyr::mutate(
                        group = dplyr::recode(
                                group,
                                mean_prev_bagged   = "bagged_flower",
                                mean_prev_unbagged = "unbagged_flower"
                        ),
                        group = factor(
                                group,
                                levels = c("unbagged_flower", "bagged_flower")
                        )
                ) |>
                
                # plot
                ggplot2::ggplot(
                        ggplot2::aes(
                                x = group,
                                y = Taxon_label,
                                fill = prevalence
                        )
                ) +
                ggplot2::geom_tile(color = "grey90") +
                ggplot2::scale_fill_viridis_c(
                        option = "magma",
                        limits = c(0, 1)
                ) +
                ggplot2::labs(
                        x = NULL,
                        y = NULL,
                        fill = "Proportion of flowers with taxon detected",
                        title = paste("Indicator", rank, "for", indicator)
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                        panel.grid = ggplot2::element_blank()
                )
}



# ---- Heat bubble visualization: insect visitation effect (bagged versus unbagged) ----
species_level_indicator_bagged_vs_unbagged_long |> 
        ggplot(aes(
                x = group,
                y = Genus,
                fill = prevalence,
                size = n_ASV
        )) +
        geom_point(shape = 21, colour = "black") +
        scale_fill_viridis_c(
                option = "viridis",
                limits = c(0, 1)
        ) +
        scale_size(range = c(2, 8)) +
        labs(
                x = NULL,
                y = "Fungal Genus",
                # Mean prevalence (presence/absence) within treatment group
                # Prevalence across flower treatments
                fill = "Proportion of flowers with taxon detected",
                # Number of ASVs agreagated into the taxon
                size = "Number of ASVs"
        ) +
        theme_minimal()

?arrange()


# ---- ISA: bagged vs unbagged (fungi associated with insect visitors) after global prevalence filtration ----
table(sample_data(phyloseq_ITS1_fungi_biosamples_prev2)$sample_type)

boxplot(sample_sums(phyloseq_ITS1_fungi_biosamples_prev2) ~ 
                sample_data(phyloseq_ITS1_fungi_biosamples_prev2)$sample_type)

# Presence/abscence ISA on prevalence filtered data (only ASVs with prev >= 3 samples)
set.seed(42)
indval_bagged_vs_unbagged_prev2_result <-
        phyloseq_ITS1_fungi_biosamples_prev2 |>
        (\(x)
         multipatt(
                 # Get ASV count matrix from phyloseq object
                 # Transpose count matrix to get: samples x ASVs
                 # Convert count matrix to presence/absence matrix
                 t(as(otu_table(x), "matrix")) > 0,
                 # Extract grouping variable values for the samples from the same phyloseq object
                 sample_data(x)$sample_type,
                 func = "IndVal.g",
                 control = how(nperm = 999)
         )
        )()

summary(indval_bagged_vs_unbagged_prev2_result)

# ---- Processin ISA results bagged Vs unbagged (after prevalence >= 3 cleaning) ----
# Extract IndVal results as dataframe
indval_bagged_vs_unbagged_prev2_df <- indval_bagged_vs_unbagged_prev2_result$sign

# Filter ISA results: ASVs significantly associated with unbagged flower group
indic_asv_bagged_vs_unbagged_prev2 <- rownames(
        indval_bagged_vs_unbagged_prev2_df[
                #indval_bagged_vs_unbagged_df$index == 2 &
                !is.na(indval_bagged_vs_unbagged_prev2_df$p.value) &
                        indval_bagged_vs_unbagged_prev2_df$p.value <= 0.01 &
                        indval_bagged_vs_unbagged_prev2_df$stat >= 0.3,
        ]
)

# Map taxonomy of indicator species 
taxonomy_indic_bagged_vs_unbagged_prev2 <- as.data.frame(taxonomy_table_ITS1)[indic_asv_bagged_vs_unbagged_prev2,]

# Build master data frame for visualization
filtered_highIndVal_bagged_vs_unbagged_prev2 <- indval_bagged_vs_unbagged_prev2_result$sign |>
        as.data.frame() |>
        dplyr::filter(!is.na(p.value), p.value <= 0.01, stat >= 0.3) |>
        (\(df) {
                # Extract otu table for selected ASVs
                pa <- otu_table(phyloseq_ITS1_fungi_biosamples_prev2)[rownames(df), ]
                # Convert OTU table to presence/absence matrix to calculate prevalence of ASV by group
                pa <- as(pa, "matrix") > 0
                
                # Sample grouping variable (aligned with otu_table column names)
                grp <- sample_data(phyloseq_ITS1_fungi_biosamples_prev2)$sample_type
                
                # Ensure ASVs are rows
                if (!taxa_are_rows(phyloseq_ITS1_fungi_biosamples_prev2)) {
                        pa <- t(pa)
                }
                
                tibble::tibble(
                        ASV = rownames(pa),
                        # Group for which IndVal was calculated
                        indicator_group = dplyr::recode(
                                as.character(df$index),
                                `1` = "bagged_flower",
                                `2` = "unbagged_flower"
                        ),
                        # IndVal strength of association of indicator ASV
                        stat = df$stat,
                        # calculate prevalence of ASVs across groups
                        prevalence_bagged_flower =
                                rowMeans(pa[, grp == "bagged_flower", drop = FALSE]),
                        prevalence_unbagged_flower =
                                rowMeans(pa[, grp == "unbagged_flower", drop = FALSE])
                )
        })() |> 
        dplyr::left_join(
                as.data.frame(tax_table(phyloseq_ITS1_fungi_biosamples_prev2)) |>
                        tibble::rownames_to_column("ASV"),
                by = "ASV"
        )

# ---- Heat map visualization for unbagged Vs bagged indicators (the best solution so far) ----
# Creating patches
p_unbagged_prev2 <- make_indicator_heatmap(
        filtered_highIndVal_bagged_vs_unbagged_prev2,
        "unbagged_flower"
)

p_bagged_prev2 <- make_indicator_heatmap(
        filtered_highIndVal_bagged_vs_unbagged_prev2,
        "bagged_flower"
)

# Assemble patches in one figure adjusting their hights to number of rows
n_unbagged_prev2 <- n_distinct(
        filtered_highIndVal_bagged_vs_unbagged_prev2 |>
                dplyr::filter(indicator_group == "unbagged_flower", !is.na(Species)) |>
                pull(Family)
)

n_bagged_prev2 <- n_distinct(
        filtered_highIndVal_bagged_vs_unbagged_prev2 |>
                dplyr::filter(indicator_group == "bagged_flower", !is.na(Species)) |>
                pull(Family)
)

# Collect patchwork
(p_unbagged_prev2 / p_bagged_prev2) +
        plot_layout(
                heights = c(n_unbagged_prev2, n_bagged_prev2),
                guides = "collect"
        ) &
        theme(legend.position = "right")



# ---- ISA: Are there any fungal indicator ASVs/Taxa for pollinated and unpollinated flowers among natural pollination (unbagged group)? ----
# ---- Sanity checks prior ISA ----
# Class imbalance sanity-check:
# 1) 38 samples with no pollination;
# 2) 15 samples with high pollination i.e. PI > 34 (50% of maximal fruit set);
# 3) 34 samples with medium-high pollination PI > 10
subset_samples(phyloseq_ITS1_fungi_biosamples, sample_type == "unbagged_flower" & pi_clsm == 0)
subset_samples(phyloseq_ITS1_fungi_biosamples, sample_type == "unbagged_flower" & pi_clsm >= 34)
subset_samples(phyloseq_ITS1_fungi_biosamples, sample_type == "unbagged_flower" & pi_clsm >= 7)
# 4) is_pollination_clsm => True Vs False = 68 Vs 38
sample_data(phyloseq_ITS1_fungi_biosamples_prev2) |> 
        as("data.frame") |> 
        filter(sample_type == "unbagged_flower") |> 
        filter(is_pollination_clsm == TRUE) |>
        filter(pi_clsm > 4) |> 
        count()

sample_data(phyloseq_ITS1_fungi_biosamples_prev2) |> 
        as("data.frame") |> 
        filter(sample_type == "unbagged_flower") |> 
        filter(is_pollination_clsm == FALSE) |> 
        count()

# Check sequencing depths in two contrast pollination phenotype groups
tibble::tibble(
        sequencing_depth = c(
                sample_sums(
                        subset_samples(
                                phyloseq_ITS1_fungi_biosamples,
                                sample_type == "unbagged_flower" & pi_clsm == 0
                        )
                ),
                sample_sums(
                        subset_samples(
                                phyloseq_ITS1_fungi_biosamples,
                                sample_type == "unbagged_flower" & pi_clsm >= 43
                        )
                )
        ),
        group = c(
                rep("unpollinated", 38),
                rep("highly_pollinated", 15)
        )
        ) |>
        ggplot(aes(x = group, y = sequencing_depth)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
        labs(
                x = NULL,
                y = "Sequencing depth (reads per sample)"
        ) +
        theme_minimal()

# ---- Approach 1: not pollinated Vs pollination PI > 10 (on unbagged flower class only) ----
# Trial ISA with bigger pollinated class (lower pollination intensity cut off => weaker signal from pollinators)
# Classes more balanced (34 pollinated Vs 38 unpollinated)
set.seed(42)
indval_unpollinated_vs_pollinated_result <-
        phyloseq_ITS1_fungi_biosamples_prev2 |>
        subset_samples(
                sample_type == "unbagged_flower" & (pi_clsm == 0 | pi_clsm >= 10)
                ) |> 
        # Pruned global singletons (abundance = 1) couldn't be indicators anyway
        (\(x) prune_taxa(taxa_sums(x) > 1, x) )() |> 
        # Prune ASVs with prevalence == 1 (found only in one sample and also couldn't be indicator)
        (\(x) prune_taxa(apply(otu_table(x) > 0, 1, sum) > 1, x))()|> 
        (\(x)
         multipatt(
                 # Get ASV count matrix from phyloseq object
                 # Transpose count matrix to get: samples x ASVs
                 # Convert count matrix to presence/absence matrix
                 t(as(otu_table(x), "matrix")) > 0,
                 # Extract grouping variable values for the samples from the same phyloseq object
                 as(sample_data(x)$pi_clsm, "logical"),
                 func = "IndVal.g",
                 control = how(nperm = 999)
         )
        )()

summary(indval_unpollinated_vs_pollinated_result)

# Subuliphorum camptosporum is a specific species of entomopathogenic (insect-killing) fungus
indval_unpollinated_vs_pollinated_result$sign["0cc6fb4cf9e727329279c03a75accc58",]


asv_list_indicators <- indval_unpollinated_vs_pollinated_result$sign |> 
        filter(!is.na(p.value) & index %in% c(1, 2) & p.value <= 0.05) |> 
        row.names()

tax_table(phyloseq_ITS1_fungi_biosamples_prev2)[asv_list_indicators, ] |>
        as.data.frame() |> 
        dplyr::pull(Species)

tax_table(phyloseq_ITS1_fungi_biosamples) |> View()

# Approach 2: Repeated subsampling of unpollinated control group to balance sample size with fixed high pollination group (PI > 34)
# Create a subset for indicator species analysis
phyloseq_pollination_contrast <-
        phyloseq_ITS1_fungi_biosamples |>
        subset_samples(
                sample_type == "unbagged_flower" &
                        (pi_clsm == 0 | pi_clsm >= 43)
        )

# Create a factor variable pollination_class in sample data for this subset
sample_data(phyloseq_pollination_contrast)$pollination_class <-
        factor(
                ifelse(sample_data(phyloseq_pollination_contrast)$pi_clsm >= 43,
                       "high_pollination",
                       "unpollinated"),
                levels = c("unpollinated", "high_pollination")
        )

# sequencing depth check
boxplot(sample_sums(phyloseq_pollination_contrast) ~ 
                sample_data(phyloseq_pollination_contrast)$pollination_class)

# ---- Approach 2: not pollinated Vs pollination PI > 43 (on umbagged flower class only) with repetitive downsampling of unpollinated class to tackle imbalance ----
# ---- Variant 1: 100 resampling iterations  ----
set.seed(42)

n_iter <- 100
n_pollinated <- 15

pollinated_ids <-
        sample_names(
                subset_samples(
                        phyloseq_ITS1_fungi_biosamples,
                        sample_type == "unbagged_flower" & pi_clsm >= 43
                )
        )

unpollinated_ids <-
        sample_names(
                subset_samples(
                        phyloseq_ITS1_fungi_biosamples,
                        sample_type == "unbagged_flower" & pi_clsm == 0
                )
        )

indval_unpollinated_vs_pollinated_results_downsample <-
        replicate(n_iter, {
                # sample() doesn't replace by default, but I stated it explicitly to avoid increasing prevalance as result of adding same samples several times during resampling
                sampled_unpollinated <-
                        sample(unpollinated_ids, n_pollinated, replace = FALSE)
                
                ps_resample <-
                        phyloseq_pollination_contrast |>
                        (\(x)
                         prune_samples(
                                 sample_names(x) %in%
                                         c(pollinated_ids, sampled_unpollinated),
                                 x
                         )
                        )()
                
                ps_resample |>
                        (\(x) prune_taxa(taxa_sums(x) > 1, x))() |>
                        (\(x)
                         prune_taxa(
                                 apply(otu_table(x) > 0, 1, sum) > 1,
                                 x
                         )
                        )() |>
                        (\(x)
                         multipatt(
                                 t(as(otu_table(x), "matrix")) > 0,
                                 sample_data(x)$pollination_class ==
                                         "high_pollination",
                                 func = "IndVal.g",
                                 control = how(nperm = 999)
                         )
                        )() |>
                        (\(res) {
                                if (is.null(res$sign)) return(NULL)
                                # save ASV IDs as column to bind them later
                                tibble::as_tibble(res$sign, rownames = "ASV_ID")
                        })()
        # the output will be a list length 100 (number of itteratios)
        }, simplify = FALSE)

# ---- Variant 2: increased number of resampling iterations to 1000 to explore the effect on IndVal and visualization ----
# There're 15,471,286,560 combinations without repetition when choosing 15 items from a set of 38.
# 15.5 * 10^9 

set.seed(42)

n_iter <- 1000
n_pollinated <- 15
# the output of this multipatt() is a list of 1000 runs (full multipatt objects)
indval_unpollinated_vs_pollinated_results_downsample <-
        replicate(n_iter, {
                
                # keep sizes consistent with current objects
                pollinated_ids <-
                        sample_names(
                                subset_samples(
                                        phyloseq_pollination_contrast,
                                        pollination_class == "high_pollination"
                                )
                        )
                
                unpollinated_ids <-
                        sample_names(
                                subset_samples(
                                        phyloseq_pollination_contrast,
                                        pollination_class == "unpollinated"
                                )
                        )
                
                n_pollinated <- length(pollinated_ids)
                
                # resample unpollinated control samples with no replacement
                # replacemnte might increase prevelance by introducing the same samples multiple times
                sampled_unpollinated <-
                        sample(unpollinated_ids, n_pollinated, replace = FALSE)
                
                ps_resample <-
                        phyloseq_pollination_contrast |>
                        (\(x)
                         prune_samples(
                                 sample_names(x) %in%
                                         c(pollinated_ids, sampled_unpollinated),
                                 x
                         )
                        )() |>
                        # drop unused factor levels after pruning samples
                        (\(x) { sample_data(x)$pollination_class <- droplevels(sample_data(x)$pollination_class); x })()
                
                ps_resample |>
                        # prune global ASV singletons (abundance)
                        (\(x) prune_taxa(taxa_sums(x) > 1, x))() |>
                        # prune ASV with prevalence 1 across the whole subset
                        (\(x)
                         prune_taxa(
                                 apply(otu_table(x) > 0, 1, sum) > 1,
                                 x
                         )
                        )() |>
                        (\(x)
                         multipatt(
                                 t(as(otu_table(x), "matrix")) > 0,
                                 sample_data(x)$pollination_class,
                                 func = "IndVal.g",
                                 control = how(nperm = 999)
                         )
                        )()
                
        }, simplify = FALSE)

class(indval_unpollinated_vs_pollinated_results_downsample)
length(indval_unpollinated_vs_pollinated_results_downsample)
summary(indval_unpollinated_vs_pollinated_results_downsample)

# Build a long table of per-iteration results from res$sign
indval_sign_long <-
        purrr::imap_dfr(
                indval_unpollinated_vs_pollinated_results_downsample,
                \(res, iter_id) {
                        if (is.null(res) || is.null(res$sign)) return(NULL)
                        
                        tibble::as_tibble(res$sign, rownames = "ASV_ID") |>
                                dplyr::mutate(
                                        iteration = as.integer(iter_id),
                                        group = dplyr::case_when(
                                                index == 1 ~ "unpollinated",
                                                index == 2 ~ "high_pollination",
                                                index == 3 ~ "both",
                                                TRUE ~ NA_character_
                                        )
                                )
                }
        )


# Summarise stability; focus on testable, group-specific signals (index 1 or 2)
indval_downsample_summary <-
        indval_sign_long |>
        # remove the non-contrast "both groups" pattern from indicator interpretation
        dplyr::filter(group %in% c("unpollinated", "high_pollination")) |>
        dplyr::group_by(ASV_ID, group) |>
        dplyr::summarise(
                # how often this ASV is returned with this group-pattern at all
                freq_reported = dplyr::n() / n_iter,
                
                # how often p-value is actually available (should be ~1 here, but keep it explicit)
                freq_testable = mean(!is.na(p.value)),
                
                # stability of statistical support across resamples
                freq_p_le_0_05 = mean(!is.na(p.value) & p.value <= 0.05),
                
                mean_stat = mean(stat, na.rm = TRUE),
                median_p = median(p.value, na.rm = TRUE),
                
                .groups = "drop"
        ) |>
        dplyr::arrange(dplyr::desc(freq_p_le_0_05), dplyr::desc(mean_stat), dplyr::desc(freq_reported))

View(indval_downsample_summary)

# Filter indicator ASVs 
ASV_ID_high_pollination <- indval_downsample_summary |> 
        filter(freq_reported >= 0.9 & group == "high_pollination" & freq_p_le_0_05 >= 0.05) |> 
        (\(x) x$ASV_ID)()

# Map taxonomy
as.data.frame(tax_table(phyloseq_pollination_contrast)[ASV_ID_high_pollination, ]) |> View()

indval_downsample_summary |> 
        filter(ASV_ID == "X0bac2651c3034bb37f8d785ac3b61017")

tax_table(phyloseq_pollination_contrast) |> View()
# ---- Visualize resampling effort effect ----
compute_stability_metrics <- function(sign_long, max_iter, step = 10) {
        
        iter_points <- seq(step, max_iter, by = step)
        
        purrr::map_dfr(iter_points, function(k) {
                
                sign_long |>
                        dplyr::filter(iteration <= k) |>
                        dplyr::filter(group %in% c("unpollinated", "high_pollination")) |>
                        dplyr::group_by(ASV_ID, group) |>
                        dplyr::summarise(
                                freq_reported = dplyr::n() / k,
                                freq_p_le_0_05 = mean(!is.na(p.value) & p.value <= 0.05),
                                median_p = median(p.value, na.rm = TRUE),
                                mean_stat = mean(stat, na.rm = TRUE),
                                .groups = "drop"
                        ) |>
                        dplyr::mutate(n_iter = k)
                
        })
}

stability_long <-
        compute_stability_metrics(
                indval_sign_long,
                max_iter = n_iter,
                step = 10
        )

top_asv <-
        indval_downsample_summary |>
        dplyr::filter(group == "high_pollination") |>
        dplyr::arrange(desc(freq_reported), desc(mean_stat)) |>
        dplyr::slice_head(n = 150) |>
        dplyr::pull(ASV_ID)

stability_long |>
        dplyr::filter(ASV_ID %in% top_asv) |>
        tidyr::pivot_longer(
                cols = c(freq_reported, freq_p_le_0_05, median_p),
                names_to = "metric",
                values_to = "value"
        ) |>
        ggplot(aes(x = n_iter, y = value, color = ASV_ID)) +
        geom_line() +
        facet_wrap(~ metric, scales = "free_y") +
        theme_minimal() +
        labs(
                x = "Number of resampling iterations",
                y = "Metric value",
                title = "Stabilisation of ISA-derived metrics under resampling"
        ) +
        theme(legend.position = "none")



# ---- Variant 3: resample unpollinated flower class 1000 times use prevalence filtered dataset as input ----
# Created a phyloseq object from prevalence filtered phyloseq
phyloseq_pollination_contrast_prev2 <-
        phyloseq_ITS1_fungi_biosamples_prev2 |>
        subset_samples(
                sample_type == "unbagged_flower" &
                        (pi_clsm == 0 | pi_clsm >= 43)
        )

# Create a factor variable pollination_class in sample data for this subset
sample_data(phyloseq_pollination_contrast_prev2)$pollination_class <-
        factor(
                ifelse(sample_data(phyloseq_pollination_contrast_prev2)$pi_clsm >= 43,
                       "high_pollination",
                       "unpollinated"),
                levels = c("unpollinated", "high_pollination")
        )
# Run multipatt() with 1000 downsampling and comparisons 
set.seed(42)

n_iter <- 1000
n_pollinated <- 15
# the output of this multipatt() is a list of 1000 runs (full multipatt objects)
indval_unpollinated_vs_pollinated_results_downsample_prev2 <-
        replicate(n_iter, {
                
                # keep sizes consistent with current objects
                pollinated_ids <-
                        sample_names(
                                subset_samples(
                                        phyloseq_pollination_contrast_prev2,
                                        pollination_class == "high_pollination"
                                )
                        )
                
                unpollinated_ids <-
                        sample_names(
                                subset_samples(
                                        phyloseq_pollination_contrast_prev2,
                                        pollination_class == "unpollinated"
                                )
                        )
                
                n_pollinated <- length(pollinated_ids)
                
                # resample unpollinated control samples with no replacement
                # replacemnte might increase prevelance by introducing the same samples multiple times
                sampled_unpollinated <-
                        sample(unpollinated_ids, n_pollinated, replace = FALSE)
                
                ps_resample <-
                        phyloseq_pollination_contrast_prev2 |>
                        (\(x)
                         prune_samples(
                                 sample_names(x) %in%
                                         c(pollinated_ids, sampled_unpollinated),
                                 x
                         )
                        )() |>
                        # drop unused factor levels after pruning samples
                        (\(x) { sample_data(x)$pollination_class <- droplevels(sample_data(x)$pollination_class); x })()
                
                ps_resample |>
                        # prune global ASV singletons (abundance)
                        (\(x) prune_taxa(taxa_sums(x) > 1, x))() |>
                        # prune ASV with prevalence 1 across the whole subset
                        (\(x)
                         prune_taxa(
                                 apply(otu_table(x) > 0, 1, sum) > 1,
                                 x
                         )
                        )() |>
                        (\(x)
                         multipatt(
                                 t(as(otu_table(x), "matrix")) > 0,
                                 sample_data(x)$pollination_class,
                                 func = "IndVal.g",
                                 control = how(nperm = 999)
                         )
                        )()
                
        }, simplify = FALSE)

class(indval_unpollinated_vs_pollinated_results_downsample_prev2)
length(indval_unpollinated_vs_pollinated_results_downsample_prev2)
summary(indval_unpollinated_vs_pollinated_results_downsample_prev2)

# Build a long table of per-iteration results from res$sign
indval_sign_long_prev2 <-
        purrr::imap_dfr(
                indval_unpollinated_vs_pollinated_results_downsample_prev2,
                \(res, iter_id) {
                        if (is.null(res) || is.null(res$sign)) return(NULL)
                        
                        tibble::as_tibble(res$sign, rownames = "ASV_ID") |>
                                dplyr::mutate(
                                        iteration = as.integer(iter_id),
                                        group = dplyr::case_when(
                                                index == 1 ~ "unpollinated",
                                                index == 2 ~ "high_pollination",
                                                index == 3 ~ "both",
                                                TRUE ~ NA_character_
                                        )
                                )
                }
        )


# Summarise stability; focus on testable, group-specific signals (index 1 or 2)
indval_downsample_prev2_summary <-
        indval_sign_long_prev2 |>
        # remove the non-contrast "both groups" pattern from indicator interpretation
        dplyr::filter(group %in% c("unpollinated", "high_pollination")) |>
        dplyr::group_by(ASV_ID, group) |>
        dplyr::summarise(
                # how often this ASV is returned with this group-pattern at all
                freq_reported = dplyr::n() / n_iter,
                
                # how often p-value is actually available (should be ~1 here, but keep it explicit)
                freq_testable = mean(!is.na(p.value)),
                
                # stability of statistical support across resamples
                freq_p_le_0_05 = mean(!is.na(p.value) & p.value <= 0.05),
                
                mean_stat = mean(stat, na.rm = TRUE),
                median_p = median(p.value, na.rm = TRUE),
                
                .groups = "drop"
        ) |>
        dplyr::arrange(dplyr::desc(freq_p_le_0_05), dplyr::desc(mean_stat), dplyr::desc(freq_reported))

# Visualize stability of indicator indentification during resampling
stability_long_prev2 <-
        compute_stability_metrics(
                indval_sign_long_prev2,
                max_iter = n_iter,
                step = 10
        )

top_asv <-
        indval_downsample_prev2_summary |>
        dplyr::filter(group == "high_pollination") |>
        dplyr::arrange(desc(freq_reported), desc(mean_stat)) |>
        dplyr::slice_head(n = 150) |>
        dplyr::pull(ASV_ID)

stability_long_prev2 |>
        dplyr::filter(ASV_ID %in% top_asv) |>
        tidyr::pivot_longer(
                cols = c(freq_reported, freq_p_le_0_05, median_p),
                names_to = "metric",
                values_to = "value"
        ) |>
        ggplot(aes(x = n_iter, y = value, color = ASV_ID)) +
        geom_line() +
        facet_wrap(~ metric, scales = "free_y") +
        theme_minimal() +
        labs(
                x = "Number of resampling iterations",
                y = "Metric value",
                title = "Stabilisation of ISA-derived metrics under resampling"
        ) +
        theme(legend.position = "none")

# Filter stable indicator ASVs
indicator_ASVs_prev2 <-
        indval_downsample_prev2_summary |>
        dplyr::filter(
                freq_reported >= 0.7,
                freq_p_le_0_05 >= 0.25,
                mean_stat >= 0.3,
                group %in% c("high_pollination", "unpollinated")
        )

# Map taxonomy for indicator ASVs
taxonomy_indicators_prev2 <-
        tax_table(phyloseq_pollination_contrast_prev2)[
                indicator_ASVs_prev2$ASV_ID,
        ] |>
        as.data.frame() |>
        tibble::rownames_to_column("ASV_ID") |>
        dplyr::left_join(
                indicator_ASVs_prev2,
                by = "ASV_ID"
        )
View(taxonomy_indicators_prev2)
taxonomy_indicators_prev2$Species


# ---- Aproach 3: strightforward comparison by is_pollination_clsm (68 pollinated Vs 38 unpollinated samples on unbagged flowers only) ----
# Though classes are imbalanced the sample size is bigger so the permutation test should be more powerfull
# Create a phyloseq object from prevalence filtered phyloseq
phyloseq_pollination_contrast_is_pollination_clsm_prev2 <-
        phyloseq_ITS1_fungi_biosamples_prev2 |>
        subset_samples(
                sample_type == "unbagged_flower" &
                        !is.na(is_pollination_clsm)
        )

ntaxa(phyloseq_pollination_contrast_is_pollination_clsm_prev2)
nsamples(phyloseq_pollination_contrast_is_pollination_clsm_prev2)
summary(sample_data(phyloseq_pollination_contrast_is_pollination_clsm_prev2)$is_pollination_clsm)

# Presence/abscence ISA on prevalence filtered data by boolean is_pollination_clsm on ubagged flowers only
set.seed(42)
indval_bagged_vs_unbagged_prev2_is_pollination_clsm_result <-
        phyloseq_pollination_contrast_is_pollination_clsm_prev2 |>
        (\(x)
         multipatt(
                 # Get ASV count matrix from phyloseq object
                 # Transpose count matrix to get: samples x ASVs
                 # Convert count matrix to presence/absence matrix
                 t(as(otu_table(x), "matrix")) > 0,
                 # Extract grouping variable values for the samples from the same phyloseq object
                 sample_data(x)$is_pollination_clsm,
                 func = "IndVal.g",
                 control = how(nperm = 999)
         )
        )()

summary(indval_bagged_vs_unbagged_prev2_is_pollination_clsm_result)

# ---- Processin ISA results bagged Vs unbagged (after prevalence >= 3 cleaning) ----
# Extract IndVal results as dataframe
indval_bagged_vs_unbagged_prev2_is_pollination_clsm_df <- indval_bagged_vs_unbagged_prev2_is_pollination_clsm_result$sign

# Filter ISA results: ASVs significantly associated with unbagged flower group
indic_asv_bagged_vs_unbagged_prev2_is_pollination_clsm <- rownames(
        indval_bagged_vs_unbagged_prev2_is_pollination_clsm_df[
                #indval_bagged_vs_unbagged_df$index == 2 &
                !is.na(indval_bagged_vs_unbagged_prev2_is_pollination_clsm_df$p.value) &
                        indval_bagged_vs_unbagged_prev2_is_pollination_clsm_df$p.value <= 0.05 &
                        indval_bagged_vs_unbagged_prev2_is_pollination_clsm_df$stat >= 0.25,
        ]
)

# Map taxonomy of indicator species 
taxonomy_indic_bagged_vs_unbagged_prev2_is_pollination_clsm <- as.data.frame(taxonomy_table_ITS1)[indic_asv_bagged_vs_unbagged_prev2_is_pollination_clsm,]

# Build master data frame for visualization
filtered_highIndVal_bagged_vs_unbagged_prev2_is_pollination_clsm <- indval_bagged_vs_unbagged_prev2_is_pollination_clsm_result$sign |>
        as.data.frame() |>
        dplyr::filter(!is.na(p.value), p.value <= 0.05, stat >= 0.25) |>
        (\(df) {
                # Extract otu table for selected ASVs
                pa <- otu_table(phyloseq_pollination_contrast_is_pollination_clsm_prev2)[rownames(df), ]
                # Convert OTU table to presence/absence matrix to calculate prevalence of ASV by group
                pa <- as(pa, "matrix") > 0
                
                # Sample grouping variable (aligned with otu_table column names)
                grp <- sample_data(phyloseq_pollination_contrast_is_pollination_clsm_prev2)$is_pollination_clsm
                
                # Ensure ASVs are rows
                if (!taxa_are_rows(phyloseq_pollination_contrast_is_pollination_clsm_prev2)) {
                        pa <- t(pa)
                }
                
                tibble::tibble(
                        ASV = rownames(pa),
                        # Group for which IndVal was calculated
                        indicator_group = dplyr::recode(
                                as.character(df$index),
                                `1` = "unpollinated",
                                `2` = "pollinated"
                        ),
                        # IndVal strength of association of indicator ASV
                        stat = df$stat,
                        # calculate prevalence of ASVs across groups
                        prevalence_unpollinated_flower =
                                rowMeans(pa[, grp == FALSE, drop = FALSE]),
                        prevalence_pollinated_flower =
                                rowMeans(pa[, grp == TRUE, drop = FALSE])
                )
        })() |> 
        dplyr::left_join(
                as.data.frame(tax_table(phyloseq_pollination_contrast_is_pollination_clsm_prev2)) |>
                        tibble::rownames_to_column("ASV"),
                by = "ASV"
        )

# ---- Heat map visualization for unbagged Vs bagged indicators (the best solution so far) ----
make_indicator_pollination_heatmap <- function(df, indicator) {
        
        df |>
                dplyr::filter(indicator_group == indicator) |>
                dplyr::filter(!is.na(Species)) |>
                dplyr::group_by(Species) |>
                dplyr::summarise(
                        n_ASV = dplyr::n(),
                        mean_prev_unpollinated   = mean(prevalence_unpollinated_flower, na.rm = TRUE),
                        mean_prev_pollinated = mean(prevalence_pollinated_flower, na.rm = TRUE),
                        ordering =
                                mean(prevalence_unpollinated_flower - prevalence_pollinated_flower, na.rm = TRUE),
                        .groups = "drop"
                ) |>
                # Cut-off for visualizing number of aggregated ASVs into the Taxon
                dplyr::filter(n_ASV >= 1) |> 
                dplyr::arrange(ordering) |>
                dplyr::mutate(
                        Species_label = paste0(Species, " (", n_ASV, " ASVs)"),
                        Species_label = factor(Species_label, levels = unique(Species_label))
                ) |>
                tidyr::pivot_longer(
                        cols = c(mean_prev_unpollinated, mean_prev_pollinated),
                        names_to = "group",
                        values_to = "prevalence"
                ) |>
                dplyr::mutate(
                        group = dplyr::recode(
                                group,
                                mean_prev_unpollinated   = "unpollinated",
                                mean_prev_pollinated = "pollinated"
                        ),
                        group = factor(group, levels = c("pollinated", "unpollinated"))
                ) |>
                ggplot(aes(x = group, y = Species_label, fill = prevalence)) +
                geom_tile(color = "grey90") +
                scale_fill_viridis_c(
                        option = "magma",
                        limits = c(0, 1)
                ) +
                guides(
                        fill = guide_colorbar(
                                title.position = "top",
                                title.hjust = 0.5
                        )
                ) +
                labs(
                        x = NULL,
                        y = NULL,
                        fill = "Proportion of flowers\nwith taxon detected",
                        title = paste("Indicator Species for", indicator, 
                                      "(IndVal.g ≥ 0.25, permutation test p < 0.05)")
                ) +
                theme_minimal() +
                theme(panel.grid = element_blank(),
                      plot.title = element_text(hjust = 0)
                )
}

# Creating patches
p_unpollinated_prev2 <- make_indicator_pollination_heatmap(
        filtered_highIndVal_bagged_vs_unbagged_prev2_is_pollination_clsm,
        "unpollinated"
)

p_pollinated_prev2 <- make_indicator_pollination_heatmap(
        filtered_highIndVal_bagged_vs_unbagged_prev2_is_pollination_clsm,
        "pollinated"
)

# Assemble patches in one figure adjusting their hights to number of rows
n_unpollinated_prev2 <- n_distinct(
        filtered_highIndVal_bagged_vs_unbagged_prev2_is_pollination_clsm |>
                dplyr::filter(indicator_group == "unpollinated", !is.na(Species)) |>
                pull(Family)
)

n_pollinated_prev2 <- n_distinct(
        filtered_highIndVal_bagged_vs_unbagged_prev2_is_pollination_clsm |>
                dplyr::filter(indicator_group == "pollinated", !is.na(Species)) |>
                pull(Family)
)

# Collect patchwork
(p_unpollinated_prev2 / p_pollinated_prev2) +
        plot_layout(
                heights = c(n_unpollinated_prev2, n_pollinated_prev2),
                guides = "collect"
        ) &
        theme(legend.position = "right")


# ---- Indicator species analysis: shortlist of stigma colonizing fungi ----
# ---- Compare fungal_colonization_score group 0 vs 5 ----
# stat: IndVal = Specifity + Fidelity
# fcs stands for fungal_colonization_score
# The goal is to compare polar phenotypes: fungal_colonising_score == 0 Vs fungal_colonising_score == 5 to identify shortlist of taxa

table(sample_data(phyloseq_ITS1_clsm_hyphae)$fungal_colonization_score)
# Phenotypes are imbalanced (0: 120 obs, 5: 20) 
# Prevalence heatmap will depend on number of samples in each FCS group

# ---- Indicator species analysis for fungal colonization groups ----
set.seed(42)
indval_fcs_0vs5_result <-
        phyloseq_ITS1_clsm_hyphae |>
        subset_samples(fungal_colonization_score %in% c("0", "5")) |>
        (\(x) prune_taxa(taxa_sums(x) > 0, x))() |>
        (\(x)
         multipatt(
                 # Get ASV count matrix from phyloseq object
                 # Transpose count matrix to get: samples x ASVs
                 # Convert count matrix to presence/absence matrix
                 t(as(otu_table(x), "matrix")) > 0,
                 # Extract grouping variable values for the samples from the same phyloseq object
                 sample_data(x)$fungal_colonization_score,
                 func = "IndVal.g",
                 control = how(nperm = 999)
         )
        )()


summary(indval_fcs_0vs5_result)

class(indval_fcs_0vs5_result)
typeof(indval_fcs_0vs5_result)
names(indval_fcs_0vs5_result)

# Specificity (A): Is this feature specific to the group?
as.data.frame(indval_fcs_0vs5_result$A)["16ccfe6feedde0ada135825c37a5057a",]
# Fidelity(B): Is this feature typical for the group?
as.data.frame(indval_fcs_0vs5_result$B)["16ccfe6feedde0ada135825c37a5057a",]
# sign
indval_fcs_0vs5_result$sign["16ccfe6feedde0ada135825c37a5057a",]


# Extract IndVal results as dataframe
# ASV id | s.0 | s.5 | index | stat | p.value
# index: 1 ASV present only in group 1 (s.0); 2 ASV present only in group 2 (s.5); 3 ASV present in both groups
indval_fcs_0vs5_df <- indval_fcs_0vs5_result$sign


# Filter ISA results: ASVs significantly associated with score = 5
indicator_asvs_fcs_5 <- rownames(
        indval_fcs_0vs5_df[
                indval_fcs_0vs5_df$index == 2 &
                        !is.na(indval_fcs_0vs5_df$p.value) &
                        indval_fcs_0vs5_df$p.value <= 0.01,
        ]
)


# Map taxonomy of indicator species 
taxonomy_indicators_fcs_5 <- as.data.frame(taxonomy_table_ITS1)[indicator_asvs_fcs_5,]

sum(indicator_asvs_fcs_5 %in% rownames(as.data.frame(taxonomy_table_ITS1)))

View(taxonomy_indicators_fcs_5)

# ---- Diagnostic commands for ISA ----

taxonomy_indicators_fcs_0 <- as.data.frame(taxonomy_table_ITS1)[indicator_asvs_fcs_0,]



# ---- Visualize indicator Genera using heatmap and median RA (for imbalanced classes prevalence heatmaps are dodgy) ----
# Select indicator ASVs for FCS 0 and 5 groups
indicator_asvs_fcs_0vs5 <-
        rownames(
                indval_fcs_0vs5_result$sign[
                        !is.na(indval_fcs_0vs5_result$sign$p.value) &
                                indval_fcs_0vs5_result$sign$p.value <= 0.05,
                ]
        )

length(indicator_asvs_fcs_0vs5)

taxonomy_indicators_fcs_0vs5 <-
        as.data.frame(taxonomy_table_ITS1) |>
        tibble::rownames_to_column("ASV") |>
        filter(ASV %in% indicator_asvs_fcs_0vs5)

genus_asv_counts <-
        taxonomy_indicators_fcs_0vs5 |>
        filter(!is.na(Genus)) |>
        count(Genus, name = "n_asv")

# Build genus-level phyloseq object
phyloseq_fcs_0vs5_genus <-
        phyloseq_ITS1_clsm_hyphae |>
        subset_samples(fungal_colonization_score %in% c("0", "5")) |>
        (\(x) prune_taxa(indicator_asvs_fcs_0vs5, x))() |>
        tax_glom(taxrank = "Genus") |>
        (\(x) prune_taxa(taxa_sums(x) > 0, x))()

# Transform to relative abundance
phyloseq_fcs_0vs5_genus_ra <-
        transform_sample_counts(
                phyloseq_fcs_0vs5_genus,
                \(x) x / sum(x)
        )

# Extract matrices
genus_ra_matrix <-
        t(as(otu_table(phyloseq_fcs_0vs5_genus_ra), "matrix"))

metadata_fcs_0vs5 <-
        data.frame(sample_data(phyloseq_fcs_0vs5_genus_ra)) |>
        mutate(fcs = factor(fungal_colonization_score, levels = c("0", "5")))

# Median RA per genus x FCS 0/5
genus_median_ra <-
        genus_ra_matrix |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "sample_id") |>
        left_join(
                metadata_fcs_0vs5 |>
                        tibble::rownames_to_column(var = "sample_id") |>
                        dplyr::select(sample_id, fcs),
                by = "sample_id"
        ) |>
        dplyr::select(-sample_id) |>
        group_by(fcs) |>
        summarise(
                across(where(is.numeric), \(x) median(x, na.rm = TRUE)),
                .groups = "drop"
        ) |>
        pivot_longer(
                cols = -fcs,
                names_to = "taxon_id",
                values_to = "median_ra"
        ) |>
        left_join(
                as.data.frame(tax_table(phyloseq_fcs_0vs5_genus)) |>
                        tibble::rownames_to_column("taxon_id") |>
                        dplyr::select(taxon_id, Genus),
                by = "taxon_id"
        ) |>
        dplyr::select(-taxon_id)

# Matrix for heatmap
heatmap_df <-
        genus_median_ra |>
        left_join(genus_asv_counts, by = "Genus") |>
        mutate(
                log10_median_ra = log10(median_ra + 1e-4)
        )

heatmap_matrix <-
        heatmap_df |>
        dplyr::select(Genus, fcs, log10_median_ra) |>
        pivot_wider(
                names_from = fcs,
                values_from = log10_median_ra
        ) |>
        column_to_rownames("Genus") |>
        as.matrix()

# Order genera by number of ASVs
genus_order <-
        heatmap_df |>
        distinct(Genus, n_asv) |>
        arrange(desc(n_asv)) |>
        pull(Genus)

heatmap_matrix <- heatmap_matrix[genus_order, ]

# Plot heatmap
col_fun <- colorRamp2(
        c(-4, -2, 0),
        c("#fff7bc", "#feb24c", "#f03b20")
)

asv_numbers_per_genus <-
        genus_asv_counts$n_asv[
                match(rownames(heatmap_matrix), genus_asv_counts$Genus)
        ]

row_ha <- rowAnnotation(
        `ASVs` = anno_text(
                asv_numbers_per_genus,
                gp = gpar(fontsize = 9),
                just = "center"
        ),
        annotation_name_side = "top",
        width = unit(1.2, "cm")
)

Heatmap(
        heatmap_matrix,
        name = "log10(median RA)",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        left_annotation = row_ha,
        row_names_side = "left",
        column_title = "Fungal colonization score",
        row_title = "Fungal genera (ordered by ASV richness)",
        rect_gp = gpar(col = "grey60", lwd = 0.6),
        heatmap_legend_param = list(
                title = "log10(median RA)"
        )
)







# ---- Indicator species analysis: fungal colonization gradient (FCS 0 to 5) ----
# Goal:
# Identify ASVs whose occurrence patterns are significantly associated
# with specific levels of fungal colonization score (FCS),
# reflecting increasing intensity of fungal colonization observed by CLSM.

# Subset phyloseq object to samples with defined fungal colonization score
phyloseq_ITS1_indval_fcs_all <- phyloseq_ITS1_clsm_hyphae |>
        subset_samples(!is.na(fungal_colonization_score)) |>
        (\(x) prune_taxa(taxa_sums(x) > 0, x))()

# Extract ASV count matrix (samples x ASVs)
asv_count_matrix_fcs_all <- phyloseq_ITS1_indval_fcs_all |>
        otu_table() |>
        as("matrix") |>
        (\(x) if (taxa_are_rows(phyloseq_ITS1_indval_fcs_all)) t(x) else x)()

# Extract grouping variable (fungal colonization score, ordered factor)
fungal_colonization_group <- phyloseq_ITS1_indval_fcs_all |>
        sample_data() |>
        (\(x) x$fungal_colonization_score)() |>
        droplevels()

# Inspect group sizes (important for interpretation)
table(fungal_colonization_group)

# Convert ASV counts to presence / absence
asv_presence_absence_fcs_all <- (asv_count_matrix_fcs_all > 0) * 1

# Indicator species analysis (IndVal.g)
indval_fcs_all_result <- multipatt(
        asv_presence_absence_fcs_all,
        fungal_colonization_group,
        func = "IndVal.g",
        control = how(nperm = 999)
)

# Overview of results
summary(indval_fcs_all_result)

# Extract IndVal result table
# Rows = ASVs
# Columns:
#   s.0 ... s.5  -> association with each FCS level
#   index       -> group (or group combination) with highest association
#   stat        -> IndVal statistic
#   p.value     -> permutation-based significance
indval_fcs_all_df <- indval_fcs_all_result$sign

dim(indval_fcs_all_df)

# Optional: count how many ASVs are associated with each FCS level
colSums(indval_fcs_all_df[, grep("^s\\.", colnames(indval_fcs_all_df))], na.rm = TRUE)

# ASVs significantly associated with at least one FCS level
indicator_asvs_fcs_all <-
        rownames(
                indval_fcs_all_result$sign[
                        !is.na(indval_fcs_all_result$sign$p.value) &
                                indval_fcs_all_result$sign$p.value <= 0.05,
                ]
        )

length(indicator_asvs_fcs_all)

# Taxonomy for ASVs associated with at least one FCS level aggregated at g__
phyloseq_ITS1_fcs_indicators_genus <-
        phyloseq_ITS1_clsm_hyphae |>
        (\(x) prune_taxa(indicator_asvs_fcs_all, x))() |>
        tax_glom(taxrank = "Genus") |>
        (\(x) prune_taxa(taxa_sums(x) > 0, x))()

taxa_names(phyloseq_ITS1_fcs_indicators_genus) <-
        as.vector(
                tax_table(phyloseq_ITS1_fcs_indicators_genus)[, "Genus"]
        )

# Build presence/abscence matrix (sample x Genus)
genus_presence_absence_matrix <-
        t(as(otu_table(phyloseq_ITS1_fcs_indicators_genus), "matrix")) > 0

# Metadata
metadata_fcs_heatmap <-
        data.frame(sample_data(phyloseq_ITS1_fcs_indicators_genus))

table(metadata_fcs_heatmap$fungal_colonization_score)

# Calculate prevalence for Genus per FCS group
genus_prevalence_by_fcs <-
        genus_presence_absence_matrix |>
        as.data.frame() |>
        dplyr::mutate(
                fungal_colonization_score =
                        metadata_fcs_heatmap$fungal_colonization_score
        ) |>
        dplyr::group_by(fungal_colonization_score) |>
        dplyr::summarise(
                across(
                        where(is.logical),
                        \(x) mean(x, na.rm = TRUE)
                )
        ) |>
        dplyr::ungroup()

#
genus_order <-
        genus_prevalence_long |>
        dplyr::group_by(Genus) |>
        dplyr::summarise(
                peak_fcs =
                        fungal_colonization_score[
                                which.max(prevalence)
                        ]
        ) |>
        dplyr::arrange(as.numeric(peak_fcs)) |>
        dplyr::pull(Genus)

genus_prevalence_long$Genus <-
        factor(
                genus_prevalence_long$Genus,
                levels = genus_order
        )


# Long format for plotting
genus_prevalence_long <-
        genus_prevalence_by_fcs |>
        pivot_longer(
                cols = -fungal_colonization_score,
                names_to = "Genus",
                values_to = "prevalence"
        )

# Plot heatmap
ggplot(
        genus_prevalence_long,
        aes(
                x = fungal_colonization_score,
                y = Genus,
                fill = prevalence
        )
) +
        geom_tile(color = "grey80", linewidth = 0.2) +
        scale_fill_viridis(
                name = "Prevalence\n(fraction of samples)",
                limits = c(0, 1)
        ) +
        labs(
                title = "Prevalence of fungal genera across flower colonization stages",
                subtitle = "Genera selected by IndVal analysis (ITS1, presence/absence)",
                x = "Fungal colonization score (FCS)",
                y = "Fungal genus"
        ) +
        theme_bw() +
        theme(
                axis.text.y = element_text(size = 7),
                axis.text.x = element_text(size = 9),
                panel.grid = element_blank()
        )





################################################################################
################ Fungi of particular interest ##################################
################################################################################

# ---- Feature tracking -----
# ---- 1) One of the major cocoa pathogens -- s__Moniliophthora_pernicosa: only one sample pkk11 (agroforest, unbagged, pollinated, with callose plugs) ----
Moniliophthora_perniciosa_asv_id <- "10cd1dfa9011ad5bcc032fa70b1ed188"

# Check that ASV exists in the phyloseq object
Moniliophthora_perniciosa_asv_id %in% taxa_names(phyloseq_ITS1_fungi_biosamples)

# Extract counts of this ASV across samples
Moniliophthora_perniciosa_asv_counts <- phyloseq_ITS1_fungi_biosamples |>
        (\(x) prune_taxa(Moniliophthora_perniciosa_asv_id, x))() |>
        otu_table() |>
        as("matrix")

# Ensure orientation: samples as rows
if (taxa_are_rows(phyloseq_ITS1_fungi_biosamples)) {
        Moniliophthora_perniciosa_asv_counts <- t(Moniliophthora_perniciosa_asv_counts)
}

# Samples where ASV is present
samples_with_Moniliophthora_perniciosa_asv <- rownames(Moniliophthora_perniciosa_asv_counts)[Moniliophthora_perniciosa_asv_counts[, 1] > 0]

length(samples_with_Moniliophthora_perniciosa_asv)
samples_with_Moniliophthora_perniciosa_asv
# pkk11



# ---- 2) Ant parasite -- s__Simplicillium_formicae: ----
simplicillium_asvs <- c(
        "0bac2651c3034bb37f8d785ac3b61017",
        "6f3f29830dc35e0e5a8236c4ba4d42cf",
        "56a92fe8dabaab6c9cc24d6fc896e7c8",
        "37eab8bdaedbfdd0bbd3e42baf6d0080",
        "e1fd5b576122825cdef43fe4be6a9699",
        "b37ea8f13596a95e9ddc65baa0d9bc93",
        "77fb76baac4bf993f0c980d174dec0a2",
        "92392e250ec445c7d452b761998fc891",
        "1b282e01a989f0faa410bbe2f002b12c"
)

simplicillium_pa <-
        phyloseq_ITS1_fungi_biosamples_prev2 |>
        (\(x) prune_taxa(simplicillium_asvs, x))() |>
        (\(x) {
                mat <- as(otu_table(x), "matrix")
                if (!taxa_are_rows(x)) mat <- t(mat)
                mat > 0
        })()

# samples with S. formicae
samples_with_simplicillium <-
        colnames(simplicillium_pa)[
                colSums(simplicillium_pa) > 0
        ]

samples_with_simplicillium

# metadata of these samples
table(
        sample_data(phyloseq_ITS1_fungi_biosamples_prev2)[
                samples_with_simplicillium,
                "sample_type"
        ]
)

table(
        sample_data(phyloseq_ITS1_fungi_biosamples_prev2)[
                samples_with_simplicillium,
                "is_pollination_clsm"
        ]
)





################################################################################

# ---- HEATMAPS ---- 
# ---- Extract OTU matrix ----
otu_mat <- as(otu_table(phyloseq_ITS1_fungi_biosamples), "matrix")

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

# ---- Mean RA per group and genus ----
mean_ra_genus <- long_df |>
        group_by(sample_type, Genus) |>
        summarise(mean_ra = mean(rel_abund), .groups = "drop")

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
                title = "Taxonomic structure of fungil communities",
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

###########

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


######

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
phyloseq_rel <- transform_sample_counts(
        phyloseq_prev2_ab0,
        function(x) x / sum(x)
)

# Sum RA of endosymbionts per sample
endo_ra <- otu_table(phyloseq_rel)[endo_asvs, , drop = FALSE] |>
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
