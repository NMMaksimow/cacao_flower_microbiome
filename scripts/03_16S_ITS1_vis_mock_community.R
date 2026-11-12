## ============================================================================
## 03 — Mock community: expected vs observed, 16S and ITS1
##
## ZymoBIOMICS Microbial Community DNA Standard D6305, added at the PCR step,
## one mock per sublibrary (14), after
## decontamination and the off-target filter. Replaces 03a (16S) and 03b (ITS1).
##
## Two figures, each with two expected bars ("Expected (DNA mass)": the
## passport DNA shares within the marker; "Expected (in silico)": amplicon
## templates from script 00) next to the 14 mocks:
##   16S             8 bacteria; mock ASVs assigned to members by exact sequence
##   ITS1 (v2 only), 2 fungi: merged-read OTU97 (main) and forward-read
##                   OTU97 (sensitivity) side by side. Saccharomyces is absent
##                   from the merged reads (its ITS1 insert is too long for the
##                   read pairs to overlap) but present in R1.
## Tables: for each fungal mock member, the share of mock reads assigned at
## order, family, genus and species, merged vs R1 (v2); expected vs observed
## per member.
##
## Input (01a, 01b, 00):
##   results/rds/ps_16S_bacteria.rds
##   results/rds/ps_ITS1_otu97_fungi.rds
##   results/rds/ps_ITS1_r1_otu97_fungi.rds          (v2 only)
##   results/tables/00_mock_expected_composition.csv, 00_mock_expected_amplicons.csv
## Output:
##   results/figures/03_16S_mock_community_exp_insilico_obs.png
##   results/figures/03_ITS1_otu97_merged_vs_r1_mock_community_exp_insilico_obs.png   (v2)
##   results/tables/03_ITS1_otu97_merged_vs_r1_mock_by_rank.csv                      (v2)
##   results/tables/03_mock_expected_insilico_vs_observed.csv
## ============================================================================

# ---- 0. Packages ---------------------------------------------------------
library(here)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(scales)
library(ggtext)    # element_markdown(): italic species names in legend labels
library(Biostrings) # readDNAStringSet(): ASV sequences, matched to the in-silico amplicons

# Data version switch: CFM_VERSION, RESULTS_DIR
source(here("scripts", "utils_paths.R"))

dir.create(file.path(RESULTS_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RESULTS_DIR, "tables"),  showWarnings = FALSE, recursive = TRUE)

# ---- 1. Expected compositions --------------------------------------------
# 16S: 8 bacterial genera, 12.5% each. Notes on SILVA taxonomy quirks:
#   - Salmonella is not annotated to genus level in SILVA: family level, g__Unassigned.
#   - Escherichia coli appears as g__Escherichia-Shigella (SILVA complex genus).
#   - SILVA 144 (v2): Lactobacillus fermentum is Limosilactobacillus, and Enterococcus
#     moved to Streptococcaceae. (Staphylococcus and Enterococcus also carry the
#     SILVA 144 mark "--other", dropped at import in utils_read_taxonomy.R.)
#   - v2 assigns 2 mock ASVs to g__Leclercia-Silvania-Enterobacter (~15k reads). One of
#     them (3e86e962...) is identical to the 7th 16S copy of Salmonella enterica in the
#     Zymo genome (script 00); it is "Other" in the genus-labelled figures below and
#     Salmonella in the sequence-matched in-silico figures (section 5).
expected_16S <- tibble(
        Family_Genus = list(
                v1 = c("f__Listeriaceae | g__Listeria", "f__Bacillaceae | g__Bacillus",
                       "f__Staphylococcaceae | g__Staphylococcus",
                       "f__Lactobacillaceae | g__Lactobacillus",
                       "f__Enterococcaceae | g__Enterococcus",
                       "f__Pseudomonadaceae | g__Pseudomonas",
                       "f__Enterobacteriaceae | g__Escherichia-Shigella",
                       "f__Enterobacteriaceae | g__Unassigned"),          # Salmonella
                v2 = c("f__Listeriaceae | g__Listeria", "f__Bacillaceae | g__Bacillus",
                       "f__Staphylococcaceae | g__Staphylococcus",
                       "f__Lactobacillaceae | g__Limosilactobacillus",
                       "f__Streptococcaceae | g__Enterococcus",
                       "f__Pseudomonadaceae | g__Pseudomonas",
                       "f__Enterobacteriaceae | g__Escherichia-Shigella",
                       "f__Enterobacteriaceae | g__Unassigned")           # Salmonella
        )[[CFM_VERSION]],
        species   = c("Listeria monocytogenes", "Bacillus subtilis", "Staphylococcus aureus",
                      "Lactobacillus fermentum", "Enterococcus faecalis",
                      "Pseudomonas aeruginosa", "Escherichia coli", "Salmonella enterica"),
        Abundance = 0.125
)

# ITS1: 2 fungal genera, 50% each among fungi.
expected_ITS1 <- tibble(
        Family_Genus = c("f__Cryptococcaceae | g__Cryptococcus",
                         "f__Saccharomycetaceae | g__Saccharomyces"),
        species      = c("Cryptococcus neoformans", "Saccharomyces cerevisiae"),
        Abundance    = 0.5
)

# ---- 2. Observed compositions --------------------------------------------
# Mock samples only; genus-level sums as relative abundance; any genus that is
# not an expected member becomes "Other".
observed_mock <- function(ps, expected) {
        ps |>
                subset_samples(sample_type == "mock_community") |>
                (\(x) prune_taxa(taxa_sums(x) > 0, x))() |>
                tax_glom(taxrank = "Genus", NArm = FALSE) |>
                transform_sample_counts(\(x) x / sum(x)) |>
                psmelt() |>
                filter(Abundance > 0) |>
                mutate(Family_Genus = paste(coalesce(Family, "f__Unassigned"),
                                            coalesce(Genus,  "g__Unassigned"), sep = " | "),
                       Family_Genus = if_else(Family_Genus %in% expected$Family_Genus,
                                              Family_Genus, "Other")) |>
                summarise(Abundance = sum(Abundance), .by = c(Sample, Family_Genus))
}

# ---- 3. Plot -------------------------------------------------------------
# Stacked bars: the expected bar(s) first (x = 1, 2, ...), then a gap marked by
# a dashed line, then the mock samples. One panel, no facets.
#   expected      members, species and expected shares; sets the legend
#   expected_bars named list of expected tables drawn as bars (default: one bar,
#                 "Expected"); the names become the x labels
plot_mock <- function(observed, expected, title, expected_bars = list(Expected = expected)) {
        samples    <- sort(unique(observed$Sample))
        n_expected <- length(expected_bars)
        fill_order <- c(expected$Family_Genus, "Other")
        positions  <- setNames(c(seq_len(n_expected), seq(n_expected + 2, length.out = length(samples))),
                               c(names(expected_bars), samples))

        # Legend labels: Family | Genus on line 1, expected species (italic) on line 2
        fill_labels <- setNames(c(paste0(expected$Family_Genus, "<br>*", expected$species, "*"),
                                  "Other"), fill_order)
        # Member colours: RColorBrewer "Paired" (needs at least 3); Other in light grey
        fill_values <- c(setNames(RColorBrewer::brewer.pal(max(3, nrow(expected)), "Paired")[seq_len(nrow(expected))],
                                  expected$Family_Genus),
                         Other = "grey80")

        expected_long <- imap(expected_bars, \(e, bar) e |> select(Family_Genus, Abundance) |> mutate(Sample = bar)) |>
                list_rbind()

        bind_rows(observed, expected_long) |>
                mutate(Family_Genus = factor(Family_Genus, levels = fill_order),
                       x_pos        = positions[Sample]) |>
                ggplot(aes(x = x_pos, y = Abundance, fill = Family_Genus)) +
                geom_col(colour = "white", linewidth = 0.3, width = 0.85) +   # bar outline and width
                geom_vline(xintercept = n_expected + 1, linetype = "dashed", # Expected | Observed separator
                           colour = "gray60", linewidth = 0.4) +
                scale_x_continuous(breaks = positions, labels = names(positions),
                                   expand = expansion(add = 0.6)) +
                scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.03))) +
                scale_fill_manual(values = fill_values, labels = fill_labels, drop = FALSE) +
                labs(title = title, y = "Relative abundance", x = NULL,
                     fill = "Family | Genus<br>*Expected*") +
                theme_bw(base_size = 10) +                                    # base text size
                theme(axis.text.x          = element_text(angle = 45, hjust = 1),
                      panel.grid.major.x   = element_blank(),
                      legend.text          = element_markdown(lineheight = 1.4),
                      legend.title         = element_markdown(),
                      legend.key.spacing.y = unit(5, "pt"))                   # space between legend keys
}

SUBTITLE <- "after decontamination and off-target sequences removal"

ps_16S       <- readRDS(file.path(RESULTS_DIR, "rds", "ps_16S_bacteria.rds"))
ps_ITS1      <- readRDS(file.path(RESULTS_DIR, "rds", "ps_ITS1_otu97_fungi.rds"))
observed_ITS1 <- observed_mock(ps_ITS1, expected_ITS1)

# The figures (section 5) show two expected bars: by DNA mass (the tables above)
# and by in-silico PCR (script 00). The one-bar versions were dropped 2026-11-12.

# ---- 4. ITS1 merged vs R1-only: assignment by rank (v2) ------------------
if (CFM_VERSION == "v2") {
        ps_ITS1_r1 <- readRDS(file.path(RESULTS_DIR, "rds", "ps_ITS1_r1_otu97_fungi.rds"))

        # Share of all mock reads assigned to each member's name at each rank
        members <- tibble(
                member  = c("Cryptococcus neoformans", "Saccharomyces cerevisiae"),
                Order   = c("o__Tremellales", "o__Saccharomycetales"),
                Family  = c("f__Cryptococcaceae", "f__Saccharomycetaceae"),
                Genus   = c("g__Cryptococcus", "g__Saccharomyces"),
                Species = c("s__Cryptococcus_neoformans", "s__Saccharomyces_cerevisiae")
        ) |>
                pivot_longer(-member, names_to = "rank", values_to = "label")

        mock_by_rank <- list(merged = ps_ITS1, r1 = ps_ITS1_r1) |>
                imap(\(ps, dataset) {
                        mock        <- subset_samples(ps, sample_type == "mock_community")
                        taxon_reads <- taxa_sums(mock)
                        tax         <- as.data.frame(tax_table(mock))
                        members |>
                                mutate(dataset = dataset,
                                       reads   = map2_dbl(rank, label, \(r, l) sum(taxon_reads[tax[[r]] %in% l])),
                                       share_of_mock_reads = reads / sum(taxon_reads))   # all mock reads
                }) |>
                list_rbind() |>
                mutate(rank = factor(rank, levels = c("Order", "Family", "Genus", "Species"))) |>
                arrange(member, rank, dataset) |>
                select(dataset, member, rank, label, reads, share_of_mock_reads)

        write_csv(mock_by_rank, file.path(RESULTS_DIR, "tables", "03_ITS1_otu97_merged_vs_r1_mock_by_rank.csv"))
        print(mock_by_rank, n = Inf)
}

# ---- 5. Expected from in-silico PCR (script 00) -------------------------
# Same figures with a second expected bar: "Expected (DNA mass)" is the bar
# above (the passport DNA shares within the marker), "Expected (in silico)" the
# share of amplicon templates computed by 00_16S_ITS1_mock_expected_insilico_pcr.R
# (DNA mass / genome size x amplifiable rRNA copies; ITS1 merged reads only
# for inserts short enough to merge). The figures above stay as they are.
EXPECTED_MASS     <- "Expected (DNA mass)"
EXPECTED_INSILICO <- "Expected (in silico)"

insilico_composition <- read_csv(file.path(RESULTS_DIR, "tables", "00_mock_expected_composition.csv"),
                                 show_col_types = FALSE) |>
        filter(basis == "default")   # bacteria: archive genomes; yeasts: Zymo manual Table 2
insilico_amplicons <- read_csv(file.path(RESULTS_DIR, "tables", "00_mock_expected_amplicons.csv"),
                               show_col_types = FALSE) |>
        filter(!off_target)

# Expected (in silico) table for one marker and read set, in the members' order
expected_insilico <- function(marker_name, read_set_name, expected) {
        e <- insilico_composition |>
                filter(marker == marker_name, read_set == read_set_name) |>
                transmute(Family_Genus = .data[[paste0("label_", CFM_VERSION)]], Abundance = expected_share)
        stopifnot(setequal(e$Family_Genus, expected$Family_Genus))
        e[match(expected$Family_Genus, e$Family_Genus), ]
}

# 16S: assign each mock ASV to a member by its exact sequence where it equals an
# in-silico amplicon. One Salmonella 16S copy differs from the others by 1 nt and
# SILVA 144 calls it Leclercia-Silvania-Enterobacter; by sequence it is Salmonella.
# Only the Family and Genus of matched ASVs change; everything else as before.
relabel_by_sequence <- function(ps, marker_name) {
        repseqs  <- readDNAStringSet(file.path(EXPORT_DIR, "CFM_16S_dada2_repseqs", "dna-sequences.fasta"))
        seq_of   <- setNames(as.character(repseqs), names(repseqs))
        amp      <- filter(insilico_amplicons, marker == marker_name)
        member   <- amp$member[match(seq_of[taxa_names(ps)], amp$insert)]
        label    <- insilico_composition$label_v2[match(member, insilico_composition$member)]
        if (CFM_VERSION == "v1") label <- insilico_composition$label_v1[match(member, insilico_composition$member)]
        matched  <- which(!is.na(label))
        tt <- tax_table(ps)
        tt[matched, "Family"] <- str_extract(label[matched], "^f__[^ ]+")
        tt[matched, "Genus"]  <- str_extract(label[matched], "g__\\S+$")
        tax_table(ps) <- tt
        cat(sprintf("16S: %d ASVs matched to a member by sequence\n", length(matched)))
        ps
}

observed_16S_seq <- observed_mock(relabel_by_sequence(ps_16S, "16S"), expected_16S)

fig_mock_16S_insilico <- plot_mock(
        observed_16S_seq, expected_16S,
        paste0("16S mock community: expected (DNA mass, in silico) vs observed (ZymoBIOMICS D6305)\n", SUBTITLE,
               "; mock ASVs assigned to members by sequence"),
        expected_bars = setNames(list(expected_16S, expected_insilico("16S", "merged", expected_16S)),
                                 c(EXPECTED_MASS, EXPECTED_INSILICO)))

ggsave(file.path(RESULTS_DIR, "figures", "03_16S_mock_community_exp_insilico_obs.png"),
       plot = fig_mock_16S_insilico, width = 9.5, height = 5, dpi = 300)

# Per member: expected (DNA mass, in silico) and observed (mean and median share
# across the mocks); log2 ratio observed / in silico (0 = as expected)
compare_expected <- function(observed, expected, marker_name, read_set_name) {
        insilico <- expected_insilico(marker_name, read_set_name, expected)
        observed |>
                filter(Family_Genus != "Other") |>
                complete(Sample, Family_Genus = expected$Family_Genus, fill = list(Abundance = 0)) |>
                summarise(observed_mean = mean(Abundance), observed_median = median(Abundance), .by = Family_Genus) |>
                right_join(tibble(Family_Genus = expected$Family_Genus, species = expected$species,
                                  expected_dna_mass = expected$Abundance, expected_insilico = insilico$Abundance),
                           by = "Family_Genus") |>
                mutate(marker = marker_name, read_set = read_set_name,
                       log2_observed_vs_insilico = log2(observed_mean / expected_insilico), .before = 1)
}

comparison <- list(compare_expected(observed_16S_seq, expected_16S, "16S", "merged"),
                   compare_expected(observed_ITS1, expected_ITS1, "ITS1", "merged"))

# ITS1 merged vs R1-only (v2): each panel with the in-silico bar of its read set
if (CFM_VERSION == "v2") {
        observed_ITS1_r1 <- observed_mock(ps_ITS1_r1, expected_ITS1)
        fig_mock_ITS1_r1_insilico <- (
                plot_mock(observed_ITS1, expected_ITS1, "Merged R1+R2, OTU97 (main)",
                          expected_bars = setNames(list(expected_ITS1, expected_insilico("ITS1", "merged", expected_ITS1)),
                                                   c(EXPECTED_MASS, EXPECTED_INSILICO))) |
                plot_mock(observed_ITS1_r1, expected_ITS1, "R1 only, OTU97 (sensitivity)",
                          expected_bars = setNames(list(expected_ITS1, expected_insilico("ITS1", "R1_only", expected_ITS1)),
                                                   c(EXPECTED_MASS, EXPECTED_INSILICO)))) +
                plot_layout(guides = "collect") +
                plot_annotation(title = paste0("ITS1 mock community, merged vs forward reads only, with in-silico expectation\n", SUBTITLE),
                                theme = theme(plot.title = element_text(size = 11)))
        ggsave(file.path(RESULTS_DIR, "figures", "03_ITS1_otu97_merged_vs_r1_mock_community_exp_insilico_obs.png"),
               plot = fig_mock_ITS1_r1_insilico, width = 15, height = 5, dpi = 300)
        comparison <- c(comparison, list(compare_expected(observed_ITS1_r1, expected_ITS1, "ITS1", "R1_only")))
}

comparison <- list_rbind(comparison)
write_csv(comparison, file.path(RESULTS_DIR, "tables", "03_mock_expected_insilico_vs_observed.csv"))
print(comparison |> mutate(across(where(is.numeric), \(x) round(x, 3))) |>
              select(marker, read_set, species, expected_dna_mass, expected_insilico, observed_mean,
                     log2_observed_vs_insilico), n = Inf, width = Inf)

cat("03 done.\n")
