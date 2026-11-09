## ============================================================================
## 03 — Mock community: expected vs observed, 16S and ITS1
##
## ZymoBIOMICS D6300 DNA mixture, one mock per sublibrary (14), after
## decontamination and the off-target filter. Replaces 03a (16S) and 03b (ITS1).
##
## Three figures:
##   16S             8 bacteria, 12.5% each
##   ITS1            2 fungi, rescaled to 50% each among fungi (each is 2% of
##                   the total DNA; bacterial DNA does not amplify with ITS1 primers)
##   ITS1 merged vs R1-only (v2 only): the same ITS1 figure for the merged-read
##                   OTU97 data (main) and the forward-read-only OTU97 data
##                   (sensitivity branch), side by side. Saccharomyces is absent
##                   from the merged reads (its ITS1 amplicon is too long for the
##                   read pairs to overlap) but present in R1.
## And one table (v2 only): for each fungal mock member, the share of mock reads
## assigned at order, family, genus and species, merged vs R1.
##
## Input (01a, 01b):
##   results/rds/ps_16S_bacteria.rds
##   results/rds/ps_ITS1_otu97_fungi.rds
##   results/rds/ps_ITS1_r1_otu97_fungi.rds          (v2 only)
## Output:
##   results/figures/03_16S_mock_community_exp_obs.png
##   results/figures/03_ITS1_otu97_mock_community_exp_obs.png
##   results/figures/03_ITS1_otu97_merged_vs_r1_mock_community_exp_obs.png   (v2)
##   results/tables/03_ITS1_otu97_merged_vs_r1_mock_by_rank.csv              (v2)
## ============================================================================

# ---- 0. Packages ---------------------------------------------------------
library(here)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(scales)
library(ggtext)    # element_markdown(): italic species names in legend labels

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
#   - v2 also assigns 2 mock ASVs to g__Leclercia-Silvania-Enterobacter
#     (s__Enterobacter_ludwigii, ~15k reads); not a D6300 member, so it falls in "Other".
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
# Stacked bars: the Expected bar at x = 1, a gap at x = 2 marked by a dashed
# line, the mock samples from x = 3. One panel, no facets.
plot_mock <- function(observed, expected, title) {
        samples    <- sort(unique(observed$Sample))
        fill_order <- c(expected$Family_Genus, "Other")
        positions  <- setNames(c(1, seq(3, length.out = length(samples))), c("Expected", samples))

        # Legend labels: Family | Genus on line 1, expected species (italic) on line 2
        fill_labels <- setNames(c(paste0(expected$Family_Genus, "<br>*", expected$species, "*"),
                                  "Other"), fill_order)
        # Member colours: RColorBrewer "Paired" (needs at least 3); Other in light grey
        fill_values <- c(setNames(RColorBrewer::brewer.pal(max(3, nrow(expected)), "Paired")[seq_len(nrow(expected))],
                                  expected$Family_Genus),
                         Other = "grey80")

        bind_rows(observed, expected |> select(Family_Genus, Abundance) |> mutate(Sample = "Expected")) |>
                mutate(Family_Genus = factor(Family_Genus, levels = fill_order),
                       x_pos        = positions[Sample]) |>
                ggplot(aes(x = x_pos, y = Abundance, fill = Family_Genus)) +
                geom_col(colour = "white", linewidth = 0.3, width = 0.85) +   # bar outline and width
                geom_vline(xintercept = 2, linetype = "dashed",              # Expected | Observed separator
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

fig_mock_16S <- plot_mock(observed_mock(ps_16S, expected_16S), expected_16S,
                          paste0("16S mock community: expected vs observed (D6300 ZymoBIOMICS)\n", SUBTITLE))
fig_mock_ITS1 <- plot_mock(observed_ITS1, expected_ITS1,
                           paste0("ITS1 mock community: expected vs observed (D6300 ZymoBIOMICS)\n", SUBTITLE))

ggsave(file.path(RESULTS_DIR, "figures", "03_16S_mock_community_exp_obs.png"),
       plot = fig_mock_16S, width = 9, height = 5, dpi = 300)
ggsave(file.path(RESULTS_DIR, "figures", "03_ITS1_otu97_mock_community_exp_obs.png"),
       plot = fig_mock_ITS1, width = 8, height = 5, dpi = 300)

# ---- 4. ITS1 merged vs R1-only (v2) --------------------------------------
if (CFM_VERSION == "v2") {
        ps_ITS1_r1 <- readRDS(file.path(RESULTS_DIR, "rds", "ps_ITS1_r1_otu97_fungi.rds"))

        # Two panels with one shared legend; panel titles name the read set
        fig_mock_ITS1_r1 <- (plot_mock(observed_ITS1, expected_ITS1, "Merged R1+R2, OTU97 (main)") |
                             plot_mock(observed_mock(ps_ITS1_r1, expected_ITS1), expected_ITS1,
                                       "R1 only, OTU97 (sensitivity)")) +
                plot_layout(guides = "collect") +
                plot_annotation(title = paste0("ITS1 mock community, merged vs forward reads only\n", SUBTITLE),
                                theme = theme(plot.title = element_text(size = 11)))

        ggsave(file.path(RESULTS_DIR, "figures", "03_ITS1_otu97_merged_vs_r1_mock_community_exp_obs.png"),
               plot = fig_mock_ITS1_r1, width = 14, height = 5, dpi = 300)   # width: two panels

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

cat("03 done.\n")
