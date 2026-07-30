# ---- Packages ----
library(here)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(scales)
library(ggtext)    # for element_markdown() — bold text in legend labels

# ---- Load intermediate objects (saved by 01b) ----
ps_ITS1_otu97_fungi    <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi.rds"))

# ---- Output directories ----
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

################################################################################
#### Mock community: expected vs observed (D6300 ZymoBIOMICS, 2 fungi) ########
################################################################################

# Expected composition for D6300 ZymoBIOMICS DNA mixture (fungal fraction only):
# 2 fungal genera, each 2% of total DNA → rescaled to 50% each among fungi,
# since bacterial DNA is not expected to amplify with ITS1 primers.
# Notes:
#   - Cryptococcus neoformans (Cryptococcaceae in UNITE) — may be detected.
#   - Saccharomyces cerevisiae — NOT detectable: ITS1 amplicon is too long for
#     paired-end read merging, illustrating a known limitation of ITS1 metabarcoding.
expected_mock_ITS1 <- tibble(
        Family_Genus = c(
                "f__Cryptococcaceae | g__Cryptococcus",
                "f__Saccharomycetaceae | g__Saccharomyces"  # not detectable: ITS1 too long
        ),
        Abundance = c(0.5, 0.5),
        Type      = "Expected",
        Sample    = "Expected"
)

# Agglomerate observed mock community reads to genus level, convert to relative abundance,
# then collapse any genus not present in the expected composition to "Other".
# Family_Genus strings from psmelt() are matched directly against expected labels
# (family names verified from tax_table(ps_ITS1_fungi)).
observed_mock_ITS1 <- ps_ITS1_otu97_fungi |>
        subset_samples(sample_type == "mock_community") |>
        (\(ps) prune_taxa(taxa_sums(ps) > 0, ps))() |>
        tax_glom(taxrank = "Genus", NArm = FALSE) |>
        transform_sample_counts(\(x) x / sum(x)) |>
        psmelt() |>
        filter(Abundance > 0) |>
        mutate(
                Family_Genus = paste(
                        if_else(is.na(Family), "f__Unassigned", Family),
                        if_else(is.na(Genus),  "g__Unassigned", Genus),
                        sep = " | "
                ),
                Type = "Observed"
        ) |>
        mutate(
                Family_Genus = if_else(
                        Family_Genus %in% expected_mock_ITS1$Family_Genus,
                        Family_Genus, "Other"
                )
        ) |>
        group_by(Sample, Family_Genus, Type) |>
        summarise(Abundance = sum(Abundance), .groups = "drop")

# Factor levels: expected genera in the same order as the expected tibble, then Other.
fill_order <- c(expected_mock_ITS1$Family_Genus, "Other")

# ── Numeric x-axis layout ────────────────────────────────────────────────────
# Expected bar sits at x = 1, a visual gap occupies x = 2 (no bar drawn),
# observed samples start at x = 3.  This puts Expected and Observed in the
# same panel without faceting while keeping a clear visual separation.
obs_samples <- sort(unique(observed_mock_ITS1$Sample))
n_obs       <- length(obs_samples)

pos_expected <- 1
pos_gap      <- 2                                      # empty — no bar here
pos_observed <- seq(pos_gap + 1, length.out = n_obs)  # 3, 4, 5, ...

pos_map <- c(setNames(pos_expected, "Expected"),
             setNames(pos_observed, obs_samples))

# sample_id is the phyloseq rowname (= Sample after psmelt rownames_to_column)
sid_map_mock <- setNames(obs_samples, obs_samples)

mock_plot_data <- bind_rows(observed_mock_ITS1, expected_mock_ITS1) |>
        mutate(
                Family_Genus = factor(Family_Genus, levels = fill_order),
                x_pos        = pos_map[Sample]
        )

# x-axis break positions and labels — use full sample_id for observed
obs_labels   <- vapply(obs_samples,
                       function(s) { l <- sid_map_mock[[s]]; if (is.na(l) || l == "") s else l },
                       character(1))
x_breaks <- c(pos_expected, pos_observed)
x_labels <- c("Expected", obs_labels)

# Legend labels: Family | Genus on line 1, expected organism (bold) on line 2.
# Replace XXX with the species name from the ZymoBIOMICS D6300 certificate.
# Note: Saccharomyces was not detectable — ITS1 amplicon too long for paired-end merging.
fill_labels_ITS1 <- c(
        "f__Cryptococcaceae | g__Cryptococcus<br>*Cryptococcus neoformans*",
        "f__Saccharomycetaceae | g__Saccharomyces<br>*Saccharomyces cerevisiae*",
        "Other"
)

fig_mock_ITS1 <- ggplot(mock_plot_data,
                         aes(x = x_pos, y = Abundance, fill = Family_Genus)) +
        geom_col(colour = "white", linewidth = 0.3, width = 0.85) +
        # Dashed vertical separator between Expected and first Observed bar
        geom_vline(xintercept = pos_gap, linetype = "dashed",
                   colour = "gray60", linewidth = 0.4) +
        scale_x_continuous(
                breaks = x_breaks,
                labels = x_labels,
                expand = expansion(add = 0.6)
        ) +
        scale_y_continuous(labels = percent_format(),
                           expand = expansion(mult = c(0, 0.03))) +
        scale_fill_manual(
                values = c(
                        setNames(
                                RColorBrewer::brewer.pal(length(fill_order) - 1, "Paired"),
                                fill_order[fill_order != "Other"]
                        ),
                        "Other" = "grey80"
                ),
                labels = setNames(fill_labels_ITS1, fill_order)
        ) +
        labs(
                title = "ITS1 mock community: expected vs observed (D6300 ZymoBIOMICS)\nafter decontamination and off-target sequences removal",
                y     = "Relative abundance",
                x     = NULL,
                fill  = "Family | Genus<br>*Expected*"
        ) +
        theme_bw(base_size = 10) +
        theme(axis.text.x          = element_text(angle = 45, hjust = 1),
              panel.grid.major.x   = element_blank(),
              legend.text          = element_markdown(lineheight = 1.4),
              legend.title = element_markdown(),
              # increase spacing between legend keys
              legend.key.spacing.y = unit(5, "pt"))

################################################################################
# Save
################################################################################

ggsave(here("results", "figures", "03b_ITS1_otu97_mock_community_exp_obs.png"),
       fig_mock_ITS1, width = 8, height = 5, dpi = 300)