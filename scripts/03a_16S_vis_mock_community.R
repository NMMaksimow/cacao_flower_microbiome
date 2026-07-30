# ---- Packages ----
library(here)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(scales)
library(ggtext)    # for element_markdown() — bold text in legend labels

# ---- Load intermediate objects (saved by 01a) ----
ps_16S_raw      <- readRDS(here("results", "rds", "ps_16S_raw.rds"))
ps_16S_decontam <- readRDS(here("results", "rds", "ps_16S_decontam.rds"))
ps_16S_bacteria <- readRDS(here("results", "rds", "ps_16S_bacteria.rds"))
# decontam_prev055_16S  <- readRDS(here("results", "rds", "decontam_prev055_16S.rds"))

# ---- Output directories ----
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

################################################################################
#### Mock community: expected vs observed (D6300 ZymoBIOMICS, 8 bacteria) ######
################################################################################

# Expected composition for D6300 ZymoBIOMICS DNA mixture: 8 bacterial genera, 12.5% each.
# Notes on SILVA taxonomy quirks:
#   - Salmonella is not annotated to genus level in SILVA → represented at Family level as g__NA.
#   - Escherichia coli appears as g__Escherichia-Shigella (SILVA complex genus).
expected_mock_16S <- tibble(
        Family_Genus = c(
                "f__Listeriaceae | g__Listeria",
                "f__Bacillaceae | g__Bacillus",
                "f__Staphylococcaceae | g__Staphylococcus",
                "f__Lactobacillaceae | g__Lactobacillus",
                "f__Enterococcaceae | g__Enterococcus",
                "f__Pseudomonadaceae | g__Pseudomonas",
                "f__Enterobacteriaceae | g__Escherichia-Shigella",
                "f__Enterobacteriaceae | g__Unassigned"   # Salmonella
        ),
        Abundance = rep(0.125, 8),
        Type      = "Expected",
        Sample    = "Expected"
)

# Agglomerate observed mock community reads to genus level, convert to relative abundance,
# then collapse any genus not present in the expected composition to "Other".
observed_mock_16S <- ps_16S_bacteria |>
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
                        Family_Genus %in% expected_mock_16S$Family_Genus,
                        Family_Genus, "Other"
                )
        ) |>
        group_by(Sample, Family_Genus, Type) |>
        summarise(Abundance = sum(Abundance), .groups = "drop")

# Factor levels: expected genera in the same order as the expected tibble, then Other.
fill_order <- c(expected_mock_16S$Family_Genus, "Other")

# ── Numeric x-axis layout ────────────────────────────────────────────────────
# Expected bar sits at x = 1, a visual gap occupies x = 2 (no bar drawn),
# observed samples start at x = 3.  This puts Expected and Observed in the
# same panel without faceting while keeping a clear visual separation.
obs_samples <- sort(unique(observed_mock_16S$Sample))
n_obs       <- length(obs_samples)

pos_expected <- 1
pos_gap      <- 2                                      # empty — no bar here
pos_observed <- seq(pos_gap + 1, length.out = n_obs)  # 3, 4, 5, ...

pos_map <- c(setNames(pos_expected, "Expected"),
             setNames(pos_observed, obs_samples))

# In the current pipeline sample_id is the phyloseq row name (= Sample).
sid_map_mock <- setNames(obs_samples, obs_samples)

mock_plot_data <- bind_rows(observed_mock_16S, expected_mock_16S) |>
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
# FIXME Replace XXX with the species name from the ZymoBIOMICS D6300 certificate.
fill_labels_16S <- c(
        "f__Listeriaceae | g__Listeria<br>*Listeria monocytogenes*",
        "f__Bacillaceae | g__Bacillus<br>*Bacillus subtilis*",
        "f__Staphylococcaceae | g__Staphylococcus<br>*Staphylococcus aureus*",
        "f__Lactobacillaceae | g__Lactobacillus<br>*Lactobacillus fermentum*",
        "f__Enterococcaceae | g__Enterococcus<br>*Enterococcus faecalis*",
        "f__Pseudomonadaceae | g__Pseudomonas<br>*Pseudomonas aeruginosa*",
        "f__Enterobacteriaceae | g__Escherichia-Shigella<br>*Escherichia coli*",
        "f__Enterobacteriaceae | g__Unassigned<br>*Salmonella enterica*",   # Salmonella
        "Other"
)

fig_mock_16S <- ggplot(mock_plot_data,
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
                labels = setNames(fill_labels_16S, fill_order)
        ) +
        labs(
                title = "16S mock community: expected vs observed (D6300 ZymoBIOMICS)\nafter decontamination and off-target sequences removal",
                y     = "Relative abundance",
                x     = NULL,
                fill  = "Family | Genus<br>*Expected*"
        ) +
        theme_bw(base_size = 10) +
        theme(axis.text.x          = element_text(angle = 45, hjust = 1),
              panel.grid.major.x   = element_blank(),
              legend.text          = element_markdown(lineheight = 1.4),
              legend.title = element_markdown(),
              legend.key.spacing.y = unit(5, "pt"))


################################################################################
# Save
################################################################################

ggsave(here("results", "figures", "03a_16S_mock_community_exp_obs.png"),
       fig_mock_16S, width = 9, height = 5, dpi = 300)