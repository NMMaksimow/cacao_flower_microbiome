# ---- Packages ----
library(here)
library(tidyverse)
library(phyloseq)
library(scales)
library(ggtext)    # element_markdown() for bold group headers in legend

# ---- Load intermediate objects (saved by 01a) ----
ps_16S_raw      <- readRDS(here("results", "rds", "ps_16S_raw.rds"))
ps_16S_decontam <- readRDS(here("results", "rds", "ps_16S_decontam.rds"))
ps_16S_bacteria <- readRDS(here("results", "rds", "ps_16S_bacteria.rds"))

# ---- Output directories ----
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

################################################################################
##########  Effect of decontamination on control samples  ######################
################################################################################

# ── Expected mock community taxa (Family | Genus, SILVA taxonomy) ─────────────
# Used to assign taxa to the "Mock community" legend group.
# Notes: Salmonella → g__Unassigned; E. coli → g__Escherichia-Shigella in SILVA.
mock_expected_16S <- c(
        "f__Listeriaceae | g__Listeria",
        "f__Bacillaceae | g__Bacillus",
        "f__Staphylococcaceae | g__Staphylococcus",
        "f__Lactobacillaceae | g__Lactobacillus",
        "f__Enterococcaceae | g__Enterococcus",
        "f__Pseudomonadaceae | g__Pseudomonas",
        "f__Enterobacteriaceae | g__Escherichia-Shigella",
        "f__Enterobacteriaceae | g__Unassigned"   # Salmonella
)

# ── Long-table helper ────────────────────────────────────────────────────────
# Melts a phyloseq object (control samples only) to a tidy table.
# Species level is omitted: 16S cannot reliably resolve species, and
# Family | Genus colour coding is sufficient to show contamination origin.
make_control_long_table <- function(physeq, label) {
        physeq |>
                subset_samples(sample_type %in% c("mock_community",
                                                  "negative_pcr_control",
                                                  "extraction_control")) |>
                (\(x) prune_taxa(taxa_sums(x) > 0, x))() |>
                psmelt() |>
                mutate(
                        stage      = label,
                        Family     = if_else(is.na(Family), "f__Unassigned", Family),
                        Genus      = if_else(is.na(Genus),  "g__Unassigned", Genus),
                        Family_Genus = paste(Family, Genus, sep = " | ")
                )
}

control_long_tbl <- bind_rows(
        make_control_long_table(ps_16S_raw,      "Before decontam"),
        make_control_long_table(ps_16S_decontam, "After decontam"),
        make_control_long_table(ps_16S_bacteria, "After off-target removal")
) |>
        mutate(stage = factor(stage, levels = c("Before decontam",
                                                "After decontam",
                                                "After off-target removal")))

# ── Union threshold: taxa with mean within-sample RA > 2 % in any condition ──
# IMPORTANT: OTUs must be aggregated to Family_Genus level BEFORE computing RA.
# Without this step, a genus split across 5 OTUs (each 2 %) would show mean
# RA of 2 % per OTU rather than 10 % total, failing the threshold unfairly
# while single-OTU contaminants at the same total RA pass.
# Threshold set at 2 % for 16S (1 % produces too many hard-to-see labels).
RA_THRESHOLD <- 0.05

# Step 1: aggregate OTUs → Family_Genus per sample (sum reads first)
control_long_agg <- control_long_tbl |>
        group_by(stage, Sample, sample_type, Family_Genus) |>
        summarise(Abundance = sum(Abundance), .groups = "drop")

# Step 2: within-sample RA on aggregated data
control_long_ra <- control_long_agg |>
        group_by(stage, Sample) |>
        mutate(sample_total = sum(Abundance),
               RA = if_else(sample_total > 0, Abundance / sample_total, 0)) |>
        ungroup()

# Step 3: union criterion — max mean RA across all (stage × sample_type) cells
union_taxa <- control_long_ra |>
        group_by(stage, sample_type, Family_Genus) |>
        summarise(mean_RA = mean(RA), .groups = "drop") |>
        group_by(Family_Genus) |>
        summarise(max_mean_RA = max(mean_RA), .groups = "drop") |>
        filter(max_mean_RA > RA_THRESHOLD) |>
        arrange(desc(max_mean_RA)) |>
        pull(Family_Genus)

message(sprintf("Union taxa (RA > %.0f%% in at least one condition): %d",
                100 * RA_THRESHOLD, length(union_taxa)))
print(union_taxa)

# ── Assign taxa to legend groups ──────────────────────────────────────────────
mock_in_union   <- intersect(mock_expected_16S, union_taxa)
contam_in_union <- setdiff(union_taxa, mock_expected_16S)

# ── Shared colour palette ─────────────────────────────────────────────────────
n_mock   <- length(mock_in_union)
n_contam <- length(contam_in_union)
n_total  <- n_mock + n_contam

if (n_total <= 12) {
        all_colors <- RColorBrewer::brewer.pal(max(n_total, 3), "Paired")
} else {
        all_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_total)
}
mock_colors   <- all_colors[seq_len(n_mock)]
contam_colors <- all_colors[seq_len(n_contam) + n_mock]

# ── Two-tier legend: header entries as separate factor levels ─────────────────
# Header strings are added as factor levels but never appear as data values →
# they show up in the legend only. Their key boxes are hidden via override.aes.
HEADER_MOCK   <- "\u25b8 Mock community"       # ▸ Mock community
HEADER_CONTAM <- "\u25b8 Lab contamination"    # ▸ Lab contamination

# Full ordered level sequence: header + taxa for each group, then Other
fill_levels <- c(HEADER_MOCK,   mock_in_union,
                 HEADER_CONTAM, contam_in_union,
                 "Other")

# Palette: "transparent" for header entries (no bars), colours for taxa
shared_palette <- c(
        setNames(rep("transparent", 2), c(HEADER_MOCK, HEADER_CONTAM)),
        setNames(mock_colors,   mock_in_union),
        setNames(contam_colors, contam_in_union),
        "Other" = "grey80"
)

# Labels: headers in bold, taxa as-is
fill_labels <- c(
        setNames(paste0("**", c(HEADER_MOCK, HEADER_CONTAM), "**"),
                 c(HEADER_MOCK, HEADER_CONTAM)),
        setNames(mock_in_union,   mock_in_union),
        setNames(contam_in_union, contam_in_union),
        "Other" = "Other"
)

# override.aes vectors (one value per fill_levels entry, in order)
header_idx  <- which(fill_levels %in% c(HEADER_MOCK, HEADER_CONTAM))
n_lev       <- length(fill_levels)
oa_fill     <- unname(shared_palette[fill_levels])
oa_fill[header_idx]   <- "transparent"
oa_colour   <- rep("grey40", n_lev)
oa_colour[header_idx] <- "transparent"

# ── Collapse non-union taxa to "Other" and re-sum ─────────────────────────────
control_long_collapsed <- control_long_tbl |>
        mutate(Family_Genus = if_else(
                Family_Genus %in% union_taxa, Family_Genus, "Other"
        )) |>
        group_by(stage, Sample, sample_type, Family_Genus) |>
        summarise(Abundance = sum(Abundance), .groups = "drop") |>
        mutate(Family_Genus = factor(Family_Genus, levels = fill_levels))

# ── Colour-coded x-axis labels by sample QC status ───────────────────────────
# Black  = reads present after off-target removal (informative control)
# Green  = reads in raw but 0 after off-target removal (worked as empty control)
# Red    = 0 reads in raw data (lab failure)
raw_sums   <- sample_sums(ps_16S_raw)
final_sums <- sample_sums(ps_16S_bacteria)

all_control_samples <- control_long_collapsed |> pull(Sample) |> unique()

# Map phyloseq row name (Sample) → full sample_id from metadata column
sid_map <- control_long_tbl |>
        select(Sample, sample_id) |>
        distinct() |>
        tibble::deframe()   # named vector: Sample → sample_id

x_labels <- setNames(
        sapply(all_control_samples, function(s) {
                label <- sid_map[[s]]
                if (is.na(label) || label == "") label <- s   # fallback to row name
                raw   <- raw_sums[s]
                final <- final_sums[s]
                if (is.na(raw) || raw == 0) {
                        sprintf("<span style='color:#D55E00'>%s</span>", label)   # red — lab failure
                } else if (!is.na(final) && final == 0) {
                        sprintf("<span style='color:#009E73'>%s</span>", label)   # green — empty control
                } else {
                        label                                                       # black — has reads
                }
        }),
        all_control_samples
)

# ── Common scale + guide (applied to all three ggplot objects) ────────────────
common_fill <- scale_fill_manual(
        values = shared_palette,
        labels = fill_labels,
        breaks = fill_levels,
        drop   = FALSE,
        name   = sprintf("Taxa with RA \u2265 %d%%", round(100 * RA_THRESHOLD))
)
common_guide <- guides(fill = guide_legend(
        ncol = 1,
        override.aes = list(fill = oa_fill, colour = oa_colour)
))
common_x_labels <- scale_x_discrete(labels = x_labels)
common_theme <- theme_bw() +
        theme(axis.text.x  = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5),
              legend.text  = ggtext::element_markdown(lineheight = 1.2))

# ── Individual before / after figures (printed, not saved) ───────────────────
plot_stage <- function(stage_label) {
        ggplot(control_long_collapsed |> filter(stage == stage_label),
               aes(x = Sample, y = Abundance, fill = Family_Genus)) +
                geom_bar(stat = "identity") +
                facet_grid(~ sample_type, scales = "free_x", space = "free_x") +
                common_fill +
                common_guide +
                common_x_labels +
                common_theme +
                labs(title = stage_label, y = "Read count", x = NULL)
}

print(plot_stage("Before decontam"))
print(plot_stage("After decontam"))
print(plot_stage("After off-target removal"))

# ── Combined before / after figure (saved) ────────────────────────────────────
fig_control_bars_16S <- ggplot(
        control_long_collapsed,
        aes(x = Sample, y = Abundance, fill = Family_Genus)
) +
        geom_bar(stat = "identity") +
        facet_grid(stage ~ sample_type, scales = "free_x", space = "free_x") +
        common_fill +
        common_guide +
        common_x_labels +
        common_theme +
        labs(y = "Read count", x = NULL)

ggsave(here("results", "figures", "02a_16S_controls_pipeline_stages.png"),
       fig_control_bars_16S, width = 16, height = 11, dpi = 300)
