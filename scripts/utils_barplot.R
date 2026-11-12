## ============================================================================
## utils_barplot — shared steps of the taxonomic barplots 06c and 06d
##
## Sourced by 06c_ITS1_vis_barplot_class.R and 06d_16S_ITS1_vis_barplot_microscopy.R,
## so that a taxon gets the same label, the same "Other" rule and the same colour
## in both. 06a/06b keep their own copies of the same rules (left unchanged).
##
## Rules (as 06a/06b):
##   - relative abundance (RA) per flower at the chosen rank;
##   - a taxon is shown by name if its mean RA over all flowers >= threshold,
##     otherwise it goes to "Other"; "Other" is summed within each flower first;
##   - unassigned at the rank -> "<prefix>__Unassigned", drawn second to last.
## ============================================================================

# Rank -> prefix used for the unassigned label
RANK_PREFIX <- c(Phylum = "p", Class = "c", Order = "o", Family = "f", Genus = "g")

# Flower samples only, RA at `rank`, one row per flower x taxon label.
# Labels that mean "not placed at this rank" (NA, and UNITE's
# "<taxon>_cls_Incertae_sedis" style names) all become one Unassigned label.
# Summing within flower and label matters: tax_glom(NArm = FALSE) keeps one NA
# row per higher taxon, and averaging those rows separately would divide the
# Unassigned share by their number.
taxa_long <- function(ps, rank) {
  unassigned <- paste0(RANK_PREFIX[[rank]], "__Unassigned")
  ps |>
    subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower")) |>
    tax_glom(taxrank = rank, NArm = FALSE) |>
    transform_sample_counts(\(x) x / sum(x)) |>
    psmelt() |>
    mutate(taxon = .data[[rank]],
           taxon = if_else(is.na(taxon) | str_detect(taxon, "Incertae_sedis$"), unassigned, taxon)) |>
    summarise(Abundance = sum(Abundance), .by = c(Sample, taxon))
}

# Taxa shown by name, in legend order: most abundant first, Unassigned second
# to last, "Other" last.
top_taxa_levels <- function(long, threshold) {
  shown <- long |>
    summarise(mean_RA = mean(Abundance), .by = taxon) |>
    filter(mean_RA >= threshold) |>
    arrange(desc(mean_RA)) |>
    pull(taxon)
  unassigned <- str_subset(shown, "__Unassigned$")
  c(setdiff(shown, unassigned), unassigned, "Other")
}

# Rare taxa to "Other", summed within each flower; returns a factor in legend order
collapse_other <- function(long, levels) {
  long |>
    mutate(taxon = if_else(taxon %in% levels, taxon, "Other")) |>
    summarise(Abundance = sum(Abundance), .by = c(Sample, taxon)) |>
    mutate(taxon = factor(taxon, levels = levels))
}

# 16S phylum colours: the 06a rule (kelly palette without white and near-black,
# in legend order), so 06d matches 06a exactly
colours_16S_phylum <- function(levels) {
  named <- setdiff(levels, "Other")
  base  <- pals::kelly()[-c(1, 2)]
  pal   <- if (length(named) <= length(base)) base[seq_along(named)] else colorRampPalette(base)(length(named))
  c(setNames(pal, named), Other = "grey80")
}

# ITS1 class colours: shades of the 06b phylum colour, so a class reads as part
# of its phylum (teal = Ascomycota, tan = Basidiomycota). Darkest = most abundant.
# Classes of any other phylum get a blue ramp. Unassigned pink and Other grey as in 06b.
PHYLUM_RAMPS <- list(
  p__Ascomycota    = c("#00594C", "#66C2A5", "#CDEBE2"),   # teal, around the 06b Ascomycota colour
  p__Basidiomycota = c("#8C510A", "#E5C494", "#F6E8C3"),   # tan, around the 06b Basidiomycota colour
  other            = c("#253494", "#8DA0CB", "#DDE3F2")    # blue, any other phylum
)

# Class levels reordered by phylum (phylum by total RA), then by RA within the
# phylum; returns the levels and their colours
its1_class_levels_colours <- function(levels, long, ps) {
  tt <- as.data.frame(tax_table(ps))
  phylum_of <- tt |> distinct(Class, Phylum) |> filter(!is.na(Class)) |> deframe()
  named <- setdiff(levels, c("Other", str_subset(levels, "__Unassigned$")))
  ra    <- long |> summarise(mean_RA = mean(Abundance), .by = taxon) |> deframe()
  info  <- tibble(taxon = named, phylum = phylum_of[named], mean_RA = ra[named]) |>
    mutate(phylum_RA = sum(mean_RA), .by = phylum) |>
    arrange(desc(phylum_RA), desc(mean_RA))
  colours <- info |>
    mutate(ramp = if_else(phylum %in% names(PHYLUM_RAMPS), phylum, "other")) |>
    group_by(ramp) |>
    mutate(colour = colorRampPalette(PHYLUM_RAMPS[[first(ramp)]])(n() + 1)[seq_len(n())]) |>
    ungroup()
  new_levels <- c(info$taxon, str_subset(levels, "__Unassigned$"), "Other")
  list(levels  = new_levels,
       colours = c(setNames(colours$colour, colours$taxon),
                   setNames("#E68FAC", str_subset(levels, "__Unassigned$")),   # pink, as 06b
                   Other = "grey80"))
}

# Theme of the 06a/06b barplots
bar_theme <- function(x_angle = 45, x_size = 9) {
  theme_bw() +
    theme(strip.placement  = "outside",
          strip.background = element_blank(),
          strip.text.x     = element_text(size = 8),
          strip.text.y     = element_text(size = 9),
          panel.spacing.x  = unit(0.15, "lines"),
          legend.position  = "right",
          axis.text.x      = element_text(angle = x_angle, hjust = 1,
                                          vjust = if (x_angle == 90) 0.5 else 1, size = x_size))
}

# ── Sampling deviations on flower-level barplots (294 set) ────────────────────
# Label colours by deviation type (Okabe-Ito: orange, blue, reddish purple)
DEVIATION_COLS <- c(no_buffer = "#E69F00", two_flowers = "#0072B2", insect_present = "#CC79A7")

# Colour the x label of each deviation flower by its deviation type and put a
# triangle of the same colour above its bar. The flower-level plots label bars
# by flower_id (sample id without the c/p prefix), which is shared by the
# bagged and unbagged rows; the triangle shows which row the deviation is in.
#   p          flower-level ggplot with x = flower_id
#   deviations one row per deviation flower: flower_id, sampling_deviation and
#              the facet columns of p (sample_type, management_type, farm_id, tree_id)
#   label_size size of the x labels (as in the plot)
mark_deviations <- function(p, deviations, label_size = 4) {
  label_colour <- \(x) coalesce(DEVIATION_COLS[deviations$sampling_deviation[match(x, deviations$flower_id)]], "black")
  p +
    # Filled triangle above the bar. The fill aesthetic already maps the taxa, so
    # the triangle fill is given directly, one colour per deviation flower
    geom_point(data = deviations, aes(x = flower_id, y = 1.045, colour = sampling_deviation),
               inherit.aes = FALSE, shape = 25, size = 1.8, stroke = 0.6,
               fill = DEVIATION_COLS[deviations$sampling_deviation]) +
    scale_colour_manual(values = DEVIATION_COLS, breaks = names(DEVIATION_COLS), name = "Sampling deviation") +
    # Legend keys: the same filled triangles
    guides(colour = guide_legend(override.aes = list(shape = 25, fill = unname(DEVIATION_COLS), size = 3))) +
    scale_x_discrete(drop = TRUE, labels = \(x) ifelse(
      label_colour(x) == "black", x,
      paste0("<span style='color:", label_colour(x), "'>**", x, "**</span>"))) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.08),
                       breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.01))) +
    theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = label_size))
}

# ── Anchoring the two dominant taxa to the edges of the box ───────────────────
# In a stacked bar only the top segment (ending at 100%) and the bottom segment
# (starting at 0%) share a common baseline across bars; the middle segments float
# and are hard to compare (Cleveland & McGill 1984; P. Schloss, Riffomonas).
# So the most abundant taxon goes to the top, the second most abundant to the
# bottom, Unassigned and Other to the middle. The first factor level is drawn at
# the top of the stack and listed first in the legend.
#   levels   shown taxa (as from top_taxa_levels() or its1_class_levels_colours())
#   ranking  named taxa by decreasing mean RA (no Unassigned, no Other)
#   group_of optional named vector taxon -> higher taxon (e.g. class -> phylum);
#            given, the taxa of one higher taxon stay together as a block
anchor_levels <- function(levels, ranking, group_of = NULL) {
  special <- c(str_subset(levels, "__Unassigned$"), "Other")
  named   <- setdiff(levels, special)
  top     <- ranking[1]
  bottom  <- ranking[2]
  if (is.null(group_of) || identical(group_of[[top]], group_of[[bottom]])) {
    return(c(top, setdiff(named, c(top, bottom)), special, bottom))
  }
  g         <- group_of[named]
  top_block <- c(top, setdiff(named[g %in% g[[top]]], top))
  bot_block <- c(setdiff(named[g %in% g[[bottom]]], bottom), bottom)
  middle    <- setdiff(named, c(top_block, bot_block))
  c(top_block, middle, special, bot_block)
}

# Class -> phylum, from the taxonomy table (for anchor_levels(group_of = ...))
class_to_phylum <- function(ps) {
  as.data.frame(tax_table(ps)) |> distinct(Class, Phylum) |> filter(!is.na(Class)) |> deframe()
}
