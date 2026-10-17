# ---- Packages ----
library(here)
library(tidyverse)
library(phyloseq)
library(ggforce)
library(ggh4x)
library(ggtext)      # element_markdown(): colour has to live inside the label
library(patchwork)

# ---- Output directories ----
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

################################################################################
# Load curated taxa database from CSV
# Columns: Domain | Kingdom | Genus | Species | GenusSyn | SpeciesSyn | SH |
#          context | role | note | DOI
# Search is genus-level; GenusSyn enables synonym matching.
################################################################################

curated_db <- read_delim(
    here("data", "curated_taxa_db.csv"),
    delim          = ";",
    trim_ws        = TRUE,
    locale         = locale(encoding = "UTF-8"),
    show_col_types = FALSE
) |>
    filter(!is.na(Genus), Genus != "") |>
    mutate(
        Genus    = str_trim(Genus),
        GenusSyn = str_trim(GenusSyn),
        role     = if_else(is.na(role) | role == "", "unknown", role)
    )

# ---- Dataset-specific genus sets ----
# ITS1: Fungi + Stramenopila (UNITE eukaryotes release covers both)
its1_genus_set <- curated_db |>
    filter(Kingdom %in% c("Fungi", "Stramenopila")) |>
    pull(Genus) |> unique()

# 16S: Bacteria + Stramenopila (off-target eukaryotes in unfiltered 16S data)
# Stramenopila intentionally in both: Phytophthora can appear in 16S as off-target
s16_genus_set <- curated_db |>
    filter(Domain == "Bacteria" | Kingdom == "Stramenopila") |>
    pull(Genus) |> unique()

# ---- Role colour palette ----
role_colors <- c(
    "pathogen"        = "#D55E00",   # vermilion  — cacao pathogens
    "biocontrol"      = "#009E73",   # green      — biocontrol agents
    "symbiont"        = "#0072B2",   # blue       — plant / phylloplane symbionts
    "insect symbiont" = "#CC79A7",   # pink       — endosymbionts of insects
    "insect pathogen" = "#E69F00",   # orange     — entomopathogenic fungi
    "unknown"         = "#999999"    # grey       — role not yet assigned in CSV
)

# Guard. "insect symbiont" (Wolbachia) and "insect pathogen" (Simplicillium) used
# to be absent from this palette, so ggplot dropped their points and both genera
# appeared as empty rows -- indistinguishable from "not detected". A missing
# colour must now stop the script instead of quietly deleting data.
missing_roles <- setdiff(unique(curated_db$role), names(role_colors))
if (length(missing_roles) > 0) {
    stop("role values without a colour in role_colors: ",
         paste(missing_roles, collapse = ", "),
         ". Add them, or points for those genera will be dropped silently.")
}

################################################################################
# Role colours for TEXT
#
# The role used to be a filled square in the left margin. It is now the colour of
# the genus name itself, and that changes what the palette has to satisfy.
# Okabe-Ito above is designed for fills: as small text on white, four of the five
# roles in use fail WCAG AA (4.5:1) -- insect pathogen is worst at 2.3:1. A large
# solid square does not care about contrast; text does.
#
# WHY THE DARKENING IS UNIFORM AND MUST STAY UNIFORM
#
# The obvious optimisation -- darken each colour by just enough to reach 4.5:1 --
# is the worst option available. Okabe-Ito stays distinguishable under colour
# blindness largely through differences in LIGHTNESS, and equalising every colour
# to the same contrast ratio equalises exactly that:
#
#                       AA as text      min pairwise dE, deuteranopia
#   original (squares)  2.3 .. 5.2:1        18.1
#   darken 0.30         4.6 .. 8.7:1        11.6   <- in use
#   per-colour fit      4.5:1 everywhere     3.6   <- destroys the palette
#
# So the whole palette is shifted down by one amount, which keeps the lightness
# ordering intact. The cost is stated plainly: deuteranopia separation drops from
# 18.1 to 11.6. The tightest pair is pathogen vs insect pathogen, and
# insect pathogen is a single genus on these figures (Simplicillium).
#
# role_colors stays the canonical fill palette: it is what 00_curated_shared.rds
# hands to anything downstream.
ROLE_TEXT_DARKEN <- 0.30
role_colors_text <- colorspace::darken(role_colors, ROLE_TEXT_DARKEN,
                                       space = "HCL")

# Wrap a row label in the colour of its role. The colour has to travel inside the
# string: axis.text.y = element_text(colour = <vector>) is a theme element and is
# recycled from the start in every facet panel, so it mismatches as soon as the
# panels hold different rows -- which is the case on all four figures here.
md_label <- function(txt, role)
    sprintf("<span style='color:%s'>%s</span>", role_colors_text[role], txt)

# Contrast printed at every run so the palette choice stays visible in the log
# rather than buried in this comment.
local({
    lum <- function(h) {
        cc <- grDevices::col2rgb(h) / 255
        cc <- ifelse(cc <= 0.03928, cc / 12.92, ((cc + 0.055) / 1.055)^2.4)
        0.2126 * cc[1, ] + 0.7152 * cc[2, ] + 0.0722 * cc[3, ]
    }
    message(sprintf("Role text palette (darken %.2f), contrast against white:\n%s",
                    ROLE_TEXT_DARKEN,
                    paste(sprintf("  %-16s %s -> %s  %4.1f:1",
                                  names(role_colors), role_colors,
                                  role_colors_text,
                                  1.05 / (lum(role_colors_text) + 0.05)),
                          collapse = "\n")))
})

# ---- Per-genus primary / secondary role ----
# Multiple species rows per genus are collapsed; primary = most frequent known role.
genus_roles <- curated_db |>
    group_by(Genus) |>
    summarise(
        Kingdom = first(Kingdom[!is.na(Kingdom) & Kingdom != ""]),
        Domain  = first(Domain[!is.na(Domain)  & Domain  != ""]),
        primary_role = {
            r <- role[role != "unknown"]
            if (length(r) == 0) "unknown"
            else names(sort(table(r), decreasing = TRUE))[1]
        },
        secondary_role = {
            r <- role[role != "unknown"]
            if (n_distinct(r) > 1) names(sort(table(r), decreasing = TRUE))[2]
            else NA_character_
        },
        .groups = "drop"
    )

# ---- Genus lookup: primary names + synonyms → role ----
genus_lookup <- bind_rows(
    genus_roles |> rename(Genus_lookup = Genus),
    curated_db |>
        filter(!is.na(GenusSyn), GenusSyn != "") |>
        distinct(Genus, GenusSyn) |>
        left_join(genus_roles, by = "Genus") |>
        rename(Genus_lookup = GenusSyn) |>
        select(Genus_lookup, Kingdom, Domain, primary_role, secondary_role)
) |>
    distinct(Genus_lookup, .keep_all = TRUE) |>
    mutate(Genus_lookup = str_trim(Genus_lookup))


# ---- Load phyloseq objects ----
ps_ITS1_raw <- readRDS(here("results", "rds", "ps_ITS1_raw.rds"))
ps_16S_raw  <- readRDS(here("results", "rds", "ps_16S_raw.rds"))

################################################################################
# search_taxa(): find curated genera in a phyloseq object.
# `genera`: optional character vector restricting which genera to search.
# Searches both primary Genus names and synonyms (GenusSyn) via genus_lookup.
################################################################################

search_taxa <- function(physeq, dataset_label, genera = NULL) {

    tt <- tax_table(physeq) |>
        as.data.frame() |>
        rownames_to_column("ASV")

    genus_col <- names(tt)[tolower(names(tt)) == "genus"]

    all_genera <- unique(str_trim(c(curated_db$Genus, curated_db$GenusSyn))) |>
        (\(x) x[!is.na(x) & x != ""])()

    # Restrict to dataset-specific genus set if provided
    if (!is.null(genera)) all_genera <- intersect(all_genera, genera)

    hits <- tt |>
        filter(
            .data[[genus_col]] %in% paste0("g__", all_genera) |
            .data[[genus_col]] %in% all_genera
        ) |>
        mutate(
            Genus_clean = str_remove(.data[[genus_col]], "^g__"),
            dataset     = dataset_label
        )

    if (nrow(hits) == 0) {
        message(sprintf("No curated genera found in %s dataset.", dataset_label))
        return(tibble())
    }

    otu <- otu_table(physeq)
    if (!taxa_are_rows(physeq)) otu <- t(otu)
    otu_df <- as.data.frame(as.matrix(otu))

    sample_meta <- sample_data(physeq) |>
        as.data.frame() |>
        (\(x) { class(x) <- "data.frame"; x })()

    bio_samples <- rownames(sample_meta)[
        sample_meta$sample_type %in% c("bagged_flower", "unbagged_flower")
    ]

    otu_bio <- otu_df[hits$ASV, bio_samples, drop = FALSE]

    hits <- hits |>
        mutate(
            n_samples_present       = rowSums(otu_bio > 0),
            prevalence              = n_samples_present / length(bio_samples),
            mean_reads_when_present = ifelse(
                n_samples_present > 0,
                rowSums(otu_bio) / n_samples_present,
                0
            )
        ) |>
        left_join(
            genus_lookup |>
                rename(Genus_clean = Genus_lookup) |>
                select(-Kingdom, -Domain),
            by = "Genus_clean"
        ) |>
        rename(role = primary_role) |>
        select(
            ASV, dataset, Genus_clean, role, secondary_role,
            n_samples_present, prevalence, mean_reads_when_present,
            any_of(c("Family", "Order", "Class", "Phylum", "Species", "SH"))
        )

    presence_long <- otu_bio |>
        rownames_to_column("ASV") |>
        pivot_longer(-ASV, names_to = "sample_id", values_to = "n_reads") |>
        filter(n_reads > 0) |>
        left_join(hits |> select(ASV, Genus_clean, role),
                  by = "ASV", relationship = "many-to-many") |>
        group_by(Genus_clean, role) |>
        summarise(
            n_samples  = n_distinct(sample_id),
            sample_ids = paste(sort(unique(sample_id)), collapse = ", "),
            .groups    = "drop"
        )

    message(sprintf("\n=== %s: curated genera detected in biological samples ===",
                    dataset_label))
    print(presence_long)

    attr(hits, "presence_long") <- presence_long
    hits
}

# ITS1: search only Fungi + Stramenopila genera
hits_ITS1 <- search_taxa(ps_ITS1_raw, "ITS1", genera = its1_genus_set)
# 16S: search only Bacteria + Stramenopila genera (off-targets)
hits_16S  <- search_taxa(ps_16S_raw,  "16S",  genera = s16_genus_set)

curated_hits <- bind_rows(hits_ITS1, hits_16S)

message("\n=== Curated taxa found in dataset ===")
print(curated_hits |> select(dataset, Genus_clean, role, n_samples_present,
                              prevalence, mean_reads_when_present))

presence_ITS1 <- attr(hits_ITS1, "presence_long") |> mutate(dataset = "ITS1")
presence_16S  <- attr(hits_16S,  "presence_long") |> mutate(dataset = "16S")
presence_all  <- bind_rows(presence_ITS1, presence_16S)

write_csv(presence_all, here("results", "tables", "00_curated_taxa_sample_lists.csv"))
write_csv(curated_hits, here("results", "tables", "00_curated_taxa_hits.csv"))

################################################################################
# Summary table: one row per dataset × genus (dataset-specific genus sets)
################################################################################

n_bio_ITS1 <- sum(sample_data(ps_ITS1_raw)$sample_type %in%
                      c("bagged_flower", "unbagged_flower"))
n_bio_16S  <- sum(sample_data(ps_16S_raw)$sample_type  %in%
                      c("bagged_flower", "unbagged_flower"))

curated_summary <- curated_hits |>
    mutate(n_bio_total = if_else(dataset == "ITS1", n_bio_ITS1, n_bio_16S)) |>
    group_by(dataset, Genus_clean, role, n_bio_total) |>
    summarise(
        n_ASVs            = n(),
        n_samples_present = max(n_samples_present),
        prevalence        = n_samples_present / n_bio_total[1],
        .groups           = "drop"
    ) |>
    mutate(detected = TRUE)

# Expand grid only over the relevant genus set per dataset
curated_summary_full <- bind_rows(
    tidyr::expand_grid(
        genus_roles |> filter(Genus %in% its1_genus_set) |>
            rename(Genus_clean = Genus, role = primary_role),
        dataset = "ITS1"
    ),
    tidyr::expand_grid(
        genus_roles |> filter(Genus %in% s16_genus_set) |>
            rename(Genus_clean = Genus, role = primary_role),
        dataset = "16S"
    )
) |>
    left_join(curated_summary, by = c("Genus_clean", "role", "dataset")) |>
    mutate(
        detected          = replace_na(detected, FALSE),
        n_ASVs            = replace_na(n_ASVs, 0L),
        n_samples_present = replace_na(n_samples_present, 0L),
        prevalence        = replace_na(prevalence, 0)
    ) |>
    arrange(dataset, role, desc(prevalence))

write_csv(curated_summary_full, here("results", "tables", "00_curated_taxa_summary.csv"))

################################################################################
# Per-sample long table (shared object for downstream use)
################################################################################

make_curated_long <- function(physeq, hits, dataset_label) {
    if (nrow(hits) == 0) return(tibble())

    otu_full <- otu_table(physeq)
    if (!taxa_are_rows(physeq)) otu_full <- t(otu_full)
    sample_totals <- colSums(as.matrix(otu_full))

    ps_sub <- prune_taxa(hits$ASV, physeq) |>
        subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))

    psmelt(ps_sub) |>
        rename(ASV = OTU) |>
        mutate(RA = Abundance / sample_totals[as.character(Sample)]) |>
        left_join(hits |> select(ASV, Genus_clean, role), by = "ASV") |>
        group_by(Sample, management_type, farm_id, tree_id, sample_type,
                 Genus_clean, role) |>
        summarise(n_reads = sum(Abundance), RA = sum(RA), .groups = "drop") |>
        mutate(dataset = dataset_label)
}

curated_long <- bind_rows(
    make_curated_long(ps_ITS1_raw, hits_ITS1, "ITS1"),
    make_curated_long(ps_16S_raw,  hits_16S,  "16S")
)

saveRDS(
    list(
        curated_long         = curated_long,
        curated_summary_full = curated_summary_full,
        curated_taxa         = genus_roles |>
                                   select(Genus, Kingdom, role = primary_role),
        role_colors          = role_colors
    ),
    here("results", "rds", "00_curated_shared.rds")
)
message("Shared object saved → results/rds/00_curated_shared.rds")

################################################################################
# Constants shared by all figures
################################################################################

MGMT_LEVELS <- c("inside_forest", "near_forest", "agroforest", "full_sun")
FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

# Treatment colours, identical to 10d, 09d and 16d so the reader carries one
# mapping across the whole project.
GROUP_COLS <- c(bagged = "#4C72B0", unbagged = "#DD8452")

# ASV-level panel: an ASV must be detected in at least this many flowers. Same
# threshold as 16d, where rarity itself is the phenomenon under study.
MIN_DETECT <- 10

################################################################################
# Detection versus comparison
#
# The dumbbell shows a DIFFERENCE of prevalences, n_bag/n - n_unbag/n. A genus
# detected in k flowers altogether therefore satisfies
#
#     |difference|  <=  k / n
#
# because the maximum is reached when every detection falls in one treatment. At
# k = 1 the difference cannot exceed 0.7 pp whatever the biology.
#
# The constraint is on the effect SCALE, not on precision. A Newcombe score
# interval puts Moniliophthora (k = 1) at -0.007 +- 0.028: a precise estimate of
# almost exactly zero. Calling such a row "not estimable" would be wrong; it is
# the absolute prevalence difference that cannot say anything about a genus this
# rare, because the interesting question there is a ratio, not a difference.
#
# So the threshold is derived rather than chosen: if the smallest difference
# worth reading is DIFF_MIN, a genus needs at least DIFF_MIN * n detections to be
# able to show it. Unlike the flip thresholds below, this one is computable
# before the treatment split has been looked at, so it cannot be calibrated
# towards a wanted list.
DIFF_MIN <- 0.10          # same value as ASV_FLIP_MIN: one notion of "readable"

# Flowers per treatment, read off the data rather than written as 147: if the
# filtering upstream ever changes, the threshold must not be computed from a
# stale number.
N_PER_GROUP <- curated_long |>
    distinct(dataset, Sample, sample_type) |>
    count(dataset, sample_type) |>
    pull(n) |>
    unique()
stopifnot(length(N_PER_GROUP) == 1)

DIFF_MIN_N <- ceiling(DIFF_MIN * N_PER_GROUP)   # 15 flowers at DIFF_MIN = 0.10

################################################################################
# Heterogeneity marks
#
# Two independent axes hide inside a genus row, and pooling either one can make a
# treatment effect appear or vanish:
#   -- ASVs of one genus can move in opposite directions (two Wolbachia ASVs
#      differing by two nucleotides sit in opposite treatments and cancel);
#   -- the same genus can favour bagged flowers on one farm and unbagged on another.
# A marked row is one that must not be read as a single story.
#
# WHY A MAGNITUDE THRESHOLD AND NOT JUST A SIGN
#
# A pure sign criterion is not the neutral choice, it is the most permissive one: on
# the farm axis it marks 19 of 27 genera, so the mark becomes background and stops
# distinguishing anything. There is no threshold-free option -- not marking is also a
# choice. The thresholds below are therefore explicit, tunable, and their effect is
# printed at every run (see the sensitivity table further down). Set either to 0 to
# recover the pure sign criterion.
#
# Scale. The SE of a difference of two proportions at p ~ 0.5 is sqrt(2*0.25/n).
# Pooled across flowers n = 147 per treatment -> SE ~ 0.058; within one farm n = 21
# -> SE ~ 0.154. Farms are therefore ~2.6x noisier, which is why the two thresholds
# differ rather than being one round number applied twice.
#
# Calibration, stated plainly: the thresholds were chosen to retain Wolbachia, whose
# effect is independently confirmed by the within-tree permutation test in 16d
# (q = 0.004) and by MaAsLin 3 in 16e (q = 0.0013). That is calibration against an
# external positive control, not selection of a threshold that produces a desired
# list. On the ASV axis a threshold of 0.20 already discards Wolbachia (+0.177 and
# -0.102), which is the hard ceiling for that axis.
ASV_FLIP_MIN  <- 0.10   # ~1.7 SE at n = 147 per treatment group
FARM_FLIP_MIN <- 0.20   # ~1.3 SE at n = 21 per farm x treatment

MIXED_MARK <- " \u2021"   # ASVs within the genus disagree
FARM_MARK  <- " *"        # farms disagree

# Strict comparisons against zero as well as the threshold: with thr = 0 a plain
# `d >= thr` would be satisfied by d == 0 on both sides and mark every genus.
flips <- function(d, thr) any(d > 0 & d >= thr) & any(d < 0 & d <= -thr)

# How many genera each threshold would mark, printed so the dependence of the figure
# on this choice is visible at every run rather than buried in the code.
flip_counts <- function(df, d_col, thresholds = c(0, 0.05, 0.10, 0.20)) {
    vapply(thresholds, function(thr)
        df |>
            group_by(dataset, Genus_clean) |>
            summarise(f = flips(.data[[d_col]], thr), .groups = "drop") |>
            pull(f) |> sum(),
        integer(1))
}

# n_seqs per genus per dataset, used in the row labels
n_seqs <- curated_hits |>
    count(dataset, Genus_clean, name = "n_seqs")

################################################################################
# Axis 1 — per-ASV detection
#
# curated_long is aggregated to genus, which is the right unit for the curated
# database (roles are assigned per genus) but hides the case that motivated this
# rework. This table keeps the ASV resolution.
################################################################################

asv_detection <- function(physeq, hits, dataset_label) {
    if (nrow(hits) == 0) return(tibble())

    # as(x, "matrix") and not as.matrix(): the latter leaves the otu_table class
    # attached, so the comparison below would return an S4 object carrying
    # logicals, which fails phyloseq's validity check further down.
    mat <- as(otu_table(physeq), "matrix")
    if (!taxa_are_rows(physeq)) mat <- t(mat)

    sdat <- sample_data(physeq) |> as.data.frame() |>
        (\(x) { class(x) <- "data.frame"; x })()
    flowers <- rownames(sdat)[sdat$sample_type %in%
                                  c("bagged_flower", "unbagged_flower")]

    det <- mat[hits$ASV, flowers, drop = FALSE] > 0
    bag <- sdat[flowers, "sample_type"] == "bagged_flower"

    tibble(
        ASV         = rownames(det),
        dataset     = dataset_label,
        Genus_clean = hits$Genus_clean[match(rownames(det), hits$ASV)],
        n_det       = rowSums(det),
        n_bagged    = rowSums(det[, bag, drop = FALSE]),
        n_unbag     = rowSums(det[, bag == FALSE, drop = FALSE])
    ) |>
        mutate(
            rate_bagged = n_bagged / sum(bag),
            rate_unbag  = n_unbag  / sum(bag == FALSE),
            diff_rate   = rate_bagged - rate_unbag
        )
}

asv_tbl <- bind_rows(
    asv_detection(ps_ITS1_raw, hits_ITS1, "ITS1"),
    asv_detection(ps_16S_raw,  hits_16S,  "16S")
)

asv_kept <- asv_tbl |> filter(n_det >= MIN_DETECT)

mixed_genera <- asv_kept |>
    group_by(dataset, Genus_clean) |>
    summarise(n_asv_kept = n(),
              asv_flip   = flips(diff_rate, ASV_FLIP_MIN),
              .groups    = "drop")

################################################################################
# Axis 2 — per-farm prevalence
#
# Computed here rather than inside the farm figure because the pooled figure needs
# the reversal flag derived from it.
################################################################################

prev_farm_raw <- curated_long |>
    group_by(dataset, Genus_clean, role, farm_id, sample_type) |>
    summarise(prevalence = mean(RA > 0), .groups = "drop") |>
    pivot_wider(names_from = sample_type, values_from = prevalence,
                values_fill = 0) |>
    rename(rate_bagged = bagged_flower, rate_unbag = unbagged_flower) |>
    mutate(d = rate_bagged - rate_unbag)

farm_flip_tbl <- prev_farm_raw |>
    group_by(dataset, Genus_clean) |>
    summarise(n_farm     = n(),
              n_farm_pos = sum(d > 0),
              n_farm_neg = sum(d < 0),
              max_abs_d  = max(abs(d)),
              farm_flip  = flips(d, FARM_FLIP_MIN),
              .groups    = "drop")

# ---- Sensitivity of both marks to the threshold ----
THR_GRID <- c(0, 0.05, 0.10, 0.20)
message(sprintf(
    "\nHeterogeneity marks, genera flagged at |delta| >= %s:\n  farm axis (* , %d genera): %s\n  ASV axis  (\u2021 , %d genera): %s\n  in use: ASV_FLIP_MIN = %.2f, FARM_FLIP_MIN = %.2f",
    paste(sprintf("%.2f", THR_GRID), collapse = " / "),
    nrow(farm_flip_tbl), paste(flip_counts(prev_farm_raw, "d", THR_GRID), collapse = " / "),
    n_distinct(paste(asv_kept$dataset, asv_kept$Genus_clean)),
    paste(flip_counts(asv_kept, "diff_rate", THR_GRID), collapse = " / "),
    ASV_FLIP_MIN, FARM_FLIP_MIN))

write_csv(
    asv_kept |>
        left_join(mixed_genera |> select(dataset, Genus_clean, asv_flip),
                  by = c("dataset", "Genus_clean")) |>
        arrange(dataset, desc(abs(diff_rate))),
    here("results", "tables", "00d_curated_taxa_asv.csv"))

################################################################################
# Genus-level prevalence, detected taxa only
#
# Undetected genera are dropped from every figure: seven of them contribute a
# third of the vertical space while saying one thing that a printed list says
# better. They stay in 00_curated_taxa_summary.csv with detected = FALSE.
################################################################################

# k in |difference| <= k / N_PER_GROUP: flowers holding the genus, both
# treatments pooled. curated_long is already aggregated to genus and carries
# psmelt's zeros, so this is a plain count of the positive rows.
n_det_tbl <- curated_long |>
    group_by(dataset, Genus_clean) |>
    summarise(n_det_genus = sum(RA > 0), .groups = "drop")

prev_genus <- curated_long |>
    group_by(dataset, Genus_clean, role, sample_type) |>
    summarise(prevalence = mean(RA > 0), .groups = "drop") |>
    pivot_wider(names_from = sample_type, values_from = prevalence,
                values_fill = 0) |>
    rename(rate_bagged = bagged_flower, rate_unbag = unbagged_flower) |>
    mutate(diff_rate = rate_bagged - rate_unbag) |>
    filter(rate_bagged > 0 | rate_unbag > 0) |>
    left_join(n_seqs,    by = c("dataset", "Genus_clean")) |>
    left_join(n_det_tbl, by = c("dataset", "Genus_clean")) |>
    left_join(mixed_genera  |> select(dataset, Genus_clean, asv_flip),
              by = c("dataset", "Genus_clean")) |>
    left_join(farm_flip_tbl |> select(dataset, Genus_clean, farm_flip,
                                      n_farm_pos, n_farm_neg, max_abs_d),
              by = c("dataset", "Genus_clean")) |>
    mutate(
        asv_flip  = replace_na(asv_flip, FALSE),
        farm_flip = replace_na(farm_flip, FALSE),
        unit      = if_else(dataset == "ITS1", "OTUs", "ASVs"),
        # can a treatment difference be read off this row at all (see DIFF_MIN)
        comparable = n_det_genus >= DIFF_MIN_N,
        block = factor(comparable, levels = c(TRUE, FALSE),
                       labels = c("treatment difference shown",
                                  "prevalence difference\nbounded by rarity")),
        # the unit differs between datasets, so a genus present in both (Hannaella,
        # Phytophthora) still gets two distinct labels and cannot collapse to one row
        #
        # Marks are suppressed off the rarity block: both of them assert something
        # about a difference. No mark currently falls there, so this changes
        # nothing today, but the rule has to hold at other threshold settings too.
        label_plain = sprintf("%s (%d %s)%s%s", Genus_clean, n_seqs, unit,
                              if_else(asv_flip  & comparable, MIXED_MARK, ""),
                              if_else(farm_flip & comparable, FARM_MARK,  "")),
        # the axis label carries its own colour; label_plain stays for joins,
        # console messages and anything that must not see HTML
        label = md_label(label_plain, role),
        dataset = factor(dataset, levels = c("ITS1", "16S"))
    )

# Reported per dataset, not pooled: a genus searched in both marker sets can be
# absent from one and present in the other (Phytophthora is missing from ITS1 but
# turns up as an off-target in 16S), and pooling the two would call it undetected.
undetected <- curated_summary_full |>
    filter(!detected) |>
    arrange(dataset, Genus_clean)
message(sprintf("Not detected, omitted from that marker's panel: %s",
                paste(sprintf("%s (%s)", undetected$Genus_clean, undetected$dataset),
                      collapse = ", ")))

# Marks carried into the summary table so the figure is checkable without a re-run
write_csv(
    curated_summary_full |>
        left_join(prev_genus |> as_tibble() |>
                      select(dataset, Genus_clean, asv_flip, farm_flip,
                             n_farm_pos, n_farm_neg, max_abs_d),
                  by = c("dataset", "Genus_clean")),
    here("results", "tables", "00_curated_taxa_summary.csv"))

# One row order for every figure, computed once here so a taxon sits on the same
# line in all four panels. Two keys: the rarity block first, and inside each
# block the only quantity that means anything there -- effect size where a
# difference is shown, number of detections where it is not.
ROW_LEVELS <- bind_rows(
    prev_genus |> filter(!comparable) |> arrange(n_det_genus),
    prev_genus |> filter(comparable)  |> arrange(diff_rate)
) |> pull(label)

prev_genus <- prev_genus |> mutate(label = factor(label, levels = ROW_LEVELS))

################################################################################
# Dumbbell: both values on ONE axis, the segment between them IS the effect.
#
# A diverging layout forces the eye to compare two lengths measured from a shared
# centre in opposite directions; this does not.
#
# The role used to sit in a reserved negative margin as a coloured square. It is
# now the colour of the genus name, so the margin is gone and the panel starts
# just left of zero.
################################################################################

# `split`: draw the treatment contrast only on rows where it can be read
#   (comparable == TRUE) and replace it with a single neutral point plus the
#   detection count on the rest. FALSE by default, which reproduces the previous
#   behaviour exactly -- 00c is filtered to comparable rows and needs no split,
#   and 00d is ASV-level and carries no such column at all.
#   data = NULL means "inherit from the plot"; a formula is applied to the
#   plot's data, so inheritance still works and the function stays reusable.
dumbbell_layers <- function(point_size = 2.6, split = FALSE) {
    contrast_data <- if (split) (~ dplyr::filter(.x, comparable)) else NULL
    rare_layers <- if (!split) list() else list(
        # pooled prevalence: the row still says "it is here, this often", but the
        # grammar changes so that no difference can be read off it. Grey is set
        # outside aes() so it stays out of scale_colour_manual and adds no third
        # key to the treatment legend.
        geom_point(data = ~ dplyr::filter(.x, !comparable),
                   aes(y = label, x = n_det_genus / (2 * N_PER_GROUP)),
                   size = point_size, colour = "grey45"),
        geom_text(data = ~ dplyr::filter(.x, !comparable),
                  aes(y = label, x = n_det_genus / (2 * N_PER_GROUP) + 0.015,
                      label = sprintf("n = %d", n_det_genus)),
                  hjust = 0, size = 2.4, colour = "grey45")
    )
    c(list(
        geom_segment(data = contrast_data,
                     aes(y = label, yend = label,
                         x = rate_unbag, xend = rate_bagged),
                     colour = "grey70", linewidth = 0.6),
        geom_point(data = contrast_data,
                   aes(y = label, x = rate_unbag, colour = "unbagged"),
                   size = point_size),
        geom_point(data = contrast_data,
                   aes(y = label, x = rate_bagged, colour = "bagged"),
                   size = point_size)),
      rare_layers,
      list(
        scale_colour_manual(values = GROUP_COLS, name = NULL),
        theme_bw(base_size = 10),
        # axis.text.y and plot.subtitle are markdown so the role colour can be
        # carried inside the strings. plot.caption deliberately is NOT: the
        # captions contain *, | and >=, which a markdown parser would mangle.
        theme(legend.position    = "bottom",
              plot.title         = element_text(size = 12, face = "bold"),
              plot.subtitle      = element_markdown(size = 8, colour = "grey40"),
              plot.caption       = element_text(size = 7, colour = "grey40"),
              axis.text.y        = element_markdown(),
              panel.grid.major.y = element_blank(),
              strip.background   = element_rect(fill = "grey92"))
    ))
}

MARK_CAPTION <- sprintf(
    paste0("\u2021 ASVs of the genus disagree in direction (|delta prevalence| >= %.2f, see 00d)",
           "  |  * farms disagree (|delta| >= %.2f, see 00c)",
           "  |  a marked row must not be read as one story"),
    ASV_FLIP_MIN, FARM_FLIP_MIN)

# Two lines on purpose: as one line this runs wider than the 9in canvas and the
# left end is clipped away.
# The fill scale is gone with the swatch, so there is no role legend any more.
# This replaces it and is better than a legend was: every role is written in its
# own colour, so the key demonstrates the encoding instead of describing it.
# Separator is a plain unicode middle dot, not an HTML entity -- gridtext parses
# only a subset of HTML and a literal character is the safer choice.
ROLES_SHOWN <- sort(unique(curated_long$role))
ROLE_LEGEND <- paste0(
    "Genus name colour = ecological role: ",
    paste(md_label(ROLES_SHOWN, ROLES_SHOWN), collapse = " \u00b7 "))

BLOCK_CAPTION <- sprintf(
    paste0("A genus found in k flowers can differ by at most k/%d between treatments, ",
           "so the contrast is drawn only where k >= %d\n",
           "below that the row reports detection alone, ",
           "with the grey point at prevalence over both treatments"),
    N_PER_GROUP, DIFF_MIN_N)

################################################################################
# Fig 1 — genus-level prevalence, bagged vs unbagged
################################################################################

p_dot <- ggplot(prev_genus) +
    dumbbell_layers(split = TRUE) +
    scale_x_continuous(
        breaks = seq(0, 1, 0.25),           # no break in the swatch margin
        limits = c(-0.02, 1.02),
        labels = scales::percent_format(accuracy = 1)
    ) +
    # block outside, dataset inside, and not the other way round: the rarity
    # block of 16S is a single row, and a rotated strip label cannot fit beside
    # one row. With the block outermost its strip spans all nine rare rows.
    ggh4x::facet_nested(block + dataset ~ ., scales = "free_y", space = "free_y") +
    labs(
        title    = "Curated ecologically relevant taxa — detection in cacao flowers",
        subtitle = paste0(
            ROLE_LEGEND, "<br>",
            "Row order differs between blocks: by difference between treatments ",
            "above, by number of detections below"),
        caption  = paste(MARK_CAPTION, BLOCK_CAPTION, sep = "\n"),
        x = "flowers with the taxon detected", y = NULL
    ) +
    theme(strip.text.y = element_text(size = 7))

ggsave(here("results", "figures", "00_curated_taxa_prevalence.png"),
       plot = p_dot, width = 9, height = 8, dpi = 300)

################################################################################
# Fig 1b — the same grammar, one column per farm
#
# Row order is inherited from Fig 1 rather than recomputed inside each farm:
# otherwise a taxon changes line between columns and cannot be tracked.
################################################################################

farm_mgmt <- curated_long |> distinct(farm_id, management_type)

prev_farm <- prev_farm_raw |>
    left_join(prev_genus |> as_tibble() |>
                  select(dataset, Genus_clean, label, comparable),
              by = c("dataset", "Genus_clean")) |>
    # drop the genera omitted from Fig 1, and the rarity block with them: this
    # figure is nothing but the contrast, so a row that cannot carry one has no
    # business here. Per-farm cells are even smaller (n = 21), which only makes
    # the point sharper.
    filter(!is.na(label), comparable) |>
    left_join(farm_mgmt, by = "farm_id") |>
    mutate(
        dataset         = factor(dataset, levels = c("ITS1", "16S")),
        management_type = factor(management_type, levels = MGMT_LEVELS),
        farm_id         = factor(farm_id, levels = FARM_LEVELS)
    )

p_dot_farm <- ggplot(prev_farm) +
    dumbbell_layers(point_size = 1.8) +
    scale_x_continuous(
        breaks = c(0, 0.5, 1),              # fewer breaks: narrow farm columns
        limits = c(-0.02, 1.02),
        labels = scales::percent_format(accuracy = 1)
    ) +
    ggh4x::facet_nested(dataset ~ management_type + farm_id,
                        scales = "free_y", space = "free_y") +
    labs(
        title    = "Curated taxa — detection per farm",
        subtitle = paste0(ROLE_LEGEND, "<br>",
                          "Row order as in the pooled figure; ",
                          "farms in the canonical order"),
        caption  = paste(
            MARK_CAPTION,
            sprintf(paste0("genera found in fewer than %d flowers are omitted here ",
                           "and reported as detection only in 00_"), DIFF_MIN_N),
            sep = "\n"),
        x = "flowers with the taxon detected within the farm", y = NULL
    ) +
    theme(axis.text.x  = element_text(size = 6),
          strip.text.x = element_text(size = 7.5))

ggsave(here("results", "figures", "00c_curated_taxa_farm.png"),
       plot = p_dot_farm, width = 20, height = 9, dpi = 300)

################################################################################
# Fig 1c — ASV level, only where the genus hides a disagreement
#
# A full ASV panel would be unreadable: ~98 ITS1 OTUs and 20 16S ASVs clear the
# detection threshold, and 65 of them belong to Hannaella alone. Only the genera
# whose ASVs point both ways are shown, and within each only the eight with the
# largest difference; the rest are summarised in the panel strip.
################################################################################

ASV_PER_GENUS <- 8

asv_mixed <- asv_kept |>
    semi_join(mixed_genera |> filter(asv_flip), by = c("dataset", "Genus_clean")) |>
    group_by(dataset, Genus_clean) |>
    mutate(rank = rank(-abs(diff_rate), ties.method = "first"),
           n_omitted_up   = sum(diff_rate > 0 & rank > ASV_PER_GENUS),
           n_omitted_down = sum(diff_rate < 0 & rank > ASV_PER_GENUS)) |>
    ungroup() |>
    filter(rank <= ASV_PER_GENUS) |>
    left_join(genus_roles |> select(Genus_clean = Genus, role = primary_role),
              by = "Genus_clean") |>
    # the role now feeds a colour inside the label string, so a join miss would
    # produce <span style='color:NA'> rather than a visibly missing swatch
    (\(d) { stopifnot(!any(is.na(d$role))); d })() |>
    mutate(
        # feature hash disambiguates ASVs sharing a genus name, as in 16d and the
        # random-forest importance figures
        label_plain = sprintf("%s [%s]", Genus_clean, substr(ASV, 1, 6)),
        label  = md_label(label_plain, role),
        facet  = sprintf("%s — %s (%s)", dataset, Genus_clean,
                         if_else(n_omitted_up + n_omitted_down > 0,
                                 sprintf("%d more up, %d more down",
                                         n_omitted_up, n_omitted_down),
                                 "all shown"))
    )
asv_mixed <- asv_mixed |>
    mutate(label = factor(label, levels = asv_mixed |> arrange(diff_rate) |>
                              pull(label) |> unique()))

p_asv <- ggplot(asv_mixed) +
    dumbbell_layers(point_size = 2.4) +
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       limits = c(-0.02, 1.02),
                       labels = scales::percent_format(accuracy = 1)) +
    facet_wrap(~ facet, scales = "free_y", ncol = 2) +
    labs(
        title    = "Where a genus hides a disagreement: ASV-level detection",
        # "at least" and not ">=": the subtitle is markdown now, and a leading
        # > would be read as a blockquote
        subtitle = paste0(ROLE_LEGEND, "<br>",
                          sprintf(paste0("genera marked \u2021; up to %d ASVs each, ",
                                         "detected in at least %d flowers"),
                                  ASV_PER_GENUS, MIN_DETECT)),
        caption  = "full ASV table in results/tables/00d_curated_taxa_asv.csv",
        x = "flowers with the ASV detected", y = NULL
    ) +
    theme(strip.text = element_text(size = 8, face = "bold"))

ggsave(here("results", "figures", "00d_curated_taxa_asv_mixed.png"),
       plot = p_asv, width = 11, height = 8, dpi = 300)

################################################################################
# Fig 2 — heatmap: how much of the taxon is there WHEN it is there
#
# The median over all flowers was zero for any taxon present in fewer than half of
# them, which left 81% of the grid grey and threw away the abundance information
# entirely. Taking the median over the flowers where the taxon was detected answers
# "how much, when present" and leaves grey to mean "absent from that farm x
# treatment altogether" -- 31% of cells. Prevalence, the other half of the question,
# is what the dumbbell figures show, so the two no longer duplicate each other.
#
# magma is perceptually uniform and colour-blind safe; direction = -1 puts high
# abundance at the dark end, which reads better on a white ground.
################################################################################

heat_tbl <- curated_long |>
    group_by(dataset, Genus_clean, farm_id, sample_type) |>
    summarise(med_RA = if (any(RA > 0)) median(RA[RA > 0]) else NA_real_,
              n_pos  = sum(RA > 0), .groups = "drop") |>
    left_join(prev_genus |> as_tibble() |>
                  select(dataset, Genus_clean, label, block),
              by = c("dataset", "Genus_clean")) |>
    filter(!is.na(label)) |>
    left_join(farm_mgmt, by = "farm_id") |>
    mutate(
        log_RA          = log10(med_RA),
        farm_id         = factor(farm_id, levels = FARM_LEVELS),
        management_type = factor(management_type, levels = MGMT_LEVELS),
        # short labels: they repeat under all seven farms
        sample_type     = factor(sample_type,
                                 levels = c("bagged_flower", "unbagged_flower"),
                                 labels = c("bgd", "unbgd")),
        dataset         = factor(dataset, levels = c("ITS1", "16S"))
    )

p_heat <- ggplot(heat_tbl, aes(x = sample_type, y = label, fill = log_RA)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    # Treatment is the x variable INSIDE each farm panel, so the two cells of one
    # farm are neighbours and the comparison needs no memory. Faceting by treatment
    # instead -- the previous layout -- put them seven columns and a strip apart,
    # and colour cannot be carried that far by eye. The upper strip restores the
    # management gradient, matching the nested strips already used in 00c.
    # all 27 genera stay: abundance-when-present is meaningful however rare the
    # taxon is. The block only groups the rows, so that the row order agrees with
    # 00_ instead of looking like an unexplained reshuffle.
    ggh4x::facet_nested(block + dataset ~ management_type + farm_id,
                        scales = "free_y", space = "free_y") +
    scale_fill_viridis_c(
        option    = "magma",
        direction = -1,                 # high abundance = dark
        na.value  = "grey90",
        name      = expression(log[10]~"median RA")
    ) +
    labs(
        title    = "Curated taxa — relative abundance where the taxon is present",
        # <br> and not "\n": the subtitle is markdown here, where a bare newline
        # opens a new paragraph and leaves a gap
        subtitle = paste0(
            ROLE_LEGEND, "<br>",
            "rows grouped as in 00_, where the lower block is too rare to carry ",
            "a treatment contrast -- which says nothing about the abundances here",
            "<br>",
            "median across the flowers in which it was detected; the two cells of ",
            "a farm sit side by side; grey = absent from that farm x treatment"),
        x = NULL, y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
        strip.background   = element_rect(fill = "grey92"),
        strip.text.x       = element_text(size = 7.5),
        panel.grid         = element_blank(),
        # tight inside a farm, a visible gap between farms, so each pair reads as a unit
        panel.spacing.x    = unit(0.25, "lines"),
        legend.position    = "right",
        axis.text.x        = element_text(size = 7),
        # p_heat does not go through dumbbell_layers(), so it needs the two
        # markdown elements of its own
        axis.text.y        = element_markdown(),
        plot.subtitle      = element_markdown()
    )

ggsave(here("results", "figures", "00b_curated_heatmap.png"),
       plot = p_heat, width = 13, height = 9, dpi = 300)

################################################################################
# Console summary
################################################################################

message("00_curated_taxa_check done")
message(sprintf("  ITS1 detected: %s",
                paste(sort(unique(hits_ITS1$Genus_clean)), collapse = ", ")))
message(sprintf("  16S  detected: %s",
                paste(sort(unique(hits_16S$Genus_clean)),  collapse = ", ")))
message(sprintf(paste0("  figures: %d genus rows -- %d with the contrast shown, ",
                       "%d detection only (k < %d); %d marked \u2021, %d marked *; ",
                       "%d ASV rows"),
                nrow(prev_genus),
                sum(prev_genus$comparable), sum(!prev_genus$comparable), DIFF_MIN_N,
                sum(prev_genus$asv_flip  & prev_genus$comparable),
                sum(prev_genus$farm_flip & prev_genus$comparable),
                nrow(asv_mixed)))
message(sprintf("  detection only: %s",
                paste(sort(prev_genus$label_plain[!prev_genus$comparable]),
                      collapse = ", ")))
message(sprintf("  heatmap: %d of %d cells grey (taxon absent from that farm x treatment)",
                sum(is.na(heat_tbl$med_RA)), nrow(heat_tbl)))
