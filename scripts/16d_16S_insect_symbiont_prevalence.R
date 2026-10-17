## ============================================================================
## 16d — Insect-symbiont prevalence: bagged vs unbagged flowers
##
## Why a separate analysis rather than another ANCOM-BC2 run.
##
## g__Wolbachia (ASV 9bd6b381) is higher in bagged flowers in all seven farms in
## 16c, yet the pooled model (16a) reports q = 0.145 AND passed_ss = FALSE. Its raw
## pattern explains both: the median is zero in BOTH groups, it is detected in 44 of
## 147 bagged and 18 of 147 unbagged flowers, and where present in bagged flowers it
## can reach 41% of the community. ANCOM-BC2 models the mean log-ratio, so for a
## taxon like this the estimate hinges on how zeros are treated -- exactly what the
## pseudo-count sensitivity check flags, and exactly what it fails.
##
## The signal is therefore about PRESENCE, not abundance, and is tested as such here:
## a within-tree permutation test on detection. That respects the same blocking the
## pooled models express as (1|tree_id), assumes no distribution, and handles a rare
## binary outcome honestly.
##
## The hypothesis under test is the user's: bags do not exclude very small insects,
## so trapped insects spend longer in contact with the flower and leave more of their
## endosymbionts behind. Wolbachia is an obligate intracellular endosymbiont and
## cannot grow on a flower surface, so its DNA must come from insect tissue -- which
## is what makes presence a meaningful marker of insect material. Testing the whole
## symbiont guild rather than Wolbachia alone is what separates "more insects under
## the bag" from "different insects under the bag".
##
## Input:
##   results/rds/ps_16S_bacteria_biosamples.rds   (post-decontam, bacteria)
##   results/rds/ps_16S_raw.rds                   (pre-decontam, for the controls)
##   qiime2/export/CFM_16S_dada2_repseqs/dna-sequences.fasta
## Output:
##   results/figures/16d_insect_symbiont_prevalence.png
##   results/tables/16d_insect_symbiont_prevalence.csv
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

# Genera with a documented insect-symbiont or insect-gut association. Reproductive
# manipulators and obligate endosymbionts first, then bee/flower-visitor gut taxa.
GUILD_GENERA <- c("Wolbachia", "Spiroplasma", "Rickettsia", "Arsenophonus",
                  "Cardinium", "Sodalis", "Buchnera",
                  "Serratia", "Lactobacillus", "Fructobacillus",
                  "Gilliamella", "Snodgrassella", "Orbus", "Apibacter", "Bartonella")

# Detection threshold. Deliberately looser than the 10% prevalence used elsewhere in
# the project: here rarity IS the phenomenon, and a 10% cut (>= 30 of 294 flowers)
# would discard the second Wolbachia ASV, whose contrast with the first carries the
# most information about the hypothesis.
MIN_DETECT <- 10

MIN_LIB    <- 500    # same library-size floor as 16a/16b
N_PERM     <- 9999   # within-tree label permutations
FDR_ALPHA  <- 0.05
SEED       <- 42

# ── 1. Libraries & data ───────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(patchwork)

# Taxonomy labels: deepest assigned rank, prefix kept, shared with 16a-18b.
source(here("scripts", "utils_tax_label.R"))

ps      <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds"))
ps_pre  <- readRDS(here("results", "rds", "ps_16S_raw.rds"))

cat(sprintf("Loaded: %d taxa x %d samples (post-decontam)\n", ntaxa(ps), nsamples(ps)))

set.seed(SEED)

# ── 2. Helper functions ───────────────────────────────────────────────────────

# Counts as taxa x samples, whichever orientation the object happens to use.
count_matrix <- function(physeq) {
        m <- as(otu_table(physeq), "matrix")
        if (taxa_are_rows(physeq) == FALSE) m <- t(m)
        m
}

# Wilson interval via prop.test(correct = FALSE); base R, no extra dependency.
wilson_ci <- function(x, n) {
        if (n == 0) return(c(NA_real_, NA_real_))
        as.numeric(prop.test(x, n, correct = FALSE)$conf.int)
}

# ── 3. Flowers, guild, detection matrix ───────────────────────────────────────

cat("\n── 3. Guild selection ───────────────────────────────────────────────────\n")

meta <- as(sample_data(ps), "data.frame")
flowers <- rownames(meta)[meta$sample_type %in% c("bagged_flower", "unbagged_flower") &
                          sample_sums(ps)[rownames(meta)] >= MIN_LIB]

meta_fl <- meta[flowers, ]
is_bag  <- meta_fl$sample_type == "bagged_flower"
tree    <- factor(meta_fl$tree_id)
counts  <- count_matrix(ps)[, flowers, drop = FALSE]

cat(sprintf("  flowers: %d (%d bagged, %d unbagged) on %d trees\n",
            length(flowers), sum(is_bag), sum(is_bag == FALSE), nlevels(tree)))

# Sequencing depth is a nuisance for ANY prevalence analysis: a shallower library
# detects fewer taxa regardless of biology. Reported here, and carried into the GLMM
# below as log10(library size).
depth <- colSums(counts)
cat(sprintf("  library size, median: bagged %.0f  unbagged %.0f  (Wilcoxon p = %.3f)\n",
            median(depth[is_bag]), median(depth[is_bag == FALSE]),
            wilcox.test(depth ~ is_bag)$p.value))

tax   <- as.data.frame(tax_table(ps))
guild <- rownames(tax)[grepl(paste(GUILD_GENERA, collapse = "|"), tax$Genus)]
det   <- counts[guild, , drop = FALSE] > 0
keep  <- rowSums(det) >= MIN_DETECT

cat(sprintf("  guild ASVs: %d matched the genus list, %d detected in >= %d flowers\n",
            length(guild), sum(keep), MIN_DETECT))

det   <- det[keep, , drop = FALSE]
guild <- rownames(det)

# ── 4. Within-tree permutation test on detection ──────────────────────────────

cat("\n── 4. Within-tree permutation test ──────────────────────────────────────\n")

# Statistic: (flowers with the taxon among bagged) - (among unbagged). Under the null
# the treatment label is exchangeable WITHIN a tree, where 3 bagged and 3 unbagged
# flowers were collected; permuting across trees would destroy the blocking that
# 16a/16b encode as (1|tree_id).
#
# One sign matrix is drawn and reused for every ASV. That is deliberate: the taxa
# share the permutation, so their null statistics keep the correlation the observed
# data have, which is the behaviour Benjamini-Hochberg expects.
sign_obs  <- ifelse(is_bag, 1, -1)
sign_null <- replicate(N_PERM, {
        s <- numeric(length(flowers))
        for (tr in levels(tree)) {
                idx    <- which(tree == tr)
                s[idx] <- sample(sign_obs[idx])   # shuffle labels inside the tree
        }
        s
})

perm_p <- function(d) {
        obs  <- sum(d * sign_obs)
        null <- as.vector(d %*% sign_null)
        # two-sided, with the 1/(N+1) floor that a Monte Carlo p-value cannot go below
        (1 + sum(abs(null) >= abs(obs))) / (N_PERM + 1)
}

res <- tibble::tibble(
        taxon     = guild,
        genus     = tax[guild, "Genus"],
        tax_label = make_tax_label(ps, guild),
        n_bagged  = as.integer(rowSums(det[, is_bag, drop = FALSE])),
        n_unbag   = as.integer(rowSums(det[, is_bag == FALSE, drop = FALSE]))
) |>
        dplyr::mutate(
                rate_bagged = n_bagged / sum(is_bag),
                rate_unbag  = n_unbag  / sum(is_bag == FALSE),
                diff_rate   = rate_bagged - rate_unbag,
                p_perm      = vapply(guild, function(i) perm_p(det[i, ]), numeric(1)),
                q_perm      = p.adjust(p_perm, method = "BH"),
                signif      = q_perm < FDR_ALPHA
        )

# Wilson intervals for the figure
ci <- t(mapply(wilson_ci, res$n_bagged, sum(is_bag)))
res$bag_lo <- ci[, 1]; res$bag_hi <- ci[, 2]
ci <- t(mapply(wilson_ci, res$n_unbag, sum(is_bag == FALSE)))
res$unb_lo <- ci[, 1]; res$unb_hi <- ci[, 2]

cat(sprintf("  tested %d ASVs, %d significant at BH q < %.2f\n",
            nrow(res), sum(res$signif), FDR_ALPHA))

# ── 5. Effect size: binomial GLMM with the same blocking ──────────────────────

cat("\n── 5. Binomial GLMM (odds ratios) ───────────────────────────────────────\n")

# Mirrors the 16a/16b design on a binary response, plus sequencing depth. Rare taxa
# make these fits fragile; failures are counted and reported rather than dropped.
glmm_or <- function(d) {
        fit <- tryCatch(
                suppressMessages(suppressWarnings(lme4::glmer(
                        d ~ sample_type + farm_id + log10(depth) + (1 | tree_id),
                        data = cbind(meta_fl, depth = depth), family = binomial))),
                error = function(e) NULL)
        if (is.null(fit)) return(NA_real_)
        cf <- lme4::fixef(fit)
        i  <- grep("^sample_typeunbagged", names(cf))
        # model contrast is unbagged vs bagged; invert so the OR reads bagged vs unbagged
        if (length(i) == 1) exp(-cf[i]) else NA_real_
}

res$or_bagged <- vapply(guild, function(i) glmm_or(det[i, ]), numeric(1))
cat(sprintf("  GLMM converged for %d of %d ASVs\n",
            sum(is.na(res$or_bagged) == FALSE), nrow(res)))

# ── 6. Genus-level aggregate ──────────────────────────────────────────────────

cat("\n── 6. Genus-level aggregate ─────────────────────────────────────────────\n")

# Shown deliberately: the two Wolbachia ASVs point in opposite directions and very
# nearly cancel once merged, so a genus-level analysis would have reported nothing.
det_genus <- rowsum(det * 1, res$genus) > 0
res_genus <- tibble::tibble(
        genus    = rownames(det_genus),
        n_asv    = as.integer(table(res$genus)[rownames(det_genus)]),
        n_bagged = as.integer(rowSums(det_genus[, is_bag, drop = FALSE])),
        n_unbag  = as.integer(rowSums(det_genus[, is_bag == FALSE, drop = FALSE]))
) |>
        dplyr::mutate(rate_bagged = n_bagged / sum(is_bag),
                      rate_unbag  = n_unbag  / sum(is_bag == FALSE),
                      p_perm      = vapply(genus,
                                           function(g) perm_p(det_genus[g, ]), numeric(1)),
                      q_perm      = p.adjust(p_perm, method = "BH"))

for (i in seq_len(nrow(res_genus)))
        cat(sprintf("  %-18s %d ASV  bagged %3d  unbagged %3d  q = %.3f\n",
                    res_genus$genus[i], res_genus$n_asv[i], res_genus$n_bagged[i],
                    res_genus$n_unbag[i], res_genus$q_perm[i]))

# ── 7. Are the two Wolbachia ASVs relatives? ──────────────────────────────────

cat("\n── 7. Wolbachia ASVs ────────────────────────────────────────────────────\n")

wol <- res$taxon[grepl("Wolbachia", res$genus)]
if (length(wol) == 2) {
        cat(sprintf("  co-occurrence: %d flowers carry both, %d carry only one\n",
                    sum(det[wol[1], ] & det[wol[2], ]),
                    sum(xor(det[wol[1], ], det[wol[2], ]))))

        d_tree <- ape::keep.tip(phy_tree(ps), wol) |> cophenetic()
        cat(sprintf("  patristic distance on the 16S tree: %.4f substitutions/site\n",
                    d_tree[1, 2]))

        fa <- here("qiime2", "export", "CFM_16S_dada2_repseqs", "dna-sequences.fasta")
        if (file.exists(fa)) {
                sq <- Biostrings::readDNAStringSet(fa)[wol]
                if (length(unique(Biostrings::width(sq))) == 1) {
                        ch <- lapply(as.character(sq), function(s) strsplit(s, "")[[1]])
                        cat(sprintf("  sequence identity: %.1f%% (%d/%d positions match)\n",
                                    100 * mean(ch[[1]] == ch[[2]]),
                                    sum(ch[[1]] == ch[[2]]), length(ch[[1]])))
                } else {
                        # unequal length: a proper alignment would be needed, not attempted here
                        cat(sprintf("  sequence lengths differ (%s) -- identity not computed\n",
                                    paste(Biostrings::width(sq), collapse = " vs ")))
                }
        }
} else {
        cat(sprintf("  %d Wolbachia ASVs passed the detection threshold\n", length(wol)))
}

# ── 8. Confounding checks ─────────────────────────────────────────────────────

cat("\n── 8. Confounding checks ────────────────────────────────────────────────\n")

# (a) Controls. Reagent contamination is the standard alternative explanation for any
# sporadic taxon; Lactobacillus and Serratia in particular turn up in kit reagents.
meta_pre <- as(sample_data(ps_pre), "data.frame")
ctrl     <- rownames(meta_pre)[meta_pre$sample_type %in%
                               c("extraction_control", "negative_pcr_control",
                                 "mock_community")]
det_ctrl <- count_matrix(ps_pre)[intersect(guild, taxa_names(ps_pre)), ctrl,
                                 drop = FALSE] > 0
res$n_controls <- rowSums(det_ctrl)[res$taxon]
cat(sprintf("  controls: %d samples; guild ASVs detected in any control: %d of %d\n",
            length(ctrl), sum(res$n_controls > 0, na.rm = TRUE), nrow(res)))

# (b) Extraction order. Bagged flowers were extracted before unbagged on every farm,
# so treatment and extraction batch are confounded everywhere EXCEPT one batch that
# happens to contain both. Twelve flowers settle nothing, but hiding the one
# unconfounded slice would be worse than showing how small it is.
shared <- names(which(table(meta_fl$extraction_batch, meta_fl$sample_type)[, 1] > 0 &
                      table(meta_fl$extraction_batch, meta_fl$sample_type)[, 2] > 0))
in_sh  <- meta_fl$extraction_batch %in% shared
cat(sprintf("  extraction batches containing both treatments: %s (n = %d flowers)\n",
            paste(shared, collapse = ", "), sum(in_sh)))
res$sh_bagged <- as.integer(rowSums(det[, in_sh & is_bag, drop = FALSE]))
res$sh_unbag  <- as.integer(rowSums(det[, in_sh & (is_bag == FALSE), drop = FALSE]))

# (c) Flowers where an insect was actually seen during microscopy.
ins <- rownames(meta_fl)[meta_fl$is_insect_present %in% c(TRUE, "TRUE")]
cat(sprintf("  flowers with a visible insect: %d (%s); guild ASVs detected in them: %d\n",
            length(ins), paste(ins, collapse = ", "),
            if (length(ins) > 0) sum(rowSums(det[, ins, drop = FALSE]) > 0) else 0L))

# ── 9. Figure ─────────────────────────────────────────────────────────────────

cat("\n── 9. Building figure ───────────────────────────────────────────────────\n")

GROUP_COLS <- c(bagged = "#4C72B0", unbagged = "#DD8452")   # as in 10d and 09d

# Dumbbell: one row per ASV, a segment from the unbagged detection rate to the bagged
# one. Rows are ordered by the difference, so the direction of the guild is readable
# at a glance rather than taxon by taxon.
dumbbell <- function(d, ylab_col, title, subtitle) {
        d <- dplyr::mutate(d, row = reorder(.data[[ylab_col]], rate_bagged - rate_unbag))
        ggplot(d) +
                geom_segment(aes(y = row, yend = row, x = rate_unbag, xend = rate_bagged),
                             colour = "grey70", linewidth = 0.6) +
                geom_point(aes(y = row, x = rate_unbag, colour = "unbagged"), size = 2.6) +
                geom_point(aes(y = row, x = rate_bagged, colour = "bagged"), size = 2.6) +
                geom_text(aes(y = row, x = pmax(rate_bagged, rate_unbag) + 0.04,
                              label = sprintf("q = %.3f", q_perm)),
                          hjust = 0, size = 2.7, colour = "grey30") +
                scale_colour_manual(values = GROUP_COLS, name = NULL) +
                # headroom on the right for the q-value labels
                scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                   expand = expansion(mult = c(0.02, 0.28))) +
                labs(title = title, subtitle = subtitle,
                     x = "flowers with the taxon detected", y = NULL) +
                theme_bw(base_size = 10) +
                theme(legend.position   = "bottom",
                      plot.title        = element_text(size = 11, face = "bold"),
                      plot.subtitle     = element_text(size = 8, colour = "grey40"),
                      panel.grid.major.y = element_blank())
}

# The feature hash disambiguates the row labels: the two Wolbachia ASVs share a genus
# name and would otherwise be drawn on top of each other, which is the whole point of
# the panel. Same convention as the random-forest importance figures.
p_asv <- dumbbell(dplyr::mutate(res, tax_label = sprintf("%s [%s]", tax_label,
                                                         substr(taxon, 1, 6))),
                  "tax_label", "Per ASV",
                  sprintf("within-tree permutation test, %d shuffles  |  BH across %d ASVs",
                          N_PERM, nrow(res)))
p_gen <- dumbbell(dplyr::mutate(res_genus, tax_label = genus), "tax_label",
                  "Aggregated by genus",
                  "the two Wolbachia ASVs run in opposite directions and cancel here")

fig <- (p_asv / p_gen) +
        plot_layout(guides = "collect", heights = c(nrow(res), nrow(res_genus))) +
        plot_annotation(
                title   = "Insect-symbiont detection in bagged vs unbagged flowers",
                caption = sprintf(paste("detection = at least one read  |  ASVs detected in >= %d of %d flowers",
                                        "|  labels permuted within tree (3 bagged + 3 unbagged each)"),
                                  MIN_DETECT, length(flowers)),
                theme = theme(plot.title   = element_text(size = 13, face = "bold"),
                              plot.caption = element_text(size = 7, colour = "grey40"))
        ) &
        theme(legend.position = "bottom")

# ── 10. Save outputs ──────────────────────────────────────────────────────────

cat("\n── 10. Saving outputs ───────────────────────────────────────────────────\n")

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

ggsave(here("results", "figures", "16d_insect_symbiont_prevalence.png"),
       plot = fig, width = 9, height = 7, dpi = 300)
cat("  Saved: 16d_insect_symbiont_prevalence.png\n")

res |>
        dplyr::select(taxon, genus, tax_label, n_bagged, n_unbag,
                      rate_bagged, rate_unbag, diff_rate,
                      p_perm, q_perm, signif, or_bagged,
                      n_controls, sh_bagged, sh_unbag) |>
        dplyr::arrange(p_perm) |>
        readr::write_csv(here("results", "tables",
                              "16d_insect_symbiont_prevalence.csv"))
cat("  Saved: 16d_insect_symbiont_prevalence.csv\n")

cat(sprintf("\n16d done. %d of %d guild ASVs differ in detection (BH q < %.2f).\n",
            sum(res$signif), nrow(res), FDR_ALPHA))
