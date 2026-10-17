## ============================================================================
## 16e — Insect-symbiont guild under MaAsLin 3: a second engine for 16d
##
## 16d tested the guild with a within-tree permutation test on detection and found
## three ASVs differing between bagged and unbagged flowers: two Wolbachia variants
## in OPPOSITE directions (44/18 and 7/22) and one Lactobacillus (3/15). That test is
## exact for this blocked design, but it is a single hand-rolled engine, and it asks
## only about presence.
##
## MaAsLin 3 fits two models per feature and reports them separately:
##   prevalence -- logistic, presence/absence, tested against zero
##   abundance  -- linear on log(TSS), among samples where the feature is present
## The second is new information: conditional on Wolbachia being there at all, is it
## more abundant in bagged flowers? 16d cannot answer that, and ANCOM-BC2 could not
## either, because the zeros dominated its log-ratio estimate.
##
## WHY THE MODEL IS FITTED ON MORE FEATURES THAN IT REPORTS
##
## MaAsLin 3 defends the abundance model against compositionality by testing each
## coefficient against the MEDIAN coefficient across features (median_comparison_-
## abundance = TRUE), not against zero. That null only means something when the
## feature set is large enough for the median to represent "most of the community
## did not change". On the seven guild ASVs alone the median would be the median of
## seven numbers, which is not a null at all.
##
## So the model is fitted on the guild PLUS the 10%-prevalence background that 16a
## and 16b used (~235 features), and Benjamini-Hochberg is then applied to the seven
## pre-specified guild ASVs only. The median is a nuisance parameter and is better
## estimated from the whole community; the family of hypotheses is the pre-specified
## guild. Both q-values are reported so the difference is visible, not hidden.
##
## The permutation test in 16d remains the primary result: with detections as low as
## 3 versus 15 across 49 blocks, a Wald test leans on asymptotics that do not hold
## here. This script is corroboration plus the abundance half.
##
## Input:
##   results/rds/ps_16S_bacteria_biosamples.rds
##   results/tables/16d_insect_symbiont_prevalence.csv   (permutation results)
## Output:
##   results/figures/16e_maaslin3_vs_permutation.png
##   results/tables/16e_maaslin3_vs_permutation.csv
##   results/tables/16e_maaslin3/                        (raw package output)
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

# Guild definition -- KEEP IN SYNC WITH 16d. Duplicated rather than shared so each
# script stays runnable on its own, which is the convention in this project.
GUILD_GENERA <- c("Wolbachia", "Spiroplasma", "Rickettsia", "Arsenophonus",
                  "Cardinium", "Sodalis", "Buchnera",
                  "Serratia", "Lactobacillus", "Fructobacillus",
                  "Gilliamella", "Snodgrassella", "Orbus", "Apibacter", "Bartonella")
MIN_DETECT  <- 10     # guild ASVs: detected in >= this many flowers (as in 16d)

BACKGROUND_PREV <- 0.10   # background for the compositional median: the same 10%
                          # prevalence filter 16a/16b applied (231 ASVs)
MIN_LIB    <- 500
FDR_ALPHA  <- 0.05

# ── 1. Libraries & data ───────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(maaslin3)

source(here("scripts", "utils_tax_label.R"))

ps   <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds"))
perm <- readr::read_csv(here("results", "tables",
                             "16d_insect_symbiont_prevalence.csv"),
                        show_col_types = FALSE)

cat(sprintf("Loaded: %d taxa x %d samples | maaslin3 %s\n",
            ntaxa(ps), nsamples(ps), as.character(packageVersion("maaslin3"))))

# ── 2. Samples and feature set ────────────────────────────────────────────────

cat("\n── 2. Feature set ───────────────────────────────────────────────────────\n")

meta <- as(sample_data(ps), "data.frame")
flowers <- rownames(meta)[meta$sample_type %in% c("bagged_flower", "unbagged_flower") &
                          sample_sums(ps)[rownames(meta)] >= MIN_LIB]

counts <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps) == FALSE) counts <- t(counts)
counts <- counts[, flowers, drop = FALSE]

n_det <- rowSums(counts > 0)
tax   <- as.data.frame(tax_table(ps))

guild <- rownames(tax)[grepl(paste(GUILD_GENERA, collapse = "|"), tax$Genus)]
guild <- guild[n_det[guild] >= MIN_DETECT]

background <- names(n_det)[n_det >= BACKGROUND_PREV * length(flowers)]
fit_set    <- union(background, guild)

cat(sprintf("  flowers: %d | guild ASVs: %d | background ASVs: %d | fitted: %d\n",
            length(flowers), length(guild), length(background), length(fit_set)))
stopifnot(all(guild %in% fit_set))

# ── 3. Run MaAsLin 3 ──────────────────────────────────────────────────────────

cat("\n── 3. Running MaAsLin 3 ─────────────────────────────────────────────────\n")

# Features as columns, samples as rows, as the package expects.
input_data <- as.data.frame(t(counts[fit_set, , drop = FALSE]))

input_meta <- meta[flowers, c("sample_type", "farm_id", "tree_id")]
# unbagged is the reference level, so a POSITIVE coefficient reads "higher in bagged"
# -- the same direction as diff_rate and or_bagged in 16d, so the tables can be laid
# side by side without mentally inverting anything.
input_meta$sample_type <- factor(input_meta$sample_type,
                                 levels = c("unbagged_flower", "bagged_flower"))

out_dir <- here("results", "tables", "16e_maaslin3")

# Sequencing depth is not added as a covariate: MaAsLin 3 normalises by total sum
# scaling, and 16d showed the two groups do not differ in library size
# (medians 10447 vs 12343, Wilcoxon p = 0.62).
invisible(maaslin3(
        input_data     = input_data,
        input_metadata = input_meta,
        output         = out_dir,
        fixed_effects  = c("sample_type", "farm_id"),
        random_effects = "tree_id",          # 3 bagged + 3 unbagged per tree
        reference      = "farm_id,ib",       # deterministic reference for 7 farms
        max_significance  = FDR_ALPHA,       # package default is 0.1
        plot_summary_plot = FALSE,           # plots would cover all ~235 features
        plot_associations = FALSE,
        cores          = 1,
        verbosity      = "WARN"
))

# ── 4. Restrict to the pre-specified guild and re-correct ─────────────────────

cat("\n── 4. Guild results ─────────────────────────────────────────────────────\n")

res <- readr::read_tsv(file.path(out_dir, "all_results.tsv"), show_col_types = FALSE) |>
        dplyr::filter(metadata == "sample_type", feature %in% guild) |>
        dplyr::select(taxon = feature, model, coef, stderr, null_hypothesis,
                      p_maaslin = pval_individual,
                      q_background = qval_individual,   # BH over all fitted features
                      n_nonzero = N_not_zero) |>
        # BH within the seven pre-specified ASVs -- the family chosen for this analysis
        dplyr::group_by(model) |>
        dplyr::mutate(q_guild = p.adjust(p_maaslin, method = "BH")) |>
        dplyr::ungroup()

stopifnot(dplyr::n_distinct(res$taxon[res$model == "prevalence"]) == length(guild))

cmp <- res |>
        dplyr::left_join(dplyr::select(perm, taxon, tax_label, n_bagged, n_unbag,
                                       q_perm, or_bagged), by = "taxon") |>
        dplyr::mutate(
                label = sprintf("%s [%s]", tax_label, substr(taxon, 1, 6)),
                # 16d's odds ratio is bagged vs unbagged; a coefficient above zero
                # means the same thing, so the two must agree in sign
                agrees = sign(coef) == sign(log(or_bagged))
        )

for (mdl in c("prevalence", "abundance")) {
        d <- dplyr::filter(cmp, model == mdl) |> dplyr::arrange(p_maaslin)
        cat(sprintf("\n  --- %s model ---\n", mdl))
        for (i in seq_len(nrow(d)))
                cat(sprintf("  %-34s coef %+6.2f  q_guild %6.4f  q_backgr %6.4f  n>0 %3d%s\n",
                            d$label[i], d$coef[i], d$q_guild[i], d$q_background[i],
                            d$n_nonzero[i],
                            if (mdl == "prevalence" && isFALSE(d$agrees[i]))
                                    "  <- DISAGREES with 16d" else ""))
}

cat(sprintf("\n  significant within the guild (q < %.2f): prevalence %d, abundance %d\n",
            FDR_ALPHA,
            sum(cmp$q_guild[cmp$model == "prevalence"] < FDR_ALPHA, na.rm = TRUE),
            sum(cmp$q_guild[cmp$model == "abundance"]  < FDR_ALPHA, na.rm = TRUE)))

# Features MaAsLin 3 could not fit are reported, not silently absent.
miss <- setdiff(paste(rep(guild, 2), rep(c("prevalence", "abundance"), each = length(guild))),
                paste(res$taxon, res$model))
if (length(miss) > 0)
        cat(sprintf("  models not returned (too few non-zero samples): %s\n",
                    paste(miss, collapse = ", ")))

# ── 5. Figure ─────────────────────────────────────────────────────────────────

cat("\n── 5. Building figure ───────────────────────────────────────────────────\n")

# Rows share one order across both panels so a taxon can be tracked between them.
row_order <- cmp |>
        dplyr::filter(model == "prevalence") |>
        dplyr::arrange(coef) |>
        dplyr::pull(label)

panel <- function(mdl, title, subtitle, show_perm) {
        d <- cmp |>
                dplyr::filter(model == mdl) |>
                dplyr::mutate(label  = factor(label, levels = row_order),
                              signif = q_guild < FDR_ALPHA)
        # headroom on the right for the q labels; symmetric so zero stays centred
        rng <- max(abs(c(d$coef - 1.96 * d$stderr, d$coef + 1.96 * d$stderr)), na.rm = TRUE)
        ggplot(d, aes(y = label, x = coef, colour = signif)) +
                geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
                geom_errorbarh(aes(xmin = coef - 1.96 * stderr,
                                   xmax = coef + 1.96 * stderr),
                               height = 0, linewidth = 0.6) +
                geom_point(size = 2.6) +
                geom_text(aes(x = rng * 1.15,
                              label = if (show_perm)
                                      sprintf("q %.3f | perm %.3f", q_guild, q_perm)
                              else sprintf("q %.3f", q_guild)),
                          hjust = 0, size = 2.7, colour = "grey30") +
                scale_colour_manual(values = c(`TRUE` = "#D6604D", `FALSE` = "grey55"),
                                    labels = c(`TRUE` = sprintf("q < %.2f", FDR_ALPHA),
                                               `FALSE` = "not significant"),
                                    breaks = c(TRUE, FALSE), drop = FALSE, name = NULL) +
                scale_x_continuous(limits = c(-rng * 1.1, rng * 2.2)) +
                labs(title = title, subtitle = subtitle,
                     x = "coefficient (bagged vs unbagged)", y = NULL) +
                theme_bw(base_size = 10) +
                theme(legend.position    = "bottom",
                      plot.title         = element_text(size = 11, face = "bold"),
                      plot.subtitle      = element_text(size = 8, colour = "grey40"),
                      panel.grid.major.y = element_blank())
}

fig <- (panel("prevalence", "Prevalence model",
              "presence / absence, logistic, tested against zero; perm = permutation q from 16d",
              TRUE) /
        panel("abundance", "Abundance model",
              "log(TSS) among flowers where the taxon is present, tested against the median coefficient",
              FALSE)) +
        plot_layout(guides = "collect") +
        plot_annotation(
                title    = "Insect-symbiont guild under MaAsLin 3",
                subtitle = sprintf(paste("fitted on %d features (guild + 10%% prevalence background)",
                                         "so the compositional median is meaningful;",
                                         "BH applied to the %d pre-specified guild ASVs"),
                                   length(fit_set), length(guild)),
                caption  = "positive coefficient = higher in bagged flowers  |  bars are 95% Wald intervals",
                theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                              plot.subtitle = element_text(size = 9, colour = "grey40"),
                              plot.caption  = element_text(size = 7, colour = "grey40"))
        ) &
        theme(legend.position = "bottom")

# ── 6. Save outputs ───────────────────────────────────────────────────────────

cat("\n── 6. Saving outputs ────────────────────────────────────────────────────\n")

ggsave(here("results", "figures", "16e_maaslin3_vs_permutation.png"),
       plot = fig, width = 9, height = 8, dpi = 300)
cat("  Saved: 16e_maaslin3_vs_permutation.png\n")

cmp |>
        dplyr::select(taxon, tax_label, model, n_bagged, n_unbag, n_nonzero,
                      coef, stderr, null_hypothesis,
                      p_maaslin, q_guild, q_background,
                      q_perm, or_bagged, agrees) |>
        dplyr::arrange(model, p_maaslin) |>
        readr::write_csv(here("results", "tables",
                              "16e_maaslin3_vs_permutation.csv"))
cat("  Saved: 16e_maaslin3_vs_permutation.csv\n")

cat("\n16e done.\n")
