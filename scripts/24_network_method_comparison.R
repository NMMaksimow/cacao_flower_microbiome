## ============================================================================
## 24 — Cross-method network comparison: SpiecEasi vs SparCC vs SPRING
##
## The three co-occurrence methods rest on different assumptions (SparCC:
## pairwise correlation with a threshold cutoff; SpiecEasi/SPRING: partial
## correlation via a graphical model that assumes n >> p, which does not hold
## here — n=21, p=88-134 genera per farm-group). Before reporting any single
## method's result, this checks whether the three actually agree, and how much
## formal statistical support each one can offer.
##
## Three comparison layers:
##   1. Hub overlap: do the methods agree on which taxa are central? (Jaccard
##      on hub taxon ID sets, per farm x treatment x marker x method pair)
##   2. Direction agreement: do the methods agree on the sign of the bagged vs
##      unbagged difference for each whole-network topology property?
##   3. Statistical validity inventory: how many farm x property tests does
##      each method support with a non-NA permutation p-value, and how many of
##      those are significant? (SPRING has none for topology — its netCompare
##      permutation test fails to converge; see 23_spring_netcompare_diag.txt)
##
## Input:
##   results/tables/{19a,19b,20a,20b,21a,21b}_network_hubs.csv
##   results/tables/{19a,19b,20a,20b,21a,21b}_network_comparison.csv
## Output:
##   results/tables/24_method_comparison_hub_jaccard.csv
##   results/tables/24_method_comparison_direction_agreement.csv
##   results/tables/24_method_comparison_significance_inventory.csv
##   results/figures/24_method_comparison_direction.png
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

FARM_LEVELS <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
SIG_ALPHA   <- 0.05

# method/marker -> table prefix. Order fixed here sets the method colour and
# facet order used throughout (SpiecEasi, SparCC, SPRING).
PREFIX_MAP <- tibble::tribble(
        ~method,      ~marker, ~prefix,
        "SpiecEasi",  "16S",   "19a",
        "SpiecEasi",  "ITS1",  "19b",
        "SparCC",     "16S",   "20a",
        "SparCC",     "ITS1",  "20b",
        "SPRING",     "16S",   "21a",
        "SPRING",     "ITS1",  "21b"
)
METHOD_LEVELS <- c("SpiecEasi", "SparCC", "SPRING")

DIR_COL_UNBAGGED <- "#D6604D"   # warm: property higher in unbagged (same as EDGE_COL_POS)
DIR_COL_BAGGED   <- "#4393C3"   # cool: property higher in bagged (same as EDGE_COL_NEG)
DIR_COL_TIED     <- "grey85"    # no difference (e.g. vertConnect = 0 in both)
DIR_COL_NA       <- "grey95"    # farm/property combination absent for this method

# ── 1. Libraries ──────────────────────────────────────────────────────────────

library(here)
library(tidyverse)

# ── 2. Load all hub and comparison tables ─────────────────────────────────────

cat("\n── 2. Loading tables ─────────────────────────────────────────────────────\n")

hubs_all <- PREFIX_MAP |>
        pmap(\(method, marker, prefix) {
                readr::read_csv(here("results", "tables", paste0(prefix, "_network_hubs.csv")),
                                show_col_types = FALSE) |>
                        mutate(method = method, marker = marker, .before = 1)
        }) |>
        list_rbind()

cmp_all <- PREFIX_MAP |>
        pmap(\(method, marker, prefix) {
                readr::read_csv(here("results", "tables", paste0(prefix, "_network_comparison.csv")),
                                show_col_types = FALSE) |>
                        mutate(method = method, marker = marker, .before = 1)
        }) |>
        list_rbind()

cat(sprintf("  hubs: %d rows across %d method x marker tables\n", nrow(hubs_all), nrow(PREFIX_MAP)))
cat(sprintf("  comparison: %d rows across %d method x marker tables\n", nrow(cmp_all), nrow(PREFIX_MAP)))

# ── 3. Layer 1 — Hub taxon Jaccard overlap between method pairs ──────────────

cat("\n── 3. Hub overlap between method pairs ───────────────────────────────────\n")

jaccard <- function(a, b) {
        if (length(a) == 0 && length(b) == 0) return(NA_real_)   # neither method found a hub here
        length(intersect(a, b)) / length(union(a, b))
}

hub_sets <- hubs_all |>
        summarise(taxa = list(taxon), .by = c(marker, farm, treatment, method))

hub_wide    <- hub_sets |> pivot_wider(names_from = method, values_from = taxa)
method_pairs <- combn(METHOD_LEVELS, 2, simplify = FALSE)

# One pass per method pair: .data[[...]] pulls each pair's two list-columns of
# taxon IDs out of hub_wide, map2_dbl runs jaccard() row by row.
hub_jaccard <- method_pairs |>
        map(\(pair) {
                hub_wide |>
                        transmute(
                                marker, farm, treatment,
                                method_a = pair[1], method_b = pair[2],
                                n_hubs_a = lengths(.data[[pair[1]]]),
                                n_hubs_b = lengths(.data[[pair[2]]]),
                                jaccard  = map2_dbl(.data[[pair[1]]], .data[[pair[2]]], jaccard)
                        )
        }) |>
        list_rbind() |>
        arrange(marker, farm, treatment, method_a, method_b)

readr::write_csv(hub_jaccard, here("results", "tables", "24_method_comparison_hub_jaccard.csv"))
cat("  Saved: 24_method_comparison_hub_jaccard.csv\n")
cat(sprintf("  mean Jaccard across all farm x treatment x marker x method-pair: %.3f\n",
            mean(hub_jaccard$jaccard, na.rm = TRUE)))

# ── 4. Layer 2 — Direction agreement on whole-network topology ───────────────

cat("\n── 4. Direction agreement (unbagged vs bagged) ───────────────────────────\n")

direction_df <- cmp_all |>
        filter(scope == "whole") |>
        mutate(
                diff      = val_unbagged - val_bagged,
                direction = case_when(
                        is.na(diff)  ~ NA_character_,
                        diff  > 0    ~ "higher in unbagged",
                        diff  < 0    ~ "higher in bagged",
                        TRUE         ~ "tied"
                )
        ) |>
        select(marker, farm, property, method, direction)

agreement_df <- direction_df |>
        filter(!is.na(direction)) |>
        summarise(
                n_methods_reporting = n(),
                n_agree_majority    = max(table(direction)),
                all_agree           = n_distinct(direction) == 1,
                consensus_direction = names(which.max(table(direction))),
                .by = c(marker, farm, property)
        ) |>
        arrange(marker, property, farm)

readr::write_csv(agreement_df, here("results", "tables", "24_method_comparison_direction_agreement.csv"))
cat("  Saved: 24_method_comparison_direction_agreement.csv\n")

agree_summary <- agreement_df |>
        summarise(pct_all_3_agree = mean(all_agree), .by = marker)
cat("  share of farm x property combinations where all 3 methods agree on direction:\n")
print(as.data.frame(agree_summary), row.names = FALSE)

# ── 5. Layer 3 — Statistical validity inventory ──────────────────────────────

cat("\n── 5. Statistical validity inventory ─────────────────────────────────────\n")

sig_inventory <- cmp_all |>
        filter(scope == "whole") |>
        summarise(
                n_tests_possible = n(),
                n_tests_with_pval = sum(!is.na(pval)),
                n_significant     = sum(!is.na(pval) & pval < SIG_ALPHA),
                .by = c(method, marker)
        ) |>
        arrange(marker, method)

readr::write_csv(sig_inventory, here("results", "tables", "24_method_comparison_significance_inventory.csv"))
cat("  Saved: 24_method_comparison_significance_inventory.csv\n")
print(as.data.frame(sig_inventory), row.names = FALSE)

# ── 6. Figure: direction agreement tile grid ──────────────────────────────────

cat("\n── 6. Building comparison figure ─────────────────────────────────────────\n")

plot_df <- direction_df |>
        mutate(
                method    = factor(method, levels = METHOD_LEVELS),
                direction = factor(direction,
                                   levels = c("higher in unbagged", "higher in bagged", "tied"))
        ) |>
        left_join(agreement_df |> select(marker, farm, property, all_agree),
                  by = c("marker", "farm", "property"))

# Tile fill = direction of the bagged/unbagged difference for that method.
# Black outline = all 3 methods agree on direction for that farm x property
# (drawn identically in every method column so the reader can scan across a
# row of facets for an unbroken black outline).
# facet_grid rows = marker, columns = method; adjust facet formula to change layout.
p_comparison <- ggplot(plot_df, aes(x = property, y = toupper(farm), fill = direction)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        geom_tile(data = filter(plot_df, all_agree), colour = "black",
                  linewidth = 0.8, fill = NA) +
        scale_fill_manual(
                values   = c("higher in unbagged" = DIR_COL_UNBAGGED,
                            "higher in bagged"    = DIR_COL_BAGGED,
                            "tied"                = DIR_COL_TIED),
                na.value = DIR_COL_NA,
                name     = "Direction"
        ) +
        facet_grid(marker ~ method) +
        labs(
                title    = "Direction of bagged vs unbagged difference across methods",
                subtitle = "Black outline = all 3 methods agree on direction",
                x = NULL, y = NULL
        ) +
        theme_minimal(base_size = 10) +
        theme(
                axis.text.x     = element_text(angle = 45, hjust = 1),   # property names
                panel.grid      = element_blank(),
                strip.text      = element_text(face = "bold"),
                legend.position = "bottom"
        )

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(
        here("results", "figures", "24_method_comparison_direction.png"),
        plot   = p_comparison,
        width  = 10,
        height = 6,
        dpi    = 300
)
cat("  Saved: 24_method_comparison_direction.png\n")

# ── 7. Recommendation rubric ──────────────────────────────────────────────────

cat("\n── 7. Method recommendation rubric ───────────────────────────────────────\n")

for (m in METHOD_LEVELS) {
        sig_m <- sig_inventory |> filter(method == m)
        jac_m <- hub_jaccard |> filter(method_a == m | method_b == m)
        cat(sprintf("\n%s:\n", m))
        cat(sprintf("  assumption fit: %s\n", switch(m,
                SparCC    = "pairwise correlation, no n>>p requirement — best fit for n=21",
                SpiecEasi = "graphical model, assumes n>>p; here p=88-134 >> n=21 (exploratory)",
                SPRING    = "graphical model, same n>>p problem as SpiecEasi (exploratory)"
        )))
        cat(sprintf("  statistical validity: %d of %d whole-topology tests have a p-value (%.0f%%), %d significant at p<%.2f\n",
                    sum(sig_m$n_tests_with_pval), sum(sig_m$n_tests_possible),
                    100 * sum(sig_m$n_tests_with_pval) / sum(sig_m$n_tests_possible),
                    sum(sig_m$n_significant), SIG_ALPHA))
        cat(sprintf("  mean hub-set Jaccard with the other two methods: %.3f\n",
                    mean(jac_m$jaccard, na.rm = TRUE)))
}

cat("\nRead the three numbers above per method before deciding what to lead with:\n")
cat("a method with few/no valid p-values cannot carry a significance claim alone;\n")
cat("a property/direction that all three methods reproduce (see the figure's\n")
cat("black outlines) is stronger evidence than any single method's p-value.\n")
