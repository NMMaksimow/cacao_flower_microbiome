## ============================================================================
## 26a — 16S bacteria CLR clustermap
##
## Pipeline:
##   1. Keep flower biosamples; prevalence filter >= 10% of samples
##      (same feature space as 08a PCoA and 10a PERMANOVA)
##   2. Multiplicative zero replacement (zCompositions::cmultRepl, CZM method)
##      — falls back to +0.5 pseudocount if cmultRepl fails; the method actually
##        used is printed to the console and in the heatmap subtitle
##   3. CLR transformation: log(x_i) - mean(log(x)) per sample
##      → Euclidean distance on CLR = Aitchison distance (aitchison_CLR in 08a)
##   4. Ward.D2 hierarchical clustering of samples and taxa
##   5. ComplexHeatmap with a sequencing-depth bar track and 5 metadata strips
##
## Input:
##   results/rds/ps_16S_bacteria_biosamples.rds
## Output:
##   results/figures/26a_16S_clr_clustermap.png
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

PREV_THRESH <- 0.10        # taxon present in >= 10% of samples (as in 08a)
ZERO_METHOD <- "CZM"       # cmultRepl: count-zero multiplicative replacement
LINK_METHOD <- "ward.D2"   # Ward minimises within-cluster variance
DIST_METHOD <- "euclidean" # Euclidean on CLR = Aitchison distance

# ── 1. Libraries ──────────────────────────────────────────────────────────────

library(here)
library(phyloseq)
library(zCompositions)   # cmultRepl() — multiplicative zero replacement
library(ComplexHeatmap)  # Heatmap(), HeatmapAnnotation()
library(circlize)        # colorRamp2() — continuous colour mapping
library(RColorBrewer)
library(grid)            # gpar()

# ── 2. Load data ──────────────────────────────────────────────────────────────
# Flower-only subset mirrors 08a, so the feature space after the prevalence
# filter is identical to the one used for PCoA and PERMANOVA.

ps <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds")) |>
    subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))

cat(sprintf("Loaded 16S bacteria biosamples: %d taxa x %d samples\n",
            ntaxa(ps), nsamples(ps)))

# ── 3. Helper: prevalence filter ──────────────────────────────────────────────
# Filter before CLR — removing rare taxa stabilises the geometric mean
# and prevents cmultRepl from being overwhelmed by near-zero features.

filter_taxa_ps <- function(ps, prev_thresh) {
    otu_mat <- as(otu_table(ps), "matrix")
    if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)   # orient: taxa x samples

    keep <- rowMeans(otu_mat > 0) >= prev_thresh    # fraction of samples with taxon > 0

    cat(sprintf("  Filter: %d / %d taxa retained (prev >= %.0f%%)\n",
                sum(keep), nrow(otu_mat), 100 * prev_thresh))

    prune_taxa(keep, ps)
}

# ── 4. Helper: CLR transformation ─────────────────────────────────────────────
# Returns a samples x taxa CLR matrix. The zero-replacement method actually used
# is stored in attr(, "zero_method") so it can be reported in the figure.

prepare_clr <- function(ps) {

    # Extract counts: orient as samples x taxa for cmultRepl
    otu_mat <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)

    clr_mat <- tryCatch({

        # Same call as 08a (aitchison_CLR). z.delete = FALSE keeps all features.
        replaced <- zCompositions::cmultRepl(
            otu_mat,
            label          = 0,
            method         = ZERO_METHOD,
            output         = "p-counts",
            suppress.print = TRUE,
            z.warning      = 0.9,    # warn only above 90% zeros
            z.delete       = FALSE   # do not drop sparse features/samples
        )

        # CLR: for each sample subtract the log geometric mean
        out <- t(apply(replaced, 1, function(x) log(x) - mean(log(x))))
        attr(out, "zero_method") <- ZERO_METHOD
        out

    }, error = function(e) {
        message(sprintf("  [WARNING] cmultRepl failed: %s", conditionMessage(e)))
        pseudo <- otu_mat + 0.5
        out    <- t(apply(pseudo, 1, function(x) log(x) - mean(log(x))))
        attr(out, "zero_method") <- "pseudocount +0.5"
        out
    })

    cat(sprintf("  CLR matrix: %d samples x %d taxa | zero replacement used: %s\n",
                nrow(clr_mat), ncol(clr_mat), attr(clr_mat, "zero_method")))
    clr_mat   # samples x taxa
}

# ── 5. Helper: safe NA-preserving character conversion ────────────────────────
# as.character(NA) returns "NA" string — ComplexHeatmap won't treat as missing.
# This helper returns real R NA so na_col = "grey80" works correctly.

safe_chr <- function(x) {
    out           <- as.character(x)
    out[is.na(x)] <- NA_character_
    out
}

# ── 6. Run pipeline ───────────────────────────────────────────────────────────

cat("\n── Filtering ────────────────────────────────────────────────────────────\n")
ps_filt <- filter_taxa_ps(ps, PREV_THRESH)

cat("\n── CLR transformation ───────────────────────────────────────────────────\n")
clr_mat <- prepare_clr(ps_filt)

# ── 7. Build annotation ───────────────────────────────────────────────────────

# Heatmap body: taxa (rows) x samples (columns)
mat  <- t(clr_mat)

# Align metadata and library size (reads before the prevalence filter) to columns
meta  <- as(sample_data(ps), "data.frame")[colnames(mat), , drop = FALSE]
depth <- sample_sums(ps)[colnames(mat)]

levels_of <- function(x) sort(unique(na.omit(safe_chr(x))))

# sample_type (2 levels: bagged_flower / unbagged_flower)
st_vals   <- levels_of(meta$sample_type)
st_colors <- setNames(brewer.pal(3, "Set1")[seq_along(st_vals)], st_vals)

# management_type (4 levels, alphabetical order)
mt_vals   <- levels_of(meta$management_type)
mt_colors <- setNames(brewer.pal(max(3, length(mt_vals)), "RdYlGn")[seq_along(mt_vals)], mt_vals)

# farm_id (7 farms)
fi_vals   <- levels_of(meta$farm_id)
fi_colors <- setNames(brewer.pal(max(3, length(fi_vals)), "Dark2")[seq_along(fi_vals)], fi_vals)

# tree_id — many values; interpolate from Paired palette
ti_vals   <- levels_of(meta$tree_id)
ti_colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(ti_vals)), ti_vals)

# is_pollination_clsm (logical TRUE / FALSE / NA)
ip_colors <- c("TRUE" = "#2ca02c", "FALSE" = "#d62728")

# Top annotation: depth bar track, then metadata strips.
# log10_depth — unbagged flowers are sequenced deeper than bagged (see 08a); if a
# dendrogram split follows this track, depth or zero replacement, not biology,
# may be driving the clustering.
# Legends are shown for every strip except tree_id (~49 levels, unreadable).
ha <- HeatmapAnnotation(
    log10_depth = anno_barplot(
        log10(depth),
        baseline   = "min",                         # bars start at the shallowest sample
        height     = unit(8, "mm"),                 # height of the depth track
        gp         = gpar(fill = "grey40", col = NA),
        border     = FALSE,
        axis_param = list(gp = gpar(fontsize = 5))
    ),
    sample_type         = safe_chr(meta$sample_type),
    management_type     = safe_chr(meta$management_type),
    farm_id             = safe_chr(meta$farm_id),
    tree_id             = safe_chr(meta$tree_id),
    is_pollination_clsm = safe_chr(meta$is_pollination_clsm),
    col = list(
        sample_type         = st_colors,
        management_type     = mt_colors,
        farm_id             = fi_colors,
        tree_id             = ti_colors,
        is_pollination_clsm = ip_colors
    ),
    na_col               = "grey80",
    show_legend          = c(log10_depth = FALSE, sample_type = TRUE,
                             management_type = TRUE, farm_id = TRUE,
                             tree_id = FALSE, is_pollination_clsm = TRUE),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_rot  = 0,                       # horizontal names; avoids overlap with the depth track
    annotation_name_gp   = gpar(fontsize = 7),
    simple_anno_size     = unit(3, "mm"),           # height of each metadata strip
    gap                  = unit(0.5, "mm")
)

# ── 8. Colour scale for heatmap body ─────────────────────────────────────────
# Symmetric diverging palette centred at 0 (CLR mean per sample = 0 by definition).
# Clipped at 1st/99th percentile to avoid compression by extreme outliers.

clr_range <- quantile(mat, c(0.01, 0.99), na.rm = TRUE)
abs_max   <- max(abs(clr_range))

col_fun <- colorRamp2(
    c(-abs_max, 0, abs_max),
    c("#313695", "white", "#a50026")   # blue–white–red diverging
)

# ── 9. Subtitle ───────────────────────────────────────────────────────────────

subtitle <- sprintf(
    "prev ≥ %.0f%%  |  zero repl: %s  |  dist: %s (Aitchison)  |  linkage: %s  |  %d taxa x %d samples",
    100 * PREV_THRESH, attr(clr_mat, "zero_method"), DIST_METHOD, LINK_METHOD,
    nrow(mat), ncol(mat)
)

# ── 10. Build heatmap ─────────────────────────────────────────────────────────

hm <- Heatmap(
    mat,
    name            = "CLR",
    column_title    = paste0("16S bacteria biosamples\n", subtitle),
    column_title_gp = gpar(fontsize = 8),

    top_annotation  = ha,

    # Ward.D2 on Euclidean = Aitchison distance for CLR data
    clustering_method_rows      = LINK_METHOD,
    clustering_method_columns   = LINK_METHOD,
    clustering_distance_rows    = DIST_METHOD,
    clustering_distance_columns = DIST_METHOD,

    col = col_fun,

    show_row_names    = FALSE,
    show_column_names = FALSE,

    row_dend_width     = unit(10, "mm"),            # width of the taxa dendrogram
    column_dend_height = unit(15, "mm"),            # height of the sample dendrogram

    heatmap_legend_param = list(
        title         = "CLR",
        direction     = "vertical",
        legend_height = unit(30, "mm")
    ),

    use_raster     = TRUE,
    raster_quality = 3
)

# ── 11. Save ──────────────────────────────────────────────────────────────────

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

out_fp <- here("results", "figures", "26a_16S_clr_clustermap.png")

# Figure size in inches; fonts are in points, so a larger canvas = smaller text
png(out_fp, width = 16, height = 12, units = "in", res = 300)
draw(
    hm,
    column_title         = "16S bacteria — CLR clustermap",
    column_title_gp      = gpar(fontsize = 13, fontface = "bold"),
    merge_legend         = TRUE,                    # CLR + annotation legends in one column
    heatmap_legend_side  = "right",
    annotation_legend_side = "right",
    padding              = unit(c(5, 15, 10, 5), "mm")   # bottom, left, top, right
)
dev.off()

cat(sprintf("\nSaved: %s\n", out_fp))
