## ============================================================================
## 26b — ITS1 fungi CLR clustermap: ASV vs OTU97 side by side
##
## Pipeline (applied independently to each dataset):
##   1. Keep flower biosamples; prevalence filter >= 10% of samples
##      (OTU97: same feature space as 08b PCoA and 10b PERMANOVA)
##   2. Multiplicative zero replacement (zCompositions::cmultRepl, CZM method)
##      — falls back to +0.5 pseudocount if cmultRepl fails; the method actually
##        used is printed to the console and in the heatmap subtitle
##   3. CLR transformation: log(x_i) - mean(log(x)) per sample
##      → Euclidean distance on CLR = Aitchison distance (aitchison_CLR in 08b)
##   4. Ward.D2 hierarchical clustering of samples and taxa
##   5. ComplexHeatmap with a sequencing-depth bar track and 5 metadata strips
##
## Layout: 2 heatmaps side by side plus a shared annotation legend column on the
## right. ASV and OTU97 have different feature sets, so they cannot be
## concatenated with + (requires identical row counts); each is drawn into its
## own grid viewport instead.
##
## Input:
##   results/rds/ps_ITS1_fungi_biosamples.rds       (ASV)
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds (OTU97)
## Output:
##   results/figures/26b_ITS1_clr_clustermap.png
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

PREV_THRESH <- 0.10        # taxon present in >= 10% of samples (as in 08b)
ZERO_METHOD <- "CZM"       # cmultRepl: count-zero multiplicative replacement
LINK_METHOD <- "ward.D2"   # linkage for hclust — Ward minimises within-cluster variance
DIST_METHOD <- "euclidean" # Euclidean on CLR = Aitchison distance

# ── 1. Libraries ──────────────────────────────────────────────────────────────

library(here)
library(phyloseq)
library(zCompositions)   # cmultRepl() — multiplicative zero replacement
library(ComplexHeatmap)  # Heatmap(), HeatmapAnnotation(), Legend(), packLegend()
library(circlize)        # colorRamp2() — continuous colour mapping for heatmap body
library(RColorBrewer)
library(grid)            # gpar(), viewport(), grid.layout()

# ── 2. Load data ──────────────────────────────────────────────────────────────
# Flower-only subset mirrors 08b, so the OTU97 feature space after the
# prevalence filter is identical to the one used for PCoA and PERMANOVA.

ps_asv   <- readRDS(here("results", "rds", "ps_ITS1_fungi_biosamples.rds")) |>
    subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))
ps_otu97 <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds")) |>
    subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))

cat(sprintf("Loaded  ASV: %d taxa x %d samples\n", ntaxa(ps_asv),   nsamples(ps_asv)))
cat(sprintf("Loaded OTU97: %d taxa x %d samples\n", ntaxa(ps_otu97), nsamples(ps_otu97)))

# The shared annotation legend assumes both datasets carry the same samples.
stopifnot(setequal(sample_names(ps_asv), sample_names(ps_otu97)))

# ── 3. Helper: prevalence filter ──────────────────────────────────────────────
# Removes rare taxa to reduce noise before CLR.
# Filtering BEFORE CLR is important: the geometric mean in CLR is sensitive
# to many near-zero features.

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

    # Extract count matrix: orient as samples x taxa for cmultRepl
    otu_mat <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)   # samples x taxa

    # Step 1: Multiplicative zero replacement (same call as 08b, aitchison_CLR).
    # cmultRepl() estimates non-zero replacements from the marginal count
    # distributions, preserving the compositional structure better than a
    # flat pseudocount. output="p-counts" returns pseudo-counts on the
    # original count scale (not proportions).
    clr_mat <- tryCatch({

        replaced <- zCompositions::cmultRepl(
            otu_mat,
            label          = 0,
            method         = ZERO_METHOD,
            output         = "p-counts",
            suppress.print = TRUE,
            z.warning      = 0.9,   # warn only above 90% zeros (not 80%)
            z.delete       = FALSE  # warn but do NOT delete sparse features/samples
        )

        # Step 2: CLR — subtract the log geometric mean of each sample
        out <- t(apply(replaced, 1, function(x) log(x) - mean(log(x))))
        attr(out, "zero_method") <- ZERO_METHOD
        out

    }, error = function(e) {
        # Fallback when data is too sparse for cmultRepl (convergence failure).
        # Simple +0.5 pseudocount before CLR — less principled; reported in subtitle.
        message(sprintf("  [WARNING] cmultRepl failed: %s", conditionMessage(e)))
        pseudo <- otu_mat + 0.5
        out    <- t(apply(pseudo, 1, function(x) log(x) - mean(log(x))))
        attr(out, "zero_method") <- "pseudocount +0.5"
        out
    })

    cat(sprintf("  CLR matrix: %d samples x %d taxa | zero replacement used: %s\n",
                nrow(clr_mat), ncol(clr_mat), attr(clr_mat, "zero_method")))
    clr_mat  # returns samples x taxa
}

# ── 5. Helper: colour-safe character conversion ────────────────────────────────
# as.character(NA) → "NA" (a string), which ComplexHeatmap won't treat as NA.
# This helper preserves real R NA so na_col takes effect.

safe_chr <- function(x) {
    out         <- as.character(x)
    out[is.na(x)] <- NA_character_
    out
}

# ── 6. Helper: annotation colour maps ─────────────────────────────────────────
# Built once and reused by both heatmaps and by the shared legend, so colours
# are guaranteed to match across panels.

build_annot_colors <- function(meta) {
    levels_of <- function(x) sort(unique(na.omit(safe_chr(x))))

    st_vals <- levels_of(meta$sample_type)
    mt_vals <- levels_of(meta$management_type)
    fi_vals <- levels_of(meta$farm_id)
    ti_vals <- levels_of(meta$tree_id)

    list(
        # sample_type: Set1 (red / blue)
        sample_type     = setNames(brewer.pal(3, "Set1")[seq_along(st_vals)], st_vals),
        # management_type: RdYlGn, levels in alphabetical order
        management_type = setNames(brewer.pal(max(3, length(mt_vals)), "RdYlGn")[seq_along(mt_vals)], mt_vals),
        # farm_id: Dark2 (8 colours, 7 farms)
        farm_id         = setNames(brewer.pal(max(3, length(fi_vals)), "Dark2")[seq_along(fi_vals)], fi_vals),
        # tree_id: Paired interpolated to the number of trees
        tree_id         = setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(ti_vals)), ti_vals),
        # is_pollination_clsm (logical)
        is_pollination_clsm = c("TRUE" = "#2ca02c", "FALSE" = "#d62728")
    )
}

NA_COL <- "grey80"   # colour for missing metadata in every strip

# ── 7. Helper: shared annotation legend ───────────────────────────────────────
# tree_id has no legend: ~49 levels would be unreadable; the strip is only meant
# to show whether flowers of one tree sit together.

build_annot_legend <- function(cols) {
    discrete_legend <- function(title, col_map) {
        Legend(title     = title,
               at        = c(names(col_map), "NA"),
               legend_gp = gpar(fill = c(unname(col_map), NA_COL)),
               title_gp  = gpar(fontsize = 9, fontface = "bold"),
               labels_gp = gpar(fontsize = 8))
    }
    packLegend(
        discrete_legend("sample_type",         cols$sample_type),
        discrete_legend("management_type",     cols$management_type),
        discrete_legend("farm_id",             cols$farm_id),
        discrete_legend("is_pollination_clsm", cols$is_pollination_clsm),
        direction = "vertical",
        gap       = unit(5, "mm")                     # vertical space between legends
    )
}

# ── 8. Helper: build one ComplexHeatmap ───────────────────────────────────────
# clr_mat     — samples x taxa CLR matrix (output of prepare_clr)
# meta, depth — sample metadata and library sizes (any order; aligned here)
# title       — string for column_title of this heatmap

build_heatmap <- function(clr_mat, meta, depth, cols, title, hm_name) {

    # Transpose: heatmap body = taxa (rows) x samples (columns)
    mat <- t(clr_mat)

    # Align metadata and library size to matrix column order
    meta  <- meta[colnames(mat), , drop = FALSE]
    depth <- depth[colnames(mat)]

    # ── Top annotation: depth bar track + metadata strips ─────────────────────
    # log10_depth — unbagged flowers are sequenced deeper than bagged (see 08b);
    # if a dendrogram split follows this track, depth or zero replacement, not
    # biology, may be driving the clustering.
    ha <- HeatmapAnnotation(
        log10_depth = anno_barplot(
            log10(depth),
            baseline   = "min",                       # bars start at the shallowest sample
            height     = unit(8, "mm"),               # height of the depth track
            gp         = gpar(fill = "grey40", col = NA),
            border     = FALSE,
            axis_param = list(gp = gpar(fontsize = 5))
        ),
        sample_type         = safe_chr(meta$sample_type),
        management_type     = safe_chr(meta$management_type),
        farm_id             = safe_chr(meta$farm_id),
        tree_id             = safe_chr(meta$tree_id),
        is_pollination_clsm = safe_chr(meta$is_pollination_clsm),
        col                  = cols,
        na_col               = NA_COL,
        show_legend          = FALSE,               # one shared legend is drawn separately
        show_annotation_name = TRUE,
        annotation_name_side = "left",              # right side collides with the CLR legend
        annotation_name_rot  = 0,                       # horizontal names; avoids overlap with the depth track
        annotation_name_gp   = gpar(fontsize = 7),
        simple_anno_size     = unit(3, "mm"),       # height of each metadata strip
        gap                  = unit(0.5, "mm")
    )

    # ── Subtitle with analysis parameters ─────────────────────────────────────
    subtitle <- sprintf(
        "prev ≥ %.0f%%  |  zero repl: %s  |  dist: %s (Aitchison)  |  linkage: %s  |  %d taxa x %d samples",
        100 * PREV_THRESH, attr(clr_mat, "zero_method"), DIST_METHOD, LINK_METHOD,
        nrow(mat), ncol(mat)
    )

    # ── Colour scale for CLR values ────────────────────────────────────────────
    # Symmetric diverging palette centred at 0 (mean CLR per sample = 0 by definition).
    # Range from 1st to 99th percentile of the matrix to avoid outlier compression.
    clr_range <- quantile(mat, c(0.01, 0.99), na.rm = TRUE)
    abs_max   <- max(abs(clr_range))   # keep scale symmetric around 0

    col_fun <- colorRamp2(
        c(-abs_max, 0, abs_max),
        c("#313695", "white", "#a50026")   # blue–white–red diverging
    )

    # ── Build Heatmap object ───────────────────────────────────────────────────
    Heatmap(
        mat,
        name            = hm_name,   # must be unique across heatmaps in the same device
        column_title    = paste0(title, "\n", subtitle),
        column_title_gp = gpar(fontsize = 8),

        top_annotation  = ha,

        # Clustering: Ward.D2 on Euclidean distance.
        # On CLR-transformed data, Euclidean distance equals the Aitchison
        # distance — the natural distance for compositional data.
        clustering_method_rows      = LINK_METHOD,
        clustering_method_columns   = LINK_METHOD,
        clustering_distance_rows    = DIST_METHOD,
        clustering_distance_columns = DIST_METHOD,

        col = col_fun,

        # No row/column labels — too many features and samples to be readable
        show_row_names    = FALSE,
        show_column_names = FALSE,

        row_dend_width     = unit(10, "mm"),        # width of the taxa dendrogram
        column_dend_height = unit(15, "mm"),        # height of the sample dendrogram

        heatmap_legend_param = list(
            title         = "CLR",
            direction     = "vertical",
            legend_height = unit(30, "mm")
        ),

        use_raster = TRUE,   # rasterise the heatmap body for faster rendering
        raster_quality = 3
    )
}

# ── 9. Run pipeline for both datasets ─────────────────────────────────────────

cat("\n── ASV dataset ──────────────────────────────────────────────────────────\n")
ps_asv_filt <- filter_taxa_ps(ps_asv, PREV_THRESH)
clr_asv     <- prepare_clr(ps_asv_filt)

cat("\n── OTU97 dataset ────────────────────────────────────────────────────────\n")
ps_otu97_filt <- filter_taxa_ps(ps_otu97, PREV_THRESH)
clr_otu97     <- prepare_clr(ps_otu97_filt)

# Metadata and colours (shared); library size per dataset, before the filter
meta <- as(sample_data(ps_otu97), "data.frame")
cols <- build_annot_colors(meta)

# ── 10. Build heatmap objects ─────────────────────────────────────────────────

hm_asv <- build_heatmap(
    clr_asv, meta, sample_sums(ps_asv), cols,
    title   = sprintf("ITS1 ASV  (%d taxa after filter)", ntaxa(ps_asv_filt)),
    hm_name = "CLR_ASV"
)

hm_otu97 <- build_heatmap(
    clr_otu97, meta, sample_sums(ps_otu97), cols,
    title   = sprintf("ITS1 OTU97  (%d taxa after filter)", ntaxa(ps_otu97_filt)),
    hm_name = "CLR_OTU97"
)

annot_legend <- build_annot_legend(cols)

# ── 11. Save ──────────────────────────────────────────────────────────────────

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

out_fp <- here("results", "figures", "26b_ITS1_clr_clustermap.png")

# Figure size in inches; fonts are in points, so a larger canvas = smaller text
png(out_fp, width = 26, height = 12, units = "in", res = 300)

grid.newpage()

# Main title centred at the top
grid.text(
    "ITS1 fungi — CLR clustermap (ASV vs OTU97)",
    x = 0.5, y = 0.98,
    gp = gpar(fontsize = 13, fontface = "bold")
)

# Grid layout: 1 row x 3 columns (ASV, OTU97, legend).
# widths controls the relative width of each column; the last one is the legend.
pushViewport(viewport(
    y = 0.47, height = 0.94,
    layout = grid.layout(nrow = 1, ncol = 3,
                         widths = unit(c(1, 1, 0.12), "null"))
))

# Left: ASV heatmap
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(hm_asv, newpage = FALSE,
     padding = unit(c(2, 5, 5, 5), "mm"))     # bottom, left, top, right
popViewport()

# Middle: OTU97 heatmap
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(hm_otu97, newpage = FALSE,
     padding = unit(c(2, 5, 5, 5), "mm"))
popViewport()

# Right: shared annotation legend
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(annot_legend, x = unit(2, "mm"), y = unit(0.5, "npc"), just = "left")
popViewport()

popViewport()
dev.off()
cat(sprintf("\nSaved: %s\n", out_fp))
