## ============================================================================
## 26c — 16S bacteria + ITS1 OTU97 CLR clustermaps: Ward.D2 vs UPGMA linkage
##
## Pipeline (applied independently to each marker):
##   1. Keep flower biosamples; prevalence filter >= 10% of samples
##      (same feature space as 08a/08b PCoA and 10a/10b PERMANOVA)
##   2. Multiplicative zero replacement (zCompositions::cmultRepl, CZM method)
##      — falls back to +0.5 pseudocount if cmultRepl fails; the method actually
##        used is printed to the console and in the heatmap subtitle
##   3. CLR transformation: log(x_i) - mean(log(x)) per sample
##      → Euclidean distance on CLR = Aitchison distance (aitchison_CLR in 08a/08b)
##   4. Hierarchical clustering of samples and taxa with two linkages:
##      top row = Ward.D2, bottom row = average (UPGMA)
##   5. ComplexHeatmap with a sequencing-depth bar track and 5 metadata strips
##
## Layout: 2 x 2 heatmaps (columns = 16S | ITS1; rows = Ward.D2 | UPGMA) plus a
## shared annotation legend column on the right. The two markers have different
## feature sets, so heatmaps cannot be concatenated with + (requires identical
## row counts); each is drawn into its own grid viewport instead.
##
## Input:
##   results/rds/ps_16S_bacteria_biosamples.rds
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output:
##   results/figures/26c_16S_ITS1_clr_clustermap.png
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

PREV_THRESH <- 0.10        # taxon present in >= 10% of samples (as in 08a/08b)
ZERO_METHOD <- "CZM"       # cmultRepl: count-zero multiplicative replacement
DIST_METHOD <- "euclidean" # Euclidean on CLR = Aitchison distance

# Linkage per figure row (top -> bottom)
LINK_METHODS <- c(ward = "ward.D2",   # minimises within-cluster variance; compact, balanced clusters
                  upgma = "average")  # mean pairwise distance; usually preserves distances best

# ── 1. Libraries ──────────────────────────────────────────────────────────────

library(here)
library(phyloseq)
library(zCompositions)   # cmultRepl()
library(ComplexHeatmap)  # Heatmap(), HeatmapAnnotation(), Legend(), packLegend()
library(circlize)        # colorRamp2()
library(RColorBrewer)
library(grid)            # gpar(), viewport(), grid.layout() etc.

# ── 2. Load data ──────────────────────────────────────────────────────────────
# Flower-only subset mirrors 08a/08b, so the feature space after the prevalence
# filter is identical to the one used for PCoA and PERMANOVA.

ps_16S <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds")) |>
    subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))
ps_ITS <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds")) |>
    subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower"))

cat(sprintf("Loaded  16S: %d taxa x %d samples\n", ntaxa(ps_16S), nsamples(ps_16S)))
cat(sprintf("Loaded ITS1: %d taxa x %d samples\n", ntaxa(ps_ITS), nsamples(ps_ITS)))

# The shared annotation legend assumes both markers carry the same samples.
stopifnot(setequal(sample_names(ps_16S), sample_names(ps_ITS)))

# ── 3. Helper: prevalence filter ──────────────────────────────────────────────

filter_taxa_ps <- function(ps, prev_thresh, label = "") {
    otu_mat <- as(otu_table(ps), "matrix")
    if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)   # orient: taxa x samples

    keep <- rowMeans(otu_mat > 0) >= prev_thresh    # fraction of samples with taxon > 0

    cat(sprintf("  [%s] %d / %d taxa retained (prev >= %.0f%%)\n",
                label, sum(keep), nrow(otu_mat), 100 * prev_thresh))

    prune_taxa(keep, ps)
}

# ── 4. Helper: CLR transformation ─────────────────────────────────────────────
# Returns a samples x taxa CLR matrix. The zero-replacement method actually used
# is stored in attr(, "zero_method") so it can be reported in the figure.

prepare_clr <- function(ps, label = "") {

    otu_mat <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)   # samples x taxa

    clr_mat <- tryCatch({
        # Same call as 08a/08b (aitchison_CLR)
        replaced <- zCompositions::cmultRepl(
            otu_mat,
            label          = 0,
            method         = ZERO_METHOD,
            output         = "p-counts",
            suppress.print = TRUE,
            z.warning      = 0.9,    # warn only above 90% zeros
            z.delete       = FALSE   # keep all features/samples
        )
        out <- t(apply(replaced, 1, function(x) log(x) - mean(log(x))))
        attr(out, "zero_method") <- ZERO_METHOD
        out

    }, error = function(e) {
        message(sprintf("  [%s][WARNING] cmultRepl failed: %s", label, conditionMessage(e)))
        pseudo <- otu_mat + 0.5
        out    <- t(apply(pseudo, 1, function(x) log(x) - mean(log(x))))
        attr(out, "zero_method") <- "pseudocount +0.5"
        out
    })

    cat(sprintf("  [%s] CLR matrix: %d samples x %d taxa | zero replacement used: %s\n",
                label, nrow(clr_mat), ncol(clr_mat), attr(clr_mat, "zero_method")))
    clr_mat
}

# ── 5. Helper: safe NA-preserving character conversion ────────────────────────

safe_chr <- function(x) {
    out           <- as.character(x)
    out[is.na(x)] <- NA_character_
    out
}

# ── 6. Helper: annotation colour maps ─────────────────────────────────────────
# Built once and reused by every heatmap and by the shared legend, so colours
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

# ── 7. Helper: top annotation (depth bar track + metadata strips) ─────────────
# depth — library size per sample (reads before the prevalence filter), aligned
#         to the heatmap columns. Unbagged flowers are sequenced deeper than
#         bagged (see 08a); if a dendrogram split follows this track, depth or
#         zero replacement, not biology, may be driving the clustering.

build_annotation <- function(meta, depth, cols) {
    HeatmapAnnotation(
        log10_depth = anno_barplot(
            log10(depth),
            baseline = "min",                         # bars start at the shallowest sample
            height   = unit(8, "mm"),                 # height of the depth track
            gp       = gpar(fill = "grey40", col = NA),
            border   = FALSE,
            axis_param = list(gp = gpar(fontsize = 5))
        ),
        sample_type         = safe_chr(meta$sample_type),
        management_type     = safe_chr(meta$management_type),
        farm_id             = safe_chr(meta$farm_id),
        tree_id             = safe_chr(meta$tree_id),
        is_pollination_clsm = safe_chr(meta$is_pollination_clsm),
        col                  = cols,
        na_col               = NA_COL,
        show_legend          = FALSE,                 # one shared legend is drawn separately
        show_annotation_name = TRUE,
        annotation_name_side = "left",                # right side collides with the CLR legend
        annotation_name_rot  = 0,                       # horizontal names; avoids overlap with the depth track
        annotation_name_gp   = gpar(fontsize = 7),
        simple_anno_size     = unit(3, "mm"),         # height of each metadata strip
        gap                  = unit(0.5, "mm")
    )
}

# ── 8. Helper: shared annotation legend ───────────────────────────────────────
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

# ── 9. Helper: build one Heatmap object ───────────────────────────────────────
# clr_mat     — samples x taxa CLR matrix (from prepare_clr)
# meta, depth — sample metadata and library sizes (any order; aligned here)
# link_method — hclust linkage used for both samples and taxa

build_heatmap <- function(clr_mat, meta, depth, cols, title, hm_name, link_method) {

    mat  <- t(clr_mat)   # taxa x samples
    meta <- meta[colnames(mat), , drop = FALSE]
    ha   <- build_annotation(meta, depth[colnames(mat)], cols)

    # Symmetric CLR colour scale clipped at 1st/99th percentile
    clr_range <- quantile(mat, c(0.01, 0.99), na.rm = TRUE)
    abs_max   <- max(abs(clr_range))
    col_fun   <- colorRamp2(c(-abs_max, 0, abs_max),
                            c("#313695", "white", "#a50026"))   # blue–white–red

    subtitle <- sprintf(
        "prev ≥ %.0f%%  |  zero repl: %s  |  dist: %s (Aitchison)  |  linkage: %s  |  %d taxa x %d samples",
        100 * PREV_THRESH, attr(clr_mat, "zero_method"), DIST_METHOD, link_method,
        nrow(mat), ncol(mat)
    )

    Heatmap(
        mat,
        name            = hm_name,                    # must be unique on the device
        column_title    = paste0(title, "\n", subtitle),
        column_title_gp = gpar(fontsize = 8),
        top_annotation  = ha,
        clustering_method_rows      = link_method,
        clustering_method_columns   = link_method,
        clustering_distance_rows    = DIST_METHOD,
        clustering_distance_columns = DIST_METHOD,
        col               = col_fun,
        show_row_names    = FALSE,
        show_column_names = FALSE,
        row_dend_width     = unit(10, "mm"),          # width of the taxa dendrogram
        column_dend_height = unit(15, "mm"),          # height of the sample dendrogram
        heatmap_legend_param = list(
            title         = "CLR",
            direction     = "vertical",
            legend_height = unit(30, "mm")
        ),
        use_raster     = TRUE,
        raster_quality = 3
    )
}

# ── 10. Run pipeline for both markers ─────────────────────────────────────────
# CLR is computed once per marker; only the linkage differs between rows.

cat("\n── 16S bacteria ─────────────────────────────────────────────────────────\n")
ps_16S_filt <- filter_taxa_ps(ps_16S, PREV_THRESH, "16S")
clr_16S     <- prepare_clr(ps_16S_filt, label = "16S")

cat("\n── ITS1 OTU97 fungi ─────────────────────────────────────────────────────\n")
ps_ITS_filt <- filter_taxa_ps(ps_ITS, PREV_THRESH, "ITS1")
clr_ITS     <- prepare_clr(ps_ITS_filt, label = "ITS1")

# Metadata, library size (before the prevalence filter) and colours
meta  <- as(sample_data(ps_16S), "data.frame")
cols  <- build_annot_colors(meta)
depth_16S <- sample_sums(ps_16S)
depth_ITS <- sample_sums(ps_ITS)

# ── 11. Build heatmap objects: one per marker x linkage ───────────────────────

hm_list <- lapply(names(LINK_METHODS), \(lk) list(
    `16S` = build_heatmap(
        clr_16S, meta, depth_16S, cols,
        title       = sprintf("16S bacteria  (%d ASVs)  —  %s", ncol(clr_16S), LINK_METHODS[[lk]]),
        hm_name     = paste0("CLR_16S_", lk),
        link_method = LINK_METHODS[[lk]]
    ),
    ITS1 = build_heatmap(
        clr_ITS, meta, depth_ITS, cols,
        title       = sprintf("ITS1 OTU97 fungi  (%d OTUs)  —  %s", ncol(clr_ITS), LINK_METHODS[[lk]]),
        hm_name     = paste0("CLR_ITS1_", lk),
        link_method = LINK_METHODS[[lk]]
    )
)) |> setNames(names(LINK_METHODS))

annot_legend <- build_annot_legend(cols)

# ── 12. Save — 2 x 2 heatmaps + legend column ─────────────────────────────────

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
out_fp <- here("results", "figures", "26c_16S_ITS1_clr_clustermap.png")

# Figure size in inches; fonts are in points, so a larger canvas = smaller text
png(out_fp, width = 26, height = 24, units = "in", res = 300)

grid.newpage()

# Main title
grid.text(
    "16S bacteria vs ITS1 fungi — CLR clustermap (top: Ward.D2, bottom: UPGMA)",
    x  = 0.5, y = 0.985,
    gp = gpar(fontsize = 14, fontface = "bold")
)

# Grid layout: 2 rows (linkages) x 3 columns (16S, ITS1, legend).
# widths controls the relative width of each column; the last one is the legend.
pushViewport(viewport(
    y = 0.48, height = 0.96,
    layout = grid.layout(nrow = 2, ncol = 3,
                         widths = unit(c(1, 1, 0.12), "null"))
))

for (i in seq_along(hm_list)) {
    for (j in seq_along(hm_list[[i]])) {
        pushViewport(viewport(layout.pos.row = i, layout.pos.col = j))
        draw(hm_list[[i]][[j]], newpage = FALSE,
             padding = unit(c(2, 5, 5, 5), "mm"))    # bottom, left, top, right
        popViewport()
    }
}

# Shared annotation legend, spanning both rows of the third column
pushViewport(viewport(layout.pos.row = 1:2, layout.pos.col = 3))
draw(annot_legend, x = unit(2, "mm"), y = unit(0.5, "npc"), just = "left")
popViewport()

popViewport()
dev.off()

cat(sprintf("\nSaved: %s\n", out_fp))
