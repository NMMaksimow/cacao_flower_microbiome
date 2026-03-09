library(phyloseq)
library(ComplexHeatmap)
library(circlize)    # colorRamp2
library(vegan)       # vegdist, decostand
library(dplyr)
library(RColorBrewer)

# ── Предполагается что phyloseq объект называется ps ──────────────────────────
# ps <- readRDS("results/rds/ps_bacteria_decontam.rds")

# ── Цветовые палитры ──────────────────────────────────────────────────────────
sample_type_colors <- c(
        biological       = "#4C72B0",
        extraction_ctrl  = "#55A868",
        negative_ctrl    = "#C44E52",
        positive_ctrl    = "#8172B2"
)
flower_colors <- c(closed = "#4C72B0", open = "#DD8452")

farm_ids  <- sort(unique(sample_data(ps)$farm_id))
farm_pal  <- setNames(brewer.pal(max(3, length(farm_ids)), "Set2")[seq_along(farm_ids)], farm_ids)


# ── Подготовка данных ─────────────────────────────────────────────────────────
ps_ra <- transform_sample_counts(ps, function(x) x / sum(x))  # TSS

otu   <- as.matrix(otu_table(ps_ra))                 # ASV × samples
if (!taxa_are_rows(ps_ra)) otu <- t(otu)             # → ASV × samples

# Top-50 по среднему обилию
top50 <- names(sort(rowMeans(otu), decreasing = TRUE))[1:50]
otu50 <- log10(otu[top50, ] + 1e-6)                  # log10 transform

# Genus-подписи из tax_table
tax   <- as.data.frame(tax_table(ps_ra))
get_genus <- function(asv) {
        g <- tax[asv, "Genus"]
        ifelse(is.na(g) | g == "", "unclassified", gsub("g__", "", g))
}
row_labels <- make.unique(sapply(top50, get_genus))

# ── Аннотации образцов (колонки) ──────────────────────────────────────────────
sdf <- as.data.frame(sample_data(ps_ra))

col_ann <- HeatmapAnnotation(
        sample_type = sdf$sample_type,
        farm        = sdf$farm_id,
        flower_type = sdf$flower_type,
        col = list(
                sample_type = sample_type_colors,
                farm        = farm_pal,
                flower_type = flower_colors
        ),
        annotation_legend_param = list(
                sample_type = list(title = "sample type", title_gp = gpar(fontsize = 10)),
                farm        = list(title = "farm",        title_gp = gpar(fontsize = 10)),
                flower_type = list(title = "flower type", title_gp = gpar(fontsize = 10))
        ),
        show_annotation_name = TRUE
)

# ── Тепловая карта ────────────────────────────────────────────────────────────
col_fun <- colorRamp2(
        c(min(otu50), median(otu50), max(otu50)),
        c("#FFFFCC", "#FD8D3C", "#800026")   # YlOrRd-style
)

ht_abund <- Heatmap(
        otu50,
        name                  = "log₁₀(RA)",
        col                   = col_fun,
        top_annotation        = col_ann,
        clustering_method_rows    = "average",
        clustering_method_columns = "average",
        clustering_distance_rows    = "euclidean",
        clustering_distance_columns = "euclidean",
        show_row_names        = TRUE,
        show_column_names     = FALSE,
        row_names_gp          = gpar(fontsize = 8),
        row_labels            = row_labels,
        column_title          = sprintf(
                "Relative abundance — %d samples × top-50 ASVs\n(log₁₀, Euclidean, average linkage)",
                ncol(otu50)
        ),
        column_title_gp       = gpar(fontsize = 11),
        heatmap_legend_param  = list(
                title    = "log₁₀(RA + 1e-6)",
                at       = c(-6, -4, -2, 0),
                direction = "horizontal"
        )
)

draw(ht_abund,
     heatmap_legend_side   = "bottom",
     annotation_legend_side = "right",
     merge_legend           = FALSE)


# ── BC матрица (на относительных обилиях) ────────────────────────────────────
otu_t   <- t(otu_table(ps_ra))          # samples × ASVs
if (taxa_are_rows(ps_ra)) otu_t <- t(otu_t)

bc_dist <- as.matrix(vegdist(otu_t, method = "bray"))

# Аннотации строк И колонок одинаковые
rc_ann <- rowAnnotation(
        sample_type = sdf$sample_type,
        farm        = sdf$farm_id,
        col         = list(sample_type = sample_type_colors, farm = farm_pal),
        show_annotation_name = TRUE,
        annotation_legend_param = list(
                sample_type = list(title = "sample type"),
                farm        = list(title = "farm")
        )
)
ca_ann <- HeatmapAnnotation(
        sample_type = sdf$sample_type,
        farm        = sdf$farm_id,
        col         = list(sample_type = sample_type_colors, farm = farm_pal),
        show_annotation_name = FALSE,
        show_legend           = FALSE
)

col_bc <- colorRamp2(c(0, 0.5, 1), c("#FFFFCC", "#FD8D3C", "#800026"))

ht_bc <- Heatmap(
        bc_dist,
        name                  = "BC distance",
        col                   = col_bc,
        top_annotation        = ca_ann,
        right_annotation      = rc_ann,
        clustering_method_rows    = "average",
        clustering_method_columns = "average",
        show_row_names        = FALSE,
        show_column_names     = FALSE,
        column_title          = sprintf(
                "Bray-Curtis dissimilarity — %d samples (average linkage)", nrow(bc_dist)
        ),
        column_title_gp       = gpar(fontsize = 11),
        heatmap_legend_param  = list(
                title     = "Bray-Curtis",
                at        = c(0, 0.25, 0.5, 0.75, 1),
                direction = "horizontal"
        )
)

draw(ht_bc,
     heatmap_legend_side    = "bottom",
     annotation_legend_side = "right")


# ── CLR трансформация через vegan::decostand ──────────────────────────────────
# Псевдосчёт +1 перед TSS, потом CLR
ps_pseudo <- transform_sample_counts(ps, function(x) x + 1)
ps_clr_ra <- transform_sample_counts(ps_pseudo, function(x) x / sum(x))

otu_clr_t <- t(as.matrix(otu_table(ps_clr_ra)))
if (taxa_are_rows(ps_clr_ra)) otu_clr_t <- t(otu_clr_t)

clr_mat     <- decostand(otu_clr_t, method = "clr")   # samples × ASVs
ait_dist    <- as.matrix(dist(clr_mat, method = "euclidean"))

col_ait <- colorRamp2(
        quantile(ait_dist[upper.tri(ait_dist)], c(0, 0.5, 1)),
        c("#FFFFCC", "#FD8D3C", "#800026")
)

ht_ait <- Heatmap(
        ait_dist,
        name                  = "Aitchison",
        col                   = col_ait,
        top_annotation        = ca_ann,       # те же аннотации что у BC
        right_annotation      = rc_ann,
        clustering_method_rows    = "average",
        clustering_method_columns = "average",
        show_row_names        = FALSE,
        show_column_names     = FALSE,
        column_title          = sprintf(
                "Aitchison distance (CLR + Euclidean) — %d samples", nrow(ait_dist)
        ),
        column_title_gp       = gpar(fontsize = 11),
        heatmap_legend_param  = list(
                title     = "Aitchison distance",
                direction = "horizontal"
        )
)

draw(ht_ait,
     heatmap_legend_side    = "bottom",
     annotation_legend_side = "right")


