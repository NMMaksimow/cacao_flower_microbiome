## ============================================================================
## 19a — Network visualization: 16S bacteria (bagged vs unbagged) — SpiecEasi-MB
##
## Standalone visualization script: loads 19a_network_results.rds and renders
## a 7 x 2 patchwork figure. Run this locally after downloading results from HPC.
## Does not require NetCoMi, SpiecEasi, or SPRING to be installed.
##
## Input:
##   results/rds/19a_network_results.rds
##   results/rds/ps_16S_bacteria_biosamples.rds
## Output:
##   results/figures/19a_network_panels_spieceasi.png
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

PREV_THRESH_NET <- 0.20    # 20% per group; used in caption only
N_PERM          <- 1000    # netCompare permutations; used in caption only
STARS_THRESH    <- 0.1     # STARS stability threshold; used in caption only
AGGREGATE_GENUS <- TRUE    # must match the compute script
MIN_LIB         <- 500
FARM_LEVELS     <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
EDGE_COL_POS    <- "#D6604D"   # warm = positive correlation
EDGE_COL_NEG    <- "#4393C3"   # cool = negative correlation
PHY_NA_COLOUR   <- "grey70"    # fallback colour: unassigned phylum, or overflow bucket
MAX_PHY_COLOURS <- 20L         # D3 category20 capacity; beyond this the tail is pooled

# ── 1. Libraries ──────────────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(ggraph)
library(igraph)
library(patchwork)
library(ggsci)

"%||%" <- function(a, b) if (!is.null(a)) a else b

ps_raw <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds"))
cat(sprintf("Loaded: %d taxa x %d samples\n", ntaxa(ps_raw), nsamples(ps_raw)))

# ── 2. Visualization helpers ──────────────────────────────────────────────────

# 16S (SILVA): falls back through Family -> Order -> Class -> Phylum before ID prefix
make_tax_label <- function(physeq, feature_ids) {
        tt <- as.data.frame(tax_table(physeq))
        vapply(feature_ids, function(id) {
                if (!id %in% rownames(tt)) return(substr(id, 1, 12))
                g   <- as.character(tt[id, "Genus"])
                f   <- as.character(tt[id, "Family"])
                o   <- as.character(tt[id, "Order"])
                cls <- as.character(tt[id, "Class"])
                p   <- as.character(tt[id, "Phylum"])
                bad <- function(x) is.na(x) || nchar(x) <= 3
                if (!bad(g)) g
                else if (!bad(f))   paste0(f,   " (fam.)")
                else if (!bad(o))   paste0(o,   " (ord.)")
                else if (!bad(cls)) paste0(cls, " (cls.)")
                else if (!bad(p))   paste0(p,   " (phy.)")
                else substr(id, 1, 12)
        }, character(1))
}

net_to_igraph <- function(net_raw, which = 1) {
        adj  <- if (which == 1) net_raw$adjaMat1 else net_raw$adjaMat2
        # NetCoMi stores adjaMat as a NON-NEGATIVE transform of the association, so
        # the sign of each edge survives only in assoMat. Taking the sign off adjaMat
        # makes every edge look positive.
        asso <- if (which == 1) net_raw$assoMat1 else net_raw$assoMat2
        g <- igraph::graph_from_adjacency_matrix(
                abs(adj), mode = "undirected", weighted = TRUE, diag = FALSE
        )
        el <- igraph::as_edgelist(g, names = FALSE)
        if (nrow(el) > 0) igraph::E(g)$sign <- asso[el]
        g
}

# ggraph panel for one farm x treatment network.
# Nodes coloured by phylum, hub taxa ringed in black and labelled.
# top_phy: phyla that have their own colour (normally all phyla in the figure).
# phy_palette: named colour vector; its names define the factor levels, so every
# panel emits the same set of levels regardless of which phyla it contains.
plot_net_panel <- function(farm_result, which_net, farm, treatment_label,
                           layout_algo = "fr", phy_palette, top_phy) {
        title <- sprintf("%s — %s", toupper(farm), treatment_label)

        if (is.null(farm_result)) {
                return(
                        ggplot() + theme_void(base_size = 8) +
                                annotate("text", x = 0.5, y = 0.5, hjust = 0.5,
                                         label = "No network", size = 2.5, colour = "grey50") +
                                labs(title = title) +
                                theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5))
                )
        }

        g <- net_to_igraph(farm_result$raw, which = which_net)

        hub_raw <- tryCatch(
                if (which_net == 1) farm_result$analyzed$hubs$hubs1
                else                farm_result$analyzed$hubs$hubs2,
                error = function(e) NULL
        )
        # NetCoMi returns hubs as an UNNAMED character vector of taxon IDs;
        # handle both forms so labels and borders are drawn either way
        hub_names <- if (is.null(hub_raw))          character(0)
                     else if (is.null(names(hub_raw))) as.character(hub_raw)
                     else                              names(hub_raw)

        tt_raw  <- as.data.frame(tax_table(ps_raw))
        phy_raw <- tt_raw[igraph::V(g)$name, "Phylum"]
        # Values outside top_phy only occur in the overflow case; otherwise every
        # phylum has its own level. NA stays NA and is drawn with na.value, so an
        # unannotated node is never silently merged into a named phylum.
        # factor() over the full palette keeps levels identical across panels.
        igraph::V(g)$phylum    <- factor(
                ifelse(is.na(phy_raw) | phy_raw %in% top_phy,
                       phy_raw, setdiff(names(phy_palette), top_phy)[1]),
                levels = names(phy_palette)
        )
        igraph::V(g)$degree    <- igraph::degree(g)
        igraph::V(g)$is_hub    <- igraph::V(g)$name %in% hub_names

        igraph::V(g)$is_labeled <- igraph::V(g)$is_hub
        igraph::V(g)$tax_label  <- ifelse(
                igraph::V(g)$is_hub,
                make_tax_label(ps_raw, igraph::V(g)$name),
                ""
        )

        if (igraph::ecount(g) == 0) {
                return(
                        ggraph(g, layout = "circle") +
                                geom_node_point(aes(size = degree, colour = phylum)) +
                                # drop = FALSE keeps unused levels so the colour mapping
                                # matches the shared legend panel exactly
                                scale_colour_manual(values = phy_palette, name = "Phylum",
                                                    na.value = PHY_NA_COLOUR, drop = FALSE) +
                                scale_size(range = c(1, 4), guide = "none") +
                                labs(title = paste0(title, " (no edges)")) +
                                theme_graph(base_size = 8) +
                                theme(plot.title = element_text(size = 7, face = "bold",
                                                                hjust = 0.5))
                )
        }

        set.seed(123)

        ggraph(g, layout = layout_algo) +
                geom_edge_link(aes(colour = sign > 0),
                               alpha = 0.55, width = 0.5, show.legend = FALSE) +
                geom_node_point(aes(size = degree, colour = phylum)) +
                geom_node_point(aes(filter = is_hub, size = degree),
                                shape = 21, colour = "black", fill = NA,
                                stroke = 1.5, show.legend = FALSE) +
                geom_node_text(aes(filter = is_labeled, label = tax_label),
                               size = 2.5, repel = TRUE, max.overlaps = 20,
                               colour = "black", fontface = "italic") +
                scale_edge_colour_manual(
                        values = c("TRUE" = EDGE_COL_POS, "FALSE" = EDGE_COL_NEG),
                        guide  = "none"
                ) +
                scale_colour_manual(
                        values   = phy_palette,
                        name     = "Phylum",
                        na.value = PHY_NA_COLOUR,
                        drop     = FALSE   # keep all levels so guides stay identical
                ) +
                scale_size(range = c(1, 4), guide = "none") +
                labs(title = title) +
                theme_graph(base_size = 8) +
                theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5))
}

# Standalone legend panel.
# patchwork's plot_layout(guides = "collect") does not reliably merge guides under
# ggplot2 4.x, so the legend is built once here as its own plot and placed beside
# the grid. Points are drawn with alpha = 0 (invisible in the panel) while
# override.aes restores full opacity in the legend keys themselves.
make_legend_panel <- function(phy_palette) {
        d <- data.frame(
                x      = 1,
                y      = seq_along(phy_palette),
                phylum = factor(names(phy_palette), levels = names(phy_palette))
        )
        ggplot(d, aes(x, y, colour = phylum)) +
                geom_point(alpha = 0) +
                scale_colour_manual(values = phy_palette, name = "Phylum", drop = FALSE) +
                guides(colour = guide_legend(
                        ncol         = 1,
                        override.aes = list(alpha = 1, size = 3)   # legend key styling
                )) +
                theme_void() +
                theme(
                        legend.position = "left",           # hug the panel grid
                        legend.title    = element_text(size = 10, face = "bold"),
                        legend.text     = element_text(size = 8),
                        legend.key.size = unit(0.5, "lines")
                )
}

make_panel_grid <- function(layout_algo) {
        p_list <- vector("list", 14)
        for (i in seq_along(FARM_LEVELS)) {
                farm <- FARM_LEVELS[i]
                p_list[[i]]     <- plot_net_panel(results[[farm]], 1, farm, "bagged",
                                                  layout_algo = layout_algo,
                                                  phy_palette = phy_cols,
                                                  top_phy     = top_phy)
                p_list[[i + 7]] <- plot_net_panel(results[[farm]], 2, farm, "unbagged",
                                                  layout_algo = layout_algo,
                                                  phy_palette = phy_cols,
                                                  top_phy     = top_phy)
        }
        # Suppress every panel's own legend; the shared one is added separately
        grid_part <- wrap_plots(p_list, nrow = 2, ncol = 7) &
                theme(legend.position = "none")

        # widths: grid gets 20 parts, legend column 1 part
        wrap_plots(grid_part, make_legend_panel(phy_cols),
                   nrow = 1, widths = c(20, 1)) +
                plot_annotation(
                        title   = title_str,
                        caption = caption_str,
                        theme   = theme(
                                plot.title   = element_text(size = 11, face = "bold"),
                                plot.caption = element_text(size = 7,  colour = "grey40")
                        )
                )
}

# ── 3. Load RDS + build palette ───────────────────────────────────────────────

cat("\n── 3. Loading results and building palette ───────────────────────────────\n")

results <- readRDS(here("results", "rds", "19a_network_results.rds"))
cat("  Loaded: 19a_network_results.rds\n")

# The palette is built from the taxa that are actually drawn as nodes, not from the
# full filtered phyloseq. Only taxa passing the per-farm prevalence filter reach a
# network, so this is what makes the legend describe exactly what the figure shows:
# every legend entry occurs in at least one panel, and every node has an entry.
net_nodes <- unique(unlist(lapply(results, function(r) {
        if (is.null(r)) return(NULL)
        union(colnames(r$raw$adjaMat1), colnames(r$raw$adjaMat2))
})))

tt_raw_all <- as.data.frame(tax_table(ps_raw))
phy_nodes  <- tt_raw_all[net_nodes, "Phylum"]
# Ranked by how many nodes carry the phylum, so the legend reads most- to least-common
phy_counts <- sort(table(phy_nodes[!is.na(phy_nodes)]), decreasing = TRUE)
n_phy      <- length(phy_counts)

if (n_phy <= MAX_PHY_COLOURS) {
        # Normal case: every phylum gets its own colour, nothing is pooled
        top_phy  <- names(phy_counts)
        phy_cols <- setNames(ggsci::pal_d3("category20")(n_phy), top_phy)
} else {
        # Defensive fallback if a future dataset exceeds the palette capacity
        top_phy  <- names(phy_counts)[seq_len(MAX_PHY_COLOURS - 1L)]
        phy_cols <- c(setNames(ggsci::pal_d3("category20")(length(top_phy)), top_phy),
                      "Other phyla" = PHY_NA_COLOUR)
}

cat(sprintf("  Nodes drawn: %d  |  phyla in figure: %d  |  pooled: %s\n",
            length(net_nodes), n_phy,
            if (n_phy <= MAX_PHY_COLOURS) "none" else
                    paste(n_phy - length(top_phy), "into 'Other phyla'")))
cat(sprintf("  Nodes with no phylum assignment: %d\n", sum(is.na(phy_nodes))))

# ── 4. Render and save ────────────────────────────────────────────────────────

cat("\n── 4. Building network panels ───────────────────────────────────────────\n")

title_str   <- "16S Bacteria — SpiecEasi-MB co-occurrence networks by farm (bagged vs unbagged)"
caption_str <- sprintf(
        "SpiecEasi-MB  |  prevalence >= %.0f%%  |  STARS threshold = %.2f  |  N = %d permutations  |  hub nodes = black border, italic label",
        100 * PREV_THRESH_NET, STARS_THRESH, N_PERM
)

fig <- make_panel_grid("fr")

dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(
        here("results", "figures", "19a_network_panels_spieceasi.png"),
        plot   = fig,
        width  = 24,
        height = 14,
        dpi    = 300
)
cat("  Saved: 19a_network_panels_spieceasi.png\n")
