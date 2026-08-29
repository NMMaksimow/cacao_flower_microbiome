## ============================================================================
## 19b — Co-occurrence networks by farm: ITS1 OTU97 fungi (bagged vs unbagged)
##
## Identical logic to 19a (16S). Differences:
##   - Input: ps_ITS1_otu97_fungi_biosamples.rds
##   - make_tax_label(): skips Incertae_sedis at every rank (UNITE taxonomy)
##   - Output prefix: 19b_
##
## Builds paired SpiecEasi-MB networks per farm (bagged = network 1,
## unbagged = network 2) using NetComi. Compares network topology between
## treatments within each farm.
##
## Prevalence filter 20% applied inside netConstruct() per group.
##
## Input:
##   results/rds/ps_ITS1_otu97_fungi_biosamples.rds
## Output:
##   results/rds/19b_network_results.rds
##   results/tables/19b_network_topology.csv
##   results/tables/19b_network_hubs.csv
##   results/tables/19b_network_comparison.csv
##
## Compute only — no figures. Render panels locally with 19b_ITS1_network_SpiecEasi_vis.R
## ============================================================================

# ── 0. Parameters ─────────────────────────────────────────────────────────────

PREV_THRESH_NET    <- 0.20
MAX_TAXA_NET       <- 30L     # cap per farm before netConstruct; keeps p/n ≤ 1.4 (SpiecEasi-MB needs p < ~2n)
MIN_LIB            <- 500
FARM_LEVELS        <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")
N_PERM             <- 200    # netCompare permutations; min resolvable p = 1/200 = 0.005
SPIEC_NLAMBDA      <- 20
SPIEC_LAMBDA_RATIO <- 1e-2
STARS_THRESH       <- 0.1     # STARS stability threshold (default 0.05 is too strict for n=21)
AGGREGATE_GENUS    <- TRUE    # aggregate OTUs to genus level before network construction
# MAX_TAXA_NET    <- 30L     # disabled: genus aggregation already reduces p sufficiently

# ── 1. Libraries & data ───────────────────────────────────────────────────────

library(here)
library(tidyverse)
library(phyloseq)
library(NetCoMi)
library(igraph)

"%||%" <- function(a, b) if (!is.null(a)) a else b

ps_raw <- readRDS(here("results", "rds", "ps_ITS1_otu97_fungi_biosamples.rds"))
cat(sprintf("Loaded: %d taxa x %d samples\n", ntaxa(ps_raw), nsamples(ps_raw)))

# ── 2. Helper functions ───────────────────────────────────────────────────────

# ITS1 (UNITE): skips Incertae_sedis at every rank; falls back Genus → Family →
# Order → Class → Phylum before using ID prefix
make_tax_label <- function(physeq, feature_ids) {
        tt <- as.data.frame(tax_table(physeq))
        vapply(feature_ids, function(id) {
                if (!id %in% rownames(tt)) return(substr(id, 1, 12))
                g   <- as.character(tt[id, "Genus"])
                f   <- as.character(tt[id, "Family"])
                o   <- as.character(tt[id, "Order"])
                cls <- as.character(tt[id, "Class"])
                p   <- as.character(tt[id, "Phylum"])
                bad <- function(x) is.na(x) || grepl("Incertae_sedis", x) || nchar(x) <= 3
                if (!bad(g)) g
                else if (!bad(f))   paste0(f,   " (fam.)")
                else if (!bad(o))   paste0(o,   " (ord.)")
                else if (!bad(cls)) paste0(cls, " (cls.)")
                else if (!bad(p))   paste0(p,   " (phy.)")
                else substr(id, 1, 12)
        }, character(1))
}
# Edge set implied by an adjacency matrix. NetCoMi puts 1 on the diagonal
# (self-association), so it has to be dropped or every taxon counts as an edge and
# every edge total comes out inflated by half the number of taxa.
edge_mask <- function(adj) {
        e <- adj != 0
        diag(e) <- FALSE
        e
}

# Topology metrics from adjacency matrices and netAnalyze global property slots.
extract_topology <- function(net_raw, net_ana, farm) {
        gp   <- net_ana$globalProps     # nComp1/2, avPath1/2, clustCoef1/2, modularity1/2
        gpl  <- net_ana$globalPropsLCC  # lccSize1/2, lccSizeRel1/2
        cent <- net_ana$centralities    # degree1, degree2 (NetCoMi 1.3.0 API)

        # adjaMat gives the edge set; its entries are non-negative by construction,
        # so the sign of each association has to come from assoMat.
        make_row <- function(adj, asso, deg_vec, s, treatment) {
                edge    <- edge_mask(adj)
                n_edges <- sum(edge) / 2
                n_pos   <- sum(edge & asso > 0) / 2
                n_neg   <- sum(edge & asso < 0) / 2
                p       <- nrow(adj)
                max_e   <- p * (p - 1) / 2
                dplyr::tibble(
                        farm         = farm,
                        treatment    = treatment,
                        n_taxa       = p,
                        n_edges      = as.integer(n_edges),
                        pos_edges    = as.integer(n_pos),
                        neg_edges    = as.integer(n_neg),
                        pos_edge_pct = if (n_edges > 0) n_pos / n_edges else NA_real_,
                        edge_density = if (max_e > 0) n_edges / max_e else NA_real_,
                        avg_degree   = mean(deg_vec, na.rm = TRUE),
                        n_comp       = gp[[paste0("nComp",      s)]],
                        lcc_size_rel = gpl[[paste0("lccSizeRel", s)]],
                        clust_coef   = gp[[paste0("clustCoef",  s)]],
                        modularity   = gp[[paste0("modularity", s)]],
                        avg_path     = gp[[paste0("avPath",     s)]]
                )
        }
        dplyr::bind_rows(
                make_row(net_raw$adjaMat1, net_raw$assoMat1, cent$degree1, "1", "bagged_flower"),
                make_row(net_raw$adjaMat2, net_raw$assoMat2, cent$degree2, "2", "unbagged_flower")
        )
}

extract_hubs <- function(net_ana, farm) {
        cent <- net_ana$centralities  # degree1, degree2, eigenv1, eigenv2 (NetCoMi 1.3.0 API)
        make_hub_rows <- function(hubs, deg, eig, treatment) {
                if (is.null(hubs) || length(hubs) == 0) return(dplyr::tibble())
                hub_ids <- if (is.null(names(hubs))) as.character(hubs) else names(hubs)
                dplyr::tibble(
                        farm        = farm,
                        treatment   = treatment,
                        taxon       = hub_ids,
                        tax_label   = make_tax_label(ps_raw, hub_ids),
                        degree      = deg[hub_ids],
                        eigenvector = eig[hub_ids]
                )
        }
        dplyr::bind_rows(
                make_hub_rows(net_ana$hubs$hubs1, cent$degree1, cent$eigenv1, "bagged_flower"),
                make_hub_rows(net_ana$hubs$hubs2, cent$degree2, cent$eigenv2, "unbagged_flower")
        )
}

# Global network properties from netCompare → tidy tibble with permutation p-values.
# Columns: farm, scope, property, val_bagged, val_unbagged, stat, pval, significant.
# scope "whole" and "LCC": stat = |val_unbagged - val_bagged|.
# scope "centrality_overlap" (Jaccard, ARI): val_bagged = similarity statistic, val_unbagged = NA.
extract_comparison <- function(net_cmp, net_ana, net_raw, farm) {
        if (is.null(net_cmp) || is.null(net_ana)) return(NULL)

        adj1  <- net_raw$adjaMat1
        adj2  <- net_raw$adjaMat2
        gp    <- net_ana$globalProps
        gpl   <- net_ana$globalPropsLCC
        pvg   <- net_cmp$pvalDiffGlobal
        pvgl  <- net_cmp$pvalDiffGlobalLCC
        lcc_p <- net_cmp$propertiesLCC  # density1/2, vertConnect1/2, edgeConnect1/2, natConnect1/2

        make_density <- function(adj) {
                p <- nrow(adj); n_e <- sum(edge_mask(adj)) / 2; max_e <- p * (p - 1) / 2
                if (max_e > 0) n_e / max_e else NA_real_
        }
        # Sign from assoMat, edge set from adjaMat (see extract_topology).
        make_pep <- function(adj, asso) {
                edge <- edge_mask(adj)
                if (sum(edge) > 0) sum(edge & asso > 0) / sum(edge) else NA_real_
        }

        # Slot layout varies: pvalDiffGlobal is a named list for SpiecEasi but a named
        # atomic vector for SparCC/SPRING, and with permTest = FALSE the p-value slots
        # may be absent altogether. Pull one scalar defensively so a missing slot
        # becomes NA instead of shortening the vector and breaking the tibble.
        pull1 <- function(x, nm) {
                v <- if (is.null(x)) NULL else x[[nm]]
                if (is.null(v) || length(v) != 1L) NA_real_ else as.numeric(v)
        }
        pulln <- function(x, nms) vapply(nms, function(n) pull1(x, n), numeric(1))

        whole <- dplyr::tibble(
                scope    = "whole",
                property = c("nComp", "avPath", "clustCoef", "modularity",
                             "vertConnect", "avDiss", "density", "PEP"),
                val_b    = c(gp$nComp1, gp$avPath1, gp$clustCoef1, gp$modularity1,
                             gp$vertConnect1, gp$avDiss1, make_density(adj1),
                             make_pep(adj1, net_raw$assoMat1)),
                val_u    = c(gp$nComp2, gp$avPath2, gp$clustCoef2, gp$modularity2,
                             gp$vertConnect2, gp$avDiss2, make_density(adj2),
                             make_pep(adj2, net_raw$assoMat2)),
                pval     = pulln(pvg, c("pvalnComp", "pvalavPath", "pvalClustCoef", "pvalModul",
                                        "pvalVertConnect", "pvalavDiss", "pvalDensity", "pvalPEP"))
        )

        lcc <- dplyr::tibble(
                scope    = "LCC",
                property = c("lccSize", "lccSizeRel", "density", "avPath", "clustCoef",
                             "modularity", "vertConnect", "edgeConnect", "natConnect"),
                val_b    = c(gpl$lccSize1, gpl$lccSizeRel1, lcc_p$density1, gpl$avPath1,
                             gpl$clustCoef1, gpl$modularity1, lcc_p$vertConnect1,
                             lcc_p$edgeConnect1, lcc_p$natConnect1),
                val_u    = c(gpl$lccSize2, gpl$lccSizeRel2, lcc_p$density2, gpl$avPath2,
                             gpl$clustCoef2, gpl$modularity2, lcc_p$vertConnect2,
                             lcc_p$edgeConnect2, lcc_p$natConnect2),
                pval     = pulln(pvgl, c("pvallccSize", "pvallccSizeRel", "pvalDensity",
                                         "pvalavPath", "pvalClustCoef", "pvalModul",
                                         "pvalVertConnect", "pvalEdgeConnect", "pvalNatConnect"))
        )

        overlap <- dplyr::tibble(
                scope    = "centrality_overlap",
                property = c("jacc_degree", "jacc_between", "jacc_close",
                             "jacc_eigen", "jacc_hub", "ARI", "ARI_LCC"),
                val_b    = c(pull1(net_cmp$jaccDeg,   "jacc"), pull1(net_cmp$jaccBetw, "jacc"),
                             pull1(net_cmp$jaccClose, "jacc"), pull1(net_cmp$jaccEigen, "jacc"),
                             pull1(net_cmp$jaccHub,   "jacc"),
                             pull1(net_cmp$randInd, "value"), pull1(net_cmp$randIndLCC, "value")),
                val_u    = NA_real_,
                pval     = c(pull1(net_cmp$jaccDeg,   "p.greater"), pull1(net_cmp$jaccBetw, "p.greater"),
                             pull1(net_cmp$jaccClose, "p.greater"), pull1(net_cmp$jaccEigen, "p.greater"),
                             pull1(net_cmp$jaccHub,   "p.greater"),
                             pull1(net_cmp$randInd, "pval"), pull1(net_cmp$randIndLCC, "pval"))
        )

        dplyr::bind_rows(whole, lcc, overlap) |>
                dplyr::transmute(
                        farm         = farm,
                        scope, property,
                        val_bagged   = val_b,
                        val_unbagged = val_u,
                        stat         = dplyr::if_else(!is.na(val_b) & !is.na(val_u),
                                                       abs(val_u - val_b), val_b),
                        pval,
                        significant  = dplyr::if_else(is.na(pval), NA, pval < 0.05)
                )
}
# ── 3. Global subset ──────────────────────────────────────────────────────────

cat("\n── 3. Global subset ─────────────────────────────────────────────────────\n")

ps_bvu <- ps_raw |>
        subset_samples(sample_type %in% c("bagged_flower", "unbagged_flower")) |>
        (\(x) prune_samples(sample_sums(x) >= MIN_LIB, x))() |>
        (\(x) prune_taxa(taxa_sums(x) > 0, x))()

cat(sprintf("  bagged + unbagged after lib filter: %d taxa x %d samples\n",
            ntaxa(ps_bvu), nsamples(ps_bvu)))

if (AGGREGATE_GENUS) {
        # NArm = TRUE discards OTUs without genus assignment (NA and Incertae_sedis).
        # Avoids collapsing ecologically distinct Incertae_sedis OTUs into one node.
        # Network represents only the genera with valid UNITE genus-level classification.
        ps_bvu <- suppressWarnings(tax_glom(ps_bvu, taxrank = "Genus", NArm = TRUE))
        cat(sprintf("  After genus aggregation (NArm=TRUE): %d genera x %d samples\n",
                    ntaxa(ps_bvu), nsamples(ps_bvu)))
}

# ── 4. Per-farm network construction ─────────────────────────────────────────

cat("\n── 4. Per-farm networks (SpiecEasi-MB via NetComi) ─────────────────────\n")

results <- vector("list", length(FARM_LEVELS))
names(results) <- FARM_LEVELS

for (farm in FARM_LEVELS) {
        cat(sprintf("\n  ── %s ──\n", toupper(farm)))

        ps_farm <- prune_samples(sample_data(ps_bvu)$farm_id == farm, ps_bvu)

        ps_bag  <- subset_samples(ps_farm, sample_type == "bagged_flower")
        ps_unb  <- subset_samples(ps_farm, sample_type == "unbagged_flower")

        n_bag <- nsamples(ps_bag)
        n_unb <- nsamples(ps_unb)
        cat(sprintf("  bagged=%d  unbagged=%d\n", n_bag, n_unb))

        if (n_bag < 5 || n_unb < 5) {
                message("  Skipping ", farm, ": too few samples")
                results[[farm]] <- NULL
                next
        }

        n_min         <- min(n_bag, n_unb)
        n_samp_thresh <- ceiling(PREV_THRESH_NET * n_min)

        net_raw <- tryCatch(
                netConstruct(
                        data        = ps_bag,
                        data2       = ps_unb,
                        filtTax     = "numbSamp",
                        filtTaxPar  = list(numbSamp = n_samp_thresh),
                        measure     = "spieceasi",
                        measurePar  = list(method           = "mb",
                                           nlambda          = SPIEC_NLAMBDA,
                                           lambda.min.ratio = SPIEC_LAMBDA_RATIO,
                                           pulsar.params    = list(thresh = STARS_THRESH)),
                        normMethod  = "none",
                        zeroMethod  = "none",
                        sparsMethod = "none",
                        verbose     = 0
                ),
                error = function(e) {
                        message("  netConstruct failed: ", conditionMessage(e))
                        NULL
                }
        )

        if (is.null(net_raw)) {
                results[[farm]] <- NULL
                next
        }

        n_taxa_net <- ncol(net_raw$adjaMat1)
        cat(sprintf("  taxa in network after filter: %d\n", n_taxa_net))

        if (n_taxa_net < 3) {
                message("  Skipping ", farm, ": <3 taxa after prevalence filter")
                results[[farm]] <- NULL
                next
        }

        net_ana <- tryCatch(
                netAnalyze(net_raw,
                           clustMethod = "cluster_fast_greedy",
                           hubPar      = "eigenvector",
                           verbose     = FALSE),
                error = function(e) {
                        message("  netAnalyze failed: ", conditionMessage(e))
                        NULL
                }
        )

        n_edges1 <- sum(edge_mask(net_raw$adjaMat1)) / 2
        n_edges2 <- sum(edge_mask(net_raw$adjaMat2)) / 2
        net_cmp <- if (!is.null(net_ana) && n_edges1 > 0 && n_edges2 > 0) {
                cat(sprintf("  netCompare: %d permutations ...\n", N_PERM))
                flush(stdout())
                t0 <- Sys.time()
                cmp <- tryCatch(
                        # cores = 1: parallel workers measured ~2x slower (script 22)
                        netCompare(net_ana, permTest = TRUE, nPerm = N_PERM,
                                   cores = 1, verbose = FALSE, seed = 42),
                        error = function(e) {
                                message("  netCompare failed: ", conditionMessage(e))
                                NULL
                        }
                )
                cat(sprintf("  netCompare: %s after %.1f min\n",
                            if (is.null(cmp)) "FAILED" else "done",
                            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
                flush(stdout())
                cmp
        } else {
                message(sprintf("  Skipping netCompare: edges bagged=%d unbagged=%d",
                                as.integer(n_edges1), as.integer(n_edges2)))
                NULL
        }

        results[[farm]] <- list(raw = net_raw, analyzed = net_ana, compared = net_cmp)
        cat(sprintf("  hub taxa (bagged): %d  |  hub taxa (unbagged): %d\n",
                    length(net_ana$hubs$hubs1 %||% character(0)),
                    length(net_ana$hubs$hubs2 %||% character(0))))
}

# ── 5. Save intermediate RDS ──────────────────────────────────────────────────

cat("\n── 5. Saving intermediate RDS ───────────────────────────────────────────\n")
dir.create(here("results", "rds"),     showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

saveRDS(results, here("results", "rds", "19b_network_results.rds"))
cat("  Saved: 19b_network_results.rds\n")

# ── 7. Topology summary table ─────────────────────────────────────────────────

cat("\n── 7. Building topology summary table ───────────────────────────────────\n")

# To rebuild the tables from an existing RDS without redoing netConstruct or
# netCompare: source Section 0 (parameters, libraries, helpers), then source from
# here. A no-op when the full script ran in the same session.
if (!exists("results")) {
        results <- readRDS(here("results", "rds", "19b_network_results.rds"))
        cat("  Loaded: 19b_network_results.rds\n")
}

topology_rows <- lapply(FARM_LEVELS, function(farm) {
        if (is.null(results[[farm]]) || is.null(results[[farm]]$analyzed)) return(NULL)
        tryCatch(extract_topology(results[[farm]]$raw, results[[farm]]$analyzed, farm),
                 error = function(e) { message("  topology extract failed for ", farm, ": ", conditionMessage(e)); NULL })
})

if (any(!sapply(topology_rows, is.null))) {
        dplyr::bind_rows(topology_rows) |>
                dplyr::arrange(farm, treatment) |>
                readr::write_csv(here("results", "tables", "19b_network_topology.csv"))
        cat("  Saved: 19b_network_topology.csv\n")
} else {
        cat("  No topology data to save.\n")
}

# ── 8. Hub taxa table ─────────────────────────────────────────────────────────

cat("\n── 8. Building hub taxa table ───────────────────────────────────────────\n")

hub_rows <- lapply(FARM_LEVELS, function(farm) {
        if (is.null(results[[farm]]) || is.null(results[[farm]]$analyzed)) return(NULL)
        tryCatch(extract_hubs(results[[farm]]$analyzed, farm),
                 error = function(e) { message("  hub extract failed for ", farm); NULL })
})

hub_combined <- dplyr::bind_rows(hub_rows)
if (nrow(hub_combined) > 0 && "treatment" %in% names(hub_combined)) {
        hub_combined |>
                dplyr::arrange(farm, treatment, dplyr::desc(eigenvector)) |>
                readr::write_csv(here("results", "tables", "19b_network_hubs.csv"))
        cat("  Saved: 19b_network_hubs.csv\n")
} else {
        cat("  No hub taxa found across all farms.\n")
}

# ── 9. Network comparison table (global properties + permutation p-values) ────

cat("\n── 9. Building network comparison table ─────────────────────────────────\n")

cmp_rows <- lapply(FARM_LEVELS, function(farm) {
        r <- results[[farm]]
        if (is.null(r) || is.null(r$analyzed) || is.null(r$compared)) return(NULL)
        tryCatch(
                extract_comparison(r$compared, r$analyzed, r$raw, farm),
                error = function(e) {
                        message("  comparison extract failed for ", farm, ": ", conditionMessage(e))
                        NULL
                }
        )
})

cmp_combined <- dplyr::bind_rows(cmp_rows)
if (nrow(cmp_combined) > 0) {
        cmp_combined |>
                dplyr::arrange(farm, scope, property) |>
                readr::write_csv(here("results", "tables", "19b_network_comparison.csv"))
        cat("  Saved: 19b_network_comparison.csv\n")
} else {
        cat("  No comparison data (netCompare skipped for all farms).\n")
}

n_success <- sum(!sapply(results, is.null))
cat(sprintf("\n19b done. %d / %d farms analysed successfully.\n",
            n_success, length(FARM_LEVELS)))
