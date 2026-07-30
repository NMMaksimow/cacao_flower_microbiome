# ---- Packages ----
library(here)
library(phyloseq)


# ---- Load intermediate objects (saved by 01b) ----
# ps_ITS1_raw      <- readRDS(here("results", "rds", "ps_ITS1_raw.rds"))
# ps_ITS1_decontam <- readRDS(here("results", "rds", "ps_ITS1_decontam.rds"))
# ps_ITS1_fungi    <- readRDS(here("results", "rds", "ps_ITS1_fungi.rds"))
# decontam_prev05_ITS1  <- readRDS(here("results", "rds", "decontam_prev05_ITS1.rds"))
ps_ITS1_fungi_biosamples    <- readRDS(here("results", "rds", "ps_ITS1_fungi_biosamples.rds"))

# Metadata as plain data.frame (needed for joins by sample_id)
# metadata <- sample_data(ps_ITS1_raw) |>
#         as.data.frame() |>
#         (\(x) { class(x) <- "data.frame"; x })()

# ---- Output directories ----
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

################################################################################
###################### Diagnostic rarefaction curves ###########################
################################################################################

# Create a matrix for building rarefaction curves and transpose it because vegan::rarecurve() requires samples in rows
asv_mat_ITS1_all <- ps_ITS1_fungi_biosamples |>
        otu_table() |>
        as("matrix") |>
        t()

# Keep only samples with reads after off-target removal (non-zero depth)
nonzero_mask  <- rowSums(asv_mat_ITS1_all) > 0
asv_mat_ITS1  <- asv_mat_ITS1_all[nonzero_mask, , drop = FALSE]

zero_samples  <- rownames(asv_mat_ITS1_all)[!nonzero_mask]
if (length(zero_samples) > 0)
        message("Samples with 0 reads after off-target removal (excluded from plot): ",
                paste(zero_samples, collapse = ", "))

min_depth <- min(rowSums(asv_mat_ITS1))

# Color coding using categorical variables from metadata
color_var <- "sample_type"
metadata_for_rarecurve <- as.data.frame(sample_data(ps_ITS1_fungi_biosamples))
metadata_for_rarecurve <- metadata_for_rarecurve[rownames(asv_mat_ITS1), , drop = FALSE]

pal <- c(bagged_flower = "#4C72B0", unbagged_flower = "#DD8452")
curve_cols <- pal[as.character(metadata_for_rarecurve[[color_var]])]

# Rarefaction curves — base R graphics, saved with png()/dev.off()
png(here("results", "figures", "05b_ITS1_rarefaction_curves.png"),
    width = 10, height = 7, units = "in", res = 300)

vegan::rarecurve(
        asv_mat_ITS1,
        step   = 500,
        col    = curve_cols,
        label  = FALSE,
        xlab   = "Sequencing depth",
        ylab   = "Observed ASVs",
        main   = paste("Rarefaction curves coloured by", color_var)
)

legend(
        "bottomright",
        legend = names(pal),
        col    = pal,
        lwd    = 2,
        bty    = "n"
)
dev.off()

# ---- Good's coverage ----
# C = 1 - F1/N; F1 = per-sample singleton taxa count, N = total reads per sample.
# Values near 1 indicate the library is well-saturated at this sequencing depth.
goods_cov <- 1 - rowSums(asv_mat_ITS1 == 1) / rowSums(asv_mat_ITS1)

message(sprintf("Good's coverage ITS1: median = %.3f  min = %.3f  max = %.3f",
                median(goods_cov), min(goods_cov), max(goods_cov)))

goods_df <- data.frame(
        sample_id      = names(goods_cov),
        goods_coverage = goods_cov,
        # use metadata_for_rarecurve already aligned to asv_mat_ITS1 rownames
        sample_type    = metadata_for_rarecurve[names(goods_cov), "sample_type"]
)
write.csv(goods_df,
          here("results", "tables", "05b_ITS1_goods_coverage.csv"),
          row.names = FALSE)

message("05b done")