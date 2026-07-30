# ---- Packages ----
library(here)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(scales)

# ---- Load intermediate objects (saved by 01a) ----
# ps_16S_raw      <- readRDS(here("results", "rds", "ps_16S_raw.rds"))
# ps_16S_decontam <- readRDS(here("results", "rds", "ps_16S_decontam.rds"))
# ps_16S_bacteria <- readRDS(here("results", "rds", "ps_16S_bacteria.rds"))
ps_16S_bacteria_biosamples <- readRDS(here("results", "rds", "ps_16S_bacteria_biosamples.rds"))
# decontam_prev055_16S  <- readRDS(here("results", "rds", "decontam_prev055_16S.rds"))

# ---- Output directories ----
dir.create(here("results", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "tables"),  showWarnings = FALSE, recursive = TRUE)

# There are two  approaches to tackle uneven sampling effort/sequencing depth:
# a) Rarefaction
# ?vegan::rarefy()
# vegan::rrarefy()
# vegan::drarefy() calculated a
# vegan::rarecurve()
# ?vegan::rareslope()

# I don't understand about this function whether it performs iterative subsampling, or do it once.
# ?phyloseq::rarefy_even_depth()
# b) Normalization
# ?vegan::decostand()

################################################################################
###################### Diagnostic rarefaction curves ###########################
################################################################################
# Create a matrix for building rarefaction curves and transpose it because vegan::rarecurve() requires samples in rows
asv_mat_16S_all <- ps_16S_bacteria_biosamples |>
        otu_table() |>
        as("matrix") |>
        t()

# Keep only samples with reads after off-target removal (non-zero depth)
nonzero_mask  <- rowSums(asv_mat_16S_all) > 0
asv_mat_16S   <- asv_mat_16S_all[nonzero_mask, , drop = FALSE]

zero_samples  <- rownames(asv_mat_16S_all)[!nonzero_mask]
if (length(zero_samples) > 0)
        message("Samples with 0 reads after off-target removal (excluded from plot): ",
                paste(zero_samples, collapse = ", "))

min_depth <- min(rowSums(asv_mat_16S))

# Color coding using categorical variables from metadata
color_var <- "sample_type"
metadata_for_rarecurve <- as.data.frame(sample_data(ps_16S_bacteria_biosamples))
metadata_for_rarecurve <- metadata_for_rarecurve[rownames(asv_mat_16S), , drop = FALSE]

pal <- c(bagged_flower = "#4C72B0", unbagged_flower = "#DD8452")
curve_cols <- pal[as.character(metadata_for_rarecurve[[color_var]])]

# Rarefaction curves — base R graphics, saved with png()/dev.off()
png(here("results", "figures", "05a_16S_rarefaction_curves.png"),
    width = 10, height = 7, units = "in", res = 300)
vegan::rarecurve(
        asv_mat_16S,
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
goods_cov <- 1 - rowSums(asv_mat_16S == 1) / rowSums(asv_mat_16S)

message(sprintf("Good's coverage 16S: median = %.3f  min = %.3f  max = %.3f",
                median(goods_cov), min(goods_cov), max(goods_cov)))

goods_df <- data.frame(
        sample_id      = names(goods_cov),
        goods_coverage = goods_cov,
        # use metadata_for_rarecurve already aligned to asv_mat_16S rownames
        sample_type    = metadata_for_rarecurve[names(goods_cov), "sample_type"]
)
write.csv(goods_df,
          here("results", "tables", "05a_16S_goods_coverage.csv"),
          row.names = FALSE)

# binwidth = 0.01 → one bin per percentage-point; adjust if distribution is too compressed
p_hist <- ggplot(goods_df, aes(x = goods_coverage, fill = sample_type)) +
        geom_histogram(binwidth = 0.01, colour = "white", linewidth = 0.2) +
        scale_fill_manual(values = pal, name = "Sample type") +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
        labs(x = "Good's coverage", y = "Number of samples") +
        theme_bw()

# angle = 20 to avoid label overlap for "bagged_flower" / "unbagged_flower"
p_box <- ggplot(goods_df, aes(x = sample_type, y = goods_coverage, fill = sample_type)) +
        geom_boxplot(show.legend = FALSE) +
        scale_fill_manual(values = pal) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        labs(x = "Sample type", y = "Good's coverage") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_goods <- p_hist | p_box

ggsave(here("results", "figures", "05a_16S_goods_coverage.png"),
       plot = p_goods, width = 10, height = 5, dpi = 300)

message("05a done")




