## ============================================================================
## 17 — Pollen intensity by farm and flower treatment
##
## All measured flowers shown as coloured dots (by management type).
## Shape encodes insect / spore detection.
## Threshold line at pi_clsm = PI_THRESH for reference.
## Two horizontal panels: Unbagged (top) / Bagged (bottom).
## Farms ordered along the forest-to-sun gradient on X axis.
##
## Input:
##   data/sample_metadata_with_microscopy.txt
## Output:
##   results/figures/17_high_pollination_by_farm.png
## ============================================================================

library(here)
library(dplyr)
library(ggplot2)

# ── Parameters ────────────────────────────────────────────────────────────────
PI_THRESH <- 40

FARM_ORDER <- c("ib", "vr", "sa", "kk", "mt", "vi", "yb")

MGMT_ORDER <- c("inside_forest", "near_forest", "agroforest", "full_sun")
MGMT_LABELS <- c(
    inside_forest = "Inside forest",
    near_forest   = "Near forest",
    agroforest    = "Agroforest",
    full_sun      = "Full sun"
)
MGMT_COLOURS <- c(
    "Inside forest" = "#1B7837",
    "Near forest"   = "#7FBC41",
    "Agroforest"    = "#E08214",
    "Full sun"      = "#B2182B"
)

MARKER_SHAPES <- c(
    "none"           = 16,   # filled circle
    "spore"          = 17,   # filled triangle
    "insect"         = 4,    # cross / ×
    "insect + spore" = 8     # asterisk / *
)

# ── Load & prepare data ───────────────────────────────────────────────────────
meta_raw <- read.delim(
    here("data", "sample_metadata_with_microscopy.txt"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

meta <- meta_raw |>
    filter(!grepl("^#", sample.id, fixed = FALSE)) |>
    filter(sample_type %in% c("bagged_flower", "unbagged_flower")) |>
    mutate(
        pi_clsm           = as.numeric(gsub(",", ".", pi_clsm)),
        is_insect_present = suppressWarnings(as.integer(is_insect_present)),
        is_spore_present  = suppressWarnings(as.integer(is_spore_present)),
        farm_id         = factor(farm_id, levels = FARM_ORDER),
        management_type = factor(management_type, levels = MGMT_ORDER,
                                 labels = MGMT_LABELS),
        sample_type     = factor(sample_type,
                                 levels = c("unbagged_flower", "bagged_flower"),
                                 labels = c("Unbagged", "Bagged")),
        bio_marker = case_when(
            !is.na(is_insect_present) & is_insect_present == 1 &
            !is.na(is_spore_present)  & is_spore_present  == 1 ~ "insect + spore",
            !is.na(is_insect_present) & is_insect_present == 1 ~ "insect",
            !is.na(is_spore_present)  & is_spore_present  == 1 ~ "spore",
            TRUE                                                ~ "none"
        ),
        bio_marker = factor(bio_marker,
                            levels = c("none", "spore", "insect", "insect + spore"))
    ) |>
    filter(!is.na(pi_clsm))

# n per farm × treatment for top labels
df_n <- meta |>
    count(farm_id, sample_type, name = "n") |>
    mutate(y_label = max(meta$pi_clsm, na.rm = TRUE) * 1.07)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(meta, aes(x = farm_id, y = pi_clsm,
                      colour = management_type,
                      shape  = bio_marker)) +
    geom_jitter(size = 2.8, alpha = 0.85, stroke = 1.2,
                width = 0.18, height = 0, seed = 42) +
    geom_hline(yintercept = PI_THRESH, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    geom_text(data = df_n,
              aes(x = farm_id, y = y_label, label = paste0("n=", n)),
              inherit.aes = FALSE,
              size = 2.6, colour = "grey40") +
    facet_grid(sample_type ~ ., switch = "y") +
    scale_colour_manual(values = MGMT_COLOURS, name = "Management type") +
    scale_shape_manual(values = MARKER_SHAPES, name = "Biological marker") +
    scale_y_continuous(
        name   = "Pollen intensity (pi_clsm, grains/stigma)",
        expand = expansion(mult = c(0.02, 0.08))
    ) +
    scale_x_discrete(name = "Farm") +
    labs(
        title    = paste0("Pollen intensity by farm  |  threshold ", PI_THRESH, " grains (dashed)"),
        subtitle = "All measured flowers  |  colour = management type  |  shape = insect / spore detection"
    ) +
    guides(
        colour = guide_legend(override.aes = list(shape = 16, size = 3)),
        shape  = guide_legend(override.aes = list(size = 3.5))
    ) +
    theme_bw(base_size = 11) +
    theme(
        strip.background   = element_rect(fill = "grey92", colour = NA),
        strip.text         = element_text(face = "bold", size = 11),
        strip.text.y.left  = element_text(angle = 0),
        strip.placement    = "outside",
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "bottom",
        plot.title         = element_text(face = "bold"),
        plot.subtitle      = element_text(colour = "grey40", size = 9)
    )

# ── Save ──────────────────────────────────────────────────────────────────────
out_path <- here("results", "figures", "17_high_pollination_by_farm.png")
ggsave(out_path, p, width = 11, height = 6, dpi = 200, bg = "white")
message("Saved: ", out_path)

print(p)
