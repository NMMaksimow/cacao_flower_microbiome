## ============================================================================
## utils_tax_label — deepest informative taxonomic rank, as a single label
##
## Sourced by the ANCOM-BC2 scripts (16a-18b).
##
## Returns one rank with its prefix kept: "g__Cladosporium", "f__Pleosporaceae",
## "d__Bacteria". No "(fam.)" / "(ord.)" suffix -- the g__/f__/o__ prefix already
## states the rank, so the suffix only repeated it.
##
## Keep in sync with deepest_rank() in notebooks/03_16S_ITS1_RFC.py: the RFC and
## ANCOM-BC2 figures must name the same taxon the same way, otherwise an
## importance bar chart cannot be read next to a volcano plot.
## ============================================================================

# Values that carry no information even when the rank field is populated.
# Counts in this project: g__uncultured 3331, f__uncultured 1100, o__uncultured
# 270, c__uncultured 93, f__Unknown_Family 84.
TAX_PLACEHOLDERS <- c("uncultured", "uncultured_bacterium", "metagenome",
                      "unidentified", "Unknown_Family",
                      "Incertae_sedis", "Incertae_Sedis")

# An unplaced rank is named after the taxon enclosing it, and the two references
# spell that differently:
#   UNITE  "<taxon>_<rank>_Incertae_sedis"  e.g. g__Fungi_gen_Incertae_sedis
#   SILVA  "<taxon>_Incertae_Sedis"         e.g. f__Rhizobiales_Incertae_Sedis
#                                           (no rank token, capital S)
# Both reduce to the enclosing taxon, which the repeat check below then drops.
INCERTAE_RE <- "_((gen|fam|ord|cls|phy)_)?Incertae_[Ss]edis$"

# Ranks from the top down. Two top ranks, because the references disagree on the
# word: SILVA calls it Domain (d__Bacteria), UNITE calls it Kingdom (k__Fungi).
# No taxonomy table carries both, so their order relative to each other never
# matters. Including them is what stops 912 16S ASVs and 350 ITS1 OTUs -- those
# classified no further than the top rank -- from falling back to a bare hex id.
RANKS_TOP_DOWN <- c("Domain", "Kingdom", "Phylum", "Class", "Order",
                    "Family", "Genus")

make_tax_label <- function(physeq, feature_ids) {
        tt    <- as.data.frame(tax_table(physeq))
        ranks <- intersect(RANKS_TOP_DOWN, colnames(tt))

        vapply(feature_ids, function(id) {
                if (!id %in% rownames(tt)) return(substr(id, 1, 12))

                vals  <- as.character(unlist(tt[id, ranks], use.names = FALSE))
                label <- NA_character_
                seen  <- character(0)   # bare names already claimed by higher ranks

                for (v in vals) {
                        # NA = rank not assigned; nchar <= 3 = bare prefix such as "g__"
                        if (is.na(v) || nchar(v) <= 3) next

                        # Split the rank prefix off BEFORE stripping the suffix. Doing it
                        # the other way round eats the prefix's own underscore when the
                        # placeholder carries no enclosing taxon: g__Incertae_Sedis would
                        # become "g_" instead of being recognised as a placeholder.
                        pref <- regmatches(v, regexpr("^[a-z]__", v))
                        bare <- sub("^[a-z]__", "", v) |> sub(pattern = INCERTAE_RE,
                                                              replacement = "")
                        v    <- paste0(if (length(pref)) pref else "", bare)

                        # A value repeating a higher rank marks placement inside that
                        # taxon, not a name at this level. Without this check 177 ITS1
                        # OTUs acquire a genus called Fungi, and f__Rhizobiales_Incertae_Sedis
                        # becomes a family named after its own order.
                        if (bare %in% TAX_PLACEHOLDERS || bare %in% seen) next

                        seen  <- c(seen, bare)
                        label <- v      # walking top-down, so the last survivor is deepest
                }

                if (is.na(label)) substr(id, 1, 12) else label
        }, character(1))
}
