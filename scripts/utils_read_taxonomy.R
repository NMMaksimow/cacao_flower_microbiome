## ============================================================================
## utils_read_taxonomy — QIIME 2 taxonomy.tsv as a phyloseq tax_table
##
## Sourced by 01a and 01b.
##
## Ranks are matched by their prefix, not by position. The references differ:
##   SILVA 138 (v1 16S)   d__ p__ c__ o__ f__ g__ s__
##   SILVA 144 (v2 16S)   d__ k__ p__ c__ o__ f__ g__ s__   (k__ added after d__)
##   UNITE (v1 ITS1)      k__ p__ c__ o__ f__ g__ s__ sh__  (species hypothesis)
##   UNITE (v2 ITS1)      k__ p__ c__ o__ f__ g__ s__
## Splitting by position would shift every SILVA 144 rank after Domain by one.
##
## Values keep their prefix ("g__Wolbachia"), as the downstream scripts expect.
## The table has only the ranks present in the file, top-down. "Unassigned" has
## no prefix and goes to the top rank of the file (Domain or Kingdom).
## ============================================================================

# Rank prefix -> column name, top-down
TAX_RANKS <- c(d = "Domain", k = "Kingdom", p = "Phylum", c = "Class",
               o = "Order", f = "Family", g = "Genus", s = "Species",
               sh = "UNITE_SH")

read_qiime2_taxonomy <- function(path) {
        tax_long <- readr::read_tsv(path, col_types = readr::cols(.default = "c")) |>
                dplyr::transmute(FeatureID = `Feature ID`, Taxon) |>
                # One row per rank: "d__Bacteria;k__Bacillati;..." -> rows
                tidyr::separate_longer_delim(Taxon, delim = ";") |>
                dplyr::mutate(Taxon = stringr::str_trim(Taxon)) |>
                # Empty fields (e.g. after a trailing ";") are missing ranks
                dplyr::filter(Taxon != "") |>
                dplyr::mutate(
                        rank = unname(TAX_RANKS[stringr::str_extract(Taxon, "^[a-z]+(?=__)")])
                )

        # Ranks present in this file, top-down; "Unassigned" goes to the first
        ranks <- TAX_RANKS[TAX_RANKS %in% tax_long$rank]
        tax_long$rank[is.na(tax_long$rank)] <- ranks[[1]]

        tax_long |>
                tidyr::pivot_wider(names_from = rank, values_from = Taxon) |>
                dplyr::select(FeatureID, dplyr::all_of(unname(ranks))) |>
                tibble::column_to_rownames("FeatureID") |>
                as.matrix() |>
                phyloseq::tax_table()
}
