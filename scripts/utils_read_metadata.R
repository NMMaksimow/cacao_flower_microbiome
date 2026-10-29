## ============================================================================
## utils_read_metadata — sample metadata from the TSV, with column types
##
## Sourced by 01a, 01b and 01_microscopy_data_analysis.Rmd.
##
## data/sample_metadata_with_microscopy.txt (exported from the lab workbook) is
## the single source of truth for sample metadata. Columns are typed by name, so
## adding, removing or reordering columns in the TSV does not shift the types.
## Columns not listed in col_types are read as text.
##
## Returns a tibble with one row per sample (336: 294 flowers + 42 controls) and
## `sample-id` renamed to sample_id. Factor conversion for phyloseq stays in the
## calling scripts.
## ============================================================================

read_sample_metadata <- function() {
        readr::read_tsv(
                here::here("data", "sample_metadata_with_microscopy.txt"),
                # The second line holds the QIIME 2 column types ("#q2:types ...");
                # everything after this string is ignored, so the line is skipped
                comment = "#q2:types",
                na = c("", "NA"),
                col_types = readr::cols(
                        .default = readr::col_character(),
                        # Batch dates
                        collection_batch = readr::col_date(),
                        extraction_batch = readr::col_date(),
                        pcr_batch_16s    = readr::col_date(),
                        qubit_batch_16s  = readr::col_date(),
                        pcr_batch_its1   = readr::col_date(),
                        qubit_batch_its1 = readr::col_date(),
                        # Failed PCR replicates repeated on another day and pooled
                        pcr_redo_16s  = readr::col_integer(),
                        pcr_redo_its1 = readr::col_integer(),
                        # Measurements
                        dna_concentration_hs_16s  = readr::col_double(),
                        dna_concentration_hs_its1 = readr::col_double(),
                        pi_dic                    = readr::col_double(),
                        pi_clsm                   = readr::col_double(),
                        callose_plug_index        = readr::col_double(),
                        fungal_colonization_score = readr::col_double(),
                        # 0/1 flags
                        is_specimen         = readr::col_logical(),
                        specimen_qc_pass    = readr::col_logical(),
                        is_pollination_dic  = readr::col_logical(),
                        is_pollination_clsm = readr::col_logical(),
                        is_germination      = readr::col_logical(),
                        is_hyphae_present   = readr::col_logical(),
                        is_spore_present    = readr::col_logical(),
                        is_insect_present   = readr::col_logical()
                )
        ) |>
                # Hyphen in sample-id is a QIIME 2 metadata requirement
                dplyr::rename(sample_id = `sample-id`) |>
                dplyr::mutate(
                        # "0.1%" -> 0.001 (fraction, as in the v1 objects)
                        dye_concentration = readr::parse_number(dye_concentration) / 100,
                        # Deviations from the sampling protocol; "none" first as reference
                        sampling_deviation = factor(
                                sampling_deviation,
                                levels = c("none", "no_buffer", "two_flowers",
                                           "insect_present")
                        )
                )
}
