# Generating feature-tables for gamma diversity analysis from publicly available
# microbiome datasets.
# 1) Human Microbiome Project (HMP)
# 2)

# curatedMetagenomicData: 93 datasets
# MetaPhlAn3 generated output (profiler) 
# - relative_abundance (species level taxonomy)
# - marker_presence
# - marker_abundance

# HUMAnN3 generated output (UniRef90 DB):
# - gene_families
# - pathway_coverage
# - pathway_abundance



# ---- Install packages ----
install.packages("curatedMetagenomicData")
install.packages("DT")

# ---- Activate packages ----
library(curatedMetagenomicData)
library(dplyr)
library(DT)

# ---- Explore available studies ----
# sampleMetadata — это tibble, не функция, без скобок
all_studies <- sampleMetadata |>
        distinct(study_name, body_site, disease, number_reads) |>
        arrange(study_name)

# Просмотр в интерактивной таблице
DT::datatable(all_studies)

# Посмотреть все доступные названия исследований
unique(sampleMetadata$study_name)

# Фильтр по конкретному исследованию
sampleMetadata |>
        filter(study_name == "ShiB_2015")

# ---- Получить данные исследования ----
# Сначала проверь что study существует, потом запрашивай данные:
# dryrun = TRUE показывает что будет загружено без скачивания
curatedMetagenomicData("ShiB_2015.relative_abundance", dryrun = TRUE)

# Загрузка через returnSamples() (современный API, возвращает TreeSummarizedExperiment)
tse <- sampleMetadata |>
        filter(study_name == "ShiB_2015") |>
        returnSamples("relative_abundance")


# ---- Sample metadata ----
# sampleMetadata is a huge tibble with metadata on different studies
# ShiB_2015 study is on oral cavity microbiome in healthy and periodontites (48 observations)
curatedMetagenomicData::sampleMetadata |> 
        dplyr::filter(study_name == "ShiB_2015") |> 
        str()

str(tse)


library(SummarizedExperiment)

# Транспонируем: baseline_gamma ожидает образцы × виды
ra_matrix <- t(assay(tse, "relative_abundance"))   # (n_samples, n_species)

# Метаданные образцов
meta <- as.data.frame(colData(tse))

dim(ra_matrix)   # посмотри сколько образцов и видов
head(meta)       # что за образцы

