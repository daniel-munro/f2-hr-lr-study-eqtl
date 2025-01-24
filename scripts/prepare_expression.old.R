library(tidyverse)

genes <- rtracklayer::import("data/reference/Rattus_norvegicus.Rnor_6.0.104.gtf.gz") |>
    as_tibble() |>
    filter(type == "gene") |>
    mutate(TSS = if_else(strand == "-", end, start)) |>
    select(seqnames, TSS, gene_id)

expr <- read_csv("data/expression/ResidualsMatrix_wMean.csv",
                 col_types = cols(`...1` = "c", .default = "d")) |>
    rename(gene_id = `...1`)

samples <- readxl::read_xlsx("data/Akil_bHRbLR_F0_F2_IDs.xlsx", range = "A1:E850",
                             col_types = "text") |>
    filter(!is.na(LibraryID)) |>
    select(LibraryID, RatNumForNIHstudy) |>
    deframe()

bed <- genes |>
    mutate(start = TSS - 1) |>
    select(`#chr` = seqnames,
           start,
           end = TSS,
           gene_id) |>
    inner_join(expr, by = "gene_id", relationship = "one-to-one") |>
    rename_with(~ samples[.x], starts_with("SL"))
