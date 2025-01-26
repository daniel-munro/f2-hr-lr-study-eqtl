library(tidyverse)

load_tensorqtl <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----c---dddd-ddd") |>
        rename(gene_id = phenotype_id) |>
        separate_wider_delim(variant_id, ":", names = c("chrom", "pos"),
                 cols_remove = FALSE) |>
        mutate(chrom = str_replace(chrom, "chr", ""),
               pos = as.integer(pos))
}

load_tensorqtl_ind <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "ci----c---dddd-cd") |>
        rename(gene_id = phenotype_id) |>
        separate_wider_delim(variant_id, ":", names = c("chrom", "pos"),
                 cols_remove = FALSE) |>
        mutate(chrom = str_replace(chrom, "chr", ""),
               pos = as.integer(pos))
}

load_afc <- function(afc_out) {
    read_tsv(afc_out, col_types = "cc--d--") |>
        rename(gene_id = pid,
               variant_id = sid) |>
        mutate(log2_aFC = signif(log2_aFC, 6))
}

genes <- read_tsv("data/reference/Rattus_norvegicus.Rnor_6.0.104.gtf.gz",
                  col_types = "--cii-c-c",
                  col_names = c("type", "start", "end", "strand", "etc"),
                  comment = "#") |>
    filter(type == "gene") |>
    mutate(gene_id = str_match(etc, 'gene_id "([^"]+)"')[, 2],
           gene_name = str_match(etc, 'gene_name "([^"]+)"')[, 2],
           tss = if_else(strand == "+", start, end)) |>
    select(gene_id, gene_name, strand, tss)
stopifnot(sum(duplicated(genes$gene_id)) == 0)

alleles <- read_tsv("data/reference/alleles.txt.gz", col_types = "ccc",
                    col_names = c("variant_id", "ref", "alt"))

alleles_bLR_bHR <- read_tsv("data/genotypes/bLR_bHR_alleles.tsv.gz", col_types = "cccddc")

inner_join(alleles, alleles_bLR_bHR, by = "variant_id", relationship = "one-to-one") |>
    with(table(ref.x == ref.y, alt.x == alt.y))

afc <- load_afc("data/eqtl/f2.aFC.txt")

top_assoc <- load_tensorqtl("data/eqtl/f2.cis_qtl.txt.gz") |>
    left_join(afc, by = c("gene_id", "variant_id"), relationship = "one-to-one") |>
    left_join(genes, by = "gene_id", relationship = "one-to-one") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
    relocate(gene_name, .after = gene_id) |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos)

write_tsv(top_assoc, "data/eqtl/f2.top_assoc.txt")

eqtls_ind <- load_tensorqtl_ind("data/eqtl/f2.cis_independent_qtl.txt.gz") |>
    left_join(afc, by = c("gene_id", "variant_id"), relationship = "one-to-one") |>
    left_join(genes, by = "gene_id", relationship = "many-to-one") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
    relocate(gene_name, .after = gene_id) |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos) |>
    left_join(select(alleles_bLR_bHR, variant_id, alt_group), by = "variant_id",
              relationship = "many-to-one") |>
    mutate(
        log2_aFC_bLR_bHR = case_when(
            alt_group == "bLR" ~ log2_aFC,
            alt_group == "bHR" ~ -log2_aFC,
            .default = NA
        ),
        upregulated_in = case_when(
            log2_aFC > 0 & alt_group == "bLR" ~ "bLR",
            log2_aFC > 0 & alt_group == "bHR" ~ "bHR",
            log2_aFC < 0 & alt_group == "bLR" ~ "bHR",
            log2_aFC < 0 & alt_group == "bHR" ~ "bLR",
            .default = "equal"
        )
    ) |>
    select(-alt_group) |>
    rename(log2_aFC_alt_ref = log2_aFC)

write_tsv(eqtls_ind, "data/eqtl/f2.eqtls_indep.txt")

