library(tidyverse)

T_SMR <- function(z_eQTL, z_GWAS) {
    (z_eQTL^2 * z_GWAS^2) / (z_eQTL^2 + z_GWAS^2)
}

p_SMR <- function(T_SMR) {
    pchisq(T_SMR, df = 1, lower.tail = FALSE)
}

eqtls <- read_tsv(
    "data/eqtl/f2.eqtls_indep.txt",
    col_types = cols(gene_id = "c", gene_name = "c", variant_id = "c",
                     ref = "c", alt = "c", slope = "d", slope_se = "d",
                     .default = "-")
) |>
    mutate(z_eQTL = slope / slope_se) |>
    select(gene_id, gene_name, variant_id, ref, alt, z_eQTL)

traits <- read_csv("data/gwas_pub/phenotype_data/data_dictionary.csv", col_types = "ccc")

gwas <- tibble(trait = traits$measure) |>
    reframe(
        read_tsv(
            str_glue("data/gwas_pub/gwas_results/u01_huda_akil_{trait}.loco.mlma.gz"),
            col_types = "-c-cc-dd-"
        ) |>
            filter(SNP %in% eqtls$variant_id),
        .by = trait
    ) |>
    mutate(z_GWAS = b / se) |>
    select(trait, variant_id = SNP, A1, A2, z_GWAS)

df <- eqtls |>
    full_join(gwas, by = "variant_id", relationship = "many-to-many") |>
    mutate(z_GWAS = if_else(alt == A1, z_GWAS, -1 * z_GWAS),
           T_SMR = T_SMR(z_eQTL, z_GWAS),
           p_SMR = p_SMR(T_SMR)) |>
    select(-ref, -alt, -A1, -A2)

write_tsv(df, "data/coloc/colocs.tsv")

sig <- df |>
    filter(p_SMR < 0.05 / n(), .by = trait)

write_tsv(sig, "data/coloc/colocs_sig.tsv")
