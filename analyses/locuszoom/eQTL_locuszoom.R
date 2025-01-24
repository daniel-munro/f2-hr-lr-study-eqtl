# Locuszoom plot for AABR07071904.1 (ENSRNOG00000052237) eQTL

# To get f2.all_cis.ENSRNOG00000052237.txt, I ran:
# python scripts/tensorqtl_all_pairs_for_gene.py \
#     data/tensorqtl/nominal/f2.cis_qtl_pairs.1.parquet \
#     ENSRNOG00000052237 \
#     data/tensorqtl/f2.all_cis.ENSRNOG00000052237.txt

# To get VCF of cis-window SNPs, I ran:
# plink2 --bfile data/tensorqtl/geno \
#     --chr chr1 \
#     --from-bp 93836295 \
#     --to-bp 95836297 \
#     --recode vcf id-paste=iid \
#     --out ENSRNOG00000052237
# bgzip ENSRNOG00000052237.vcf
# tabix ENSRNOG00000052237.vcf.gz

library(tidyverse)

load_geno <- function(chrom, start, end, samples) {
    filename <- str_glue("analysis/locuszoom/ENSRNOG00000052237.vcf.gz")
    rng <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end))
    gt <- VariantAnnotation::readGT(filename, param = rng)
    geno <- apply(gt, 2, function(x) c("0/0" = 0, "0/1" = 1, "1/1" = 2)[x])
    rownames(geno) <- rownames(gt)
    t(geno)[samples, ]
}

LD_with_top <- function(variant_id, chrom, pos, top_id, samples) {
    geno <- load_geno(chrom[1], min(pos) - 1, max(pos) + 1, samples)
    map_dbl(variant_id, ~ if (sd(geno[, .x]) == 0) {0} else {cor(geno[, .x], geno[, top_id]) ^ 2})
}

gene <- "ENSRNOG00000052237"

# tss <- read_tsv("data/expression/expression.iqn.bed",
#                       col_types = cols(gene_id = "c", end = "i", .default = "-")) |>
#     filter(gene_id == gene) |>
#     pull(end)
tss <- 94836296

# top_snp <- read_tsv("data/eqtl/f2.eqtls_indep.txt", col_types = "cc-cii----d-------") |>
#     filter(gene_id == gene)
top_snp <- "chr1:95003407"

# log10_threshold <- read_tsv("data/eqtl/f2.top_assoc.txt", col_types = "cc-cii---------d-") |>
#     filter(gene_id == gene)
log10_threshold <- -log10(0.000472019)

plekhf1_snp <- "chr1:94610123"

# Get samples to subset genotypes when calculating LD:
samples <- read_lines("data/samples.txt")

# Load all cis-window p-values for this gene:
pvals <- read_tsv(str_glue("analysis/locuszoom/f2.all_cis.{gene}.txt"),
                  col_types = "-ci---d--") |>
    separate(variant_id, c("chrom", "pos"), sep=":", convert = TRUE,
             remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", "") |> as.integer()) |>
    mutate(top = pval_nominal == min(pval_nominal),
           LD = LD_with_top(variant_id, chrom, pos,
                            variant_id[top][ceiling(sum(top) / 2)],
                            samples)) |>
    mutate(log10p = -log10(pval_nominal),
           pos = pos / 1e6)

pvals |>
    ggplot(aes(x = pos, y = log10p, color = LD)) +
    geom_hline(yintercept = log10_threshold, lty = "12", linewidth = 0.5, color = "#555555") +
    geom_vline(xintercept = tss / 1e6, color = "#888888") +
    annotate("text", x = tss / 1e6, y = 10, label = "'TSS '*(italic('AABR07071904.1'))", parse = TRUE, hjust = -0.03) +
    geom_point(size = 1) +
    geom_point(data = filter(pvals, variant_id == top_snp), color = "purple", shape = 18, size = 3) +
    # annotate("text", x = pvals$pos[pvals$variant_id == top_snp], y = pvals$pos[pvals$variant_id == top_snp]) +
    geom_text(data = filter(pvals, variant_id == top_snp), label = top_snp, hjust = 0.2, vjust = -0.8, color = "black") +
    geom_point(data = filter(pvals, variant_id == plekhf1_snp), color = "black", shape = 15, size = 2) +
    geom_text(data = filter(pvals, variant_id == plekhf1_snp),
              label = "italic('Plekhf1')*'missense variant'", parse = TRUE, color = "black", hjust = 1, vjust = -0.5) +
    scale_color_gradientn(colors = c("#000066", "#000066", "#377eb8", "#377eb8", "#208020", "#208020", "#ff7f00", "#ff7f00", "#e41a1c", "#e41a1c"),
                          values = c(0, 0.1999, 0.2, 0.3999, 0.4, 0.5999, 0.6, 0.7999, 0.8, 1),
                          breaks = c(0.2, 0.4, 0.6, 0.8)) +
    expand_limits(y = 55) +
    theme_classic() +
    theme(
        strip.text = element_text(color = "black"),
        legend.position = c(0.9, 0.7),
    ) +
    xlab("Position on chr1 (Mb)") +
    ylab(expression(-log[10]*"(p-value)")) +
    labs(color = expression("LD "*(r^2)))

ggsave("analysis/locuszoom/AABR07071904.1_locuszoom.png", width = 5, height = 3, device = png)
