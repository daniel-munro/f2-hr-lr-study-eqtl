## Locuszoom plot for AABR07071904.1 (ENSRNOG00000052237) eQTL

## Get eQTL sumstats for all cis-window SNPs
python scripts/tensorqtl_all_pairs_for_gene.py \
    data/tensorqtl/nominal/f2.cis_qtl_pairs.1.parquet \
    ENSRNOG00000052237 \
    data/analysis/f2.all_cis.ENSRNOG00000052237.txt

## Get VCF of cis-window SNPs
plink2 --bfile data/tensorqtl/geno \
    --chr chr1 \
    --from-bp 93836295 \
    --to-bp 95836297 \
    --recode vcf id-paste=iid \
    --out ENSRNOG00000052237
bgzip data/analysis/ENSRNOG00000052237.vcf
tabix data/analysis/ENSRNOG00000052237.vcf.gz
