######################################
## Prepare genotypes and covariates ##
######################################

awk '{print $1, $1}' data/samples.txt > tmp.txt
plink2 --make-bed \
    --bfile data/genotype_data/u01_huda_akil_genotypes \
    --fa ~/ratgtex/ref_rn6/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    --ref-from-fa force \
    --keep tmp.txt \
    --maf 0.01 \
    --mac 2 \
    --max-alleles 2 \
    --out data/tensorqtl/geno
rm tmp.txt

## Prune genotypes to compute covariate PCs

plink2 \
    --bfile data/tensorqtl/geno \
    --geno 0.00 \
    --maf 0.05 \
    --indep-pairwise 200 100 0.1 \
    --out data/covar/geno
plink2 \
    --bfile data/tensorqtl/geno \
    --extract data/covar/geno.prune.in \
    --export vcf bgz id-paste=iid \
    --out data/covar/geno

Rscript ~/ratgtex/scripts/covariates.R data/covar/geno.vcf.gz data/tensorqtl/expression.iqn.bed.gz 5 5 data/covar/covar.txt

##################
## eQTL mapping ##
##################

module load cuda
python3 ~/ratgtex/scripts/run_tensorqtl.py \
    data/tensorqtl/geno \
    data/tensorqtl/expression.iqn.bed.gz \
    data/tensorqtl/f2.cis_qtl.txt.gz \
    --covariates data/covar/covar.txt \
    --mode cis

module load cuda
python3 ~/ratgtex/scripts/run_tensorqtl.py \
    data/tensorqtl/geno \
    data/tensorqtl/expression.iqn.bed.gz \
    data/tensorqtl/f2.cis_independent_qtl.txt.gz \
    --covariates data/covar/covar.txt \
    --cis_output data/tensorqtl/f2.cis_qtl.txt.gz \
    --mode cis_independent

###########################################
## Get all significant SNPs per cis-eQTL ##
###########################################

mkdir -p data/tensorqtl/nominal
module load cuda
python3 -m tensorqtl \
    data/tensorqtl/geno \
    data/tensorqtl/expression.iqn.bed.gz \
    f2 \
    --covariates data/covar/covar.txt \
    --output_dir data/tensorqtl/nominal \
    --mode cis_nominal

python3 ~/ratgtex/scripts/tensorqtl_all_signif.py \
    data/tensorqtl/f2.cis_qtl.txt.gz \
    data/tensorqtl/nominal/f2 \
    data/tensorqtl/f2.cis_qtl_signif.txt.gz \
    --fdr 0.05

#########
## aFC ##
#########

## Get VCF for aFC.py

plink2 \
    --bfile data/tensorqtl/geno \
    --recode vcf id-paste=iid \
    --out data/afc/geno
bgzip data/afc/geno.vcf
tabix -f data/afc/geno.vcf.gz

python3 ~/ratgtex/tools/aFC/aFC.py \
    --vcf data/afc/geno.vcf.gz \
    --pheno data/afc/expression.log2.bed.gz \
    --qtl <(python3 ~/ratgtex/scripts/prepare_qtl_for_afc.py data/tensorqtl/f2.cis_qtl.txt.gz data/tensorqtl/f2.cis_independent_qtl.txt.gz) \
    --cov data/covar/covar.txt \
    --log_xform 1 \
    --output data/afc/f2.aFC.txt

####################
## Process output ##
####################

Rscript process_eqtl_files.R
