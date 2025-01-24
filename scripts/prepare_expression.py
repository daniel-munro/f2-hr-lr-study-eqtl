from gtfparse import read_gtf
import numpy as np
import pandas as pd
import rnaseqnorm

cpm = pd.read_csv('data/expression/ResidualsMatrix_wMean.csv', index_col=0)
cpm.index.name = 'gene_id'
iqn = rnaseqnorm.normalize_quantiles(cpm)
iqn = rnaseqnorm.inverse_normal_transform(iqn)

counts = pd.read_csv('data/expression/GSE225746_Elaine_es103rn_FCgene_4R.txt.gz', sep='\t', index_col=0)
counts.index.name = 'gene_id'
counts = counts.loc[cpm.index, cpm.columns]
log2 = np.log2(counts + 1)

samples = pd.read_csv('data/Akil_bHRbLR_F0_F2_IDs.csv', dtype=str)
# keep samples where LibraryID is not missing
samples = samples.loc[~samples['LibraryID'].isnull(), :]
samples = samples.set_index('LibraryID')

# Replace LibraryID column headers with numeric IDs used for genotypes
log2.columns = samples.loc[log2.columns, 'RatNumForNIHstudy']
iqn.columns = samples.loc[iqn.columns, 'RatNumForNIHstudy']
# Keep only samples with genotype information
samples_geno = pd.read_csv('data/samples_geno.txt', header=None, dtype=str)[0].values
samples_geno = [x for x in samples_geno if x in log2.columns]
print(f'Keeping {len(samples_geno)} samples with genotype data')
log2 = log2.loc[:, samples_geno]
iqn = iqn.loc[:, samples_geno]

anno = read_gtf('data/reference/Rattus_norvegicus.Rnor_6.0.104.gtf.gz')
anno = anno.loc[anno['feature'] == 'gene', :]
anno['tss'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
anno['start'] = anno['tss'] - 1  # BED coordinates are 0-based
anno['end'] = anno['tss']
anno = anno.loc[anno['seqname'].isin([str(i) for i in range(1, 21)]), :]
anno['#chr'] = anno['seqname'].astype(int)  # for sorting
anno = anno.sort_values(['#chr', 'start'])
anno = anno[['#chr', 'start', 'end', 'gene_id']]

log2 = anno.merge(log2.reset_index(), on='gene_id', how='inner')
iqn = anno.merge(iqn.reset_index(), on='gene_id', how='inner')

log2.to_csv(f'data/expression/expression.log2.bed', sep='\t', index=False, float_format='%g')
iqn.to_csv(f'data/expression/expression.iqn.bed', sep='\t', index=False, float_format='%g')
