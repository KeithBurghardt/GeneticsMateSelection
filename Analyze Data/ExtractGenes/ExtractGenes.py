import pandas as pd
import numpy as np

# find all significant slices

# if oe == True then we are extracting data from the European ethnicity cohort, otherwise we are looking at the full dataset
oe=True
one_eth = ''
if oe:
    one_eth = '_one_ethnic'
   
# find all significant slices, which are recorded in a separate file. These can be easily extracted from data on p-values
significant_slices_file = 'slices_all'+one_eth+'_cog27.csv'

slices=pd.read_csv(significant_slices_file)

unique_chrom = np.unique([cfile.split('_')[0] for cfile in slices['slice'].values])
# find all genes in chromosome
print('library(biomaRt)\nmart <- useMart(\"ensembl\")\nmart <- useDataset(\"hsapiens_gene_ensembl\", mart)\nattributes <-c(\"ensembl_gene_id\",\"start_position\",\"end_position\",\"strand\",\"hgnc_symbol\",\"chromosome_name\",\"entrezgene_id\",\"ucsc\",\"band\")\nfilters <- c(\"chromosome_name\")')
for chrom in unique_chrom:
    gene_file = '/Volumes/Keith_Burghardt_EHD/GeneticClustering/genes_'+str(chrom)+'.csv'
    print('values <- list(chromosome=\"'+chrom+'\")')
    print('all.genes <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)')
    print("write.csv(all.genes, \""+gene_file+"\", row.names=FALSE)")
gene_slices = {'slice_chromosome':[],'slice_bp1':[],'slice_bp2':[]}
all_genes=[]
for cfile in slices['slice'].values:#["10_111238239-112439771"]:#slices['slice'].values:#all_cov_files:##['/Users/keithab/Dropbox/GeneticClustering/CommonSNPs/common_married_SNPs_9_92830743-92998733.csv']:#common_snp_files:#['/Users/keithab/Dropbox/GeneticClustering/CommonSNPs/common_married_SNPs_6_169671532-170606719.csv']:#common_snp_files:
    # common_married_SNPs_9_92830743-92998733.csv
    #print(cfile)
    chrom = cfile.split('_')[0]
    b1,b2 = cfile.split('_')[1].split('-')
    chrom = int(float(chrom))
    b1=int(float(b1))
    b2=int(float(b2))
    print([chrom,b1,b2])
    gene_file = '/Volumes/Keith_Burghardt_EHD/GeneticClustering/genes_'+str(chrom)+'.csv'
    if os.path.exists(gene_file):
        gene_data = pd.read_csv(gene_file)
        gene_data['start_position'] = gene_data['start_position'].astype(int)
        gene_data['end_position'] = gene_data['end_position'].astype(int)
        genes = gene_data.loc[(gene_data['start_position']<=b2) | (gene_data['end_position']>=b1)].dropna().drop_duplicates()
        all_genes.append(genes)
       
        gene_slices ['slice_chromosome']+=[chrom]*len(genes)
        gene_slices ['slice_bp1']+=[b1]*len(genes)
        gene_slices ['slice_bp2']+=[b2]*len(genes)
all_genes = pd.concat(all_genes).reset_index()

all_genes = pd.concat([all_genes,pd.DataFrame(gene_slices)],axis=1)
print(all_genes.columns)
keep_cols = list(all_genes.columns)
keep_cols.remove('index')
all_genes[keep_cols].to_csv('AllGenes'+one_eth+'.csv',index=False)
