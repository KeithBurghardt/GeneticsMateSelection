import pandas as pd
import numpy as np
## Code one can copy-paste into R to get all genes within each chromosome
# If oe == True, we only look at European ethnicity cohort
oe = True
one_eth = ''
if oe:
    one_eth = '_one_ethnic'
# file that keeps all slices with p-value below the Bonferroni correction
slices=pd.read_csv('slices_all'+one_eth+'_cog27.csv')
# file is of the form: 
# slice
# 10_122598339-123215850
# 10_111238239-112439771
# 4_66556084-69246945

unique_chrom = np.unique([cfile.split('_')[0] for cfile in slices['slice'].values])
print('library(biomaRt)\nmart <- useMart(\"ensembl\")\nmart <- useDataset(\"hsapiens_gene_ensembl\", mart)\nattributes <-c(\"ensembl_gene_id\",\"start_position\",\"end_position\",\"strand\",\"hgnc_symbol\",\"chromosome_name\",\"entrezgene_id\",\"ucsc\",\"band\")\nfilters <- c(\"chromosome_name\")')


for chrom in unique_chrom:
    gene_file = '/Volumes/Keith_Burghardt_EHD/GeneticClustering/genes_'+str(chrom)+'.csv'
    print('values <- list(chromosome=\"'+chrom+'\")')
    print('all.genes <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)')
    print("write.csv(all.genes, \""+gene_file+"\", row.names=FALSE)")
    
# extracting genes in each slice given genes in each chromosome
gene_slices = {'slice_chromosome':[],'slice_bp1':[],'slice_bp2':[]}
all_genes=[]
# for each fil
for cfile in slices['slice'].values:
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
