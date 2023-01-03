import pandas as pd
import numpy as np

# This file extracts features associated with SNP shared between married respondences who share a given slice

one_eth = ''
if oe:
    one_eth = '_one_ethnic'
all_common_snp_files = glob('/Users/keithab/Dropbox/GeneticClustering/CommonSNPs_abs_cog27_all/common_married_*'+one_eth+'*_renamed.csv')
if oe:
    all_common_snp_files = glob('/Users/keithab/Dropbox/GeneticClustering/CommonSNPs_abs_cog27_all/common_married_*'+one_eth+'*_renamed.csv')
print(all_common_snp_files)
all_gwas_snp_files = glob('/Users/keithab/Dropbox/GeneticClustering/SNP_GWAS_abs_cog27_all/*.tsv')
#for name in data['snp_name'].values: print(name) 
genes = []
print(len(all_common_snp_files))


all_sig_slices = pd.DataFrame(sig_slices_all,columns=['covar','chrom','bp1','bp2','mean_p'])
print(all_sig_slices)
unique_features = []

for cfile in all_common_snp_files:
    print(cfile)
    if oe:
        chrom = cfile.replace('_cog27','').split('_')[-5]
        b1,b2 = cfile.replace('_cog27','').split('_')[-4].replace('.csv','').split('-')
    else:
        chrom = cfile.replace('_cog27','').split('_')[-3]
        b1,b2 = cfile.replace('_cog27','').split('_')[-2].replace('.csv','').split('-')
    chrom = int(float(chrom))
    b1=int(float(b1))
    b2=int(float(b2))
    non_sig = len(all_sig_slices.loc[(all_sig_slices['chrom']==chrom) & (all_sig_slices['bp1']==b1) & (all_sig_slices['bp2']==b2),]) > 0
    if not non_sig: continue
    chrom = str(chrom)
    b1=str(b1)
    b2=str(b2)
    # files from ebi.ac.uk on traits associated with SNPs in a given slice, e.g., https://www.ebi.ac.uk/gwas/regions/chr10:122598339-123215850
    gwas_file = 'gwas-association-downloaded_2022-12-05-chromosomeName_ '+chrom+' AND chromosomePosition_[ '+b1+' TO '+b2+' ].tsv'
    if not os.path.exists(gwas_file):
        
        new_slice = chrom+':'+b1+'-'+b2
        continue
    common_data = pd.read_csv(cfile)
    gwas_data = pd.read_csv(gwas_file,sep='\t')
    if len(common_data) > 0 and len(gwas_data) > 0:
        common_gwas_snps = pd.merge(common_data,gwas_data,left_on='snp_name',right_on='SNPS',how='inner')
    else:
        print('error')

        continue
    risk = []
    common_gwas_snps['slice_bp1']=[b1]*len(common_gwas_snps)
    common_gwas_snps['slice_bp2']=[b2]*len(common_gwas_snps)    
gwas_snps_data = pd.concat(gwas_snps_data)

gwas_snps_data.to_csv('AllSNPs'+one_eth+'.csv',index=False)
