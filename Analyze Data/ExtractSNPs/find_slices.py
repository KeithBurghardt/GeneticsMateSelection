import pandas as pd
import os,random

# file to read in
oe = True
one_eth = ''
if oe:
    one_eth = '_one_ethnic'
for file in ['slices_all'+one_eth+'_cog27.csv']:#,'slices_none_cog27.csv']:
    slices = list(pd.read_csv(file)['slice'].values)
    random.shuffle(slices)
    for slice in slices:
        print(slice)
        covar = file.split('_')[1].split('.')[0]
        directory = 'common_snps_'+covar+'_cog27/'
        outfile_common = directory + 'common_married_SNPs_'+slice+one_eth+'_cog27.csv'
        if not os.path.exists(outfile_common):
            command = "python3 study_slices.py --slice "+str(slice)+' --outfile_common '+outfile_common + ' --one_ethnicity '+str(oe)
            print(command)
            os.system(command)
        # we extracted all SNPs whose names come from our local MAP file
        # we instead want to find their more universal name from https://genome.ucsc.edu/
        SNP_names = 'snp_rename_cog27/'+slice+'.csv'
        sep = ','
        with open(SNP_names) as f:
            if '\t' in f.readline():
                sep = '\t'
        SNP_rename_data = pd.read_csv(SNP_names,sep=sep)
        common_snps = pd.read_csv(outfile_common)
        common_snps['bp'] = common_snps['bp'].astype(int)
        SNP_rename_data['chromEnd'] = SNP_rename_data['chromEnd'].astype(int)
        rename_dict = {base:SNP_rename_data.loc[SNP_rename_data['chromEnd']==base,'name'].values[0] for base in common_snps['bp'].values if base in SNP_rename_data['chromEnd'].values}
        common_snps['snp_name'] = [rename_dict[pos] if pos in rename_dict.keys() else common_snps.loc[common_snps['bp']==pos,'snp_name'].values[0] for pos in common_snps['bp'].values]
        # names of known SNPs
        common_snps.to_csv(outfile_common[:-4]+'_renamed.csv',index=False)
# output file name indicates the slice position
# output file has the following columns:
# "Unnamed: 0": (legacy, ignore)
# "index": (legacy, ignore)
# chrom: chromosome SNP is located
# snp_name: SNP name
# cM: centiMorgan position
# bp: basepair position
# base: common allele shared in married cohort
