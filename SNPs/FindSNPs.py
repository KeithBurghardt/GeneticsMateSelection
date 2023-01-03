one_eth = ''
if oe:
    one_eth = '_one_ethnic'
all_common_snp_files = glob('/Users/keithab/Dropbox/GeneticClustering/CommonSNPs_abs_cog27_all/common_married_*'+one_eth+'*_renamed.csv')
if oe:
    all_common_snp_files = glob('/Users/keithab/Dropbox/GeneticClustering/CommonSNPs_abs_cog27_all/common_married_*'+one_eth+'*_renamed.csv')
print(all_common_snp_files)
all_gwas_snp_files = glob('/Users/keithab/Dropbox/GeneticClustering/SNP_GWAS_abs_cog27_all/*.tsv')
#data = pd.concat([pd.read_csv(file) for file in files])
#for name in data['snp_name'].values: print(name) 
genes = []
print(len(all_common_snp_files))

race = []
social = []
beauty = []
total = []
health = []
all_sig_slices = pd.DataFrame(sig_slices_all,columns=['covar','chrom','bp1','bp2','mean_p'])
print(all_sig_slices)
unique_features = []
#all_cov_files = ['/Users/keithab/Dropbox/GeneticClustering/CommonSNPs_abs/common_married_SNPs_'+sig_slice+'.csv' for sig_slice in ['3_73308295-74010362', '6_169671532-170606719', '9_113480169-114654601','9_115953701-115960547', '9_115960547-116299953']]
all_SNP_data={'chromosome':[],'slice_bp1':[],'slice_bp2':[],'bp':[],'snp':[],'common_allele':[],'mapped_gene':[],'prefex':[],'trait':[]}
gwas_snps_data = []
for cfile in all_common_snp_files:#all_cov_files:#['/Users/keithab/Dropbox/GeneticClustering/CommonSNPs/common_married_SNPs_9_92830743-92998733.csv']:#common_snp_files:#['/Users/keithab/Dropbox/GeneticClustering/CommonSNPs/common_married_SNPs_6_169671532-170606719.csv']:#common_snp_files:
    # common_married_SNPs_9_92830743-92998733.csv
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
    print([chrom,b1,b2])
    non_sig = len(all_sig_slices.loc[(all_sig_slices['chrom']==chrom) & (all_sig_slices['bp1']==b1) & (all_sig_slices['bp2']==b2),]) > 0
    if not non_sig: continue
    chrom = str(chrom)
    b1=str(b1)
    b2=str(b2)
    gwas_file = '/Users/keithab/Dropbox/GeneticClustering/SNP_GWAS_abs_cog27_all/gwas-association-downloaded_2022-12-05-chromosomeName_ '+chrom+' AND chromosomePosition_[ '+b1+' TO '+b2+' ].tsv'
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
    base_gwas = common_gwas_snps[['STRONGEST SNP-RISK ALLELE','base','DISEASE/TRAIT','MAPPED_GENE','SNPS']].values.astype(str)
    base_gwas[:,0] = [b.split('-')[-1] for b in base_gwas[:,0]]
    common_gwas_snps['slice_bp1']=[b1]*len(common_gwas_snps)
    common_gwas_snps['slice_bp2']=[b2]*len(common_gwas_snps)    
    gwas_snps_data.append(common_gwas_snps)
    for risk_base,common_base,trait,gene,snp in base_gwas:
        prefex = 'Lower '
        if risk_base == common_base:
            prefex = 'Higher '
        
        output = snp+' & '+common_base+' & '+gene+' & '+prefex + trait+'\\\\'

        if output not in unique_features and os.path.exists(gwas_file):
            print(output)
            unique_features.append(output)
        elif not os.path.exists(gwas_file):
            print('FILE NOT FOUND')
        total += [snp]
        if 'melanin' in trait:
            race += [snp]
        if 'cogn' in trait.lower() or 'educat' in trait.lower() or 'brain' in trait.lower() or 'intelligence' in trait.lower() or 'income' in trait.lower() or 'MTAG' in trait:
            social += [snp]
        if 'bald' in trait.lower() or 'testo' in trait.lower(): 
            beauty += [snp]
    genes += list(np.unique(base_gwas[:,3].astype(str)))
gwas_snps_data = pd.concat(gwas_snps_data)
gwas_snps_data.to_csv('AllSNPs'+one_eth+'.csv',index=False)
