import numpy as np
import pandas as pd
import pickle as pk
import os,time,sys, ast
from glob import glob
import random

def main():
    demog_data = pd.read_csv('/project/burghard_687/genetic_data/hrs_pcs15506_cog27.csv')
    eur_descent = set(demog_data.loc[demog_data['ethnicity']==1,'FID'].tolist())
    # data directory
    dump_dir = '/project/burghard_687/genetic_data/'
    # argument: --slice chrXX_bp1-bp2
    args = sys.argv[1:]
    args = {key:val for key,val in zip(args[::2],args[1::2])}
    # married pairs
    married_pairs=np.load(dump_dir+'married_pairs.npy')
    # slice from argument
    slice = args['--slice']
    # oe == True, we only look at European ethnicity cohort
    oe = ast.literal_eval(args['--one_ethnicity'])
    one_eth = ''
    if oe:
        one_eth = '_one_ethnic'

    outfile_pairs = 'ibd_pairs_'+args['--slice']+one_eth+'.csv'
    if '--outfile_pairs' in args.keys():
        outfile_pairs = args['--outfile_pairs']
    outfile_common = 'common_married_SNPs_'+args['--slice']+one_eth+'.csv'
    if '--outfile_common' in args.keys():
        outfile_common = args['--outfile_common']
    # chromosome
    chrom = slice.split('_')[0].replace('chr','')
    # bp range of slice
    slice_range = np.array(slice.split('_')[1].split('-')).astype(int)
    
    # ilash file associated with that chromosome
    # we use this to find who shared SNPs are within the slice range
    ilash_file = dump_dir+'ILASH_output/dist_hrs_' + chrom
    ibd_pairs = []
    with open(ilash_file,'r') as f:
        for line in f:
            # line tab-separated with values: family_id1,family_hap1,family_id2,family_hap2,chrom,bp1,bp2,snp1,snp2,cM,p
            # family_idX: Subject ID for person X
            # family_hapX: haploid of person X
            # chrom: chromosome
            # bp1,bp2: basepair positions of IBD segments
            # snp1,snp2: snps associated with basepair ends
            # cM: centiMorgan distance of segment
            # p: metric for segment quality
            # for each line in the IBD file...
            line = line.replace('\n','').split('\t')
            # associated info each line
            family_id1,family_hap1,family_id2,family_hap2,chrom,bp1,bp2,snp1,snp2,cM,p = line
            bp1 = float(bp1)
            bp2 = float(bp2)
            # look for segments that fully include the slice
            if bp1 <= slice_range[0] and bp2 >= slice_range[1]:
                # if pairs are marriages
                # first find all marriages with family_ids
                all_fam1_marriages = set([m1 for m1,m2 in married_pairs if m2 == float(family_id1)] + [m2 for m1,m2 in married_pairs if m1 == float(family_id1)])
                # if this is a marriage pair
                if float(family_id2) in all_fam1_marriages:#float(family_id1) in all_marriages and float(family_id2) in all_marriages:
                    if oe and float(family_id1) in eur_descent and float(family_id2) in eur_descent:
                        ibd_pairs.append([family_hap1,family_hap2])
                    elif not oe:
                        ibd_pairs.append([family_hap1,family_hap2])
    # find snps in common with ALL users within this slice
    # MAP file for this chromosome
    map = pd.read_csv('/project/arpawong_181/HRS_AsMa/keith/phased/chr'+chrom+'_imputed.map',sep='\t',header=None)
    map.columns = ['chrom','snp_name','cM', 'bp']
    # SNPs within the slice
    snp_data = map.loc[(map['bp'] >=slice_range[0])  & (map['bp'] <=slice_range[1]),]
    print(snp_data)
    # what are the line positions in the file associated with the SNPs
    snp_locations = snp_data.index.values
    # PED file associated with the chromosome
    ped_file = '/project/arpawong_181/HRS_AsMa/keith/phased/chr'+chrom+'.ped'
    # all subjects we want to analyze
    subjects = np.array(ibd_pairs).flatten()
    subjects = set([u.split('_')[0] for u in subjects])
    ped_data = {}
    # collect subset of much larger PED file associated with the subjects of interest
    with open(ped_file,'r') as p:
         for n,line in enumerate(p):
             # ignore rest of line without subject data
             subj = line[:100].split('\t')[0]
             if subj in subjects:
                 all_gen_line = np.array(line.split('\t')[6:]).astype(str)

                 # data for haplotype 1, haplotype2
                 ped_data[subj] = [all_gen_line[2*snp_locations],all_gen_line[2*snp_locations+1]]
    # check SNPs in common across all pairs
    common_snps = []
    # look at SNPs within slice       
    snp_data = snp_data.reset_index()
    for n,row in snp_data.iterrows():
        # if all elements are the same, set will be length 1
        candidate_snp = set()
        # look at SNPs shared (controlling for haplotypes)
        for p1,p2 in ibd_pairs:
            subj1 = p1.split('_')[0]
            hap1 = int(float(p1.split('_')[1]))
            subj2 = p2.split('_')[0]
            hap2 = int(float(p2.split('_')[1]))
            snp1 = ped_data[subj1][hap1][n]
            snp2 = ped_data[subj2][hap2][n]
            # this is a sanity check, SNPs must be the same
            if snp1 != snp2: print('ERROR'); break
            candidate_snp.add(snp1)
        # now that we know all the SNPs that users share, check if they are all the same
        # if all subjects only share 1 SNP, this may be important
        if len(candidate_snp) == 1: 
            row['base']= list(candidate_snp)[0]
            common_snps.append(list(row.values))
    cols = list(snp_data.columns) + ['base'] # base: G, A, T, or C
    
    organized_data = pd.DataFrame(common_snps,columns=cols)#concat(common_snps)
    organized_data.to_csv(outfile_common)
    # 'common_married_SNPs_'+args['--slice']+'.csv')
    # save marriages within slice to determine covariates
    pd.DataFrame(ibd_pairs,columns=['sample1','sample2']).to_csv(outfile_pairs,index=False)
if __name__ == "__main__":
    main()




