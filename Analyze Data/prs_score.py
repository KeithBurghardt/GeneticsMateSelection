import pandas as pd
import pickle as pk
import numpy as np
from glob import glob
import time,os,sys,ast,random,re
from scipy.stats import spearmanr

num_married = np.load('married_slices/married_slices_chr=1_comm=False_heldout.npy').shape[0]
num_rand = (num_married-1)*num_married-2
outcome_test = np.array(num_married*[1] + num_rand * [0])

slice_features = []

def file2slice(file):
    match = re.search(r'\d+_\d+-\d+', file)
    slice = None
    if match:
        slice = match.group()
    return slice


def find_clusters(slice,p_value_data,coef_columns,p_value,slice2file):
    # get slice file
    slice_file = slice2file[slice]
    # if we can't find it, ignore
    if slice_file is None:
        return None,None
    # all cluster ID names in the file (faster than Pandas)
    # (e.g., 1, 0, ... NOT newclust_1, newclust_2...)
    cluster_ids = [line.replace('\n','') for line in open(slice_file)][1:]
    # common clusters, where we remove newclust_ from the coef columns
    common_clusters = set(cluster_ids).intersection(set([c.replace('newclust','') for c in coef_columns]))
    # common cluster with appropriate names
    # we INCLUDE -1, because we cannot control for the new clusters, like "1" without including "-1"
    # Ex: clusters -1,0,1
    #     coefficiens are one-hot encoded: [1,0,0], [0,1,0], [0,0,1]
    #     if you contain neither -1,0,1 then your coefficient should be 0, but not a zero coefficient if
    #     you have -1, OR 0, OR 1
    #     -1 vs 0 might be insigificant, but obscure clusters are significant because they are not the common clusters
    #     therefore to set -1 coefficients to zero discard the many datapoints we need to include!
    common_clusters = ['newclust' + c for c in common_clusters] # if 'minus1' not in c]
    # if no common clusters, ignore
    if len(common_clusters) == 0:
        return None,None
    # else find p_value ata
    # if no significant p_values, ignore
    if len(p_value_data[common_clusters].values[p_value_data[common_clusters].values<p_value]) == 0:
         return None,None
    # finally convert clusters to pandas dataframe
    cluster_ids = pd.DataFrame({'cluster':cluster_ids}) 
    # create dummy indices
    clusters = pd.get_dummies(cluster_ids.astype(str), prefix=['cluster'])
    clusters.columns = [c.replace('cluster','newclust').replace('-1','minus1').replace('_','') for c in clusters.columns]
    return clusters,common_clusters

def add_prs_array(coef_data,p_value_data,slice_file,prs_score_array,weight,p_value,slice2file):
    # clusters that are in common across HRS, ELSA
    slice = file2slice(slice_file)
    if slice is None:
        return prs_score_array
    clusters,common_clusters = find_clusters(slice,p_value_data,coef_data.columns,p_value,slice2file)
    if clusters is None:
        return prs_score_array
    # we only consider common clusters
    p_value_data = p_value_data[common_clusters]
    for cluster_id in common_clusters:
        # do not add cluster PRS unless it is statistically significant
        if p_value_data[cluster_id].values[0] < p_value:
            # score if significant
            # defaults to 0/1 due to one-hot encoding
            score = clusters[cluster_id].values.astype(float)
            # if we add weighted score
            # we will ignore scores < 0!! - this is because linear regression incorrectly labels random values as ~ -3, less negative values are insignificant
            if coef_data[cluster_id].values[0] < 0:
                continue
            if weight:
                score = (clusters[cluster_id].values)*(coef_data[cluster_id].values[0])
            prs_score_array += score.astype(float)
    return prs_score_array


# slices
index=int(float(sys.argv[2]))
one_ethnic=True
one_eth = ''
if one_ethnic:
    one_eth = '_one_ethnic'
for weight in [True,False]:
    for p_value in [sys.argv[1]]:
        p_value = float(p_value)
        print(p_value)
        prs_array = np.array([0] * len(outcome_test)).astype(float)
        slice_features = set(slice_features)
        train80=False
        slice_prs = True
        covars = ['place','ethnicity','relig','edu','height','cog27']#'bmi'
        feature_sets = ['none','+'.join(covars)]
        absolute_features=True
        sig_slice_data = {'weight':[],'pval':[],'feature':[],'num_slices':[]}
        clust_indices = [0,0]
        for feature_set in feature_sets:
            outfile = 'slice_prs_p-value='+str(-np.log10(p_value))+'_'+feature_set+'_weight='+str(weight)+'_all_abs'+one_eth+'_'+str(clust_indices[0])+'-'+str(clust_indices[1])+'_multicluster_index='+str(index)+'.csv'
            if feature_set!= 'none':
                feature_set = 'place'
            coefs = pd.read_csv('all_slices_'+feature_set+'_coef_hrs'+one_eth+'.csv')
            pvals = pd.read_csv('all_slices_'+feature_set+'_prob_hrs'+one_eth+'.csv')
            prs_array = np.array([0] * len(outcome_test)).astype(float)
            if True:
                zip_file = 'only_zipped_slices'
                slice_files = glob(zip_file.replace('.zip','')+'/elsa_slice_*.csv')
                print('Num slice files: ',len(slice_files))
                slice2file = {}
                snp_coefs = coefs.groupby('snp')
                snp_pvals = pvals.groupby('snp')
                all_snps = set(pvals['snp'].tolist())
                for ii,slice_file in enumerate(slice_files[index*10000:(index+1)*10000]):
                    if ii % 10000 == 0:
                        print(round(ii/len(slice_files),2),'%')
                    slice = file2slice(slice_file)
                    if 'chr'+slice in all_snps:
                        per_snp_array = np.array([0] * len(outcome_test)).astype(float)
                        coef_data = snp_coefs.get_group('chr'+slice)#coefs.loc[coefs['snp']==slice,]
                        p_value_data = snp_pvals.get_group('chr'+slice)#pvals.loc[pvals['snp']==slice,]
                        # only look at "valid" coefficients (>0)
                        valid_c = [c for c in pvals.columns if 'clust' in c]
                        valid_c = [c for c in valid_c if coef_data[c].values[0] > 0]
                        if len(valid_c) == 0:
                            continue
                        # among valid coefficients, if the p-values are all above the threshold, we skip
                        if np.min(p_value_data[valid_c].values) > p_value:
                            continue
                        # otherwise, record file, find weights
                        slice2file[file2slice(slice_file)] = slice_file
                        per_snp_array = add_prs_array(coef_data,p_value_data,slice_file,per_snp_array,weight,p_value,slice2file)
                        if np.sum(np.abs(per_snp_array)) > 0:
                            print(np.min(per_snp_array),np.max(per_snp_array))#
                            if np.max(per_snp_array) > 1.5:
                                print(spearmanr(per_snp_array,outcome_test))
                            prs_array += per_snp_array 
                print('SUM OF PRS: ',np.sum(prs_array))
            pd.DataFrame({'prs':prs_array,'marriage':outcome_test}).to_csv(outfile,index=False)




