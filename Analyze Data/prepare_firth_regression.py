import numpy as np
import pandas as pd
import pickle as pk
import os,time,sys
from glob import glob
import random
from itertools import permutations,combinations
from more_itertools import consume
from multiprocessing import Pool
import warnings
import firth


def find_sex(df_pca,IDs):
    sex={}
    for ii in IDs:
        s=df_pca.loc[df_pca['IID']==ii,'sex'].values
        if len(s)>0:
            s = s[0]
        else:
            s = None
        sex[ii]=s
    return sex

def sex_dict():
    directory='/project/arpawong_181/HRS_AsMa/keith/'
    df_pca = pd.read_csv(directory+'hrs_pcs15506_cog27.csv')

    IDs=df_pca['IID'].values
    IDs = IDs[~np.isnan(IDs)]
    sex = find_sex(df_pca,IDs)
    return sex

def collect_common_slices(married_pairs,rand_pairs,community_ibd,married_only):
    directory= '/project/burghard_687/genetic_data/ILASH_output/'
    bp_ranges=pd.read_csv(directory+'Slices_Per_Chromosome.csv')
    #../ILASH_output/Slices_Per_Chromosome.csv')
    slice_features = []
    all_rand_slice = []
    all_married_slice = []
    sex = sex_dict()
    married_nodes = set(list(np.unique(np.array(married_pairs).flatten())))
    married_pairs_tuples = set([(p1,p2) for p1,p2 in married_pairs])
    rand_nodes = set(list(np.unique(np.array(rand_pairs).flatten())))
    all_nodes = married_nodes.union(rand_nodes)
    s_features = 0
    for chrom in range(1,23):
        all_slices = sorted(list(set(list(bp_ranges.loc[bp_ranges['chr']==chrom,'slice'].values))))
        # allow everything at the edge to be counted
        all_slices[-1] += 10
        slice_bins = [[b1,b2] for b1,b2 in zip(all_slices[:-1],all_slices[1:]) if b1<b2]
        rand_slice = np.zeros((len(rand_pairs),len(slice_bins)))#
        married_slice = np.zeros((len(married_pairs),len(slice_bins)))#
        
        print('CHROMOSOME: ',chrom)
        rand_file = 'slices/rand_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
        married_file = 'slices/married_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
        unique_slices = set()
        if os.path.exists(rand_file) and os.path.exists(married_file):
            
            rand_slice = np.load(rand_file)
            married_slice = np.load(married_file)
            keep_pos = []
            for slice in range(len(rand_slice[0])):
                common_pairs_rand = np.sum(rand_slice[:,slice])
                common_pairs_married = np.sum(married_slice[:,slice])
                # if more than one person was married with this slice, keep
                if common_pairs_rand > 0 and common_pairs_married > 0:
                    keep_pos.append(slice)
            slice_bins = np.array(slice_bins)
            slice_bins = slice_bins[keep_pos]
            #slice_features+=['chr'+str(chrom)+'_'+str(b1)+'-'+str(b2) for b1,b2 in slice_bins]
            s_features +=len(keep_pos)
            print([len(slice_features),s_features])
            slice_features+=['chr'+str(chrom)+'_'+str(b1)+'-'+str(b2) for b1,b2 in slice_bins]
            all_married_slice.append(np.array(married_slice)[:,keep_pos])
            all_rand_slice.append(np.array(rand_slice)[:,keep_pos])
            continue
        slice_features+=['chr'+str(chrom)+'_'+str(b1)+'-'+str(b2) for b1,b2 in slice_bins]
        file = directory+'dist_hrs_'+str(chrom)#'hrs123_qc1all'+str(chrom)
        if community_ibd:
            marriage_comm_file='../parsed_data/communities/dist_hrs_'+str(chrom)+'_marriage_communities_ilash-bins_mcl=False_married_only='+str(married_only)+'.pkl'
            if not os.path.exists(marriage_comm_file): continue
            communities = pk.load(open(marriage_comm_file,'rb'))
            for nn,[community,marriage_pairs,bp1,bp2] in enumerate(zip(communities['community'], communities['pairs'], communities['bp_min'],communities['bp_max'])):
                if len(community) <=1: continue
                community = set(community)
                mp_tuples = set([tuple(pair) for pair in marriage_pairs])
                bp_pos = slice_bins.index([bp1,bp2])
                if nn % 1000 == 0:
                    print(round(nn/len(communities['bp_min'])*100,2))

                # check if data are within communities = 1, else 0                
                if not community.isdisjoint(rand_nodes):
                    for ii,[p1,p2] in enumerate(rand_pairs):
                        if p1 in community and p2 in community:
                            rand_slice[ii][bp_pos] = 1

                if not mp_tuples.isdisjoint(married_pairs_tuples):
                    for jj,[p1,p2] in enumerate(married_pairs):
                        if tuple([p1,p2]) in mp_tuples:
                            married_slice[jj][bp_pos] = 1
        else:
            data=pd.read_csv(file,sep='\t',header=None)
            data.columns=['family_id1','sample_id1','family_id2','sample_id2','chr','bp_ind1','bp_ind2','snp_id1','snp_id2','cM','kmer_agree']
            data = data[['family_id1','family_id2','bp_ind1','bp_ind2','cM']]
            data= data.dropna()
            data['bp_ind1'] = data['bp_ind1'].astype(int)
            data['bp_ind2'] = data['bp_ind2'].astype(int)
            data['family_id1'] = data['family_id1'].astype(int)
            data['family_id2'] = data['family_id2'].astype(int)
            # remove duplicate data
            data = data.loc[(data['cM']<200)&(data['family_id1']!=data['family_id2']),]
            sorted_data=data.sort_values(by='bp_ind1')
            prev_b=-1
            for nn,[bp1,bp2] in enumerate(slice_bins):
                if nn % 10 == 0:
                    print(round(nn/len(slice_bins)*100,2))
                if prev_b>np.mean([bp1,bp2]):
                    print('ERROR',[bp1,bp2],prev_b)
                prev_b=np.mean([bp1,bp2])
                # split by data within this range
                # if the far ends of the ilash output are within this bin
                edges = data.loc[(data['bp_ind2']>=bp1)&(data['bp_ind1']<bp2),['family_id1','family_id2']].values
                # no self loops
                edges = np.array([[e1,e2] for e1,e2 in edges if e1 != e2])
                edges = set([(e1,e2) if e1<e2 else (e2,e1) for e1,e2 in edges])
                # only look at edges that may be in rand/married edges
                edges = [(e1,e2) for e1,e2 in edges if e1 in sex.keys() and e2 in sex.keys() and e1 in all_nodes and e2 in all_nodes]
                edges = set([(e1,e2) for e1,e2 in edges if sex[e1] != sex[e2]])
                if len(edges) == 0: continue
                # check if data in ilash edges
                bp_pos = slice_bins.index([bp1,bp2])
                # check if data are within communities = 1, else 0
                for ii,[p1,p2] in enumerate(rand_pairs):
                    # if random pairs within common slices
                    if not set([(p1,p2)]).isdisjoint(edges):
                        rand_slice[ii][bp_pos] = 1
                        #print(np.sum([rand_slice[kk][bp_pos] for kk in range(len(rand_slice))]))
                for jj,[p1,p2] in enumerate(married_pairs):
                    # if married pairs within common slices
                    if not set([(p1,p2)]).isdisjoint(edges):
                        married_slice[jj][bp_pos] = 1
        all_married_slice.append(np.array(married_slice))
        all_rand_slice.append(np.array(rand_slice))
        
        np.save('slices/rand_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy',rand_slice)
        np.save('slices/married_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy',married_slice)
    all_married_slice = np.concatenate(all_married_slice,axis=1)
    all_rand_slice = np.concatenate(all_rand_slice,axis=1)
    return all_married_slice,all_rand_slice,slice_features


def collect_data(married_pairs,rand_pairs,community_ibd,married_only,train80=False,output_fold=-1):
    dump_dir = '/project/burghard_687/genetic_data/'#'/project/burghard_687/genetic_data/'
    X_file = dump_dir+'slice_features_comm='+str(community_ibd)+'_cog27.npy'
    outcome_file = dump_dir+'slice_outcome_cog27.npy'
    feature_file = dump_dir+'slice_feature_cols_cog27.npy'
    fold=4
    # no cog
    X_file = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(fold)+'.npy'
    outcome_file = dump_dir+'slice_outcome_train80_fold='+str(fold)+'.npy'
    feature_file = dump_dir+'slice_feature_cols_train80_fold='+str(fold)+'.npy'
    X_file_train = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(fold)+'.npy'
    X_file_test = dump_dir+'slice_features_comm='+str(community_ibd)+'_test20_fold='+str(fold)+'.npy'
    outcome_file_test = dump_dir+'slice_outcome_test20_fold='+str(fold)+'.npy'
    outcome_file_train = dump_dir+'slice_outcome_train80_fold='+str(fold)+'.npy'
    # with cog
    X_file = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(fold)+'_cog27.npy'
    outcome_file = dump_dir+'slice_outcome_train80_fold='+str(fold)+'_cog27.npy'
    feature_file = dump_dir+'slice_feature_cols_train80_fold='+str(fold)+'_cog27.npy'
    X_file_train = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(fold)+'_cog27.npy'
    X_file_test = dump_dir+'slice_features_comm='+str(community_ibd)+'_test20_fold='+str(fold)+'_cog27.npy'
    outcome_file_test = dump_dir+'slice_outcome_test20_fold='+str(fold)+'_cog27.npy'
    outcome_file_train = dump_dir+'slice_outcome_train80_fold='+str(fold)+'_cog27.npy'


    if not (os.path.exists(X_file) and os.path.exists(outcome_file) and os.path.exists(feature_file)): 
        X_rand=np.load(dump_dir+'X_rand_cog27.npy')
        X_married=np.load(dump_dir+'X_married_cog27.npy')
        feature_cols_orig = np.load(dump_dir+'feature_cols_cog27.npy')        
        # create more features: 1 if links are the same between pairs, 0 otherwise
        #if True:#
        if not os.path.exists(dump_dir+'slice_features_cog27.npy'):
            all_married_slice,all_rand_slice,slice_features = collect_common_slices(married_pairs,rand_pairs,community_ibd,married_only)
            #PRS 
            prs_married = np.array([[0,0,0,0]]*len(married_pairs))
            prs_rand = np.array([[0,0,0,0]]*len(rand_pairs))
            X_rand = np.concatenate([X_rand,np.array(all_rand_slice),np.array(prs_rand)],axis=1)
            X_married = np.concatenate([X_married,np.array(all_married_slice),np.array(prs_married)],axis=1)
            np.save(dump_dir+'X_rand_full_cog27.npy',X_rand)
            np.save(dump_dir+'X_married_full_cog27.npy',X_married)
            np.save(dump_dir+'slice_features_cog27.npy',slice_features)
        else:
            print('LOADING LARGEST DATASETS...')
            X_rand = np.load(dump_dir+'X_rand_full_cog27.npy')
            X_married = np.load(dump_dir+'X_married_full_cog27.npy')
            slice_features = np.load(dump_dir+'slice_features_cog27.npy')
            print('DATASET LOADED...')

        outcome_folds = []
        n_folds = 5
        folds = list(range(n_folds))
        train_m_len = int((n_folds-1)/(n_folds)*len(X_married))
        fold_m_len = int(1/(n_folds)*len(X_married))
        # find all pairs who are strictly within married fold
        for vf in folds:
            if not os.path.exists(dump_dir+'slice_features_comm=False_train80_fold='+str(vf)+'_3_cog27.npy'):
                #'X_fold='+str(vf)+'.npy'):

                X_file = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(vf)+'_cog27.npy'
                outcome_file = dump_dir+'slice_outcome_train80_fold='+str(vf)+'_cog27.npy'
                feature_file = dump_dir+'slice_feature_cols_train80_fold='+str(vf)+'_cog27.npy'
                X_file_train = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(vf)+'_cog27.npy'
                X_file_test = dump_dir+'slice_features_comm='+str(community_ibd)+'_test20_fold='+str(vf)+'_cog27.npy'
                outcome_file_test = dump_dir+'slice_outcome_test20_fold='+str(vf)+'_cog27.npy'
                outcome_file_train = dump_dir+'slice_outcome_train80_fold='+str(vf)+'_cog27.npy'


                #fold_r_len = int(1/(n_folds)*len(X_rand))

                married_fold = X_married[vf*fold_m_len:(vf+1)*fold_m_len]
                married_subjects_fold = set(list(np.array(married_pairs[vf*fold_m_len:(vf+1)*fold_m_len]).flatten()))
                good_data_pos = [pos for pos,line in enumerate(married_fold) if None not in line]
                married_fold = married_fold[good_data_pos]

                rand_subjects_test = [pos for pos,[p1,p2] in enumerate(rand_pairs) if p1 in married_subjects_fold and p2 in married_subjects_fold]
                rand_fold = X_rand[rand_subjects_test]
                #good_data_pos = [pos for pos,line in enumerate(rand_fold) if None not in line]
                #rand_fold = rand_fold[good_data_pos]


                print('rand fold made...')
                X_test = np.concatenate([married_fold,rand_fold],axis=0)
                #X_test = X_collected_fold[good_data_pos]

                fold_m_len = len(married_fold)
                fold_r_len = len(rand_fold)
                outcome_fold = np.array([1]*fold_m_len+[0]*fold_r_len)
                #outcome_test = outcome_fold[good_data_pos]
                np.save(X_file_test,X_test)
                np.save(outcome_file_test,outcome_fold)
                del X_test

                # train data: look at all links disconnected from test
                # all married pairs not in the fold
                if vf == 0:
                    married_train = X_married[(vf+1)*fold_m_len:]
                elif vf == folds[-1]:
                    married_train = X_married[:vf*fold_m_len]
                else:
                    married_train = np.concatenate([X_married[:vf*fold_m_len], X_married[(vf+1)*fold_m_len:]],axis=0)
                good_data_pos = [pos for pos,line in enumerate(married_train) if None not in line]
                married_train = married_train[good_data_pos]
                #set_rand_subjects_test = set(rand_subjects_test)
                # all random pairs such that neither male nor female is in test set (no overlapping edges)
                rand_subjects_train = [pos for pos,[p1,p2] in enumerate(rand_pairs)  if p1 not in married_subjects_fold and p2 not in married_subjects_fold]#if pos not in set_rand_subjects_test]
                print('rand_train')
                rand_train = X_rand[rand_subjects_train[:330000]]
                #good_data_pos = [pos for pos,line in enumerate(rand_train) if None not in line]
                #print('rand_train good data')
                #rand_train = rand_train[good_data_pos]
                print('concatenate')
                X_train = np.concatenate([married_train,rand_train],axis=0)
                train_m_len = len(married_train)
                train_r_len = len(rand_train)
                outcome_train = np.array([1]*train_m_len+[0]*train_r_len)
                np.save(outcome_file_train,outcome_train)
                
                np.save(X_file_train,X_train)
                del rand_train
                del X_train
                print('rand_train 2')

                rand_train = X_rand[rand_subjects_train[330000:660000]]
                #good_data_pos = [pos for pos,line in enumerate(rand_train) if None not in line]
                #rand_train = rand_train[good_data_pos]
                train_r_len = len(rand_train)
                outcome_train = np.array([0]*train_r_len)
                np.save(outcome_file_train[:-4]+'_2_cog27.npy',outcome_train)
                np.save(X_file_train[:-4]+'_2_cog27.npy',rand_train)
                del	rand_train
                #del X_train
                print('rand_train 3')

                rand_train = X_rand[rand_subjects_train[660000:]]
                #good_data_pos = [pos for pos,line in enumerate(rand_train) if None not in line]
                #rand_train = rand_train[good_data_pos]
                train_r_len = len(rand_train)
                outcome_train = np.array([0]*train_r_len)
                np.save(outcome_file_train[:-4]+'_3_cog27.npy',outcome_train)
                np.save(X_file_train[:-4]+'_3_cog27.npy',rand_train)
                del	rand_train
                #del X_train
                
                feature_cols = np.array(list(np.array(feature_cols_orig).astype(str)) + list(np.array(slice_features).astype(str))+['PRS_dummy1','PRS_dummy2','PRS_dummy3','PRS_dummy4'])
                np.save(feature_file,feature_cols)


    outcome_file = dump_dir+'slice_outcome_train80_fold='+str(output_fold)+'_cog27.npy'
    feature_file = dump_dir+'slice_feature_cols_train80_fold='+str(output_fold)+'_cog27.npy'
    X_file_train = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(output_fold)+'_cog27.npy'
    X_file_test = dump_dir+'slice_features_comm='+str(community_ibd)+'_test20_fold='+str(output_fold)+'_cog27.npy'
    outcome_file_test = dump_dir+'slice_outcome_test20_fold='+str(output_fold)+'_cog27.npy'
    outcome_file_train = dump_dir+'slice_outcome_train80_fold='+str(output_fold)+'_cog27.npy'
    print(X_file_test)
    #del X_rand
    X_test = np.load(X_file_test,allow_pickle = True)
    X_train = np.concatenate([np.load(X_file_train,allow_pickle = True),np.load(X_file_train[:-4]+'_2_cog27.npy',allow_pickle = True),np.load(X_file_train[:-4]+'_3_cog27.npy',allow_pickle = True)],axis=0)
    outcome_test = np.load(outcome_file_test,allow_pickle = True)
    outcome_train = np.concatenate([np.load(outcome_file_train,allow_pickle = True),np.load(outcome_file_train[:-4]+'_2_cog27.npy',allow_pickle = True),np.load(outcome_file_train[:-4]+'_3_cog27.npy',allow_pickle = True)],axis=0)
    feature_cols = np.load(feature_file)#,allow_pickle = True)
    return [X_train,X_test],[outcome_train,outcome_test],feature_cols


def run_r(slice):
    #dump_dir = ''
    dump_dir = '/project/burghard_687/genetic_data/'
    slice,train80 = slice
    # make R command line argument
    n_folds = 5
    for fold in range(n_folds):
        if not os.path.exists(dump_dir+"slices80/test_"+slice+"_test20_fold="+str(fold)+"_cog27.csv"): continue
        if not os.path.exists(dump_dir+"slices80/train_"+slice+"_train80_fold="+str(fold)+"_cog27.csv"): continue
        if np.sum(pd.read_csv(dump_dir+"slices80/train_"+slice+"_train80_fold="+str(fold)+"_cog27.csv")["slice"].values) == 0: continue
        print('fold: ',fold)    
        command = "Rscript firth_fit.r --slice "+slice +" --train80 "+str(train80) + " --fold "+str(fold)
        print(command)
        # run command
        os.system(command)    
        # ignore folds if we train on all data
        if not train80: return

def firth_r_GWAS(num_threads):
        train80 = True
        dump_dir = '/project/burghard_687/genetic_data/'
        married_only = True
        community_ibd=False
        load_data =False
        if not os.path.exists(dump_dir + "train_X_non_slice_features_train80_fold=0_cog27.csv") or not os.path.exists(dump_dir+"train_outcome_train80_fold=0_cog27.csv"):
            load_data = True
        n_folds = 5
        folds = list(range(n_folds))[::-1]
        for fold in folds:
            if not os.path.exists(dump_dir+"train_outcome_train80_fold="+str(fold)+"_cog27.csv"):
                load_data = True

        if load_data:
            rand_pairs=np.load(dump_dir+'random_pairs_cog27.npy')
            married_pairs=np.load(dump_dir+'married_pairs_cog27.npy')    
            if train80:
                for fold in folds:
                    print('Collecting data')
                    [X_train,X_test],[outcome_train,outcome_test],feature_cols = collect_data(married_pairs,rand_pairs,community_ibd,married_only,train80,fold)
                    pd.DataFrame(feature_cols,columns=["features"]).to_csv(dump_dir+"features_cog27.csv",index=False)

                    slice_features = [f for f in feature_cols if 'chr' in f]#['chr5_30918721-31744422']#[f for f in feature_cols if 'chr' in f]
                    non_slice_features_pos = [pos for pos,f in enumerate(feature_cols) if f not in slice_features]
                    pd.DataFrame(X_train[:,non_slice_features_pos],columns = feature_cols[non_slice_features_pos]).to_csv(dump_dir+"train_X_non_slice_features_train80_fold="+str(fold)+"_cog27.csv",index=False)
                    pd.DataFrame(X_test[:,non_slice_features_pos],columns = feature_cols[non_slice_features_pos]).to_csv(dump_dir+"test_X_non_slice_features_test20_fold="+str(fold)+"_cog27.csv",index=False)
                    pd.DataFrame(outcome_train,columns = ['marriage']).to_csv(dump_dir+"train_outcome_train80_fold="+str(fold)+"_cog27.csv",index=False)
                    pd.DataFrame(outcome_test,columns = ['marriage']).to_csv(dump_dir+"test_outcome_test20_fold="+str(fold)+"_cog27.csv",index=False)
                    del outcome_train,outcome_test
                    # create a file for each slice to cut down on memory usage
                    for slice in slice_features:
                        slice_pos = np.min(np.nonzero(feature_cols == slice)[0])
                        slice_covar = X_train[:,slice_pos]
                        train_slice = pd.DataFrame(slice_covar,columns=["slice"])
                        train_slice.to_csv(dump_dir+"slices80/train_"+slice+"_train80_fold="+str(fold)+"_cog27.csv",index=False)
                        slice_covar = X_test[:,slice_pos]
                        test_slice = pd.DataFrame(slice_covar,columns=["slice"])
                        test_slice.to_csv(dump_dir+"slices80/test_"+slice+"_test20_fold="+str(fold)+"_cog27.csv",index=False)
                        del slice_pos, slice_covar,test_slice,train_slice
                    del X_train,X_test
                print('DONE Saving data')
            else:
                X,outcome,feature_cols = collect_data(married_pairs,rand_pairs,community_ibd,married_only,train80)
                pd.DataFrame(feature_cols,columns=["features"]).to_csv(dump_dir+"features_cog27.csv",index=False)
                slice_features = [f for f in feature_cols if 'chr' in f]
                non_slice_features_pos = [pos for pos,f in enumerate(feature_cols) if f not in slice_features]
                pd.DataFrame(X[:,non_slice_features_pos],columns = feature_cols[non_slice_features_pos]).to_csv(dump_dir+"train_X_non_slice_features_cog27.csv",index=False)
                pd.DataFrame(outcome_train[:,non_slice_features_pos],columns = ['marriage']).to_csv(dump_dir+"train_outcome_cog27.csv",index=False)

                # create a file for each slice to cut down on memory usage
                for slice in slice_features:
                    slice_pos = np.min(np.nonzero(feature_cols == slice)[0])
                    slice_covar = X[:,slice_pos]
                    pd.DataFrame(slice_covar,columns=["slice"]).to_csv(dump_dir+"slices/train_"+slice+"_cog27.csv",index=False)
        feature_cols = pd.read_csv(dump_dir+"features_cog27.csv").values.flatten()

        directory = '/project/burghard_687/genetic_data/'
        slice_sets = []
        covar_combos  = ["ethnicity_only","pca_distance","edu_only","income","num_marriage","num_div","cog27","relig","height","none","place","place+*cog27_","place+*nopca_"]
        for cov in covar_combos:
            for fold in range(5):
                slices = []
                files = list(glob(directory+'slices_abs_cog27/firth_'+cov+'*prob_*train80_fold='+str(fold)+'_abs_no_cog-bmi_cog27.csv'))# slices_abs_prob
                files += list(glob(directory+'slices_abs_cog27_2/firth_'+cov+'*prob_*train80_fold='+str(fold)+'_abs_no_cog-bmi_cog27.csv'))# slices_abs_prob
                for f in files:
                    chrom_pos = [n for n,c in enumerate(f.split('_')) if 'chr' in c][0]
                    slice_vals = f.split('_')[chrom_pos:chrom_pos+2]
                    slice_vals = '_'.join(slice_vals)
                    slices.append(slice_vals)        
                print([len(slices),len(files)])
                slice_sets.append(set(slices))
        # remove intersection of all folds
        shared_slices = slice_sets[0]
        for slice_set in slice_sets:
            shared_slices = shared_slices.intersection(slice_set)
        print(len(shared_slices))
        all_slices = set([f for f in feature_cols if 'chr' in f])        
        print(len(all_slices))

        missing_slices = all_slices - shared_slices
        slice_features = [[f,train80] for f in missing_slices]
       
        random.shuffle(slice_features)
        input = slice_features
        pool = Pool(num_threads)
        
        pool.map(run_r,slice_features)


def main():

    num_threads = 1
    firth_r_GWAS(num_threads)
if __name__ == "__main__":
    main()
