import os,time,sys,  random
from glob import glob
from itertools import permutations,combinations
#from more_itertools import consume
from multiprocessing import Pool
import warnings
import numpy as np
import pandas as pd
import pickle as pk

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
    # columns
    # FID,IID,sex,EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,EV11,EV12,EV13,EV14,EV15,EV16,EV17,EV18,EV19,EV20
    IDs=df_pca['IID'].values
    IDs = IDs[~np.isnan(IDs)]
    sex = find_sex(df_pca,IDs)
    return sex

def collect_common_slices(married_pairs,rand_pairs,community_ibd,married_only):
    directory= '/project/arpawong_181/burghard_687/genetic_data/ILASH_output/'
    directory= '/project/arpawong_181/HRS_AsMa/keith/phased2/'
    bp_ranges=pd.read_csv(directory+'Slices_Per_Chromosome.csv')
    slice_features = []
    all_rand_slice = np.array([])
    all_married_slice = []
    sex = sex_dict()
    sex_nodes = set(list(sex.keys()))

    married_nodes = set(list(np.unique(np.array(married_pairs).flatten())))
    married_pairs_tuples = set([(p1,p2) for p1,p2 in married_pairs])
    rand_nodes = set(list(np.unique(np.array(rand_pairs).flatten())))
    all_nodes = married_nodes.union(rand_nodes)
    rand_check = set([(p1,p2) for p1,p2 in rand_pairs])
    rand_pos = {(p1,p2):ii for ii,[p1,p2] in enumerate(rand_pairs)}
    s_features = 0
    collect = False
    for chrom in range(1,23):
        rand_file = '/project/arpawong_181/burghard_687/genetic_data/slices/rand_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
        married_file = '/project/arpawong_181/burghard_687/genetic_data/slices/married_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
        if not (os.path.exists(rand_file) and os.path.exists(married_file)):
             collect = False
    for chrom in range(1,23):
        all_slices = sorted(list(set(list(bp_ranges.loc[bp_ranges['chr']==chrom,'slice'].values))))
        # allow everything at the edge to be counted
        all_slices[-1] += 10
        slice_bins = [[b1,b2] for b1,b2 in zip(all_slices[:-1],all_slices[1:]) if b1<b2]
        #print(np.mean([np.abs(b1-b2) for b1,b2 in slice_bins]))
        #print(len(slice_bins))
        print('CHROMOSOME: ',chrom)
        rand_file = '/project/arpawong_181/burghard_687/genetic_data/slices/rand_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
        married_file = '/project/arpawong_181/burghard_687/genetic_data/slices/married_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
        unique_slices = set()
        if os.path.exists(rand_file) and os.path.exists(married_file):
            try:
                rand_slice = np.load(rand_file)
                print('remove: ',rand_file)
                pd.DataFrame({'not loaded':[]}).to_csv('loaded_'+str(chrom)+'.csv')
            except:
                continue
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
            #slice_features+=['chr'+str(chrom)+'_'+str(b1)+'-'+str(b2) for b1,b2 in slice_bins]
            #s_features +=len(keep_pos)
            slice_features+=['chr'+str(chrom)+'_'+str(b1)+'-'+str(b2) for b1,b2 in slice_bins]
            #all_married_slice.append(np.array(married_slice)[:,keep_pos])
            #all_rand_slice.append(np.array(rand_slice)[:,keep_pos])
            #if collect:
            #     if chrom == 1:
            #         all_rand_slice = np.array(rand_slice)[:,keep_pos]
            #     else:
            #         all_rand_slice = np.concatenate([all_rand_slice,np.array(rand_slice)[:,keep_pos]],axis=1)
            continue
        rand_slice = np.zeros((len(rand_pairs),len(slice_bins)),dtype=np.uint8)#
        married_slice = np.zeros((len(married_pairs),len(slice_bins)),dtype=np.uint8)#        
        file = directory+'dist_HRS_'+str(chrom)#'hrs123_qc1all'+str(chrom)
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
                #print(nn)
                if nn % 10 == 0:
                    print(round(nn/len(slice_bins)*100,2))
                if prev_b>np.mean([bp1,bp2]):
                    print('ERROR',[bp1,bp2],prev_b)
                prev_b=np.mean([bp1,bp2])
                # split by data within this range
                # if the far ends of the ilash output are within this bin
                edges = data.loc[(data['bp_ind2']>=bp1)&(data['bp_ind1']<bp2),['family_id1','family_id2']].values
                #print('time-'+str(nn))

                # no self loops
                #np.array([(e1,e2) for e1,e2 in edges if e1 != e2])
                edges = set([tuple(sorted([e1,e2])) for e1,e2 in edges if e1 != e2])
                # only look at edges that may be in rand/married edges
                edges = [(e1,e2) for e1,e2 in edges if e1 in sex_nodes and e2 in sex_nodes and e1 in all_nodes and e2 in all_nodes]
                edges = set([(e1,e2) for e1,e2 in edges if sex[e1] != sex[e2]])
                if len(edges) == 0: continue
                # check if data in ilash edges
                bp_pos = nn#slice_bins.index([bp1,bp2])
                # check if data are within communities = 1, else 0
                for [p1,p2] in edges:
                    #for ii,[p1,p2] in enumerate(rand_pairs):
                    # if random pairs within common slices
                    if not set([(p1,p2)]).isdisjoint(rand_check):
                        ii = rand_pos[(p1,p2)]
                        rand_slice[ii][bp_pos] = 1
                        #print(np.sum([rand_slice[kk][bp_pos] for kk in range(len(rand_slice))]))
                for jj,[p1,p2] in enumerate(married_pairs):
                    # if married pairs within common slices
                    if not set([(p1,p2)]).isdisjoint(edges):
                        married_slice[jj][bp_pos] = 1
            keep_pos = []
            for slice in range(len(rand_slice[0])):
                common_pairs_rand = np.sum(rand_slice[:,slice])
                common_pairs_married = np.sum(married_slice[:,slice])
                # if more than one person was married with this slice, keep
                if common_pairs_rand > 0 and common_pairs_married > 0:
                    keep_pos.append(slice)
            slice_bins = np.array(slice_bins)
            slice_bins = slice_bins[keep_pos]
            s_features +=len(keep_pos)
            slice_features+=['chr'+str(chrom)+'_'+str(b1)+'-'+str(b2) for b1,b2 in slice_bins]
            married_slice = married_slice[:,keep_pos]
            rand_slice = rand_slice[:,keep_pos]
        #rand_slice = rand_slice.astype(np.uint8)
        #married_slice = married_slice.astype(np.uint8)
        #all_married_slice.append(np.array(married_slice))
        #all_rand_slice.append(np.array(rand_slice))
        #if collect:
        #     if chrom == 1:
        #         all_rand_slice = np.array(rand_slice)[:,keep_pos]
        #     else:
        #         all_rand_slice = np.concatenate([all_rand_slice,np.array(rand_slice)[:,keep_pos]],axis=1)
        
        np.save('/project/arpawong_181/burghard_687/genetic_data/slices/rand_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy',rand_slice)
        np.save('/project/arpawong_181/burghard_687/genetic_data/slices/married_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy',married_slice)
    #all_married_slice = np.concatenate(all_married_slice,axis=1)
    #all_rand_slice = np.concatenate(all_rand_slice,axis=1)
    dump_dir = '/project/arpawong_181/burghard_687/genetic_data/'#'/project/burghard_687/genetic_data/'
    np.save(dump_dir+'slice_features_cog27.npy',slice_features)
    return #all_married_slice,all_rand_slice,slice_features

def collect_prs_data(pairs):
    dump_dir = ''#dump_dir = '/project/arpawong_181/burghard_687/genetic_data/'
    PRS_data = pd.read_csv(dump_dir+'pgs12871_4keith.csv')
    diffs=[]
    for p1,p2 in pairs:
        prs1 = PRS_data.loc[PRS_data['FID']==p1,[c for c in PRS_data.columns if 'pgs_' in c]].values
        prs2= PRS_data.loc[PRS_data['FID']==p2,[c for c in PRS_data.columns if 'pgs_' in c]].values
        if len(prs1) > 0 and len(prs2) > 0:
            # absolute difference in polygenic risk score; N.B., not absolute values
            diff = (prs1-prs2)[0]
        else:
            diff = [None]*4
        diffs.append(list(diff)) 
    diffs = np.array(diffs)
    return diffs 

def collect_data(married_pairs,rand_pairs,community_ibd,married_only,train80=False,output_fold=-1):
    dump_dir = '/project/arpawong_181/burghard_687/genetic_data/'#'/project/burghard_687/genetic_data/'
    #X_file = dump_dir+'slice_features_comm='+str(community_ibd)+'_cog27.npy'
    #outcome_file = dump_dir+'slice_outcome_cog27.npy'
    #feature_file = dump_dir+'slice_feature_cols_cog27.npy'
    fold=-1
    X_file = dump_dir+'slice_features_comm='+str(community_ibd)+'_all_data_fold='+str(fold)+'_cog27.npy'
    outcome_file = dump_dir+'slice_outcome_all_data_fold='+str(fold)+'_cog27.npy'
    feature_file = dump_dir+'slice_feature_cols_all_data_fold='+str(fold)+'_cog27.npy'
    X_file_train = dump_dir+'slice_features_comm='+str(community_ibd)+'_all_data_fold='+str(fold)+'_cog27.npy'
    outcome_file_train = dump_dir+'slice_outcome_all_data_fold='+str(fold)+'_cog27.npy'

    if train80:
        fold=4
        #X_file = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(fold)+'.npy'
        #outcome_file = dump_dir+'slice_outcome_train80_fold='+str(fold)+'.npy'
        #feature_file = dump_dir+'slice_feature_cols_train80_fold='+str(fold)+'.npy'
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

    if True:#
        #if not (os.path.exists(X_file) and os.path.exists(outcome_file) and os.path.exists(feature_file)): 
        X_rand=np.load(dump_dir+'X_rand_cog27.npy')
        X_married=np.load(dump_dir+'X_married_cog27.npy')
        #set_rp = set([(p1,p2) for p1,p2 in rand_pairs])
        #set_mp = set([(p1,p2) for p1,p2 in married_pairs])
        feature_cols_orig = np.load(dump_dir+'feature_cols_cog27.npy')        
        # create more features: 1 if links are the same between pairs, 0 otherwise
        if True:#
            #if not os.path.exists(dump_dir+'slice_features_cog27.npy'):
            #all_married_slice,all_rand_slice,slice_features = 
            collect_common_slices(married_pairs,rand_pairs,community_ibd,married_only)
            #PRS 
            #np.save(dump_dir+'slice_features_cog27.npy',slice_features)
            # END IT HERE
            return



            prs_married = np.array([[0,0,0,0]]*len(married_pairs))#collect_prs_data(married_pairs)
            X_married = np.concatenate([X_married,np.array(all_married_slice),np.array(prs_married)],axis=1)
            np.save(dump_dir+'X_married_full_cog27.npy',X_married)

            #all_rand_slice = np.concatenate([np.load('/project/burghard_687/genetic_data/slices/rand_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy') for chrom in range(1,23)],axis=1)
            #X_rand = np.concatenate([X_rand,np.array(all_rand_slice),np.array(prs_rand)],axis=1)
            #prs_rand = np.array([[0,0,0,0]]*len(rand_pairs))#collect_prs_data(rand_pairs)
            #X_rand = np.concatenate([X_rand,np.array(all_rand_slice),np.array(prs_rand)],axis=1)
            #np.save(dump_dir+'X_rand_full_cog27.npy',X_rand)
            #np.save(dump_dir+'X_rand_full_heldout_highmem_0.npy',X_rand[:int(len(X_rand)/2)])
            #np.save(dump_dir+'X_rand_full_heldout_highmem_1.npy',X_rand[int(len(X_rand)/2):])
            total_len = len(X_rand)
            num_bins = 10
            for i in range(num_bins):
                #if False:
                all_rand_slice = []
                for chrom in range(1,23):
                    print(chrom)
                    rand_file = '/project/arpawong_181/burghard_687/genetic_data/slices/rand_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
                    married_file = '/project/arpawong_181/burghard_687/genetic_data/slices/married_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
                    rand_slice = np.load(rand_file)
                    married_slice = np.load(married_file)
                    keep_pos = []
                    for slice in range(len(rand_slice[0])):
                        common_pairs_rand = np.sum(rand_slice[:,slice])
                        common_pairs_married = np.sum(married_slice[:,slice])
                        # if more than one person was married with this slice, keep
                        if common_pairs_rand > 0 and common_pairs_married > 0:
                            keep_pos.append(slice)
                    #rand_chrom = np.load('/project/arpawong_181/burghard_687/genetic_data/slices/rand_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_heldout.npy')
                    rand_slice = rand_slice[i*int(total_len/num_bins):(i+1)*int(total_len/num_bins),keep_pos]
                    all_rand_slice.append(rand_slice)
                    del rand_slice
                    del married_slice
                all_rand_slice = np.concatenate(all_rand_slice,axis=1)
                #np.concatenate([np.load('/project/arpawong_181/burghard_687/genetic_data/slices/rand_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_heldout.npy')[i*int(total_len/num_bins):(i+1)*int(total_len/num_bins)] for chrom in range(1,23)],axis=1)
                X_rand = np.concatenate([X_rand[i*int(total_len/num_bins):(i+1)*int(total_len/num_bins)],np.array(all_rand_slice),np.array(prs_rand)[i*int(total_len/num_bins):(i+1)*int(total_len/num_bins)]],axis=1)
                np.save(dump_dir+'X_rand_full_cog27_highmem_'+str(i)+'.npy',X_rand[:int(len(X_rand)/2)])
            #prs_rand = np.array([[0,0,0,0]]*len(rand_pairs))#collect_prs_data(rand_pairs)
            #X_rand = np.concatenate([X_rand,np.array(all_rand_slice),np.array(prs_rand)],axis=1)
            #np.save(dump_dir+'X_rand_full_cog27.npy',X_rand)

        else:

            print('LOADING LARGEST DATASETS...')
            #X_rand = np.load(dump_dir+'X_rand_full_cog27.npy')
            X_married = np.load(dump_dir+'X_married_full_cog27.npy')
            slice_features = np.load(dump_dir+'slice_features_cog27.npy')
            print('DATASET LOADED...')

        outcome_folds = []
        if train80:
            n_folds = 5
            folds = list(range(n_folds))
            train_m_len = int((n_folds-1)/(n_folds)*len(X_married))
            fold_m_len = int(1/(n_folds)*len(X_married)) 
            ##train_r_len = int((n_folds-1)/(n_folds)*len(X_rand))
            # find all pairs who are strictly within married fold
            #X_folds = []
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

            
            #for fold in folds:
            #    X_file = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(fold)+'.npy'
            #    outcome_file = dump_dir+'slice_outcome_train80_fold='+str(fold)+'.npy'
            #    feature_file = dump_dir+'slice_feature_cols_train80_fold='+str(fold)+'.npy'
            #    X_file_train = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(fold)+'.npy'
            #    X_file_test = dump_dir+'slice_features_comm='+str(community_ibd)+'_test20_fold='+str(fold)+'.npy'
            #    outcome_file_test = dump_dir+'slice_outcome_test20_fold='+str(fold)+'.npy'
            #    outcome_file_train = dump_dir+'slice_outcome_train80_fold='+str(fold)+'.npy'
            #    
            #    X_test = X_folds[fold]
            #    outcome_test = outcome_folds[fold]
            #    train_folds = [vf for vf in folds if vf != fold]
            #    X_train = np.concatenate([X_folds[vf] for vf in folds if vf != fold],axis=0)
            #    outcome_train = np.concatenate([outcome_folds[vf] for vf in folds if vf != fold],axis=0).flatten()
            #    feature_cols = np.array(list(np.array(feature_cols_orig).astype(str)) + list(np.array(slice_features).astype(str))+['PRS_dummy1','PRS_dummy2','PRS_dummy3','PRS_dummy4'])
            #    good_data_pos = [pos for pos,line in enumerate(X_train) if None not in line]
            #    X_train = X_train[good_data_pos]
            #    outcome_train = outcome_train[good_data_pos]
            #    
            #    good_data_pos = [pos for pos,line in enumerate(X_test) if None not in line]
            #    X_test = X_test[good_data_pos]
            #    outcome_test = outcome_test[good_data_pos]
            #    
            #    np.save(X_file_train,X_train)
            #    np.save(outcome_file_train,outcome_train) 
            #    np.save(X_file_test,X_test)
            #    np.save(outcome_file_test,outcome_test) 
            #    np.save(feature_file,feature_cols)

        else:
            n_folds = -1
            folds = list(range(n_folds))
            if n_folds == -1:
                folds = [-1]
            train_m_len = len(X_married)
            fold_m_len = len(X_married)
            for vf in folds:
                if True:

                    X_file = dump_dir+'slice_features_comm='+str(community_ibd)+'_all_data_fold='+str(vf)+'_cog27.npy'
                    outcome_file = dump_dir+'slice_outcome_all_data_fold='+str(vf)+'_cog27.npy'
                    feature_file = dump_dir+'slice_feature_cols_all_data_fold='+str(vf)+'_cog27.npy'
                    X_file_train = dump_dir+'slice_features_comm='+str(community_ibd)+'_all_data_fold='+str(vf)+'_cog27.npy'
                    outcome_file_train = dump_dir+'slice_outcome_all_data_fold='+str(vf)+'_cog27.npy'


                    married_fold = X_married[:]#random_indices]
                    married_subjects_fold = set(list(np.array(married_pairs[:]).flatten()))
                    good_data_pos = [pos for pos,line in enumerate(married_fold) if None not in line]
                    married_fold = married_fold[good_data_pos]


                    # train data: look at all links disconnected from test
                    # all married pairs not in the fold
                    married_train = married_fold[:]
                    # all random pairs such that neither male nor female is in test set (no overlapping edges)
                    rand_subjects_train = [pos for pos,[p1,p2] in enumerate(rand_pairs)]#  if p1 in married_train and p2 in married_train]#if pos not in set_rand_subjects_test]
                    print('rand_train')
                    for i in range(max_iter):
                        outcome_file_train_i = outcome_file_train.replace('.npy','_'+str(i)+'.npy')
                        X_file_train_i = X_file_train.replace('.npy','_'+str(i)+'.npy')
                        rand_train = []
                        for i in range(num_bins):
                            X_rand = np.load(dump_dir+'X_rand_full_cog27_highmem_'+str(i)+'.npy')
                            # map to correct positions
                            X_rand_pos = {int(total_len/num_bins)*i+j:j for j in range(len(X_rand))}
                            valid_pos = set(list(X_rand_pos.keys()))
                            train_pos = [X_rand_pos[pos] for pos in rand_subjects_train[i*330000:(i+1)*330000] if pos in valid_pos]
                            if len(train_pos)>0:
                                rand_train.append(X_rand[train_pos])
                        rand_train = np.concatenate(rand_train,axis=0)#X_rand[rand_subjects_train[i*330000:(i+1)*330000]]
                        print('concatenate')
                        if i == 0:
                            train_m_len = len(married_train)
                            train_r_len = len(rand_train)
                            X_train = np.concatenate([married_train,rand_train],axis=0)
                            outcome_train = np.array([1]*train_m_len+[0]*train_r_len)
                        elif i != max_iter-1:
                            train_r_len = len(rand_train)
                            outcome_train = np.array([0]*train_r_len)
                            X_train = rand_train[:]
                        else:
                            rand_train = X_rand[rand_subjects_train[i*330000:]]
                            train_r_len = len(rand_train)
                            outcome_train = np.array([0]*train_r_len)
                            X_train = rand_train[:]
                        np.save(outcome_file_train_i,outcome_train)
                        np.save(X_file_train_i,X_train)
                        del rand_train
                        del X_train

                    feature_cols = np.array(list(np.array(feature_cols_orig).astype(str)) + list(np.array(slice_features).astype(str))+['PRS_dummy1','PRS_dummy2','PRS_dummy3','PRS_dummy4'])
                    np.save(feature_file,feature_cols)
                    print('SAVING', outcome_file_train_i, ' ',X_file_train_i, ' ',feature_file)

            #n_folds = -1
            #folds = list(range(n_folds))
            #if n_folds == -1:
            #    folds = [-1]
            #train_m_len = len(X_married)
            #fold_m_len = len(X_married)
            #for vf in folds:
            #    if not os.path.exists(dump_dir+'slice_features_comm=False_all_data_fold='+str(vf)+'_cog27.npy'):
            #        #'X_fold='+str(vf)+'.npy'):
            #        X_file = dump_dir+'slice_features_comm='+str(community_ibd)+'_all_data_fold='+str(vf)+'_cog27.npy'
            #        outcome_file = dump_dir+'slice_outcome_all_data_fold='+str(vf)+'_cog27.npy'
            #        feature_file = dump_dir+'slice_feature_cols_all_data_fold='+str(vf)+'_cog27.npy'
            #        X_file_train = dump_dir+'slice_features_comm='+str(community_ibd)+'_all_data_fold='+str(vf)+'_cog27.npy'
            #        outcome_file_train = dump_dir+'slice_outcome_all_data_fold='+str(vf)+'_cog27.npy'
            #
            #        married_fold = X_married[:]#random_indices]
            #        married_subjects_fold = set(list(np.array(married_pairs[:]).flatten()))
            #        good_data_pos = [pos for pos,line in enumerate(married_fold) if None not in line]
            #        married_fold = married_fold[good_data_pos]
            #
            #        # train data: look at all links disconnected from test
            #        # all married pairs not in the fold
            #        married_train = married_fold[:]
            #        # all random pairs such that neither male nor female is in test set (no overlapping edges)
            #        rand_subjects_train = [pos for pos,[p1,p2] in enumerate(rand_pairs)]#  if p1 in married_train and p2 in married_train]#if pos not in set_rand_subjects_test]
            #        print('rand_train')
            #        rand_train = X_rand[rand_subjects_train[:330000]]
            #        print('concatenate')
            #        X_train = np.concatenate([married_train,rand_train],axis=0)
            #        train_m_len = len(married_train)
            #        train_r_len = len(rand_train)
            #        outcome_train = np.array([1]*train_m_len+[0]*train_r_len)
            #        np.save(outcome_file_train,outcome_train) 
            #       
            #        np.save(X_file_train,X_train)
            #        del rand_train
            #        del X_train
            #        print('rand_train 2')
            #        rand_train = X_rand[rand_subjects_train[330000:660000]]
            #        train_r_len = len(rand_train)
            #        outcome_train = np.array([0]*train_r_len)
            #        np.save(outcome_file_train[:-4]+'_2_cog27.npy',outcome_train) 
            #        np.save(X_file_train[:-4]+'_2_cog27.npy',rand_train)
            #        del	rand_train
            #        #del X_train
            #        print('rand_train 3')
            #        rand_train = X_rand[rand_subjects_train[660000:]]
            #        train_r_len = len(rand_train)
            #        outcome_train = np.array([0]*train_r_len)
            #        np.save(outcome_file_train[:-4]+'_3_cog27.npy',outcome_train) 
            #        np.save(X_file_train[:-4]+'_3_cog27.npy',rand_train)
            #        del	rand_train
            #       
            #        feature_cols = np.array(list(np.array(feature_cols_orig).astype(str)) + list(np.array(slice_features).astype(str))+['PRS_dummy1','PRS_dummy2','PRS_dummy3','PRS_dummy4'])
            #        np.save(feature_file,feature_cols)




            #X = np.concatenate([X_married,X_rand],axis=0)
            #outcome = np.array([1]*len(X_married)+[0]*len(X_rand))
            #feature_cols = np.array(list(feature_cols) + slice_features+['PRS_dummy1','PRS_dummy2','PRS_dummy3','PRS_dummy4'])
            #good_data_pos = [pos for pos,line in enumerate(X) if None not in line]
            #X = X[good_data_pos]
            #outcome = outcome[good_data_pos]       
            #np.save(X_file,X)
            #np.save(outcome_file,outcome) 
            #np.save(feature_file,feature_cols)

    if train80:
        outcome_file = dump_dir+'slice_outcome_train80_fold='+str(output_fold)+'_cog27.npy'
        feature_file = dump_dir+'slice_feature_cols_train80_fold='+str(output_fold)+'_cog27.npy'
        X_file_train = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(output_fold)+'_cog27.npy'
        X_file_test = dump_dir+'slice_features_comm='+str(community_ibd)+'_test20_fold='+str(output_fold)+'_cog27.npy'
        outcome_file_test = dump_dir+'slice_outcome_test20_fold='+str(output_fold)+'_cog27.npy'
        outcome_file_train = dump_dir+'slice_outcome_train80_fold='+str(output_fold)+'_cog27.npy'
        #del X_rand
        X_test = np.load(X_file_test,allow_pickle = True)
        X_train = np.concatenate([np.load(X_file_train,allow_pickle = True),np.load(X_file_train[:-4]+'_2_cog27.npy',allow_pickle = True),np.load(X_file_train[:-4]+'_3_cog27.npy',allow_pickle = True)],axis=0)
        outcome_test = np.load(outcome_file_test,allow_pickle = True)
        outcome_train = np.concatenate([np.load(outcome_file_train,allow_pickle = True),np.load(outcome_file_train[:-4]+'_2_cog27.npy',allow_pickle = True),np.load(outcome_file_train[:-4]+'_3_cog27.npy',allow_pickle = True)],axis=0)
        feature_cols = np.load(feature_file)#,allow_pickle = True)
        return [X_train,X_test],[outcome_train,outcome_test],feature_cols
    else:
        outcome_file = dump_dir+'slice_outcome_all_data_fold='+str(output_fold)+'_cog27.npy'
        feature_file = dump_dir+'slice_feature_cols_all_data_fold='+str(output_fold)+'_cog27.npy'
        X_file_train = dump_dir+'slice_features_comm='+str(community_ibd)+'_all_data_fold='+str(output_fold)+'_cog27.npy'
        outcome_file_train = dump_dir+'slice_outcome_all_data_fold='+str(output_fold)+'_cog27.npy'
        X_train = np.concatenate([np.load(X_file_train.replace('.npy','_'+str(i)+'.npy'),allow_pickle = True) for i in range(max_iter)],axis=0)
        outcome_train = np.concatenate([np.load(outcome_file_train.replace('.npy','_'+str(i)+'.npy'),allow_pickle = True) for i in range(max_iter)],axis=0)

        feature_cols = np.load(feature_file)#,allow_pickle = True)
        random_indices = list(range(len(outcome_train)))
        if output_fold != -1:
            random_indices = np.random.choice(random_indices,len(random_indices),replace=True)
        return [X_train[random_indices],[]],[outcome_train[random_indices],[]],feature_cols

        #outcome_file = dump_dir+'slice_outcome_all_data_fold='+str(output_fold)+'_cog27.npy'
        #feature_file = dump_dir+'slice_feature_cols_all_data_fold='+str(output_fold)+'_cog27.npy'
        #X_file_train = dump_dir+'slice_features_comm='+str(community_ibd)+'_all_data_fold='+str(output_fold)+'_cog27.npy'
        #outcome_file_train = dump_dir+'slice_outcome_all_data_fold='+str(output_fold)+'_cog27.npy'
        #del X_rand
        #X_train = np.concatenate([np.load(X_file_train),np.load(X_file_train[:-4]+'_2_cog27.npy'),np.load(X_file_train[:-4]+'_3_cog27.npy')],axis=0)
        #outcome_train = np.concatenate([np.load(outcome_file_train,allow_pickle = True),np.load(outcome_file_train[:-4]+'_2_cog27.npy',allow_pickle = True),np.load(outcome_file_train[:-4]+'_3_cog27.npy',allow_pickle = True)],axis=0)
        #feature_cols = np.load(feature_file)#,allow_pickle = True)
        #random_indices = list(range(len(outcome_train)))
        #if output_fold != -1:
        #    random_indices = np.random.choice(random_indices,len(random_indices),replace=True)
        #return [X_train[random_indices],[]],[outcome_train[random_indices],[]],feature_cols

        #X_file = dump_dir+'slice_features_comm='+str(community_ibd)+'_train80_fold='+str(output_fold)+'_cog27.npy'
        #X = np.load(X_file,allow_pickle = True)
        #outcome = np.load(outcome_file,allow_pickle = True)
        #feature_cols = np.load(feature_file,allow_pickle = True)
        #return X,outcome,feature_cols

def shuffle_2arrays(a, b):
    rng_state = np.random.get_state()
    np.random.shuffle(a)
    np.random.set_state(rng_state)
    np.random.shuffle(b)
def find_pvals(community_ibd,firth_fit,last_slice):
    file = 'ModelCoefs/ibd_pvals_firth='+str(firth_fit)+'.pkl'
    dump_dir = ''#'/project/arpawong_181/burghard_687/genetic_data/'
    if not os.path.exists(file):
        print('LOOKING FOR PVALS')	
        married_only = True
        rand_pairs=np.load(dump_dir+'random_pairs_cog27.npy')
        married_pairs=np.load(dump_dir+'married_pairs_cog27.npy')
        #X,outcome,feature_cols = 
        collect_data(married_pairs,rand_pairs,community_ibd,married_only)
        feature_cols = np.load(dump_dir+'slice_features_cog27.npy')
        non_slice_cols = feature_cols[list(feature_cols).index('pca_distance')]
        X_rand=np.load(dump_dir+'X_rand_cog27.npy')
        X_married=np.load(dump_dir+'X_married_cog27.npy')
        non_slice_data = np.concatenate([X_married,X_rand],axis=0)
        outcome = [1]*len(X_married)+[0]*len(X_rand)
        #X = non_slice_data
        # remove PRS dummy variables
        #X = X[:,:len(feature_cols)]
        print('TRAIN DATA SAVED')
        feature_cols_pos = {c:ii for ii, c in enumerate(feature_cols)}
        pos_class = np.array([pos for pos,o in enumerate(outcome) if o == 1])
        neg_class = np.array([pos for pos,o in enumerate(outcome) if o == 0])
        #slice_features = list(pd.read_csv('select_slices.csv')['slice'].values)#[f for f in feature_cols if 'chr' in f]
        slice_features = [f for f in feature_cols if 'chr' in f]
        bp_ranges=pd.read_csv(directory+'Slices_Per_Chromosome.csv')
        for chrom in range(1,23):
             all_slices = sorted(list(set(list(bp_ranges.loc[bp_ranges['chr']==chrom,'slice'].values))))
             all_slices[-1] += 10
             slice_bins = [[b1,b2] for b1,b2 in zip(all_slices[:-1],all_slices[1:]) if b1<b2]
             rand_file = '/project/arpawong_181/burghard_687/genetic_data/slices/rand_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
             married_file = '/project/arpawong_181/burghard_687/genetic_data/slices/married_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
             #married_slices = np.load(married_file)
             #rand_slices = np.load(rand_file)
             slices_features = np.concatenate([np.load(married_file),np.load(rand_file)],axis=0) 
        #if last_slice is not None:
        #    slice_len = 100
        #    if (last_slice+1)*slice_len > len(slice_features):
        #        slice_features = slice_features[last_slice*slice_len:]
        #    else:
        #        slice_features = slice_features[last_slice*slice_len:(last_slice+1)*slice_len]
        #print(len(slice_features))
        #covar_features = []
        #covar_features += ['RABPLACE_1.0','RABPLACE_2.0','RABPLACE_3.0','RABPLACE_4.0','RABPLACE_5.0','RABPLACE_6.0','RABPLACE_7.0','RABPLACE_8.0','RABPLACE_9.0','RABPLACE_10.0','RABPLACE_11.0']
        #covar_features += ['ethnicity_1.0','ethnicity_2.0','ethnicity_3.0','ethnicity_4.0','pca_distance']
        #covar_features += ['RAEDYRS','RARELIG_1.0','RARELIG_2.0','RARELIG_3.0','RARELIG_4.0','RARELIG_5.0','iearn','mdiv','mrct','cog27']
        #covar_features += ['bmi']
        #covar_features += ['height']
        #coef_data = {}
        coef_data['covars'] = covar_features
        coef_data['slices']={}
        old_params = [0]*(len(covar_features)+2)
        if last_slice is not None:
            file = 'ModelCoefs/ibd_pvals_firth='+str(firth_fit)+'_'+slice_features[-1]+'.pkl'

        for ii,slice in enumerate(slice_features):
            print([ii,slice])
            if ii % 10 == 0:
                print(round(ii/len(slice_features)*100,2))
            features = covar_features +  [slice]
            features_pos = [feature_cols_pos[f] for f in features]
            X2 = X[:,features_pos]
            
            np.random.shuffle(pos_class)
            np.random.shuffle(neg_class)
            # train: undersample majority class
            num_train0 = len(neg_class)#int(fract_sampled*len(neg_class))
            if firth_fit: num_train0 = len(pos_class) * 10
            num_train1 = len(pos_class)#int(fract_sampled*len(pos_class))
            
            # equal number of 0,1 class data
            train_X = np.concatenate([X2[pos_class[:num_train1]],X2[neg_class[:num_train0]]],axis=0)                
            train_outcome = np.concatenate([outcome[pos_class[:num_train1]],outcome[neg_class[:num_train0]]],axis=0)
            # positive class << negative class, so there is no overlap in data between training and testing
            #shuffle_2arrays(train_X,train_outcome)
            if np.sum(train_outcome) == len(train_outcome): continue
            print('saved data!')
            bad_data = False
            for feature in zip(train_X.T):
                if np.std(feature) == 0: bad_data = True
            if bad_data: continue
            try:
                exog = sm.add_constant(train_X)
                with warnings.catch_warnings():
                    if firth_fit:
                        coef,waldp,lrtp=firth.fit_pvalues(train_X,train_outcome)
                        pvals = waldp
                        err = [[0,0]]
                    else:
                        warnings.filterwarnings("error")
                        reg = Logit(train_outcome, exog).fit(maxiter=200)#,start_params=[0]+[0]*len(covar_features)+[2])#,full_output=False)#disp=False,full_output=False)#LogisticRegression(penalty='none',solver='newton-cg',max_iter=200)
                        err = reg.cov_params()
                        coef = reg.params#.values
                        slice_coef = coef[-1]
                        slice_err = err[-1,-1]
                        z= np.abs(slice_coef/slice_err)
                        pvals=[2*(1-norm.cdf(z))]
            except:
                print(ii)
                coef= [None]
                pvals = [None]
                err = [[0,0]]
            if coef[-1] is not None:
                # standard error = sqrt(variance) for slice
                logp = -np.log10(pvals)
            else:
                logp=-999

            coef_data['slices'][slice] = {}
            coef_data['slices'][slice]['coef'] = coef
            coef_data['slices'][slice]['err'] = err
            coef_data['slices'][slice]['logp'] = logp
        pk.dump(coef_data,open(file,'wb'))
    else:
        coef_data = pk.load(open(file,'rb'))    
    return coef_data
def run_r(slice_features):
    #dump_dir = ''
    dump_dir = '/project/arpawong_181/burghard_687/genetic_data/'
    slice,train80,one_ethnic = slice_features
    print('SLICE: ',slice)
    #slice = 'chr1_160623094-161442861'
    # make R command line argument
    n_folds = 5

    sample = 'full'
    folder = '_abs'
    one_eth = ''
    u_one_eth = ''

    if one_ethnic:
        one_eth = 'one_ethnic'
        u_one_eth = '_one_ethnic'
    for fold in [0]:#range(1,10):#range(n_folds):
        print('CHECK FILE NAMES ARE CONSISTENT IN PY AND R CODE FOR COEFS, PVALUES, ETC')
        print(fold)
        #if not os.path.exists(dump_dir+"slices80/test_"+slice+"_test20_fold="+str(fold)+"_cog27.csv"): continue
        #if not os.path.exists(dump_dir+"slices80/train_"+slice+"_train80_fold="+str(fold)+"_cog27.csv"): continue
        #if np.sum(pd.read_csv(dump_dir+"slices80/train_"+slice+"_train80_fold="+str(fold)+"_cog27.csv")["slice"].values) == 0: continue
        directory = '/project/arpawong_181/burghard_687/genetic_data/'
        slice_sets = []
        all_covars = ["place","ethnicity","relig","edu","height","cog27"]
        covar_combos  = ["ethnicity_only","pca_distance","edu_only","income","num_marriage","num_div","cog27","relig","height","none","place",'+'.join(all_covars),'+'.join(all_covars)+"+nopca"]
        all_files_found = True
        if fold > 0:
            covar_combos = ["none",'+'.join(all_covars)]
        #for cov in covar_combos:   
        #    #prob_file = dump_dir+"slices"+folder+u_one_eth+"_cog27_"+str(fold)+"/firth_"+cov+"_prob_"+slice+"_"+sample+"_train80_fold="+str(fold)+folder+"_one_ethnic="+str(one_ethnic)+"_no_cog-bmi_cog27.csv"
        #    #if not train80:
        #    #    prob_file = dump_dir+'slices_all'+folder+u_one_eth+'_cog27_'+str(fold)+'/firth_'+cov+'_prob_'+slice+'_'+sample+folder+'_no_cog-bmi_cog27_fold='+str(fold)+'.csv'
        #    #if not os.path.exists(prob_file): 
        #    #    all_files_found = False
        #    #    print(prob_file)
        #    #    break
        if True:#not all_files_found:
            print('fold: ',fold)    
            command = "Rscript firth_fit_allslices_test.r --slice "+slice +" --train80 "+str(train80) + " --fold "+str(fold)+" --one_ethnic "+str(one_ethnic)
            print(command)
            # run command
            try:
                os.system(command)    
            except:
                return
        # ignore folds if we train on all data
        #if not train80: return

def firth_r_GWAS(num_threads,train80,one_ethnic=False):
        dump_dir = '/project/arpawong_181/burghard_687/genetic_data/'
        married_only = True
        community_ibd=False
        load_data = False
        #if not os.path.exists(dump_dir + "train_X_non_slice_features_train80_fold=0_cog27.csv") or not os.path.exists(dump_dir+"train_outcome_train80_fold=0_cog27.csv"):
        #    load_data = True
        n_folds = 5
        #folds = list(range(n_folds))[::-1]
        #for fold in folds:
        #    if not os.path.exists(dump_dir+"train_outcome_train80_fold="+str(fold)+"_cog27.csv"):
        #        load_data = True
        if not train80:
            fold=-1
            if not os.path.exists(dump_dir+"train_outcome_all_data_fold="+str(fold)+"_cog27.csv"):
                load_data = True
        print('LOAD DATA = ',load_data)
        #if not os.path.exists(dump_dir+"features_cog27.csv"):
        #    rand_pairs=np.load(dump_dir+'random_pairs_cog27.npy')
        #    married_pairs=np.load(dump_dir+'married_pairs_cog27.npy')    
        #    fold = 0
        #    [X_train,X_test],[outcome_train,outcome_test],feature_cols = collect_data(married_pairs,rand_pairs,community_ibd,married_only,train80,fold)
        #    pd.DataFrame(feature_cols,columns=["features"]).to_csv(dump_dir+"features_cog27.csv",index=False)

        if load_data:
            rand_pairs=np.load(dump_dir+'random_pairs_cog27.npy')
            married_pairs=np.load(dump_dir+'married_pairs_cog27.npy')    
            if train80:
                for fold in folds:
                    print('Collecting data')
                    #[X_train,X_test],[outcome_train,outcome_test],feature_cols = 
                    if not os.path.exists(dump_dir+"features_cog27.csv"):
                        collect_data(married_pairs,rand_pairs,community_ibd,married_only,train80,fold)
                        #collect_data(married_pairs,rand_pairs,community_ibd,married_only)
                    feature_cols = np.load(dump_dir+'slice_features_cog27.npy')
                    non_slice_cols = feature_cols[list(feature_cols).index('pca_distance')]
                    X_rand=np.load(dump_dir+'X_rand_cog27.npy')
                    X_married=np.load(dump_dir+'X_married_cog27.npy')
                    non_slice_data = np.concatenate([X_married,X_rand],axis=0)
                    outcome = [1]*len(X_married)+[0]*len(X_rand)
                    X = non_slice_data
                    pd.DataFrame(feature_cols,columns=["features"]).to_csv(dump_dir+"features_cog27.csv",index=False)
                    slice_features = [f for f in feature_cols if 'chr' in f]#['chr5_30918721-31744422']#[f for f in feature_cols if 'chr' in f]
                    #if not os.path.exists(dump_dir+"test_outcome_test20_fold="+str(fold)+".csv"):
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
                for fold in [-1]:
                    #[X_train,X_test],[outcome_train,outcome_test],feature_cols = 
                    collect_data(married_pairs,rand_pairs,community_ibd,married_only,train80,fold)
                    feature_cols = np.load(dump_dir+'feature_cols_cog27.npy')
                    # ignore slice data and dummy variables
                    non_slice_cols = feature_cols[:list(feature_cols).index('pca_distance')+1]
                    X_rand=np.load(dump_dir+'X_rand_cog27.npy')
                    X_married=np.load(dump_dir+'X_married_cog27.npy')
                    non_slice_data = np.concatenate([X_married,X_rand],axis=0)
                    outcome = [1]*len(X_married)+[0]*len(X_rand)
                    #X = non_slice_data
                    slice_features = [f for f in feature_cols if 'chr' in f]
                    #if not os.path.exists(dump_dir+"test_outcome_test20_fold="+str(fold)+".csv"):
                    non_slice_features_pos = [pos for pos,f in enumerate(feature_cols) if f not in slice_features]   
                    pd.DataFrame(feature_cols,columns=["features"]).to_csv(dump_dir+"features_cog27.csv",index=False)
                    pd.DataFrame(non_slice_data,columns=non_slice_cols).to_csv(dump_dir+"train_X_non_slice_features_all_data_fold="+str(fold)+"_cog27.csv",index=False)
                    #pd.DataFrame(X_train[:,non_slice_features_pos],columns = feature_cols[non_slice_features_pos]).to_csv(dump_dir+"train_X_non_slice_features_all_data_fold="+str(fold)+"_cog27.csv",index=False)
                    pd.DataFrame(outcome,columns = ['marriage']).to_csv(dump_dir+"train_outcome_all_data_fold="+str(fold)+"_cog27.csv",index=False)
                    #slice_features = [f for f in feature_cols if 'chr' in f]
                    #bp_ranges=pd.read_csv(directory+'Slices_Per_Chromosome.csv')
                    bp_ranges = pd.read_csv('/project/arpawong_181/HRS_AsMa/keith/phased2/Slices_Per_Chromosome.csv')
                    bp_ranges_heldout=pd.read_csv('/project/arpawong_181/HRS_AsMa/keith/ELSA/Slices_Per_Chromosome_heldout.csv')
                    del non_slice_data, X_rand,X_married
                    slice_features = np.load(dump_dir+'slice_features_cog27.npy')
                    for chrom in range(1,23):
                        slice_bins = slice_features[['chr'+str(chrom)+'_' in f for f in slice_features]]
                        rand_file = '/project/arpawong_181/burghard_687/genetic_data/slices/rand_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
                        married_file = '/project/arpawong_181/burghard_687/genetic_data/slices/married_slices_chr='+str(chrom)+'_comm='+str(community_ibd)+'_cog27.npy'
                        slices_rand_features = np.load(rand_file)
                        slices_married_features = np.load(married_file)
                        #np.concatenate([np.load(married_file),np.load(rand_file)],axis=0)
                        for ii,slice in enumerate(slice_bins):
                            #slice = 'chr'+str(chrom)+'_'+str(b1)+'-'+str(b2)
                            # create a file for each slice to cut down on memory usage
                            #for slice in slice_features:
                            slice_file = dump_dir+"slices_all/train_"+slice+"_all_data_fold="+str(fold)+"_cog27.csv"
                            # if file exists AND it was made in September 2024, ignore
                            #if os.path.exists(slice_file) and os.path.getmtime(slice_file)>1725148800: continue
                            #slice_pos = np.min(np.nonzero(feature_cols == slice)[0])
                            slice_covar = np.concatenate([slices_married_features[:,ii],slices_rand_features[:,ii]],axis=0)#slices_features[:,ii]
                            train_slice = pd.DataFrame(slice_covar,columns=["slice"])
                            train_slice.to_csv(slice_file,index=False)
                        #del slice_pos, slice_covar,train_slice
                    #del X_train



                #X,outcome,feature_cols = collect_data(married_pairs,rand_pairs,community_ibd,married_only,train80)
                #[X_train,X_test],[outcome_train,outcome_test],feature_cols = collect_data(married_pairs,rand_pairs,community_ibd,married_only,train80)
                #pd.DataFrame(feature_cols,columns=["features"]).to_csv(dump_dir+"features_cog27.csv",index=False)
                #slice_features = [f for f in feature_cols if 'chr' in f]
                #non_slice_features_pos = [pos for pos,f in enumerate(feature_cols) if f not in slice_features]
                #pd.DataFrame(X[:,non_slice_features_pos],columns = feature_cols[non_slice_features_pos]).to_csv(dump_dir+"train_X_non_slice_features_cog27.csv",index=False)
                #pd.DataFrame(outcome_train[:,non_slice_features_pos],columns = ['marriage']).to_csv(dump_dir+"train_outcome_cog27.csv",index=False)

                # create a file for each slice to cut down on memory usage
                #for slice in slice_features:
                #    slice_pos = np.min(np.nonzero(feature_cols == slice)[0])
                #    slice_covar = X[:,slice_pos]
                #    pd.DataFrame(slice_covar,columns=["slice"]).to_csv(dump_dir+"slices/train_"+slice+"_cog27.csv",index=False)
        feature_cols = pd.read_csv(dump_dir+"features_cog27.csv").values.flatten()
        #slice_features = [[f,train80] for f in feature_cols if 'chr' in f]
        #slice_features = [[f,train80] for f in ['chr20_46656919-47251687','chr2_10522103-11152503','chr3_103963497-105933211','chr9_115960547-116299953','chr2_88178451-96561651','chr9_80049060-80896890','chr14_32928634-33362930','chr11_61685950-63778845','chr10_126787435-127315437','chr7_35771187-36527098','chr2_27623768-29785417','chr12_131022595-131501738','chr1_217450293-218893464','chr2_165834757-168073830','chr5_116181428-117017440','chr15_35747745-36647375','chr19_17641391-17877896','chr5_145559650-146513049','chr1_236666575-237283740','chr8_135678208-136390657','chr20_49580834-49930425','chr21_42320428-42748780','chr6_82453501-85460711','chr21_16724803-16731472','chr6_41415737-41677149','chr18_6471939-6790732','chr12_116673596-117684780','chr22_27763167-28023875','chr15_70024982-70365762','chr20_416954-691433']]
        


        directory = '/project/arpawong_181/burghard_687/genetic_data/'
        slice_sets = []
        covar_combos  = ["ethnicity_only","pca_distance","edu_only","income","num_marriage","num_div","cog27","relig","height","none","place","place+*cog27_","place+*nopca_"]
        #for cov in covar_combos:
        #    files = list(glob(directory+'slices_abs_one_ethnic_cog27/firth_'+cov+'*prob_*train80_fold='+str(fold)+'_one_ethnic='+str(one_ethnic)+'_abs_no_cog-bmi_cog27.csv'))# slices_abs_prob
        #    for fold in range(5):
        #        slices = [] # firth_edu_only_coefs_chr22_27763167-28023875_full_train80_fold=3_abs_no_cog-bmi_cog27.csv
        #        #files += list(glob(directory+'slices_abs_one_ethnic_cog27_2/firth_'+cov+'*prob_*train80_fold='+str(fold)+'_one_ethnic='+str(one_ethnic)+'_abs_no_cog-bmi_cog27.csv'))# slices_abs_prob
        #        files2 = [f for f in files if 'fold='+str(fold) in f]
        #        for f in files2:
        #            chrom_pos = [n for n,c in enumerate(f.split('_')) if 'chr' in c][0]
        #            slice_vals = f.split('_')[chrom_pos:chrom_pos+2]
        #            slice_vals = '_'.join(slice_vals)
        #            slices.append(slice_vals)        
        #        slice_sets.append(set(slices))
        # remove intersection of all folds
        #shared_slices = slice_sets[0]
        #for slice_set in slice_sets:
        #    shared_slices = shared_slices.intersection(slice_set)
        #all_slices = set([f for f in feature_cols if 'chr' in f])
        slice_feature_cols = list(np.load(dump_dir+'slice_features_cog27.npy'))        
        all_slices = set(slice_feature_cols)
        if not os.path.exists(dump_dir+'slice_features_cog27.csv'):
            pd.DataFrame(slice_feature_cols,columns=["slice"]).to_csv(dump_dir+"slice_features_cog27.csv",index=False)            
        missing_slices = all_slices# - shared_slices
        #print([len(all_slices),len(missing_slices)])
        slice_features = [[f,train80,one_ethnic] for f in missing_slices]
       
        #input = slice_features
        #pool = Pool(num_threads)
        #command = "Rscript firth_fit.r --slice ... --train80 "+str(train80) + " --fold "+str(0)+" --one_ethnic "+str(one_ethnic)
        #print(command)
        print(slice_features[:10])
        #pool.map(
        # /project/arpawong_181/burghard_687/genetic_data/slices_combined/
        slice_features = [['chr'+(sf.split('/')[-1]).replace('hrs_slice_','').replace('.csv',''),train80,one_ethnic]  for sf in list(glob(dump_dir+'slices_combined/hrs_slice_*'))]
        random.shuffle(slice_features)

        for sf in slice_features:
            run_r(sf)

def parallel_GWAS(input):
        dump_dir = ''#dump_dir = '/project/arpawong_181/burghard_687/genetic_data/'
        num_boot,fract_sampled,community_ibd,logit_fit,covar,firth_fit,last_slice = input
        print(covar)
        married_only = True
        rand_pairs=np.load(dump_dir+'random_pairs_cog27.npy')
        married_pairs=np.load(dump_dir+'married_pairs_cog27.npy')
        X,outcome,feature_cols = collect_data(married_pairs,rand_pairs,community_ibd,married_only)
        # remove PRS dummy variables
        X = X[:,:len(feature_cols)]

        feature_cols_pos = {c:ii for ii, c in enumerate(feature_cols)}
        pos_class = np.array([pos for pos,o in enumerate(outcome) if o == 1])
        neg_class = np.array([pos for pos,o in enumerate(outcome) if o == 0])
        coef_pvals = find_pvals(community_ibd,firth_fit,last_slice)
        slice_features = [f for f in feature_cols if 'chr' in f]
        if num_boot > 1:
            slice_features = [f for f in slice_features if coef_pvals['slices'][f]['logp'] > 0] # if p-value < 10^-2.5
        covar_features = []
        if 'place' in covar:
            covar_features += ['RABPLACE_1.0','RABPLACE_2.0','RABPLACE_3.0','RABPLACE_4.0','RABPLACE_5.0','RABPLACE_6.0','RABPLACE_7.0','RABPLACE_8.0','RABPLACE_9.0','RABPLACE_10.0','RABPLACE_11.0']
        if 'ethnicity' in covar:
            covar_features += ['ethnicity_1.0','ethnicity_2.0','ethnicity_3.0','ethnicity_4.0','pca_distance']
        if 'edu' in covar:
            covar_features += ['RAEDYRS','iearn','mdiv','mrct','cog27']
        if 'relig' in covar:
            covar_features += ['RARELIG_1.0','RARELIG_2.0','RARELIG_3.0','RARELIG_4.0','RARELIG_5.0']
        if 'bmi' in covar:
            covar_features += ['bmi']
        if 'height' in covar:
            covar_features += ['height']
        if 'PRS' in covar:
            covar_features += ['pgs_cog','pgs_bmi','pgs_height','pgs_educ']
        for boot in range(1,num_boot+1):
            coef_data = {}
            coef_data['covars'] = covar_features
            coef_data['slices']={}
            for ii,slice in enumerate(slice_features):
                #print(ii)
                if ii % 50 == 0:
                    print(round(ii/len(slice_features)*100,2))
                features = covar_features + [slice]
                features_pos = [feature_cols_pos[f] for f in features]
                X2 = X[:,features_pos]
                np.random.shuffle(pos_class)
                np.random.shuffle(neg_class)
                # train: undersample majority class
                num_train0 = len(neg_class)
                if firth_fit: num_train0 = 10*len(pos_class)
                num_train1 = len(pos_class)#len(pos_class)#int(fract_sampled*len(pos_class))
                #if num_boot == 1:
                #    num_train0 = int(len(neg_class))
                #    num_train1 = int(len(pos_class))

                # equal number of 0,1 class data
                train_X = np.concatenate([X2[pos_class[:num_train1]],X2[neg_class[:num_train0]]],axis=0)                
                train_outcome = np.concatenate([outcome[pos_class[:num_train1]],outcome[neg_class[:num_train0]]],axis=0)
                # positive class << negative class, so there is no overlap in data between training and testing
                #shuffle_2arrays(train_X,train_outcome)
                if num_boot > 1:
                    # keep the fraction of slices 
                    pos = np.random.choice(list(range(len(train_X))),size=int(fract_sampled*len(train_X))).astype(int)
                    train_X = train_X[pos]
                    train_outcome = train_outcome[pos]
                if np.sum(train_outcome) == 0: continue
                if np.sum(train_X[:,-1]) == 0: continue
                if logit_fit:
                    # if a feature does not vary, ignore
                    bad_data = False
                    for feature in zip(train_X.T):
                        if np.std(feature) == 0: bad_data = True
                    if bad_data: continue
                    exog = sm.add_constant(train_X)
                    try:
                        with warnings.catch_warnings():
                            if firth_fit:
                                coef,waldp,lrtp=firth.fit_pvalues(train_X,train_outcome)
                                pvals = waldp
                                err = [[0,0]]
                            else:
                                warnings.filterwarnings("error")
                                reg = Logit(train_outcome, exog,start_params=[-2]+[0]*len(covar_features)+[2]).fit(maxiter=200)#,full_output=False)#disp=False,full_output=False)#LogisticRegression(penalty='none',solver='newton-cg',max_iter=200)
                                err = reg.cov_params()
                                coef = reg.params#.values
                                slice_coef = coef[-1]
                                slice_err = err[-1,-1]
                                z= np.abs(slice_coef/slice_err)
                                pvals=[2*(1-norm.cdf(z))]
                    except:
                        continue
                    # ignore if error is unusually high (probably bad fit)
                    if err[-1,-1] > 1000000: continue
                    coef_data['slices'][slice] = {}
                    coef_data['slices'][slice]['coef'] = coef
                    coef_data['slices'][slice]['err'] = err
                    coef_data['slices'][slice]['logp'] = -np.log10(pvals)
                else:
                    continue
                    fract_sampled_shap = 0.8
                    num_train0 = int(fract_sampled_shap*len(neg_class))
                    num_train1 = int(fract_sampled_shap*len(pos_class))

                    # equal number of 0,1 class data
                    train_X = np.concatenate([X2[pos_class[:num_train1]],X2[neg_class[:num_train0]]],axis=0)                
                    train_outcome = np.concatenate([outcome[pos_class[:num_train1]],outcome[neg_class[:num_train0]]],axis=0)
                    # positive class << negative class, so there is no overlap in data between training and testing
                    #shuffle_2arrays(train_X,train_outcome)
                    test_X = np.concatenate([X2[pos_class[num_train1:]],X2[neg_class[num_train0:]]],axis=0)                
                    test_outcome = np.concatenate([outcome[pos_class[num_train1:]],outcome[neg_class[num_train0:]]],axis=0)
                    pos = np.random.choice(list(range(len(test_X))),size=int(fract_sampled*len(test_X))).astype(int)
                    test_X = test_X[pos]
                    test_outcome = test_outcome[pos]

                    # positive class << negative class, so there is no overlap in data between training and testing
                    #shuffle_2arrays(test_X,test_outcome)

                    random_forest = RandomForestClassifier(n_estimators=100)                    
                    random_forest.fit(train_X, train_outcome)
                    explainer = shap.TreeExplainer(random_forest)
                    shap_values = explainer.shap_values(test_X)
                    slice_shap = shap_values[:,-1]
                    coef_data['slices'][slice] = {}
                    coef_data['slices'][slice]['coef'] = np.mean(slice_shap)
                    coef_data['slices'][slice]['err'] = np.std(slice_shap)/np.sqrt(len(slice_shap))  
                    #shap.summary_plot(shap_values, train_X)
                    #shap_values = shap.TreeExplainer(random_forest).shap_values(test_X)
                    #shap.summary_plot(shap_values, train_X,plot_type="layered_violin",show=False)
                    #plt.savefig('ModelFigs/SHAP_'+covar+'_comm='+str(community_ibd)+'_sampled='+str(fract_sampled)+'_boot='+str(boot)+'.png',dpi=150)
                    #plt.close()
                    print('data saved')

            #if boot > 0 and fract_sampled == 1.0: break
            rand_boot = np.random.randint(0,10000*boot)
            if num_boot == 1:
                rand_boot = -1
            pk.dump(coef_data,open('ModelCoefs/Model_coefs_Logit='+str(logit_fit)+'_'+covar+'_comm='+str(community_ibd)+'_sampled='+str(fract_sampled)+'_bootID='+str(rand_boot)+'_firth='+str(firth_fit)+'.pkl','wb'))

def slice_GWAS(community_ibd,logit_fit,num_threads,firth_fit,last_slice):
    covars = ['place','ethnicity','edu','relig','bmi','height']#,'PRS']


    # this will create all permutations of covars
    covar_perms = []
    for i in range(1,len(covars)+1):
        covar_perms += list(combinations(covars,i))
    covar_perms = ['none','+'.join(covars)]+['+'.join(list(ll)) for ll in list(combinations(covars,len(covars)-1))]#['+'.join(covars)] + ['+'.join(list(ll)) for ll in list(combinations(covars,len(covars)-1))]#['none']+['+'.join(list(ll)) for ll in covar_perms]
    random.shuffle(covar_perms)
    #batch_size = int(len(covar_perms)/num_threads)+1
    #'place+ethnicity','place+ethnicity+edu','place+ethnicity+edu+bmi','place+ethnicity+edu+bmi+height']
    covar_vars = []
    num_boot = 1
    if firth_fit:
        num_boot = 1
    fract_sampled = 1.0

    inputs = [[num_boot,fract_sampled,community_ibd,logit_fit,covar,firth_fit,last_slice] for covar in covar_perms]
    #inputs_split = [covar_perms[batch_size*i:batch_size*(i+1)] if batch_size*(i+1) < len(covar_perms) else covar_perms[batch_size*i:] for i in range(num_threads)]
    if num_threads == -1:
        num_threads = len(covar_perms)
    if num_threads > 4:
        num_threads = 4

    #pool = Pool(num_threads)
    out = consume(map(parallel_GWAS,inputs))


def find_slice_pvals(data):
    slice_pvals = {}
    for key in list(data['slices'].keys()):
        chrom = float(key.split('_')[0].replace('chr',''))
        coef = data['slices'][key]['coef'][-1]
        err = data['slices'][key]['err'][-1,-1]
        z= np.abs(coef/err)
        p=2*(1-norm.cdf(z))
        slice_pvals[key]=p
    return slice_pvals

def parallel_predict(input):
        #married_pairs,rand_pairs,community_ibd,sum_slices,slice_only,num_common,data,covar = input
        married_pairs,rand_pairs,community_ibd,sum_slices,data,covar = input
        X,outcome,feature_cols = collect_data(married_pairs,rand_pairs,community_ibd,True)
        slice_only = [pos for pos,f in enumerate(feature_cols) if 'chr' in f]
        num_common = np.array([[np.sum(line[slice_only])] for line in X])

        features_used = []
        if 'place' in covar:
            features_used += ['RABPLACE_1.0','RABPLACE_2.0','RABPLACE_3.0','RABPLACE_4.0','RABPLACE_5.0','RABPLACE_6.0','RABPLACE_7.0','RABPLACE_8.0','RABPLACE_9.0','RABPLACE_10.0','RABPLACE_11.0']
        if 'ethnicity' in covar:
            features_used += ['ethnicity_1.0','ethnicity_2.0','ethnicity_3.0','ethnicity_4.0','pca_distance']
        if 'edu' in covar:
            features_used += ['RAEDYRS','RARELIG_1.0','RARELIG_2.0','RARELIG_3.0','RARELIG_4.0','RARELIG_5.0','cog27','iearn','mdiv','mrct']
        if 'bmi' in covar:
            features_used += ['bmi']
        if 'height' in covar:
            features_used += ['height']
        if 'PRS' in covar:
            features_used += ['pgs_cog','pgs_bmi','pgs_height','pgs_educ']
        covar_only_pos = [pos for pos,f in enumerate(feature_cols) if 'chr' not in f and f in features_used]
        print(features_used)
        for p_val in [0]:#,10**-3,10**-5]:#[1.1,10**-1,10**-2,10**-3,10**-4,10**-6,10**-7,0]:
            print(p_val)
            slice_features = [key for key in feature_cols if 'chr' in key and key in data['slices'].keys()]
            slice_pvals = find_slice_pvals(data)
            slice_features = [key for key in slice_features if slice_pvals[key]<p_val]
            features_used += slice_features
            slice_features_pos = [list(feature_cols).index(key) for key in slice_features]
            X = X[:,covar_only_pos + slice_features_pos]

            if sum_slices:
                features_used  += ['num_slices']
                X = np.concatenate((X,num_common),axis=1)
            pos_class = np.array([pos for pos,o in enumerate(outcome) if o == 1])
            neg_class = np.array([pos for pos,o in enumerate(outcome) if o == 0])
            num_boot = 20
            model_outcomes = {'boot':[],'Metrics':[],'None':[],'Logistic_L2':[]}#,'Logistic_L1':[],'Logistic_Elastinet':[]}
            model_coefs=[]
            for n in range(num_boot):
                        print(round(n/num_boot*100,3))
                        np.random.shuffle(pos_class)
                        np.random.shuffle(neg_class)
                        # train: undersample majority class
                        num_train0 = int(0.7*len(neg_class))
                        num_train1 = int(0.7*len(pos_class))
                        # equal number of 0,1 class data
                        train_X = np.concatenate([X[pos_class[:num_train1]],X[neg_class[:num_train1]]],axis=0)
                        train_outcome = np.concatenate([outcome[pos_class[:num_train1]],outcome[neg_class[:num_train1]]],axis=0)
                        # positive class << negative class, so there is no overlap in data between training and testing
                        #if equal_class:
                        #    test_X = np.concatenate([X[pos_class[num_train1:]],X[neg_class[-len(pos_class[num_train1:]):]]],axis=0)#X[neg_class[num_train0:]]],axis=0)
                        #    test_outcome = np.concatenate([outcome[pos_class[num_train1:]],outcome[neg_class[-len(pos_class[num_train1:]):]]],axis=0)#outcome[neg_class[num_train0:]]],axis=0)
                        #else:
                        test_X = np.concatenate([X[pos_class[num_train1:]],X[neg_class[num_train0:]]],axis=0)
                        test_outcome = np.concatenate([outcome[pos_class[num_train1:]],outcome[neg_class[num_train0:]]],axis=0)
                        if len(test_X) == 0 or len(train_X) == 0: continue
                        shuffle_2arrays(train_X,train_outcome)
                        shuffle_2arrays(test_X,test_outcome)
                        models = [LogisticRegression(penalty='none',solver='newton-cg',max_iter=200),LogisticRegressionCV(cv=5,penalty='l2',max_iter=200)]
                        metrics = [accuracy_score,f1_score,brier_score_loss,roc_auc_score]
                        metric_names = ['accuracy','accuracy_null','f1','f1_null','brier','brier_null','roc_auc','roc_auc_null']
                        model_outcomes['boot']+=[n]*len(metric_names)
                        model_outcomes['Metrics']+=metric_names
                        for ii,[model,key] in enumerate(zip(models,list(model_outcomes.keys())[2:])):
                            #print(ii/(len(list(model_outcomes.keys()))-1))
                            reg = model.fit(train_X, train_outcome)
                            if key == 'None':
                                model_coefs+=np.array(reg.coef_).tolist()

                            pred = model.predict(test_X)
                            pred_random = np.array([0]*len(pred))# assume NO marriages#np.copy(pred)
                            class_pos = list(model.classes_).index(1)
                            prob_pred = model.predict_proba(test_X)[:,class_pos]
                            #prob_pred_random = np.array([0]*len())#np.copy(prob_pred)
                            #np.random.shuffle(pred_random)
                            #np.random.shuffle(prob_pred_random)
                            metric_scores = []
                            if np.sum(test_outcome) != 0 and np.sum(test_outcome) != len(test_outcome):
                                for metric,name in zip(metrics,metric_names[::2]):
                                    model_score = metric(test_outcome,pred)
                                    null_model = metric(test_outcome,pred_random)
                                    if name == 'roc_auc':
                                        model_score = metric(test_outcome,prob_pred)
                                        null_model = None#metric(test_outcome,prob_pred_random)
                                    metric_scores+=[model_score,null_model]
                            model_outcomes[key]+=metric_scores
            print('MODEL COEFS')
            print(features_used)
            pd.DataFrame(model_outcomes).to_csv('ModelResults/ModelResults_comm='+str(community_ibd)+'_slice_p-val='+str(-np.log10(p_val))+'_sum_slice='+str(sum_slices)+'_covars='+covar+'_numboot='+str(num_boot)+'.csv',index=False)
            pd.DataFrame(data=model_coefs,columns=features_used).to_csv('ModelCoefs/ModelCoefs_Logistic_comm='+str(community_ibd)+'_slice_p-val='+str(-np.log10(p_val))+'_sum_slice='+str(sum_slices)+'_covars='+covar+'_numboot='+str(num_boot)+'.csv',index=False)
            #pd.DataFrame(data=features_used,columns=['features']).to_csv('ModelResults/ModelFeatures_comm='+str(community_ibd)+'_slice_p-val='+str(-np.log10(p_val))+'_sum_slice='+str(sum_slices)+'_covars='+covar+'.csv',index=False)

            print('FINISHED ANALYSIS')

def prs_slice_features(community_ibd,sum_slices,num_threads=-1):
    dump_dir = ''#dump_dir = '/project/arpawong_181/burghard_687/genetic_data/'
    rand_pairs=np.load(dump_dir+'random_pairs_cog27.npy')
    married_pairs=np.load(dump_dir+'married_pairs_cog27.npy')
    data = pk.load(open('ModelCoefs/Model_coefs_place+ethnicity+edu+bmi+height_comm='+str(community_ibd)+'.pkl','rb'))
    print('Downloaded!')

    covars = ['place','ethnicity','edu','bmi','height','PRS']    
    # this will create all permutations of covars
    covar_perms = []
    for i in range(1,len(covars)+1):
        covar_perms += list(combinations(covars,i))
    covar_perms =  ['+'.join(list(ll)) for ll in covar_perms]
    covar_perms = covar_perms[::-1]
    print(covar_perms)
    if num_threads == -1:
        num_threads = len(covar_perms)  
    if num_threads > 5:
        num_threads = 5
    pool = Pool(num_threads)
    inputs = [[married_pairs,rand_pairs,community_ibd,sum_slices,data,covar] for covar in covar_perms]
    #pool = Pool(num_threads)
    #out = consume(pool.map(parallel_predict,inputs))
    for ii,input in enumerate(inputs):
        #print(input)
        print(ii)
        parallel_predict(input)
    #out = zip(*pool.
    #for input in inputs:
    #    parallel_predict(input)

        
def main():
    train80=False#True
    args = sys.argv[1:]
    args = {key:val for key,val in zip(args[::2],args[1::2])}
    # args is a list of the command line args
    married_only = True
    equal_class = True
    only_covar = True
    min_feature = True
    per_chrom = True
    one_ethnic = True
    num_threads = 1
    logit_fit = True
    firth_fit =True
    last_slice = None
    if '--slice' in args.keys():
        last_slice = int(float(args['--slice']))
    for community_ibd in [False]:#[True,False]:#[True]:
        if firth_fit:
            firth_r_GWAS(num_threads,train80,one_ethnic);continue
        #slice_GWAS(community_ibd,logit_fit,num_threads,firth_fit,last_slice); continue
        #for sum_slices in [True,False]:
        #    prs_slice_features(community_ibd,sum_slices,num_threads); continue
if __name__ == "__main__":
    main()

