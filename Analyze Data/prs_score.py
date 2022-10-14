import pandas as pd
import numpy as np
from glob import glob
import time,os


dump_dir = '/project/burghard_687/genetic_data/'

# collect coefficients for all slices with p-value > threshold
train80=True
slice_prs = True
covars = ['place','ethnicity','relig','edu','height','cog27']#'bmi'
feature_sets = ['+'.join(covars)+'+nopca','+'.join(covars),'none']#['none','+'.join(covars),'+'.join(covars)+'+nopca','height','relig','edu_only','ethnicity_only','pca_distance','place','income','num_div','num_marriage']#,'place+ethnicity+relig','none']


for absolute_features in [True]:#[True,False]:
    for weight in [True,False]:
        for p_value in [10**-2,10**-2.5,10**-3,10**-3.5,10**-4,10**-4.5,10**-5]:
            for feature_set in feature_sets:
                print([p_value,feature_set])
                if slice_prs:
                    if train80:
                        n_folds = 5
                        for fold in range(n_folds):
                            outcome_file_train = dump_dir+'slice_outcome_train80_fold='+str(fold)+'_cog27.npy'
                            outcome_file_test = dump_dir+'slice_outcome_test20_fold='+str(fold)+'_cog27.npy'
                            outcome_test = np.load(outcome_file_test)
                            
                            if absolute_features:
                                slice_prob_files_glob = glob('/project/burghard_687/genetic_data/slices_abs_cog27*/firth_'+feature_set+'_prob_*_full_train80_fold='+str(fold)+'_abs_no_cog-bmi_cog27.csv')
                            else:
                                slice_prob_files_glob = glob('slices_nonabs/firth_'+feature_set+'_prob_*_full_train80_fold='+str(fold)+'_nonabs_no_cog-bmi_cog27.csv')
                            slice_prob_files = {'_'.join(f.split(feature_set)[1].split('_')[2:4]):f for f in slice_prob_files_glob if os.path.exists(f.replace('_prob_','_coefs_'))}

                            sig_slices = [slice for slice in slice_prob_files.keys() if pd.read_csv(slice_prob_files[slice])['slice'].values[0] < p_value]
                            if absolute_features:
                                sig_slice_coefs = {slice:pd.read_csv(list(glob('/project/burghard_687/genetic_data/slices_abs_cog27*/firth_'+feature_set+'_coefs_'+slice+'_full_train80_fold='+str(fold)+'_abs_no_cog-bmi_cog27.csv'))[0])['slice'].values[0] for slice in sig_slices}
                            else:
                                sig_slice_coefs = {slice:pd.read_csv('slices_nonabs/firth_'+feature_set+'_coefs_'+slice+'_full_train80_fold='+str(fold)+'_nonabs_no_cog-bmi_cog27.csv')['slice'].values[0] for slice in sig_slices}
                            outcome_train = np.load(outcome_file_train)
                            prs = np.zeros(outcome_test.shape)
                            for slice in sig_slices:
                                # what user pairs have these slices in common??
                                slice_covar_train = pd.read_csv(dump_dir+'slices80/train_'+slice+'_train80_fold='+str(fold)+'_cog27.csv')['slice']
                                weighted_coef = sig_slice_coefs[slice]*slice_covar_train.values
                                if not weight:
                                    weighted_coef = slice_covar_train.values

                                slice_covar_train = pd.read_csv(dump_dir+'slices80/train_'+slice+'_train80_fold='+str(fold)+'_cog27.csv')['slice']
                                slice_covar_test = pd.read_csv(dump_dir+'slices80/test_'+slice+'_test20_fold='+str(fold)+'_cog27.csv')['slice']
                                weighted_coef = sig_slice_coefs[slice]*slice_covar_test.values
                                if not weight:
                                    weighted_coef = slice_covar_test.values
                                try:
                                    prs += weighted_coef
                                except:
                                    print([weight,p_value,fold,slice,dump_dir+'slices80/test_'+slice+'_test20_fold='+str(fold)+'_cog27.csv'])

                            # compare to marriage probability
                            outcome = list(zip([list(prs),list(outcome_test)]))
                            outcome = np.array([list(prs),list(outcome_test)]).T
                            print('Mean PRS: ',outcome[:,0].mean())
                            if absolute_features:
                                pd.DataFrame(outcome,columns = ['prs','marriage']).to_csv('prs_abs/slice_prs_p-value='+str(-np.log10(p_value))+'_'+feature_set+'_weight='+str(weight)+'_train80='+str(train80)+'_fold='+str(fold)+'_abs_no_bmi_cog_cog27.csv',index=False)
                            else:
                                pd.DataFrame(outcome,columns = ['prs','marriage']).to_csv('prs_nonabs/slice_prs_p-value='+str(-np.log10(p_value))+'_'+feature_set+'_weight='+str(weight)+'_train80='+str(train80)+'_fold='+str(fold)+'_nonabs_no_bmi_cog_cog27.csv',index=False)

                    else:
                        outcome_file_train = dump_dir+'slice_outcome_cog27.npy'
                        slice_prob_files_glob = glob('slices/firth_'+feature_set+'_prob_*_full_cog27.csv')
                        slice_prob_files = {'_'.join(f.split('_')[-3:-1]):f for f in slice_prob_files_glob}
                        sig_slices = [slice for slice in slice_prob_files.keys() if pd.read_csv(slice_prob_files[slice])['slice'].values[0] < p_value]
                        sig_slice_coefs = {slice:pd.read_csv('slices/firth_'+feature_set+'_coefs_'+slice+'_full_cog27.csv')['slice'].values[0] for slice in sig_slices}

                        outcome_train = np.load(outcome_file_train)
                        prs = np.zeros(outcome_train.shape)
                        for slice in sig_slices:
                            print(slice)
                            # what user pairs have these slices in common??
                            slice_covar_train = pd.read_csv(dump_dir+"slices/train_"+slice+"_cog27.csv")['slice']
                            weighted_coef = sig_slice_coefs[slice]*slice_covar_train.values
                            if not weight:
                                weighted_coef = slice_covar_train.values

                            prs += weighted_coef
                        # compare to marriage probability
                        outcome = np.array([list(prs),list(outcome_train)]).T
                        pd.DataFrame(outcome,columns = ['prs','marriage']).to_csv('slice_prs_p-value='+str(-np.log10(p_value))+'_'+feature_set+'_weight='+str(weight)+'_train80='+str(train80)+'_cog27.csv',index=False)

                else: # snps
                    slice_prob_files_glob = glob('snps/firth_'+feature_set+'_prob_*_full_train80_cog27.csv')
                    slice_prob_files = {'_'.join(f.split('_')[-3:-1]):f for f in slice_prob_files_glob}
                    sig_slices = [slice for slice in slice_prob_files.keys() if pd.read_csv(slice_prob_files[slice])['slice'].values[0] < p_value]
                    sig_slice_coefs = {slice:pd.read_csv('slices/firth_'+feature_set+'_coefs_'+slice+'_full_train80_cog27.csv')['slice'].values[0] for slice in sig_slices}
                    outcome_file_train = dump_dir+'slice_outcome_train80_cog27.npy'
                    outcome_file_test = dump_dir+'slice_outcome_test20_cog27.npy'
                    outcome_test = np.load(outcome_file_test)
                    outcome_train = np.load(outcome_file_train)
                    prs = np.zeros(outcome_test.shape)
                    for slice in sig_slices:
                        # what user pairs have these slices in common??
                        slice_covar_train = pd.read_csv(dump_dir+"slices/train_"+slice+"_train80_cog27.csv")['slice']
                        slice_covar_test = pd.read_csv(dump_dir+"slices/test_"+slice+"_test20_cog27.csv")['slice']
                        weighted_coef = sig_slice_coefs[slice]*slice_covar_train.values
                        if not weight:
                             weighted_coef = slice_covar_train.values
                        prs += weighted_coef
                    # compare to marriage probability
                    if train80:
                        outcome = np.array([list(prs),list(outcome_test)]).T

                    else:
                        outcome = np.array([list(prs),list(outcome_train)]).T
                    pd.DataFrame(outcome,columns = ['prs','marriage']).to_csv('snp_prs_p-value='+str(-np.log10(p_value))+'_'+feature_set+'_weight='+str(weight)+'_train80='+str(train80)+'_cog27.csv',index=False)

