from scipy.stats import spearmanr
import pandas as pd
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
#import karateclub
# docs for karateclub's BigCLAM: https://karateclub.readthedocs.io/en/latest/_modules/karateclub/community_detection/overlapping/bigclam.html#BigClam.fit
import pickle as pk
from glob import glob
import random
from random import shuffle
#import community
from typing import Dict
#from karateclub.estimator import Estimator
from scipy.stats import binom_test
from sklearn.metrics import normalized_mutual_info_score, adjusted_mutual_info_score
from networkx.algorithms.community import greedy_modularity_communities
from networkx.algorithms.core import k_core
#from cdlib import algorithms
#import infomap
#import markov_clustering as mc
import os,time,copy
from scipy.stats import binom
from joblib import dump, load
from sklearn.neighbors import KernelDensity
import seaborn as sns

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
def find_unique(IDs,first_spouse):
    unique_pairs=[]
    for ii in IDs:
        s1=first_spouse[ii]
        if s1 is not None:
            if float(s1) > ii:
                unique_pairs.append([ii,s1])
    return unique_pairs

def find_pca(IDs,df_pca,pca_num):
    PCA={}
    for ii in IDs:
        pca = df_pca.loc[df_pca['IID']==ii,'EV'+str(pca_num)].values
        if len(pca) > 0:
            pca=pca[0]
        else:
            pca = None
        PCA[ii] = pca
    return PCA

def find_gendered_pairs(unique_pairs,sex):
    gendered_pairs = []
    for ii,jj in unique_pairs:
        if sex[ii] != sex[jj] and sex[ii] is not None and sex[jj] is not None:
            if sex[ii]==1:
                gendered_pairs.append([ii,jj])
            else:
                gendered_pairs.append([jj,ii])
    return gendered_pairs
def find_spouses(df_pca,df_links,IDs):
    first_spouse={}
    second_spouse={}
    for ii in IDs:
        s1=[None]
        s2=[None]
        if ii in df_links['IID']:
            s1 = df_links.loc[df_links['IID']==ii,'RASPIID1'].values
            s2 = df_links.loc[df_links['IID']==ii,'RASPIID2'].values
        else:
            if ii in df_links['RASPIID1']:
                s1 = df_links.loc[df_links['RASPIID1']==ii,'IID'].values
            if ii in df_links['RASPIID2']:
                s2 = df_links.loc[df_links['RASPIID2']==ii,'IID'].values

        if len(s1) > 0:
            s1 = s1[0]
        else:
            s1 = None
        if len(s2) > 0:
            s2= s2[0]
        else:
            s2 = None
        first_spouse[ii] = s1
        second_spouse[ii] = s2
    return first_spouse,second_spouse


def find_marriages():
    directory='/project/burghard_687/genetic_data/'#'/project/arpawong_181/HRS_AsMa/keith/'
    df_pca = pd.read_csv(directory+'hrs_pcs15506_cog27.csv')
    df_links = pd.read_csv(directory+'hrs_spouse_IIDs.csv')
    # columns
    # FID,IID,sex,EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,EV11,EV12,EV13,EV14,EV15,EV16,EV17,EV18,EV19,EV20
    IDs=np.unique(list(df_pca['IID'].values)+list(df_links['IID'].values) + list(df_links['RASPIID1'].values)+  list(df_links['RASPIID2'].values))
    IDs = IDs[~np.isnan(IDs)]
    sex = find_sex(df_pca,IDs)
    first_spouse, second_spouse = find_spouses(df_pca,df_links,IDs)
    unique_pairs = find_unique(IDs,first_spouse)
    gendered_pairs=find_gendered_pairs(unique_pairs,sex)
    return gendered_pairs

def null_marriages(num_marriages,num_nodes,community_size,num_pairs):
    # probability of at least X marriages appearing in a random community of size n
    G = nx.Graph([(i+num_nodes,i+2*num_nodes) for i in range(num_marriages)])
    # single nodes:
    single_nodes = num_nodes - 2*num_marriages
    G.add_nodes_from(list(range(single_nodes)))
    num_samples = 10000
    num_larger = 0
    for samp in range(num_samples):
        if samp % 1000 == 0:
            print(round(samp/num_samples*100,1),'%')
        random_nodes = random.sample(list(G.nodes),community_size)
        H = G.subgraph(random_nodes)
        num_dyads = len([cc for cc in nx.connected_components(H) if len(cc)==2])
        if num_dyads >= num_pairs:
            num_larger += 1
    return num_larger/num_samples


def cov_matrix(cov_df):
    # covariate matrix:
    cov_mat = [[col for col in cov_df.columns if col not in ['FID', 'IID']]]
    for col1 in cov_df.columns:
        for col2 in cov_df.columns:
            cov_mat.append([])
            if col1 not in ['FID', 'IID'] and col2 not in ['FID', 'IID']:
                #new_cov_df = cov_df[[col1,col2]].dropna().replace('M',-np.inf).replace('R',-np.inf).replace('F',-np.inf)
                new_cov_df = cov_df.loc[(cov_df[col1]>-np.inf)&(cov_df[col2]>-np.inf),]
                cov_mat[-1].append(spearmanr(new_cov_df[col1].values,new_cov_df[col2].values))
    
    np.save('cov_matrix.npy',np.array(cov_mat))



def pca_2dplot(IDs,gendered_pairs,df_pca,sex,only_couples,markers,out_file,plot_unmarried=False):
    gendered_pairs=np.array(gendered_pairs)
    PCA1 = find_pca(IDs,df_pca,1)
    PCA2 = find_pca(IDs,df_pca,2)
    #PCA3 = find_pca(IDs,df_pca,2)
    classes = np.unique(df_pca['ethnicity'].values)#list(np.unique([sex[key] for key in sex.keys() if sex[key] is not None]))#~np.isnan(sex[key])]))
    non_married_IDs = np.array([ID for ID in df_pca['IID'].values if ID not in IDs])
    PCA1_unmarried = find_pca(non_married_IDs,df_pca,1)
    PCA2_unmarried = find_pca(non_married_IDs,df_pca,2)

    print(df_pca['ethnicity'].values[:10])

    class_per_person = [df_pca.loc[df_pca['IID']==ii,'ethnicity'].values if ii in df_pca['IID'].values else [None] for ii in list(IDs)]
    class_per_person_unmarried = [df_pca.loc[df_pca['IID']==ii,'ethnicity'].values if ii in df_pca['IID'].values else [None] for ii in list(non_married_IDs)]
    
    class_per_person = [c[0] for c in class_per_person]
    class_per_person_unmarried = [c[0] for c in class_per_person_unmarried]
    #classes.reverse()
    print(classes)
    PCA12 = [np.array([[PCA1[iid],PCA2[iid]] for ii,iid in enumerate(IDs) if PCA1[iid] is not None and PCA2[iid] is not None and class_per_person[ii]==c]) for c in classes]#sex[ii]==c]) for c in classes]
    PCA12_unmarried = [np.array([[PCA1_unmarried[iid],PCA2_unmarried[iid]] for ii,iid in enumerate(non_married_IDs) if PCA1_unmarried[iid] is not None and PCA2_unmarried[iid] is not None and class_per_person_unmarried[ii]==c]) for c in classes]
    #male_PCA12 = np.array([[PCA1[ii],PCA2[ii]] for ii in gendered_pairs[:,0] if PCA1[ii] is not None and PCA2[ii] is not None])
    #female_PCA12 = np.array([[PCA1[ii],PCA2[ii]] for ii in gendered_pairs[:,1] if PCA1[ii] is not None and PCA2[ii] is not None])
    partner_PCA_pairs = [np.array([[PCA1[ii],PCA2[ii]],[PCA1[jj],PCA2[jj]]]) for ii,jj in gendered_pairs if PCA1[ii] is not None and PCA1[jj] is not None]


    #figure=plt.figure()
    size=5
    mew_val=0.5
    fig, ax = plt.subplots(1, 1,figsize=(5,5))
    for ppca in partner_PCA_pairs:
        plt.plot(ppca[:,0],ppca[:,1],'k-',alpha=0.2,linewidth=0.5,label='_none_')

    if markers:
        #if only_couples:
        #    plt.plot(male_PCA12[:,0],male_PCA12[:,1],'k.',markersize=size)
        #    plt.plot(female_PCA12[:,0],female_PCA12[:,1],'rx',mew=mew_val,markersize=size+2)
        #else:
        if True:
            plot_markers=['md','ko','cs','g^']
            if len(classes) == 2:
                plot_markers=['k.','rx']
            for jj,c in enumerate(classes):
                if plot_unmarried:
                    plt.plot(PCA12_unmarried[jj][:,0],PCA12_unmarried[jj][:,1],plot_markers[jj],markersize=size,alpha=0.2,markerfacecolor='none',label=['_NH White_','_NH Black_','_NH Other_','_Hispanic_'][jj])
                if 'x' in plot_markers[jj]:
                    plt.plot(PCA12[jj][:,0],PCA12[jj][:,1],plot_markers[jj],mew=mew_val,markersize=size+2,alpha=0.3)
                else:
                    plt.plot(PCA12[jj][:,0],PCA12[jj][:,1],plot_markers[jj],markersize=size,alpha=[0.5,0.2][int(float(plot_unmarried))],markerfacecolor='none',label=['_NH White_','_NH Black_','_NH Other_','_Hispanic_'][jj])#label=['NH White','NH Black','NH Other','Hispanic'][jj])#alpha=0.3,label='_none_')
                plt.plot([-10],[-10],plot_markers[jj],markersize=size,alpha=1.0,markerfacecolor='none',label=['NH White','NH Black','NH Other','Hispanic'][jj])
             
        if len(classes) == 2:
            plt.legend(['Male','Female'])
        else:
            plt.legend()
            #plt.legend(['Hispanic','NH Other','NH Black','NH White'])
    plt.xlim([-0.2,0.7])
    plt.ylim([-0.05,0.55])
    plt.xlabel('PCA 1')
    plt.ylabel('PCA 2')
    plt.savefig(out_file)#'PCA_spouse_loc_only_couples='+str(only_couples)+'_markers='+str(markers)+'.pdf')
    plt.close()


def bootstrap_mean(array):
    boot_num=10**4
    boot_mean=[]
    for ii in range(boot_num):
        resampled=np.random.choice(array,replace=True,size=len(array))
        boot_mean.append(np.mean(resampled))
    return [np.mean(boot_mean),np.std(boot_mean)]
def plot_dist(PCA_embed,gendered_pairs,col,cat):
    num_samples = 1000
    #if num_samples > len(PCA_embed.keys())/2:
    #    num_samples = int(len(PCA_embed.keys())/2)
    # create random pairs
    #random_dist = []
    if len(gendered_pairs) ==0:
        return [None]*5
    all_m = list(np.array(gendered_pairs)[:,0])
    all_f = list(np.array(gendered_pairs)[:,1])
    all_keys = set(list(PCA_embed.keys()))
    rand_means = []
    for ii in range(num_samples):
        all_m_copy = all_m[:]
        random.shuffle(all_m_copy)
        all_f_copy = all_f[:]
        random.shuffle(all_f_copy)

        random_dist = []
        for m,f in zip(all_m_copy,all_f_copy):
            if m in all_keys and f in all_keys:
                #rand_keys = list(PCA_embed.keys())#np.unique(np.array(gendered_pairs).flatten())
                #key1 = random.sample(rand_keys,num_samples)
                #key2 = random.sample(rand_keys,num_samples)
                dist = np.sqrt(np.dot(PCA_embed[m]-PCA_embed[f],PCA_embed[m]-PCA_embed[f]))
                #random_dist=np.array([np.sqrt(np.dot(PCA_embed[k1]-PCA_embed[k2],PCA_embed[k1]-PCA_embed[k2])) for k1,k2 in zip(key1,key2) if k1 != k2])
                if ~np.isnan(dist):
                    random_dist.append(dist)
        rand_means.append(np.mean(random_dist))
    #random_dist = random_dist[~np.isnan(random_dist)]
    marriage_dist=np.array([np.sqrt(np.dot(PCA_embed[key1]-PCA_embed[key2],PCA_embed[key1]-PCA_embed[key2])) for key1,key2 in gendered_pairs if key1 in PCA_embed.keys() and key2 in PCA_embed.keys()])
    marriage_dist = marriage_dist[~np.isnan(marriage_dist)]
    marriage_mean = np.mean(marriage_dist)
    # we ignore this for now
    marriage_std = 0
    #marriage_mean,marriage_std = bootstrap_mean(marriage_dist)
    random_mean,[random_q1,random_q2] = [np.mean(rand_means),np.quantile(rand_means,[0.025,0.975])]#bootstrap_mean(random_dist)#marriage_dist)
    if len(PCA_embed.keys()) < 100 or len(marriage_dist) < 100:
        print('TOO LITTLE DATA: ',[col,cat])
        return [random_mean,random_q1,random_q2,marriage_mean,marriage_std]
    random_dist = np.array(random_dist)[:,np.newaxis]
    marriage_dist = np.array(marriage_dist)[:,np.newaxis]
    # kernel density estimation plot
    if False:
        random_kde = KernelDensity(kernel='gaussian',bandwidth=0.002).fit(random_dist)
        marriage_kde = KernelDensity(kernel='gaussian',bandwidth=0.002).fit(marriage_dist)
        X_plot=np.linspace(-0.1,0.5,5000)[:,np.newaxis]
        bins=np.linspace(-1,1,10)
        random_log_dens = random_kde.score_samples(X_plot)
        marriage_log_dens = marriage_kde.score_samples(X_plot)
        fig, ax = plt.subplots(1,1,sharex=True,sharey=True)
        fig.subplots_adjust(hspace=0.1,wspace=0.05)
        ax.plot(X_plot,np.exp(random_log_dens),'k-')
        ax.plot(X_plot,np.exp(marriage_log_dens),'r--')
        ax.fill(X_plot,np.exp(random_log_dens),fc='k',alpha=0.3)
        ax.fill(X_plot,np.exp(marriage_log_dens),fc='r',alpha=0.3)
        ax.set_ylim(0,40)
        ax.set_xlim(-0.01,0.2)#np.max(random_dist))
        ax.set_xlabel('PCA Distance')
        ax.set_ylabel('PDF')
        plt.legend(['Random','Marriages'])
        plt.savefig('dist_PDF_allpcas_'+col+'_class='+str(cat)+'.pdf')
        plt.close()
    return [random_mean,random_q1,random_q2,marriage_mean,marriage_std]


def pop_random(lst):
    idx = random.randrange(0, len(lst))
    return lst.pop(idx)

def bootstrap_null_p(cov_df,pairs,true_mean):
    n = 10**3
    pairs = np.array(pairs)
    bootstrap_null = []
    if len(pairs) == 0:
        return None
    new_cov_df = cov_df[['FID','EV1','EV2']].dropna()

    for ii in range(n):
        if ii % 100 == 0:
            print(ii)
        random_ids1 = np.random.choice(pairs[:,0], size=len(pairs), replace=True)
        random_ids2 = np.random.choice(pairs[:,1], size=len(pairs), replace=True)
        random_pairs = [[r1,r2] for r1,r2 in zip(random_ids1,random_ids2)]
        xy_distribution = [[cov_df.loc[cov_df['FID']==id1,['EV1','EV2']].values,cov_df.loc[cov_df['FID']==id2,['EV1','EV2']].values] for id1,id2 in random_pairs]
        dist_distribution = [np.linalg.norm(xy1[0]-xy2[0]) for xy1,xy2 in xy_distribution if len(xy1) >0 and len(xy2) > 0]
        mean = np.mean(dist_distribution)
        bootstrap_null.append(mean)
    bootstrap_null = np.array(bootstrap_null)
    p = len(bootstrap_null[bootstrap_null < true_mean])/len(bootstrap_null)
    return p

def find_mean_dist(cov_df,gendered_pairs_cat):
    xy_distribution = [[cov_df.loc[cov_df['FID']==id1,['EV1','EV2']].values,cov_df.loc[cov_df['FID']==id2,['EV1','EV2']].values] for id1,id2 in gendered_pairs_cat]#id1 in cov_df['FID'] and id2 in cov_df['FID']]
    dist_distribution = [np.linalg.norm(xy1[0]-xy2[0]) for xy1,xy2 in xy_distribution if len(xy1) >0 and len(xy2) > 0]
    mean = np.mean(dist_distribution)
    return mean

def combine_pca_dem(directory,max_pca):
    # demographics, ignore pca
    file = directory+'hrs_pcs15506_cog27.csv'#covariates.csv'
    df_dem = pd.read_csv(file).dropna().replace('M',-np.inf).replace('R',-np.inf).replace('D',-np.inf).replace('S',-np.inf).replace('N',-np.inf).astype(float)
    #df_dem = pd.read_csv(directory+'hrs_pcs15506_cog27.csv')#race.csv')#'hrs_pcs15506.csv')
    dem_cols = [col for col in df_dem.columns if 'EV' not in col]
    df_dem = df_dem[dem_cols]

    # pcas
    df_pca = pd.read_csv(directory+'phased_data/all_chr_gds_evecs.csv',header=None)
    # rename pca cols
    df_pca.columns=['EV'+str(i) for i in range(1,max_pca+1)]

    # eigenvalues
    df_evals = pd.read_csv(directory+'phased_data/all_chr_gds_evals.csv',header=None)
    df_evals.columns=['evals']
    # de-normalize by multiplying pca components by sqrt(evals)
    df_evals = np.sqrt(df_evals['evals'].values[:max_pca])
    for i in range(1,33):
        df_pca['EV'+str(i)]=df_pca['EV'+str(i)].values*df_evals[i-1]

    # map pca values to ids
    df_ids=pd.read_csv(directory+'phased_data/all_ids.csv',header=None)
    df_ids.columns=['ids']
    df_pca['IID']=df_ids['ids'].values

    # merge PCAs and Demographics
    df_pca = pd.merge(df_dem,df_pca,on='IID')
    print(df_pca.columns)
    return df_pca

def assort_pca(cov_df,gendered_pairs,PCA_correl):
    assortativity_pca = {'pca':[],'assortativity':[],'dist_marriage':[],'p_null':[],'n_total':[],'n_dist':[]}
    #if covs == -1:
    dist_marriage = find_mean_dist(cov_df,gendered_pairs)
    
    p_null = bootstrap_null_p(cov_df,gendered_pairs,dist_marriage)
    for key in PCA_correl.keys():
        value_pairs = PCA_correl[key]
        value_pairs = np.array([[v1[0],v2[0]] for v1,v2 in value_pairs if len(v1) == 1 and len(v2) == 1])
        plt.scatter(value_pairs[:,0],value_pairs[:,1],marker='o',c='k')
        plt.xlabel(key+' Male')
        plt.ylabel(key+' Female')
        plt.savefig('PCA_scatter_all-marriages_pca='+str(key)+'.pdf')
        plt.close()
        assortativity,p = spearmanr(value_pairs[:,0],value_pairs[:,1])
        assortativity_pca['pca'].append(key)
        assortativity_pca['assortativity'].append(assortativity)
        assortativity_pca['dist_marriage'].append(dist_marriage)
        assortativity_pca['p_null'].append(p_null)
        assortativity_pca['n_total'].append(len(gendered_pairs))
        assortativity_pca['n_dist'].append(len(gendered_pairs))
    pd.DataFrame(data=assortativity_pca).to_csv('PCA_assortativity.csv',index=False)


def distance_plots(cov_df,cov_cols,category_cols):
        max_pca_dist = 10
        distance_comparisons = {'feature':[],'category':[],'marriage_mean':[],'marriage_std':[],'random_mean':[],'random_q1':[],'random_q2':[]}
        IDs = cov_df['FID'].values
        gendered_pairs = find_marriages()
        sex = find_sex(cov_df,IDs)
        PCA_embed={}
        for ID in IDs:
            pca_values = cov_df.loc[cov_df['IID']==ID,['EV'+str(i) for i in range(1,max_pca_dist)]].values.flatten()
            if len(pca_values) > 0:
                PCA_embed[int(ID)] = pca_values
        random_mean,random_q1,random_q2,marriage_mean,marriage_std = plot_dist(PCA_embed,gendered_pairs,'all','all')
        distance_comparisons['feature'].append('all')
        distance_comparisons['category'].append('all')
        distance_comparisons['marriage_mean'].append(marriage_mean)
        distance_comparisons['marriage_std'].append(marriage_std)
        distance_comparisons['random_mean'].append(random_mean)
        distance_comparisons['random_q1'].append(random_q1)
        distance_comparisons['random_q2'].append(random_q2)

        for col in cov_cols:
            print(col)

            if col not in category_cols:
                classes = ['above_median','below_median']
                median = np.median(cov_df[col].dropna().values)                
            else:
                classes = np.sort(cov_df[col].dropna().drop_duplicates().values)#list(np.unique([sex[key] for key in sex.keys() if sex[key] is not None]))
                classes = classes[(~np.isnan(classes))&(classes>-np.inf)]
                #classes = [c for c in classes if ~np.isnan(c) and c is not None]
                #classes.reverse()
            cat_list = []
            #print(classes)
            for cat in classes:
                if 'median' in str(cat):
                    IDs = cov_df.loc[cov_df[col]>median,'IID'].values
                    if cat=='below_median':
                        IDs = cov_df.loc[cov_df[col]<=median,'IID'].values
                else:
                    IDs = cov_df.loc[cov_df[col]==cat,'IID'].values
                #print(len(IDs))
                gendered_pairs_class = [[p1,p2] for p1,p2 in gendered_pairs if p1 in IDs and p2 in IDs]
                PCA_embed={}
                for ID in IDs:
                    pca_values = cov_df.loc[cov_df['IID']==ID,['EV'+str(i) for i in range(1,max_pca_dist)]].values[0].flatten()
                    if len(pca_values) > 0:
                        PCA_embed[int(ID)] = pca_values
                random_mean,random_q1,random_q2,marriage_mean,marriage_std = plot_dist(PCA_embed,gendered_pairs_class,col,cat)
                distance_comparisons['feature'].append(col)
                distance_comparisons['category'].append(cat)
                distance_comparisons['marriage_mean'].append(marriage_mean)
                distance_comparisons['marriage_std'].append(marriage_std)
                distance_comparisons['random_mean'].append(random_mean)
                distance_comparisons['random_q1'].append(random_q1)
                distance_comparisons['random_q2'].append(random_q2)
                
        pd.DataFrame(distance_comparisons).to_csv('distance_comparisons.csv',index=False)

def random_cat_p(value_pairs,assortativity_cat):
    categories = sorted(list(np.unique(value_pairs.flatten())))
    assort_boot_all = []
    for i in range(10000):
        # randomize value_pairs
        cat_m = list(value_pairs[:,0])
        random.shuffle(cat_m)    
        cat_m = np.array(cat_m).flatten()
        cat_f = list(value_pairs[:,1])
        random.shuffle(cat_f)    
        cat_f = np.array(cat_f).flatten()

        random_pairs = np.array([cat_m,cat_f]).T
        confusion_matrix=np.array([[len(random_pairs[(random_pairs[:,0]==c1)&(random_pairs[:,1]==c2)])/len(random_pairs) for c1 in categories] for c2 in categories])
        e2 = np.sum(confusion_matrix**2)
        assortativity_boot = (np.trace(confusion_matrix)-e2)/(1-e2)
        assort_boot_all.append(assortativity_boot)
    assort_boot_all = np.array(assort_boot_all)
    p = len(assort_boot_all[assort_boot_all>assortativity_cat])/len(assort_boot_all)
    return p

def assort_pca_plots(cov_df,gendered_pairs,category_cols,col):
    new_cov_df = cov_df.loc[(cov_df[col]>-np.inf),]
    value_pairs = np.array([[new_cov_df.loc[new_cov_df['FID']==id1,col].values,new_cov_df.loc[new_cov_df['FID']==id2,col].values] for id1,id2 in gendered_pairs])
    value_pairs = np.array([[v1[0],v2[0]] for v1,v2 in value_pairs if len(v1) == 1 and len(v2) == 1])
    assortativity_pca_dist_cat = {'category':[],'pca':[],'assortativity':[],'dist_marriage':[],'p_null':[],'n_category':[]}
        
    # create confusion matrix, 
    if col in category_cols:
        categories = sorted(list(np.unique(value_pairs.flatten())))
        confusion_matrix=np.array([[len(value_pairs[(value_pairs[:,0]==c1)&(value_pairs[:,1]==c2)])/len(value_pairs) for c1 in categories] for c2 in categories])
        np.save(col+'_confusion-matrix.npy',confusion_matrix)
        e2 = np.sum(confusion_matrix**2)
        assortativity_cat = (np.trace(confusion_matrix)-e2)/(1-e2)
        p = random_cat_p(value_pairs,assortativity_cat)
        print(col,' p-value: ',p)
    else:

        assortativity_cat,p = spearmanr(value_pairs[:,0],value_pairs[:,1])
        categories = ['above_median','below_median']
        median = np.median(new_cov_df[col].dropna().values)
        print(col,' p-value: ',p)

        #print('p-value: ',p)
    return assortativity_cat

def plot_covariates(cov_df,category_cols,gendered_pairs,col):
    figure=plt.figure()
    labelfonts = {'fontname':'Arial','fontsize':12}
    if col not in category_cols:
        cov_plotxy = [[cov_df.loc[cov_df['IID']==x,col].values[0],cov_df.loc[cov_df['IID']==y,col].values[0]] for x,y in gendered_pairs if x in cov_df['IID'].values and y in cov_df['IID'].values]
        cov_plotxy = np.array(cov_plotxy)

        if col == 'mdiv' or col=='mrct':
            plt.plot(cov_plotxy[:,0],cov_plotxy[:,1],'ks',alpha=0.2,markersize=8)
        else:
            plt.plot(cov_plotxy[:,0],cov_plotxy[:,1],'ks',alpha=0.2)
        if col == 'bmi':
            plt.xlim([10,60])
        elif col=='iearn':
            plt.axes().set_xscale('log')
            plt.axes().set_yscale('log')
    else:
        categories={1:'M',2:'F'}
        if col=='ethnicity':
            categories={1:'White',2:'Black',3:'Other',4:'Hispanic'}
        elif col == 'RABPLACE':
            categories={1:'New England',2:'Mid Atlantic',3:'EN Central',4:'WN Central',5:'S Atlantic',6:'ES Central',7:'WS Central',8:'Mountain',9:'Pacific',10:'US/NA Division',11:'Not US'}
        elif col == 'RARELIG':
            categories={1:'Protestant',2:'Catholic',3:'Jewish',4:'None',5:'Other'}

        cov_plotxy = [[cov_df.loc[cov_df['IID']==x,col].values[0],cov_df.loc[cov_df['IID']==y,col].values[0]] for x,y in gendered_pairs if x in cov_df['IID'].values and y in cov_df['IID'].values]
        cov_plotxy = np.array(cov_plotxy)
        
        cats = np.unique(cov_plotxy.flatten())
        cats = cats[cats>-np.inf].astype(int)
        cov_mat = np.array([[len(cov_plotxy[(cov_plotxy[:,1]==cr)&(cov_plotxy[:,0]==cc)])  for cc in cats] for cr in cats])
        green = sns.light_palette("seagreen", reverse=False, as_cmap=True)
        green.set_under('gray')#'tomato')
        ax = sns.heatmap(cov_mat,cmap=green,vmin=0.01,cbar_kws={'extend': 'min'})
        x=list(np.arange(0.5,len(cats),1))
        
        #if feature in category_features:
        my_xticks = ['']
        print([col,categories])
        if len(x)>1:
            my_xticks=[categories[val] for val in cats]
        if len(x)>2:
            plt.yticks(x, my_xticks,rotation=0,**labelfonts)
            plt.xticks(x, my_xticks,rotation=90,**labelfonts)
        else:
            plt.yticks(x, my_xticks,**labelfonts)
            plt.xticks(x, my_xticks,**labelfonts)

    plt.xlabel('Male '+col.replace('_',' ').capitalize(),**labelfonts)
    plt.ylabel('Female '+col.replace('_',' ').capitalize(),**labelfonts)
    plt.tight_layout()
    plt.savefig(col+'_mvsf.pdf')
    plt.close()


def main():


    pca_dist=True#False#True
    plot_covar=True#False#True
    assort=True
    max_pca = 32
    directory = '/project/burghard_687/genetic_data/'#'/project/arpawong_181/HRS_AsMa/keith/'
    cov_df = combine_pca_dem(directory,max_pca)
    boot_num = 100
    cov_cols = [col for col in cov_df.columns if col not in ['FID', 'IID','nhw','nhb','nho','hisp','RAHISPAN','RARACEM'] and 'EV' not in col]
    category_cols = ['RABPLACE','RARELIG','ethnicity','sex']
    # PCA correlations
    gendered_pairs = find_marriages()
    cov_matrix(cov_df[[c for c in cov_cols if c not in category_cols]])
    df_pca = cov_df.copy(deep=True)
    df_pca['IID'] = cov_df['FID'].values
    PCA_correl = {'EV'+str(i):np.array([[cov_df.loc[cov_df['FID']==id1,'EV'+str(i)].values,cov_df.loc[cov_df['FID']==id2,'EV'+str(i)].values] for id1,id2 in gendered_pairs]) for i in range(1,21)}
    IDs = list(set([f for m,f in gendered_pairs]+[m for m,f in gendered_pairs]))
    sex = find_sex(df_pca,IDs)
    for only_couples in [True,False]:
        plot_unmarried = not only_couples
        markers = True
        out_file = 'PCA_Embedding_couples='+str(only_couples)+'.pdf'
        pca_2dplot(IDs,gendered_pairs,df_pca,sex,only_couples,markers,out_file,plot_unmarried)
    if pca_dist:
        distance_plots(cov_df,cov_cols,category_cols)
        print('distances plotted')
        #assort_pca(cov_df,gendered_pairs,PCA_correl)

    cov_assortativity={'covariate':[],'assortativity':[]}
    for cov_feature_num,col in enumerate(cov_cols):
        if plot_covar:
            plot_covariates(cov_df,category_cols,gendered_pairs,col)
        if assort:
            assortativity = assort_pca_plots(cov_df,gendered_pairs,category_cols,col)
            cov_assortativity['covariate'].append(col)
            cov_assortativity['assortativity'].append(assortativity)
    if len(cov_assortativity['covariate']) == len(cov_cols):
        pd.DataFrame(cov_assortativity).to_csv('covariate_assortativity_cog27.csv',index=False)



if __name__ == "__main__":
    main()


