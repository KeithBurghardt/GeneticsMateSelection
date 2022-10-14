import pandas as pd 
import numpy as np
import os,random
# causal modeling method
#from econml.metalearners import XLearner
import scipy
from scipy.stats import spearmanr


def is_invertible(a):
    return a.shape[0] == a.shape[1] and np.linalg.matrix_rank(a) == a.shape[0]


def pca_distance(cov_df,pairs):
    max_pca = 20#int(max([float(c.replace('EV','')) for c in cov_df.columns if 'EV' in c]))
    pcas = ['EV'+str(i) for i in range(1,max_pca+1)]
    # L1 distance
    k = 1
    all_pair_dist = []
    for p1,p2 in pairs:
        f1 = cov_df.loc[cov_df['IID']==p1,pcas].values[0]
        f2 = cov_df.loc[cov_df['IID']==p2,pcas].values[0]
        pair_dist = np.sum(np.abs(f1-f2)**k)**(1/k)
        all_pair_dist.append(pair_dist)
    return all_pair_dist



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


def find_unique(IDs,spouse):
    unique_pairs=[]
    for ii in IDs:
        s=spouse[ii]
        if s is not None:
            if float(s) > ii:
                unique_pairs.append((ii,s))
    return unique_pairs

def find_gendered_pairs(unique_pairs,sex):
    gendered_pairs = []
    for ii,jj in unique_pairs:
        if sex[ii] != sex[jj] and sex[ii] is not None and sex[jj] is not None:
            if sex[ii]==1:
                gendered_pairs.append([ii,jj])
            else:
                gendered_pairs.append([jj,ii])
    return gendered_pairs

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

def find_marriages(cov_df):
    directory='/project/arpawong_181/HRS_AsMa/keith/'
    directory = '/project/burghard_687/genetic_data/'
    df_pca = pd.read_csv(directory+'hrs_pcs15506_cog27.csv')
    df_links = pd.read_csv(directory+'hrs_spouse_IIDs.csv')
    cov_df = cov_df.dropna()
    # columns
    IDs=np.unique(list(df_pca['IID'].values)+list(df_links['IID'].values) + list(df_links['RASPIID1'].values)+  list(df_links['RASPIID2'].values))
    IDs = IDs[~np.isnan(IDs)]
    sex = find_sex(df_pca,IDs)
    first_spouse, second_spouse = find_spouses(df_pca,df_links,IDs)
    unique_pairs = list(set(find_unique(IDs,first_spouse)+find_unique(IDs,second_spouse)))
    gendered_pairs=find_gendered_pairs(unique_pairs,sex)
    gendered_pairs = [tuple([p1,p2]) for p1,p2 in gendered_pairs if p1 in cov_df['IID'].values and p2 in cov_df['IID'].values]
    #making all pairs unique
    gendered_pairs = list(set(gendered_pairs))
    #ever_married = list(df_links[['IID','RASPIID1']].dropna().values)+list(df_links[['IID','RASPIID2']].dropna().values)
    #IDs = set(list(cov_df['IID'].values))
    #married_pairs = np.array([[p1,p2] for p1,p2 in ever_married if p1 in IDs and p2 in IDs])
    return gendered_pairs

def combine_pca_dem(directory,max_pca,ignore_pca=True):
    # demographics, ignore pca
    file = '/project/arpawong_181/HRS_AsMa/keith/'+'hrs_pcs15506_cog27.csv'#covariates.csv'
    # remove missing data
    df_dem = pd.read_csv(file)
    df_dem = df_dem.drop(columns='cogtot')
    df_dem = df_dem.drop(columns='bmi')
    df_dem = df_dem.dropna().replace('M',np.nan).replace('R',np.nan).replace('D',np.nan).replace('S',np.nan).replace('N',np.nan).astype(float)
    # demographic features
    dem_cols = [col for col in df_dem.columns if 'EV' not in col]
    # ignore PCA components... (see next line)
    #if ignore_pca:
    df_dem = df_dem[dem_cols]

    # pcas are collected here instead
    df_pca = pd.read_csv(directory+'phased_data/all_chr_gds_evecs.csv',header=None)
    # rename pca cols
    df_pca.columns=['EV'+str(i) for i in range(1,max_pca+1)]
    # eigenvalues
    df_evals = pd.read_csv(directory+'phased_data/all_chr_gds_evals.csv',header=None)
    df_evals.columns=['evals']
    # de-normalize by multiplying pca components by sqrt(evals)
    df_evals = np.sqrt(df_evals['evals'].values[:max_pca])
    for i in range(1,max_pca+1):
        df_pca['EV'+str(i)]=df_pca['EV'+str(i)].values*df_evals[i-1]

    # map pca values to ids
    df_ids=pd.read_csv(directory+'phased_data/all_ids.csv',header=None)
    df_ids.columns=['ids']
    df_pca['IID']=df_ids['ids'].values

    # merge PCAs and Demographics
    df_pca = pd.merge(df_dem,df_pca,on='IID')

    return df_pca

def pair_features(cov_df,cov_cols,category_cols,pairs):
    X = []
    feature_cols = []
    for col in cov_cols:
        print(col)
        # all data is heterosexual so this is not needed
        if col == 'sex': continue
        # binarize categories
        if col in category_cols:
            unique_cats = sorted(list(cov_df[col].drop_duplicates().values))
            if len(unique_cats)>2:
                for f in unique_cats:
                    feature_cols.append(col+'_'+str(f))
        else:
            feature_cols.append(col)
        features=[]
        for p1,p2 in pairs:
            f1 = cov_df.loc[cov_df['IID']==p1,col].values[0]
            f2 = cov_df.loc[cov_df['IID']==p2,col].values[0]
            if col in category_cols:
                if f1==f2:
                    features.append(f1)
                else:
                    features.append(-1)	
            else:
                if col == 'iearn':
                    f1 = np.log10(f1+1)
                    f2 = np.log10(f2+1)
                features.append(f1-f2)
        # convert categories into binary variables
        if col in category_cols:
            for f in unique_cats:
                X.append([int(f==feature) for feature in features])
        else:
            X.append(features)
    features = pca_distance(cov_df,pairs)
    X.append(features)
    feature_cols.append('pca_distance')
    X = np.transpose(np.array(X))
    return feature_cols,X


def sort_unique_sublist(list):
    return np.unique(np.array([tuple(sorted(li)) for li in list]))





def bootstrap_mean(l,alpha=0.05):
    n_boot = 1000
    if n_boot*len(l) < 10**8:
        # if data is not too big, create more bootrapped data
        n_boot = 10000
    
    means = np.mean(np.random.choice(l,replace=True,size=(n_boot,len(l))),axis=0)
    # 2 sided p-value
    p_val = np.min([len(means[means>=0])/len(means),len(means[means<=0])/len(means)])*2
    quantile=np.quantile(means,[alpha/2,1-alpha/2])
    return p_val,quantile
def standard_errors(l):
    mean = np.mean(l)
    standard_error = np.std(l)/np.sqrt(len(l))
    z_score = np.abs(mean/standard_error)
    p_val = scipy.stats.norm.sf(abs(z_score))*2
    # 95% CI
    quantile=[mean-1.96*standard_error,mean+1.96*standard_error]#np.quantile(means,[alpha/2,1-alpha/2])
    return [p_val,quantile]



def feature_cov_mat(cov_df,cov_cols,category_cols):
    
    X = []
    feature_cols = []
    IDs = cov_df['IID'].values
    for col in cov_cols:
        # all data is heterosexual so this is not needed
        if col == 'sex': continue
        # binarize categories
        if col in category_cols:
            unique_cats = sorted(list(cov_df[col].drop_duplicates().values))
            if len(unique_cats)>2:
                for f in unique_cats:
                    feature_cols.append(col+'_'+str(f))
        else:
            feature_cols.append(col)
        features=[]
        
        features = cov_df[col].values
        # convert categories into binary variables
        if col in category_cols:
            for f in unique_cats:
                X.append([int(f==feature) for feature in features])
        else:
            X.append(features)
    X = np.transpose(np.array(X))
    return feature_cols,X

def create_correl_mat(cov_df,cov_cols,category_cols,correl_file):
    # look at correlations between each col
    feature_cols,X = feature_cov_mat(cov_df,cov_cols,category_cols)
    correl_mat = np.array([[0.0]*len(feature_cols)]*len(feature_cols))
    correl_mat = []
    for ii,f1 in enumerate(feature_cols):
        line = []
        for jj,f2 in enumerate(feature_cols):
            s,p = spearmanr(X[:,ii].flatten(),X[:,jj].flatten())
            line.append(s)
        correl_mat.append(line)
    correl_mat = np.array(correl_mat)
    correl_mat2 = {f:correl_mat[:,ii] for ii,f in enumerate(feature_cols)}
    
    pd.DataFrame(correl_mat2).to_csv(correl_file,index=False)


def main():

    max_pca = 32
    directory = '/project/burghard_687/genetic_data/'
    X_rand_file = directory + 'X_rand_cog27.npy'
    X_married_file = directory + 'X_married_cog27.npy'
    feature_cols_file = directory + 'feature_cols_cog27.npy'
    rp_file = directory + 'random_pairs_cog27.npy'
    mp_file = directory+'married_pairs_cog27.npy'
    
    cov_df = combine_pca_dem(directory,max_pca)
    married_pairs = find_marriages(cov_df)
    random.shuffle(married_pairs)
    print('marriages: ',len(married_pairs))
    np.save(mp_file,married_pairs)
    df_links = pd.read_csv(directory+'hrs_spouse_IIDs.csv')  
    IDs = cov_df['IID'].values
    cov_cols = [col for col in cov_df.columns if col not in ['FID', 'IID','nhw','nhb','nho','hisp','RAHISPAN','RARACEM']]
    for pca_dim in range(21,40):
        if 'EV'+str(pca_dim) in cov_cols:
            cov_cols.remove('EV'+str(pca_dim))
    print('cols made')
    # remove dirty data
    # columns that are too correlated
    category_cols = ['RABPLACE','RARELIG','ethnicity','sex']
    cov_df = cov_df[['IID']+cov_cols].dropna()
    correl_file = '../parsed_data/FeatureCorrelations_cog27.csv'
    if not os.path.exists(correl_file):
        create_correl_mat(cov_df,cov_cols,category_cols,correl_file)
    num_nodes = len(cov_df['IID'].dropna().drop_duplicates())
    tot_num_pairs = int(num_nodes*(num_nodes-1)/2)
    np_married_pairs = np.array([[p1,p2] for p1,p2 in married_pairs])
    if not os.path.exists(rp_file):
        random_pairs = [tuple([p1,p2]) for p1 in np_married_pairs[:,0] for p2 in np_married_pairs[:,1] if tuple([p1,p2]) not in married_pairs]#np.load(rp_file).astype(int)
        random_pairs = list(set(random_pairs))
        random.shuffle(random_pairs)
        np.save(rp_file,random_pairs)
    else:
        random_pairs = np.load(rp_file)
    if not os.path.exists(X_married_file):
        feature_cols,X_married = pair_features(cov_df,cov_cols,category_cols,married_pairs)
        print(['X_Married',len(X_married)])
        np.save(feature_cols_file,feature_cols)
        np.save(X_married_file,X_married)
    if not os.path.exists(X_rand_file):
        print('creating new X array')
        feature_cols,X_rand = pair_features(cov_df,cov_cols,category_cols,random_pairs)
        print('saving data...')
        np.save(X_rand_file,X_rand)
        np.save(feature_cols_file,feature_cols)
        print(feature_cols)

if __name__ == "__main__":
    main()


