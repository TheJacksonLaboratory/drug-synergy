import numpy as np
import pandas as pd
from scipy.stats import pearsonr

def encodeNames(df):
    
    '''A.k.a. One-hot encode
    '''

    dfs = [pd.get_dummies(df.reset_index()[col]) for col in ['MODEL', 'DRUG1', 'DRUG2']]
    
    return dfs[0].join((dfs[1] + dfs[2]).fillna(0))

def ro_normalized(df):
    
    '''Takes a dataframe where index has 3 levels (MODEL, DRUG1, DRUG2),
    and there are two columns, first with ground truth scores, second with predicted scores.
    
    Output is the average pearson correlation coefficient for each drug pair
    weighted by the number of models for which a pair was measured against.
    
    See definitions in: https://www.synapse.org/#!Synapse:syn4231880/wiki/235660
    '''
    
    s = df.dropna().groupby(level=[1,2]).agg(list)

    def pearsonCorr(u, v):
       
        return ((v - v.mean())*(u - u.mean())).sum() / np.sqrt(((v - v.mean())**2).sum() * ((u - u.mean())**2).sum())

    ro = s.apply(lambda se: (np.sqrt(len(se[0])-1)*pearsonr(se[0], se[1])[0]) if len(se[0])>=2 else np.nan, axis=1).sum() /\
         s.apply(lambda se: (np.sqrt(len(se[0])-1)*max(pearsonr(se[0], se[1])[0], 1.)) if len(se[0])>=2 else np.nan, axis=1).sum()

    return ro

def loadAZchallengeIndex(se_synergy, dirAZpharmacol):
    
    def loadOne(fname):
        
        df = pd.read_csv(fname, delimiter=',', index_col=0)[['COMPOUND_A', 'COMPOUND_B']]
        
        df.loc[:, 'COMPOUND_A'] = df['COMPOUND_A'].str.lower()
        df.loc[:, 'COMPOUND_B'] = df['COMPOUND_B'].str.lower()

        df = df.set_index(['COMPOUND_A', 'COMPOUND_B'], append=True)
        df.index.names = ['MODEL', 'DRUG1', 'DRUG2']

        df_ = df.reorder_levels([0,2,1])
        df_.index.names = df.index.names
        df = pd.concat([df, df_], axis=0)
        
        index = df.index.unique().sort_values()

        print(index.shape[0], fname)
        
        return index

    se = se_synergy.copy()
    se.iloc[:] = np.nan

    se.loc[se.index.intersection(loadOne(dirAZpharmacol + 'ch1_train_combination_and_monotherapy.csv'))] = 'ch1_train'
    se.loc[se.index.intersection(loadOne(dirAZpharmacol + 'ch1_test_monotherapy.csv'))] = 'ch1_test'
    se.loc[se.index.intersection(loadOne(dirAZpharmacol + 'ch1_lb.csv'))] = 'ch1_leaderboard'
    se.loc[se.index.intersection(loadOne(dirAZpharmacol + 'ch2_test_monotherapy.csv'))] = 'ch2_test'
    se.loc[se.index.intersection(loadOne(dirAZpharmacol + 'ch2_lb.csv'))] = 'ch2_leaderboard'
    
    vc = se.value_counts()
    print(vc.sum())
    print(vc)
    
    return se