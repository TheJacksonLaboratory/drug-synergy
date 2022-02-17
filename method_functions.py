import os
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, cophenet


def getCDAdrugDistance(df_drugs_celllines, method='pearson'):

    '''Excludes pairs of points with missing data
    method: {'pearson', 'kendall', 'spearman', 'cosine'}
    '''
    
    # Project drug response vectors (n drugs, m cell lines) AUCs to cos⁡𝛼
    # space --> A
    def cosA(u, v):
        wh = (~np.isnan(u)) & (~np.isnan(v))
        return (u[wh] * v[wh]).sum() / np.sqrt((u[wh] * u[wh]).sum() * (v[wh] * v[wh]).sum())

    if method == 'cosine':
        method = cosA
    elif method not in ['pearson', 'kendall', 'spearman']:
        raise NotImplementedError
        
    A = df_drugs_celllines.corr(method=method)

    # Calculate Euclidean distance of A --> M
    def euclideanExcludeMissing(u, v):
        wh = ~np.isnan(u * v)
        return np.sqrt(((u[wh] - v[wh]) ** 2).sum())

    M = A.copy()
    M[:] = squareform(pdist(M.values, metric=euclideanExcludeMissing))

    # Hierarchically cluster M, with Ward linkage --> Linkage Z (for later use
    # in plots)
    Z = linkage(squareform(np.nan_to_num(M.values, nan=M.max().max())), method='ward', optimal_ordering=True)

    # Cophenetic distance --> C
    dfC = A.copy()
    dfC[:] = squareform(cophenet(Z))
    
    return dfC, Z


def getDistanceAndSensitivityData(df, method='pearson', sensitivity_metric='LNIC50', drug='DRUG'):
    
    """ Sensitivity data and CDA drug distance
    dfC, dfS, Z = getDistanceAndSensitivityData(df_drug_sensitivity_GDSC1)
    """
    
    dfS = df[sensitivity_metric].groupby(level=[0,1]).mean().unstack(drug).sort_index()
    dfS.columns = dfS.columns.str.lower().str.strip()
    
    dfC, Z = getCDAdrugDistance(dfS, method=method)

    return dfC, dfS, Z


def prepDfTa(df_drug_sensitivity, se_models_mutations, se_drug_targets, dname):
    
    """
    Drug targeting mutated genes
    """
    cacheFile = 'data/dfTa_%s.pklz' % dname
    if not os.path.isfile(cacheFile):
        models = df_drug_sensitivity.index.get_level_values('MODEL').unique().intersection(se_models_mutations.index)
        drugs = df_drug_sensitivity.index.get_level_values('DRUG').unique().intersection(se_drug_targets.index)
        dfT = pd.DataFrame(index=models, columns=drugs, data=False, dtype=bool)
        for imodel, model in enumerate(dfT.index):
            if imodel % 50 == 0:
                print(imodel, end=' ')
            for drug in dfT.columns:
                dfT.loc[model, drug] = np.isin(se_models_mutations[model], se_drug_targets[drug]).any()

        dfTa = pd.DataFrame(index=dfT.index, data=(dfT.values[:,None,:] + dfT.values[:,:,None]).reshape(dfT.shape[0],-1).astype(bool),
                            columns=pd.MultiIndex.from_product([dfT.columns, dfT.columns]))
        dfTa.columns = pd.MultiIndex.from_arrays([dfTa.columns.get_level_values(0).str.lower().str.strip(), dfTa.columns.get_level_values(1).str.lower().str.strip()])
        if not os.path.exists(os.path.dirname(cacheFile)):
            os.makedirs(os.path.dirname(cacheFile))
        dfTa.to_pickle(cacheFile)
        print('Saving to cache...')
    else:
        print('Loading from cache...')
        dfTa = pd.read_pickle(cacheFile)

    print('%s cell lines, %s drug pairs' % dfTa.shape)    
    
    return dfTa


def filter_synergy_pairs_list(df, dfTa, pref='Synergy pairs'):
    
    """
    Get pairs that overlap with sensitivity-mutations-targets data
    
    Usage:
    dfKS = filter_synergy_pairs_list(df_drug_synergy_Narayan, dfTa)
    """
        
    dfKS = df.copy().reset_index(['DRUG1', 'DRUG2']).apply(lambda s: (s[0], s[1]), axis=1)
    dfKS = dfKS.iloc[np.where((~dfTa.reindex(dfKS.values, axis=1).iloc[0].isna()).values)]
    dfKS = dfKS.loc[dfKS.index.intersection(dfTa.index)]
    
    dfKS = dfKS.apply(pd.Series).set_index([0, 1], append=True)
    dfKS.index.names = ['MODEL', 'DRUG1', 'DRUG2']
    dfKS['SYNERGY_SCORE'] = 1.0
    dfKS = dfKS['SYNERGY_SCORE']
    dfKS = dfKS.index.unique()
    print('%s in:\t' % pref, df.shape[0])
    print('%s out:\t' % pref, dfKS.shape[0])
    
    return df.loc[~df.index.duplicated()].reindex()


def prepDfTas(dfTa, dfKS, se_tissue_annotation, tissue=None):
    
    """
    Reshape dfTa, add tissue and synergy information
    """

    addtissue = lambda ind: pd.MultiIndex.from_frame(ind.to_frame().assign(TISSUE=se_tissue_annotation.reindex(ind.get_level_values('MODEL')).fillna('Unknown').values).astype('category'))
    
    dfKS = filter_synergy_pairs_list(dfKS, dfTa)
    
    synergeticIndex = dfKS.index
    synergeticIndex.names = ['MODEL', 'DRUG1', 'DRUG2']
    synergeticIndex = addtissue(synergeticIndex)

    dfTas = dfTa.copy()
    dfTas.index = addtissue(dfTas.index)
    dfTas = dfTas.stack().stack().reorder_levels([0,2,3,1])
    dfTas.index.names = ['MODEL', 'DRUG1', 'DRUG2', 'TISSUE']
    print('Data size:', dfTas.shape[0])

    chsize = 10 ** 7
    lsize = len(dfTas.index)
    whfound = []
    for i in range(int(np.ceil(lsize / chsize))):
        indices = i * chsize, min((i + 1) * chsize, lsize)
        found = indices[0] + np.where(dfTas.index[indices[0]:indices[1]].droplevel('TISSUE').isin(synergeticIndex.droplevel('TISSUE')))[0]
        print(i, end=' ')
        whfound.append(found)
    whfound = np.hstack(whfound)
    dfTas_synergy = dfTas.iloc[whfound]
    if not tissue is None:
        dfTas_synergy = dfTas_synergy.xs(tissue, level='TISSUE', drop_level=False)
    print()
    print('With known synergy:', dfTas_synergy.shape[0])
    
    dfTas_synergy = dfTas_synergy.to_frame()
    dfTas_synergy.columns = ['Tijk']
    dfTas_synergy = dfTas_synergy.reset_index('TISSUE')
    dfTas_synergy['SYNERGY_SCORE'] = dfKS.reindex(dfTas_synergy.index).values
    dfTas_synergy = dfTas_synergy.set_index('TISSUE', append=True)
        
    if not tissue is None:
        dfTas = dfTas.xs(tissue, level='TISSUE', drop_level=False)
    else:
        np.random.seed(0)
        dfTas = dfTas.sample(n=min(10 ** 7, dfTas.shape[0]))

    dfTas = dfTas.iloc[np.where(~dfTas.index.isin(synergeticIndex))[0]]
    dfTas_no_synergy = dfTas
    print('With unknown synergy:', dfTas_no_synergy.shape[0])
    
    dfTas_no_synergy = dfTas_no_synergy.to_frame()
    dfTas_no_synergy.columns = ['Tijk']
    dfTas_no_synergy['SYNERGY_SCORE'] = np.nan
    
    dfTas = pd.concat([dfTas_synergy, dfTas_no_synergy], axis=0)

    return dfTas[['SYNERGY_SCORE', 'Tijk']]


def split_train_test_validate_predict(df, factor=1 / 2, random_state=None):
    
    print()
    print('All pairs:\t\t', df.shape[0])
    print('Pairs with missing:\t', ((df['Cij'].isna()) | (df['Sik'].isna()) | (df['Sjk'].isna())).sum())
    
    df_temp = df.loc[(~df['Cij'].isna()) & (~df['Sik'].isna()) & (~df['Sjk'].isna()) & (~df['SYNERGY_SCORE'].isna())]
    print('All with known pairs:\t', df_temp.shape[0])
    
    df_train_test = df_temp.sample(n=int(df_temp.shape[0] * factor), random_state=random_state)
    print('Training-testing pairs:\t', df_train_test.shape[0])
    
    df_validate = df_temp.loc[df_temp.index.difference(df_train_test.index)]
    print('Validation pairs:\t', df_validate.shape[0])
    
    df_predict = df.loc[(~df['Cij'].isna()) & (~df['Sik'].isna()) & (~df['Sjk'].isna()) & (df['SYNERGY_SCORE'].isna())]
    print('Prediction pairs:\t', df_predict.shape[0])
    
    return df_train_test, df_validate, df_predict


def centerDf(df):
    
    ''' similar to the StandardScaler:
    from sklearn.preprocessing import StandardScaler
    StandardScaler().fit(df).transform(df)
    '''
    
    dft = df.copy()
    dft /= dft.std(axis=0)
    dft -= dft.mean(axis=0)

    return dft

