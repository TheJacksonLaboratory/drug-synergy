import os
import pandas as pd
import numpy as np
import pickle
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LogisticRegression, LinearRegression

from .general_functions import encodeNames


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


def getDistanceAndSensitivityData(df, method='pearson', sensitivity_metric='LNIC50', drug='DRUG', lower=True):
    
    """ Sensitivity data and CDA drug distance
    dfC, dfS, Z = getDistanceAndSensitivityData(df_drug_sensitivity_GDSC1)
    """
    
    dfS = df[sensitivity_metric].groupby(level=[0,1]).mean().unstack(drug).sort_index()
    
    if lower:
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


def makeCDAformattedData(tissue, se_drug_synergy, se_tissue_annotation, df_drug_sensitivity, se_models_mutations, se_drug_targets, synergy_cutoff, dataset, sensitivityMethod='cosine', sensitivity_metric='LNIC50', returnMore=False):
    
    dfC, dfS, Z = getDistanceAndSensitivityData(df_drug_sensitivity, method=sensitivityMethod, sensitivity_metric=sensitivity_metric)
    dfTa = prepDfTa(df_drug_sensitivity, se_models_mutations, se_drug_targets, dataset)
    dfTas = prepDfTas(dfTa, se_drug_synergy, se_tissue_annotation, tissue=tissue)

    dfTas['Cij'] = dfC.stack().reindex(pd.MultiIndex.from_frame(dfTas.index.to_frame()[['DRUG1', 'DRUG2']])).values
    dfTas['Sik'] = dfS.stack().reindex(pd.MultiIndex.from_frame(dfTas.index.to_frame()[['MODEL', 'DRUG1']])).values
    dfTas['Sjk'] = dfS.stack().reindex(pd.MultiIndex.from_frame(dfTas.index.to_frame()[['MODEL', 'DRUG2']])).values

    if returnMore:
        return dfTas, dfC, dfS, Z
    
    return dfTas


def prepareFromDCforCDA(df_DC, study, tissue, se_models_mutations, sample=None, random_state=None, sensitivity_measure='AUC', synergy='SYNERGY_ZIP', returnMore=False):
    o = se_models_mutations.index.intersection(df_DC.index.to_frame()['MODEL'].unique())
    df_temp = df_DC.loc[df_DC.index.to_frame()['MODEL'].isin(o)].xs(tissue, level='TISSUE', drop_level=False).xs(study, level='STUDY')
    
    if not sample is None:
        df_temp = df_temp.sample(min(sample, df_temp.shape[0]), random_state=random_state)
    
    df_temp = df_temp[[synergy, '%s_DRUG1' % sensitivity_measure, '%s_DRUG2' % sensitivity_measure, 'DRUG1_TARGETS', 'DRUG2_TARGETS']]
    df_temp['MUTATIONS'] = se_models_mutations.loc[df_temp.index.get_level_values('MODEL')].values
    df_temp = df_temp.rename({synergy: 'SYNERGY_SCORE', '%s_DRUG1' % sensitivity_measure: 'Sik', '%s_DRUG2' % sensitivity_measure: 'Sjk'}, axis=1)
    df_temp['Tijk'] = df_temp.apply(lambda s: len(np.intersect1d(s['DRUG1_TARGETS'] + s['DRUG2_TARGETS'], s['MUTATIONS']))>0, axis=1)
    
    assert df_temp.index.names==['MODEL', 'DRUG1', 'DRUG2', 'TISSUE'], 'Incorrect index levels names'
    
    df_temp_i = df_temp.reorder_levels([0,2,1,3])
    df_temp_i.index.names = df_temp.index.names
    
    temp = df_temp_i['Sik'].copy()
    df_temp_i['Sik'] = df_temp_i['Sjk'].copy()
    df_temp_i['Sjk'] = temp
    
    temp = df_temp_i['DRUG1_TARGETS'].copy()
    df_temp_i['DRUG1_TARGETS'] = df_temp_i['DRUG2_TARGETS'].copy()
    df_temp_i['DRUG2_TARGETS'] = temp
    
    df_temp_i = pd.concat([df_temp, df_temp_i])
    df_temp = df_temp_i.loc[~df_temp_i.index.duplicated()]
    
    se = pd.concat([df_temp.reset_index().set_index(['MODEL', 'DRUG1'])['Sik'], df_temp.reset_index().set_index(['MODEL', 'DRUG2'])['Sjk']], axis=0)
    se.index.names = ['MODEL', 'DRUG']
    df_drug_sensitivity = se.groupby(level=[0, 1]).agg(np.median).to_frame().rename({0: 'AUC'}, axis=1)
    dfC, dfS, Z = getDistanceAndSensitivityData(df_drug_sensitivity, method='cosine', sensitivity_metric='AUC', lower=False)

    df_temp['Cij'] = dfC.stack().reindex(pd.MultiIndex.from_frame(df_temp.index.to_frame()[['DRUG1', 'DRUG2']])).values
    
    df_temp = df_temp[['SYNERGY_SCORE', 'Tijk', 'Cij', 'Sik', 'Sjk']]

    if returnMore:
        return df_temp, dfC, dfS, Z
    
    return df_temp


def prepareFromAZforCDA(dir, se_tissue_annotation, se_models_mutations, se_drug_targets, fname='oi_combinations_synergy_scores_final.txt', sensitivity_measure='LNIC50'):
    
    '''{'IC50', 'LNIC50', 'H'}
    '''

    df = pd.read_csv(dir + fname, delimiter='\t', index_col=[0,1,2])
    
    temp = df.index.to_frame()
    temp.columns = ['MODEL', 'DRUG1', 'DRUG2']
    temp['TISSUE'] = se_tissue_annotation.reindex(temp['MODEL'].unique()).fillna('NA').loc[temp['MODEL']].values

    df.index = pd.MultiIndex.from_frame(temp)
       
    df = df.rename({'IC50_A': 'IC50_DRUG1', 'IC50_B': 'IC50_DRUG2', 'H_A': 'H_DRUG1', 'H_B': 'H_DRUG2'}, axis=1)
    df = df.drop(['MAX_CONC_A', 'MAX_CONC_B', 'Einf_A', 'Einf_B', 'QA'], axis=1)

    df.loc[:, 'LNIC50_DRUG1'] = np.log(df.loc[:, 'IC50_DRUG1']).values
    df.loc[:, 'LNIC50_DRUG2'] = np.log(df.loc[:, 'IC50_DRUG2']).values
    
    df['MUTATIONS'] = se_models_mutations.reindex(df.index.get_level_values('MODEL').unique()).apply(lambda s: np.array([]) if s is np.nan else s).loc[df.index.get_level_values('MODEL')].values    
    
    df = df.rename({'%s_DRUG1' % sensitivity_measure: 'Sik', '%s_DRUG2' % sensitivity_measure: 'Sjk'}, axis=1)
    
    
    df['DRUG1_TARGETS'] = se_drug_targets.reindex(df.index.get_level_values('DRUG1').unique()).apply(lambda s: np.array([]) if s is np.nan else s).loc[df.index.get_level_values('DRUG1')].values
    df['DRUG2_TARGETS'] = se_drug_targets.reindex(df.index.get_level_values('DRUG2').unique()).apply(lambda s: np.array([]) if s is np.nan else s).loc[df.index.get_level_values('DRUG2')].values
    
    df['Tijk'] = df.apply(lambda s: len(np.intersect1d(np.union1d(s['DRUG1_TARGETS'], s['DRUG2_TARGETS']), s['MUTATIONS']))>0, axis=1)
    
    
    assert df.index.names==['MODEL', 'DRUG1', 'DRUG2', 'TISSUE'], 'Incorrect index levels names'
    
    df_temp_i = df.reorder_levels([0,2,1,3])
    df_temp_i.index.names = df.index.names
    
    temp = df_temp_i['Sik'].copy()
    df_temp_i['Sik'] = df_temp_i['Sjk'].copy()
    df_temp_i['Sjk'] = temp
    
    temp = df_temp_i['DRUG1_TARGETS'].copy()
    df_temp_i['DRUG1_TARGETS'] = df_temp_i['DRUG2_TARGETS'].copy()
    df_temp_i['DRUG2_TARGETS'] = temp
    
    df_temp_i = pd.concat([df, df_temp_i])
    df = df_temp_i.loc[~df_temp_i.index.duplicated()]
    
    se = pd.concat([df.reset_index().set_index(['MODEL', 'DRUG1'])['Sik'], df.reset_index().set_index(['MODEL', 'DRUG2'])['Sjk']], axis=0)
    se.index.names = ['MODEL', 'DRUG']
    df_drug_sensitivity = se.groupby(level=[0, 1]).agg(np.median).to_frame().rename({0: sensitivity_measure}, axis=1)
    dfC, dfS, Z = getDistanceAndSensitivityData(df_drug_sensitivity, method='cosine', sensitivity_metric=sensitivity_measure, lower=False)

    df['Cij'] = dfC.stack().reindex(pd.MultiIndex.from_frame(df.index.to_frame()[['DRUG1', 'DRUG2']])).values

    df = df[['SYNERGY_SCORE', 'Tijk', 'Cij', 'Sik', 'Sjk']]

    return df.sort_index()


def trainOneTestAnother(df_one, df_another, n=10, n_sample=0.5):

    res_temp = dict()
    for i in range(n):
        print(i)

        dfTas_copy = df_one.copy()
        dfTas1 = dfTas_copy.loc[(~dfTas_copy['Cij'].isna()) & (~dfTas_copy['Sik'].isna()) & (~dfTas_copy['Sjk'].isna()) & (~dfTas_copy['SYNERGY_SCORE'].isna())]
        if not n_sample is None:
            dfTas1 = dfTas1.sample(int(n_sample*dfTas1.shape[0]), random_state=i)    
        print(dfTas1.shape)

        dfTas_copy = df_another.copy()
        dfTas2 = dfTas_copy.loc[(~dfTas_copy['Cij'].isna()) & (~dfTas_copy['Sik'].isna()) & (~dfTas_copy['Sjk'].isna()) & (~dfTas_copy['SYNERGY_SCORE'].isna())]
        if not n_sample is None:
            dfTas2 = dfTas2.sample(int(n_sample*dfTas2.shape[0]), random_state=i)    
        print(dfTas2.shape)

        print('\n', 'M:')
        res_temp[('EN', i)] = testCase(dfTas1, dfTas2, encode=True, useAllFeatures=False, cv=None)
        print('\n', 'ACDA:')
        res_temp[('ACDA', i)] = testCase(dfTas1, dfTas2, encode=False, useAllFeatures=False, cv=None)
        print('\n', 'CDA:')
        res_temp[('CDA', i)] = testCase(dfTas1, dfTas2, encode=False, useAllFeatures=False, cv=None, clf=LinearRegression())
        print('\n', 'EN-ACDA:')
        res_temp[('EN-ACDA', i)] = testCase(dfTas1, dfTas2, encode=True, useAllFeatures=True, cv=None)

    df_temp = pd.DataFrame(res_temp).T['val']
    df_temp2 = pd.concat([df_temp.groupby(level=0).mean(), df_temp.groupby(level=0).sem()], axis=1)
    
    return df_temp2, df_temp


def fit_validate_predict(inData, inSynergy, extData=None, extSynergy=None, cv=None, max_iter=10**4, clf=None, **kwargs):
    
    if clf is None:
        clf = RandomForestRegressor(random_state=0, max_depth=40, n_estimators=250)
        
    clf.fit(inData, inSynergy.astype(float))
    
    res = dict()
    pr = clf.predict(inData)
    res.update({'self': np.round(pearsonr(inSynergy.astype(float), pr)[0], 5)})
    print('On self:\t\t', res['self'])
    
    if not cv is None:
        def getPartitions(df, se, n_splits=4, shuffle=True, seed=None):
            np.random.seed(seed)
            df_ = pd.concat([df, se], axis=1)
            splitter = KFold(n_splits=n_splits, shuffle=shuffle)
            sp = [(df_.T.iloc[:, i].values, df_.T.iloc[:, j].values) for (i, j) in list(splitter.split(df_))]
            return [((u[:-1, :], v[:-1, :]), (u[-1:, :], v[-1:, :])) for u, v in sp]

        partitions = getPartitions(pd.DataFrame(inData), pd.Series(inSynergy), n_splits=cv, shuffle=True, seed=41) 

        scores = []
        for ((df_train, df_test), (se_train, se_test)) in partitions:
            clf.fit(df_train.T, se_train[0].astype(float))
            predicted = clf.predict(df_test.T)
            score = np.round(pearsonr(se_test[0], predicted)[0], 3)
            scores.append(score)

        res.update({'cv_mean': np.mean(scores)})
        res.update({'cv_std': np.std(scores)})

        print('Cross-validation:\t mean=%.3f std=%.3f' % (np.mean(scores), np.std(scores)))
        
        clf.fit(inData, inSynergy.astype(float))
    else:
        res.update({'cv_mean': np.nan})
        res.update({'cv_std': np.nan})
    
    if (not extData is None) and (not extSynergy is None):
        predicted = clf.predict(extData)
        res.update({'val': np.round(pearsonr(extSynergy.astype(float), predicted)[0], 3)})
        print('Validation data:\t', res['val'])
        return predicted, res
    
    return clf


def tryExcept(func):

    def internal(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as exception:
            print('\tError message: %s\n' % (exception))

    internal.__name__ = func.__name__
    internal.__doc__ = func.__doc__

    return internal


@tryExcept
def testCase(df_train_test, df_validate, CDA_features=['Tijk', 'Cij', 'Sik', 'Sjk'], encode=False, useAllFeatures=False, cv=10, clf=None):
    
    if encode:
        dfe_train_test = encodeNames(df_train_test)
        dfe_train_test = dfe_train_test[dfe_train_test.columns.intersection(encodeNames(df_validate).columns)]
        dfe_validate = encodeNames(df_validate).reindex(dfe_train_test.columns, axis=1).fillna(0.).astype(int)
        
        df_train_test.loc[:, dfe_train_test.columns] = dfe_train_test.values
        df_validate.loc[:, dfe_validate.columns] = dfe_validate.values

        if useAllFeatures:
            # Exclude the predictee
            features = df_train_test.columns[~df_train_test.columns.isin(['SYNERGY_SCORE'])]
        else:
            # Exclude the predictee and CDA features, keep encoded drug and models
            features = df_train_test.columns[~df_train_test.columns.isin(['SYNERGY_SCORE'] + CDA_features)]
    else:
        features = CDA_features
        
    r, res = fit_validate_predict(df_train_test[features].values, df_train_test['SYNERGY_SCORE'].values,
                         extData=df_validate[features].values, extSynergy=df_validate['SYNERGY_SCORE'].values,
                         cv=cv, clf=clf)

    df_ro = df_validate['SYNERGY_SCORE'].to_frame()
    df_ro['Predicted'] = r
    
    try:
        print("ro normalized:\t\t %.3f" % ro_normalized(df_ro))
        res.update({'ro': ro_normalized(df_ro)})
    except:
        res.update({'ro': np.nan})

    return res


def MonteCarloCrossValidation(dfTas, n=10, sample_non_synergy=False, sample_non_synergy_size=100, clf_for_CDA=LogisticRegression(), deidentify=False):
    
    '''For continuous predictee: clf_for_CDA=LinearRegression()
    '''

    res_temp = dict()
    
    for i in range(n):
        print(i)
        dfTas_copy = dfTas.copy()
        
        if sample_non_synergy:
            dfTas_copy.loc[dfTas_copy[dfTas_copy['SYNERGY_SCORE'].isna()].sample(sample_non_synergy_size, random_state=i).index, 'SYNERGY_SCORE'] = 0.
            
        df_train_test, df_validate, df_predict = split_train_test_validate_predict(dfTas_copy, factor=2/3, random_state=i)
        
        if deidentify:
            f = df_validate.index.to_frame()
            f['DRUG1'] += '_D'
            f['DRUG2'] += '_D'
            df_validate.index = pd.MultiIndex.from_frame(f)

        print('\n', 'ACDA:')
        res_temp[('ACDA', i)] = testCase(df_train_test, df_validate, encode=False, useAllFeatures=False, cv=None)
        
        print('\n', 'CDA:')
        res_temp[('CDA', i)] = testCase(df_train_test, df_validate, encode=False, useAllFeatures=False, cv=None, clf=LinearRegression())

        print('\n', 'EN:')
        res_temp[('EN', i)] = testCase(df_train_test, df_validate, encode=True, useAllFeatures=False, cv=None)
        
        print('\n', 'EN-ACDA:')
        res_temp[('EN-ACDA', i)] = testCase(df_train_test, df_validate, encode=True, useAllFeatures=True, cv=None)

    df_temp1 = pd.DataFrame(res_temp).T.sort_index()
    df_temp2 = pd.concat([df_temp1.groupby(level=0).mean().round(3),
                          df_temp1.groupby(level=0).sem().round(3)],
                         keys=['mean', 'sem'], axis=1).xs('val', level=1, axis=1)
    
    return df_temp2, df_temp1


def downsampleRun(ra, dfTas=None, pref='temp', mid='temp', rep=10, basedir='output/', cv=None):
    
    '''df_sub10000 = downsampleRun(range(500, 10000, 500), dfTas=dfTas_AZ_breast, pref='', mid='AZ_breast')
    '''
    
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    for n in list(ra):
        for i in range(rep):
            fname = basedir + '%s%s_%s_%s.pklz' % (pref, mid, i, n)

            if not os.path.isfile(fname):
                print(n, i)
                res = dict()
                df_train_test, df_validate, df_predict = split_train_test_validate_predict(dfTas, factor=2/3, random_state=i)

                df_train_test_sample = df_train_test.sample(n, random_state=i)

                print('\n', 'EN:')
                res.update({('EN', n, i): testCase(df_train_test_sample, df_validate, encode=True, useAllFeatures=False, cv=None)})

                print('\n', 'ACDA:')
                res.update({('ACDA', n, i): testCase(df_train_test_sample, df_validate, encode=False, useAllFeatures=False, cv=None)})

                print('\n', 'CDA:')
                res.update({('CDA', n, i): testCase(df_train_test_sample, df_validate, encode=False, useAllFeatures=False, cv=None, clf=LinearRegression())})

                print('\n', 'EN-ACDA:')
                res.update({('EN-ACDA', n, i): testCase(df_train_test_sample, df_validate, encode=True, useAllFeatures=True, cv=None)})

                with open(fname, 'wb') as outfile:
                    pickle.dump(res, outfile)
                
    res = dict()
    for n in list(ra):
        for i in range(rep):
            fname = basedir + '%s%s_%s_%s.pklz' % (pref, mid, i, n)

            if os.path.isfile(fname):
                with open(fname, 'rb') as outfile:
                    temp = pickle.load(outfile)
                    res.update(temp)

    df = pd.DataFrame(res).T

    df.to_csv(basedir + 'downSample_%s%s_%s_%s.csv' % (pref, mid, i, n))
                
    return df
