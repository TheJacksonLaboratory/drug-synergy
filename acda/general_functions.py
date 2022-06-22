import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import pymysql

def encodeNames(df):
    
    '''A.k.a. One-hot encoding
    ['MODEL', 'DRUG1', 'DRUG2'] should be present either in index levels or in the columns.
    The idea is from the AstraZeneca DREAM challenge second-best winning method for drug synergy prediction.
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


def fetchMySQL(s, host="host", database="database", user="user", password="password"):
    
    Database = pymysql.connect(host=host, database=database, user=user, password=password)

    if Database==None:
        print("Could not establish connection.")

        return

    DatabaseCursor = Database.cursor()

    try:
        DatabaseCursor.execute(s)
        Table = DatabaseCursor.fetchall()

    except Exception as exception:
        print(exception)
        print ("Error: unable to fetch data")

    Database.close()

    return Table


def getGeneToProteinNameAssociation(host="genome-mysql.cse.ucsc.edu", database="uniProt", user="genomep", password="password"):

    '''
     Tables in uniProt:
     ['accToKeyword', 'accToTaxon', 'author', 'bigFiles', 'citation', 'citationRc', 'citationRp', 'comment', 
     'commentType', 'commentVal', 'commonName', 'description', 'displayId', 'extDb', 'extDbRef', 'feature', 
     'featureClass', 'featureId', 'featureType', 'gene', 'geneLogic', 'history', 'info', 'keyword', 'organelle', 
     'otherAcc', 'pathogenHost', 'protein', 'proteinEvidence', 'proteinEvidenceType', 'rcType', 
     'rcVal', 'reference', 'referenceAuthors', 'tableDescriptions', 'tableList', 'taxon', 'varAcc', 'varProtein']

     To see all columns in a table:
     fetchMySQL("SHOW COLUMNS FROM my_table;")

     Usage:
        assoc = getGeneToProteinNameAssociation()
    '''

    print('Query to genome-mysql.cse.ucsc.edu ...')

    SQLquery = "SELECT gene.val, description.val FROM accToTaxon \
                LEFT JOIN gene ON accToTaxon.acc=gene.acc \
                LEFT JOIN description ON accToTaxon.acc=description.acc \
                WHERE accToTaxon.taxon=9606;"

    df = pd.DataFrame(fetchMySQL(SQLquery, host=host, database=database, user=user, password=password))
    df = df.set_index(0)[1].apply(lambda s: s.split(';')[0].split('RecName: ')[-1].split('Full=')[1])
    df = df[~df.index.isna()].dropna()
    df = df.reset_index().set_index(1)[0].groupby(level=0).unique()

    df.index.name = None
    df.name = None
    
    return df


def convertProteinNamesToGenes(se, dfgp):
    
    seConv = se.drop_duplicates()
    seConv.index = seConv.values.copy()

    preG = seConv.str.split(';').apply(lambda l: [s.strip() for s in l]).apply(lambda l: dfgp.reindex(np.unique(l)).dropna().values)
    preG = preG.apply(lambda l: [v for sub in l for v in sub]).apply(lambda l: np.unique(l).tolist())
    seConv[:] = preG.values

    se[:] = seConv.loc[se.values].values
    
    return


def centerDf(df):
    
    ''' similar to the StandardScaler:
    from sklearn.preprocessing import StandardScaler
    StandardScaler().fit(df).transform(df)
    '''
    
    dft = df.copy()
    dft /= dft.std(axis=0)
    dft -= dft.mean(axis=0)

    return dft
