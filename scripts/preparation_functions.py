import pandas as pd
import numpy as np

def prepGDSC(whichGDSC, dir, fname='{}_fitted_dose_response_25Feb20.csv.gz'):

    # CDA authors used the first generation of GDSC (MGH part)
    # Download drug screening data from:
    # https://www.cancerrxgene.org/downloads/bulk_download
    # GDSC1: Screened in 96 or 384-well plates with Cyto60 or Resazurin
    # endpoint
    # GDSC2 (recommended where available): 1536-well plates with CellTiterGlo
    # endpoint

    df = pd.read_csv(dir + fname.format(whichGDSC), index_col=0, header=0)
    df.loc[:, 'CELL_LINE_NAME'] = df.loc[:, 'CELL_LINE_NAME'].replace({'LS-1034': 'LS1034'})
    df_ = df[['CELL_LINE_NAME', 'DRUG_NAME', 'LN_IC50', 'AUC']]
    df_.columns = ['MODEL', 'DRUG', 'LNIC50', 'AUC']
    df_ = df_.set_index(['MODEL', 'DRUG'])

    return df_.groupby(level=[0, 1]).mean().sort_index()


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


def prepAZsensitivity(dir, fname='oi_combinations_synergy_scores_final.txt'):

    df_AZ = pd.read_csv(dir + fname, delimiter='\t', index_col=[0,1,2])

    def prepAZsensitivity_(df_AZ, a, b):
        df_ = df_AZ.reset_index('COMPOUND_%s' % b)[['IC50_%s' % a, 'H_%s' % a, 'Einf_%s' % a]]
        df_.columns = ['LNIC50', 'H', 'EINF']
        df_.index.names = ['MODEL', 'DRUG']
        return df_

    df = pd.concat([prepAZsensitivity_(df_AZ, 'A', 'B'), prepAZsensitivity_(df_AZ, 'B', 'A')], axis=0)
    df = df.groupby(level=[0, 1]).apply(lambda s: pd.Series(index=s.columns, data=np.median(s.values, axis=0)))
    #df = df.groupby(level=[0, 1]).mean()
    df.loc[:, 'LNIC50'] = np.log(df.loc[:, 'LNIC50']).values

    return df.sort_index()


def prepAZsynergy(dir, fname='oi_combinations_synergy_scores_final.txt'):
    
    df_ = pd.read_csv(dir + fname, delimiter='\t', index_col=[0,1,2])['SYNERGY_SCORE']
    
    df_.index.names = ['MODEL', 'DRUG1', 'DRUG2']
    df_ = df_.groupby(level=[0, 1, 2]).mean()
    df_ = df_.reset_index()
        
    # Lowercase drug names
    df_.loc[:, 'DRUG1'] = df_.loc[:, 'DRUG1'].str.lower().str.strip()
    df_.loc[:, 'DRUG2'] = df_.loc[:, 'DRUG2'].str.lower().str.strip()
    
    df_ = df_.set_index(['MODEL', 'DRUG1', 'DRUG2'])['SYNERGY_SCORE']
    
    dfw_ = df_.reorder_levels([0,2,1])
    dfw_.index.names = df_.index.names
    df_ = pd.concat([df_, dfw_], axis=0)
    df_ = df_.groupby(level=[0, 1, 2]).mean().sort_index()
    
    df_ = df_.loc[(df_.index.get_level_values('DRUG1') != df_.index.get_level_values('DRUG2'))]
    
    return df_


def prepNarayanSynergy(fname):
    
    df_ = pd.read_excel(fname, index_col=0, header=0)[['Drug1', 'Drug2', 'Drug3']]

    # Filter for dual synergy pairs only
    df_ = df_[df_['Drug3'].isna()][['Drug1', 'Drug2']]
    df_.columns = ['DRUG1', 'DRUG2']

    # Lowercase drug names
    df_.loc[:, 'DRUG1'] = df_.loc[:, 'DRUG1'].str.lower().str.strip()
    df_.loc[:, 'DRUG2'] = df_.loc[:, 'DRUG2'].str.lower().str.strip()

    # Curation to match GDSC cell lines names with the synergy data (by SD)
    cellLinesNamesCorrection = {'BxPC3': 'BxPC-3', 'FTC133': 'FTC-133', 'SKNDZ': 'SK-N-DZ', 'G361': 'G-361', 'HCT116': 'HCT-116', 'HL60': 'HL-60',
                                'HS578T': 'Hs-578-T', 'Hs578T': 'Hs-578-T', 'OVCAR3': 'OVCAR-3', 'CAOV3': 'Caov-3', 'HT1080': 'HT-1080', 'THP1': 'THP-1',
                                'SKNSH': 'SK-N-SH', 'SKOV-3': 'SK-OV-3', 'SKOV3': 'SK-OV-3', 'T-24': 'T24', 'THP-1': 'THP1', 'U266': 'U-266', 'U2OS': 'U-2-OS'}
    df_.index = pd.Index(pd.Series(df_.index).str.strip().replace(cellLinesNamesCorrection).values)
    df_.index.name = 'MODEL'
    
    df_ = df_.set_index(['DRUG1', 'DRUG2'], append=True)
    df_['SYNERGY_SCORE'] = 1
    df_ = df_['SYNERGY_SCORE'].astype(float)
    
    dfw_ = df_.reorder_levels([0,2,1])
    dfw_.index.names = df_.index.names
    df_ = pd.concat([df_, dfw_], axis=0)
    df_ = df_.groupby(level=[0, 1, 2]).mean().sort_index()
    
    return df_


def prepMutationsAZ(fname):
    
    # FATHMM.prediction {'PASSENGER/OTHER': 4436, 'CANCER': 658}
    
    se_models_mutations = pd.read_csv(fname, index_col=0)['cell_line_name'].dropna()
    se_models_mutations = se_models_mutations.reset_index().set_index('cell_line_name').groupby(level=0).apply(np.unique).sort_index()
    se_models_mutations.index.name = 'MODEL'
    se_models_mutations.name = 'MUTATIONS'
    
    return se_models_mutations


def prepMutationsGDSC(fname):
    
    df = pd.read_csv(fname).set_index('model_name')[['gene_symbol', 'cancer_driver']]
    
    # Drivers are only 4886 of 1796526 records
    se1 = df['gene_symbol'].groupby(level=0).unique() #.apply(len).sort_values(ascending=False)
    se2 = df['gene_symbol'][df['cancer_driver']].groupby(level=0).unique()
    df = pd.concat([se1, se2], axis=1, sort=False, keys=['Any Genes', 'Cancer Drivers'])
    se_models_mutations = df['Any Genes']
    se_models_mutations.index.name = 'MODEL'
    se_models_mutations.name = 'MUTATIONS'
    
    return se_models_mutations


def getTissueAnnotationGDSC(fname):
    
    # There are 43 cell lines for which we have mutation information, but no
    # tissue information.
    # pd.Series(dfTa.index.difference(df.index)).to_csv('missing_tissue_annotations.csv',
    # index=False))

    se_tissue_annotation_GDSC = pd.read_csv(fname).set_index('model_name')['tissue']
    se_tissue_annotation_GDSC = se_tissue_annotation_GDSC.loc[~se_tissue_annotation_GDSC.index.duplicated()].sort_index()
    se_tissue_annotation_GDSC.index.name = 'MODEL'
    se_tissue_annotation_GDSC.name = 'TISSUE'
    
    return se_tissue_annotation_GDSC


def getTissueAnnotationAZ(fname):

    se_tissue_annotation_AZ = pd.read_csv(fname, index_col=0)['Tissue..General.'] #['Disease.Area']
    se_tissue_annotation_AZ.name = 'TISSUE'
    se_tissue_annotation_AZ.index.name = 'MODEL'
    
    return se_tissue_annotation_AZ


def prepDrugTargetsGDSC(fname):
    
    # Curated by Grace from:
    # ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/screened_compunds_rel_8.2.csv
    
    se_ = pd.read_excel(fname).set_index('DRUG_NAME')['Known Target Genes']
    se_ = se_[~se_.index.duplicated()].sort_index()
    se_ = se_.str.split(', ').dropna()
    se_.name = 'TARGETS'
    se_.index.name = 'DRUG'
    
    return se_


def prepGenesSanger(fname):
    
    se_ = pd.read_csv(fname, index_col=0)['hgnc_symbol'].sort_values()
    se_.index = se_.values
    
    return se_


def prepDrugTargetsAZ(fname, fnameGenes):
    
    # https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20191101.csv
    
    se_ = pd.read_csv(fname, index_col=0)['Target(Official Symbol)']
    se_.name = 'TARGETS'
    se_.index.name = 'DRUG'
    
    se_drug_targets_AZ = se_
    
    se_genes_sanger = prepGenesSanger(fnameGenes)

    func = lambda pat: se_genes_sanger[se_genes_sanger.str.match('^%s.*[0-9]' % pat.strip('*'))].values.tolist()
    t = se_drug_targets_AZ.str.replace(' ', '').str.split(',').apply(pd.Series).stack()
    # print({v: func(v) for v in
    # t[t.str.contains(r'\*')].drop_duplicates().values})
    wt = {'AKT*': ['AKT1', 'AKT1S1', 'AKT2', 'AKT3', 'AKT3-IT1'], 
         'PIK3C*': ['PIK3C2A', 'PIK3C2B', 'PIK3C2G', 'PIK3C3', 'PIK3CD-AS1', 'PIK3CD-AS2', 'PIK3CDP1'], 
         'SGK*': ['SGK1', 'SGK2', 'SGK3'], 
         'IGF*R': ['IGF1R', 'IGF2R'], 
         'TOP*': ['TOP1', 'TOP1MT', 'TOP1P1', 'TOP1P2', 'TOP2A', 'TOP2B', 'TOP3A', 'TOP3B', 'TOP3BP1', 'TOPAZ1', 'TOPBP1'], 
         'ERBB*': ['ERBB2', 'ERBB3', 'ERBB4'], 
         'HDAC*': ['HDAC1', 'HDAC10', 'HDAC11', 'HDAC11-AS1', 'HDAC1P1', 'HDAC1P2', 'HDAC2', 'HDAC2-AS2', 'HDAC3', 'HDAC4', 'HDAC5', 'HDAC6', 'HDAC7', 'HDAC8', 'HDAC9'], 
         'IKBK*': ['IKBKAP-IT1', 'IKBKGP1'], 
         'MAP2K*': ['MAP2K1', 'MAP2K1P1', 'MAP2K2', 'MAP2K3', 'MAP2K4', 'MAP2K4P1', 'MAP2K5', 'MAP2K6', 'MAP2K7'], 
         'PIP5K1*': ['PIP5K1P1', 'PIP5K1P2'], 
         'PRKC*': ['PRKCA-AS1', 'PRKCQ-AS1', 'PRKCZ-AS1'], 
         'PTK2*': ['PTK2', 'PTK2B'], 
         'TOP2*': ['TOP2A', 'TOP2B'], 
         'WNT*': ['WNT1', 'WNT10A', 'WNT10B', 'WNT11', 'WNT16', 'WNT2', 'WNT2B', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT5A-AS1', 
                  'WNT5B', 'WNT6', 'WNT7A', 'WNT7B', 'WNT8A', 'WNT8B', 'WNT9A', 'WNT9B']}

    wt = {k:','.join(v) for k,v in wt.items()}
    se_ = se_drug_targets_AZ.str.replace(' ', '').str.split(',').apply(pd.Series).stack().replace(wt).unstack()
    se_ = pd.Series(index=se_.index, data=se_.apply(lambda s: ','.join(s.dropna().values), axis=1).values)
    se_ = se_.str.split(',')
    se_.name = 'TARGETS'

    return se_


def prepareSubChallenge1(df, index_train_test, index_validate, encode=True):
    
    df_temp = df.loc[(~df['Cij'].isna()) & (~df['Sik'].isna()) & (~df['Sjk'].isna()) & (~df['SYNERGY_SCORE'].isna())]
    #print('All with known pairs:\t', df_temp.shape[0])
    
    df_train_test = df_temp.loc[df_temp.reset_index('TISSUE').index.isin(index_train_test)]
    print('Training-testing pairs:\t', df_train_test.shape[0])
    
    df_validate = df_temp.loc[df_temp.reset_index('TISSUE').index.isin(index_validate)]
    print('Validation pairs:\t', df_validate.shape[0])
    
    if encode:
        dfe_train_test = encodeNames(df_train_test)
        dfe_train_test = dfe_train_test[dfe_train_test.columns.intersection(encodeNames(df_validate).columns)]
        dfe_validate = encodeNames(df_validate).reindex(dfe_train_test.columns, axis=1).fillna(0.).astype(int)
            
        df_train_test.loc[:, dfe_train_test.columns] = dfe_train_test.values
        df_validate.loc[:, dfe_validate.columns] = dfe_validate.values
        
    return df_train_test, df_validate

