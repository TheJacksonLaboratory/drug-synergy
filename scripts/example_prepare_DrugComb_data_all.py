from acda.general_functions import *

keepMonotherapyData = False

dir = '../data/'

if keepMonotherapyData:
    tempfname = dir + 'cacheDrugCombWithAllMono.pklz'
else:
    tempfname = dir + 'cacheDrugComb.pklz'

if not os.path.exists(dir):
    os.makedirs(dir)

if not os.path.isfile(tempfname):
    # Data was downloaded from "https://drugcomb.fimm.fi/jing/summary_v_1_5.csv"
    df_DC = pd.read_csv(dir + 'drugcomb_data_v1.5.csv.gz', index_col=0)
    df_DC = df_DC.set_index(['study_name', 'cell_line_name', 'drug_row', 'drug_col', 'tissue_name'])
    df_DC.index.names = ['STUDY', 'MODEL', 'DRUG1', 'DRUG2', 'TISSUE']

    # Remove entries which are duplicates
    df_DC = df_DC.sort_index(level='DRUG2', ascending=True)
    df_DC = df_DC.loc[~df_DC.index.duplicated(keep='first')]
    df_DC = df_DC.sort_index()

    # Remove entries where DRUG1 is equal to DRUG2
    df_DC = df_DC.loc[df_DC.index.to_frame()['DRUG1'] != df_DC.index.to_frame()['DRUG2']]

    # Remove entries with no combinations measures, i.e. monotherapy experiments
    if not keepMonotherapyData:
        df_DC = df_DC.loc[pd.MultiIndex.from_frame(df_DC.index.to_frame().dropna())]

    dfgp = getGeneToProteinNameAssociation()
    convertProteinNamesToGenes(df_DC['drug_row_target_name'], dfgp)
    convertProteinNamesToGenes(df_DC['drug_col_target_name'], dfgp)

    # Keep subset of the columns, rename selected, below is the list of "not selected" columns:
    # ['conc_row_unit', 'conc_col_unit', 'css_row', 'css_col', 'css_ri', 'S_sum', 'S_mean', 
    # 'S_max', 'drug_row_clinical_phase', 'drug_col_clinical_phase']
    df_DC = df_DC[['ic50_row', 'ic50_col', 'ri_row', 'ri_col', 'synergy_zip', 'synergy_loewe', 
                   'synergy_hsa', 'synergy_bliss', 'drug_row_target_name', 'drug_col_target_name']]
    df_DC.columns = ['IC50_DRUG1', 'IC50_DRUG2', 'AUC_DRUG1', 'AUC_DRUG2', 'SYNERGY_ZIP', 
                     'SYNERGY_LOEWE', 'SYNERGY_HSA', 'SYNERGY_BLISS', 'DRUG1_TARGETS', 'DRUG2_TARGETS']
    
    df_DC.to_pickle(tempfname)
else:
    df_DC = pd.read_pickle(tempfname)
    
print(df_DC)
