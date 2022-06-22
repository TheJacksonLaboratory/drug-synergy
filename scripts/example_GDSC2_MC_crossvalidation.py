from acda.method_functions import *

sdatadir = '../docs/examples/data/'

df_drug_sensitivity_GDSC2 = pd.read_csv(
    sdatadir + 'GDSC2_drug_sensitivity.csv.gz').set_index(['MODEL', 'DRUG'])

se_drug_synergy_CDA = pd.read_csv(
    sdatadir + 'CDA_synergy_pairs.csv.gz').set_index(['MODEL', 'DRUG1', 'DRUG2'])['SYNERGY_SCORE']

se_tissue_annotation_GDSC2 = pd.read_csv(
    sdatadir + 'GDSC2_tissue_annotation.csv.gz').set_index(['MODEL'])['TISSUE']

se_drug_targets_GDSC2 = pd.read_csv(
    sdatadir + 'GDSC2_drug_targets.csv.gz').set_index(['DRUG'])['TARGETS']

se_models_mutations_GDSC2 = pd.read_csv(
    sdatadir + 'GDSC2_model_mutations.csv.gz').set_index(['MODEL'])['MUTATIONS']

dfTas_GDSC2_breast = makeCDAformattedData('Breast', 
                                          se_drug_synergy_CDA, 
                                          se_tissue_annotation_GDSC2, 
                                          df_drug_sensitivity_GDSC2, 
                                          se_models_mutations_GDSC2, 
                                          se_drug_targets_GDSC2, 
                                          'GDSC2', 
                                          sensitivity_metric='LNIC50')
print(dfTas_GDSC2_breast)

df_GDSC2_breast = MonteCarloCrossValidation(dfTas_GDSC2_breast, sample_non_synergy=True)[0]
print(df_GDSC2_breast)