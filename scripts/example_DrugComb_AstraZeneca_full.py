from acda.method_functions import *

sdatadir = '../docs/examples/data/'

df_full = pd.read_csv(
    sdatadir + 'DrugComb_ASTRAZENECA_breast.csv.gz').set_index(['MODEL', 'DRUG1', 'DRUG2', 'TISSUE'])

df, dfC, dfS, Z = prepareFromDCfull(df_full, 1000, returnMore=True, random_state=0)

df_DC_AZ_breast = MonteCarloCrossValidation(df, n=3)[0]
print(df_DC_AZ_breast)
