from acda.method_functions import *
from acda.plot_functions import *

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

dfTas_GDSC2_breast, dfC_G2, dfS_G2, Z_G2 = makeCDAformattedData('Breast', 
                                                                se_drug_synergy_CDA, 
                                                                se_tissue_annotation_GDSC2, 
                                                                df_drug_sensitivity_GDSC2, 
                                                                se_models_mutations_GDSC2, 
                                                                se_drug_targets_GDSC2, 
                                                                'GDSC2', 
                                                                sensitivity_metric='LNIC50', 
                                                                returnMore=True)

se_predicted = pd.concat([sample_train_predicted(dfTas_GDSC2_breast, i) for i in range(3)], 
                         axis=1).unstack(0).groupby(level=1, axis=1).agg(np.nanmean).stack(
                             ).reorder_levels([3, 0, 1, 2]).sort_index()
print(se_predicted)
se_predicted.to_csv('predicted.csv')

fig = plotHeatmapPredictedSynergy(dfC_G2, Z_G2, se_predicted[se_predicted>=0.95].index.droplevel(-1).values)
fig.savefig('heatmap.png', dpi=300)

temp = dfTas_GDSC2_breast['SYNERGY_SCORE'].droplevel(['MODEL', 'TISSUE'])
fig = plotDendrogramWithKnownPairs(Z_G2, dfC_G2, temp[temp==1].index.unique())
fig.savefig('dendrogram.png', dpi=300)
