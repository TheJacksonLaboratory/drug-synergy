**Examples**
============

Example of loading and running a full-format dataset from DrugComb, AstraZeneca study, breast tissue. The script takes a few minutes per iteration (n=3 iterations are suggested for testing).

.. literalinclude:: ../../scripts/example_DrugComb_AstraZeneca_full.py

+---------+-------+-------+
|         |  mean |   sem |
+---------+-------+-------+
|    ACDA | 0.867 | 0.028 |
+---------+-------+-------+
|     CDA | 0.521 | 0.039 |
+---------+-------+-------+
|      EN | 0.829 | 0.020 |
+---------+-------+-------+
| EN-ACDA | 0.894 | 0.021 |
+---------+-------+-------+


Example of Monte Carlo cross-validation on an alternative format dataset: GDSC2 subset of breast tissue with the CDA synergy pairs.

.. literalinclude:: ../../scripts/example_GDSC2_MC_crossvalidation.py


Example of generating predictions on the alternative format dataset and visualizing them on a heatmap.

.. literalinclude:: ../../scripts/example_GDSC2_heatmap_predictions.py


Load data from DrugComb. The file "drugcomb_data_v1.5.csv" is 1.3Gb (compressed is only 170Mb) and can be downloaded from https://drugcomb.fimm.fi/jing/summary_v_1_5.csv.

.. literalinclude:: ../../scripts/example_prepare_DrugComb_data_all.py

