.. _input-data:

**Input data format**
=====================

Input is expected as spreadsheets of comma-separated values ``csv``.

> *NB:* In the underlying software, MODEL, DRUG and GENE identifiers are used to match and query the tables described above.
Therefore these identifiers must be curated so that they are consistent across the tables.


Known synergy pairs
-------------------
The first three columns contain model identifier, and drugs identifiers. The last column is the synergy score, or, if unavailable, a binary value 1/0 to indicate synergy or no synergy between the pair of drugs.

+-------+------------+--------------+---------------+
| MODEL | DRUG1      | DRUG2        | SCORE         |
+=======+============+==============+===============+
| 5637  | cisplatin  | sunitinib    | 1             |
+-------+------------+--------------+---------------+
| 5637  | sunitinib  | cisplatin    | 1             |
+-------+------------+--------------+---------------+
| 8505C | bortezomib | docetaxel    | 1             |
+-------+------------+--------------+---------------+
| 8505C | docetaxel  | bortezomib   | 1             |
+-------+------------+--------------+---------------+
| A172  | cisplatin  | temozolomide | 1             |
+-------+------------+--------------+---------------+
| ...   | ...        | ...          | ...           |
+-------+------------+--------------+---------------+


Sensitivity
-----------
A table of sensitivity measures; examples include LNIC50, AUC etc.

+-------+---------------------+-----------+----------+-----+
| MODEL | DRUG                | LNIC50    | AUC      | ... |
+=======+=====================+===========+==========+=====+
| 201T  | 5-Fluorouracil      | 3.738474  | 0.873656 | ... |
+-------+---------------------+-----------+----------+-----+
| 201T  | A-83-01             | 5.577933  | 0.975815 | ... |
+-------+---------------------+-----------+----------+-----+
| 201T  | ACY-1215            | 2.008393  | 0.746972 | ... |
+-------+---------------------+-----------+----------+-----+
| 201T  | AGI-6780            | 2.031296  | 0.977958 | ... |
+-------+---------------------+-----------+----------+-----+
| 201T  | AICA Ribonucleotide | 10.006271 | 0.972249 | ... |
+-------+---------------------+-----------+----------+-----+
| ...   | ...                 | ...       | ...      | ... |
+-------+---------------------+-----------+----------+-----+


Model tissue annotation
-----------------------
We recommend to perform analysis by tissue. Therefore we require tissue annoation for each model.

+---------+------------------------+
| MODEL   | TISSUE                 |
+=========+========================+
| 1181N1  | Central Nervous System |
+---------+------------------------+
| 1205Lu  | Skin                   |
+---------+------------------------+
| 1273-99 | Soft Tissue            |
+---------+------------------------+
| 1321N1  | Central Nervous System |
+---------+------------------------+
| 143B    | Bone                   |
+---------+------------------------+
| ...     | ...                    |
+---------+------------------------+

Model mutations
---------------
List of mutated genes for each model. Note that the truncated raw csv file looks as below:

```
MODEL,MUTATED GENES
201T,"IL16, TEKT4P2, TP53, NDC1, MROH7, INTS11"
22RV1,"SLC2A13, MUC19, SLC15A4, RAB40C, RHOT2"
```

+----------+-----------------------------------------------+
| MODEL    | GENES                                         |
+==========+===============================================+
| 201T     | IL16, TEKT4P2, TP53, NDC1, MROH7, INTS11, ... |
+----------+-----------------------------------------------+
| 22RV1    | SLC2A13, MUC19, SLC15A4, RAB40C, RHOT2, ...   |
+----------+-----------------------------------------------+
| 23132-87 | MYH7, EMILIN1, HECW2, CCDC39, PPARGC1B, ...   |
+----------+-----------------------------------------------+
| 42-MG-BA | LRRC74A, PILRB, CUX1, RB1, STAG2, TP53, ...   |
+----------+-----------------------------------------------+
| 451Lu    | POTEG, SPEF1, LMBRD1, BRAF, ROBO2, TP53, ...  |
+----------+-----------------------------------------------+
| ...      | ...                                           |
+----------+-----------------------------------------------+

Drug targets
------------
List of genes targeted by a drug for each drug.

+-------------------+----------------------------+
| DRUG              | GENE                       |
+===================+============================+
| (5Z)-7-Oxozeaenol | MAP3K7                     |
+-------------------+----------------------------+
| 5-Fluorouracil    | Antimetabolite (DNA & RNA) |
+-------------------+----------------------------+
| A-443654          | AKT1, AKT2, AKT3           |
+-------------------+----------------------------+
| A-770041          | LCK, FYN                   |
+-------------------+----------------------------+
| A-83-01           | TGFB1                      |
+-------------------+----------------------------+
| ...               | ...                        |
+-------------------+----------------------------+
