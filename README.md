## Overview

This repository stores drug synergy prediction codebase for Augmented Cancer Drug Atlas (ACDA) and the analysis jupyter notebooks which use the ACDA code. We augmented the drug synergy prediction modeling approach CDA described in Narayan et al. by applying a Random Forest Regression and optimization via cross-validation hyperparameter tuning. For ease of sharing and use we implemented it as a python package. See documentation at https://acda.readthedocs.io.


## Installation and Dependencies

The main prerequisite o install the Python package is python >3.8 environment. To install run:

```
pip install acda
```

The dependencies are installed automatically with the command above. See [setup.py](https://github.com/TheJacksonLaboratory/drug-synergy/blob/main/setup.py) for the very basic dependency list.

## Methods Description

Drug Synergy prediction is a complex problem typically approached with Machine Learning techniques by using molecular and pharmacological data. A recently published method Cancer Drug Atlas, CDA (Narayan et al. 2020), for drug synergy prediction in cell line models uses drug target information, knowledge of genes mutated in each model, and models' monotherapy drug sensitivity. The approach builds a logistic regression model to predict a binary synergy outcome. Here, we improved the CDA drug synergy prediction modeling approach by applying a CART-based model instead of a linear regression.

The datails of ACDA methods and benchmarking results are described in the [preprint](https://www.biorxiv.org/).


## Contents

| Directory/File  | Description                |
|-----------------|----------------------------|
| acda/           | Python source code of the package                  |
| docs/           | Source code of the documentation                   |
| scripts/        | Scripts and additional functions which are used with the package |
| ChangeLog.md    | File details changes implemented in new releases   |
| LICENSE         | The license agreement statement                    |
| [Other misc files]             | ...                                                |


## Data

The datailed examples and data requirements are in the documentation at [documentation](https://acda.readthedocs.io).

This package makes use of the data outlined below:
 + Sanger molecular data: https://cellmodelpassports.sanger.ac.uk/downloads
 + GDSC pharmacology data: https://www.cancerrxgene.org/downloads/bulk_download
 + DrugComb sensitivity and synergy data: https://drugcomb.org/download/


## References:

Narayan, R. S., Molenaar, P., Teng, J., Cornelissen, F. M. G., Roelofs, I., Menezes, R., Dik, R., et al. (2020). *A cancer drug atlas enables synergistic targeting of independent drug vulnerabilities*, **Nature Communications, 11/1: 2935. Nature Publishing Group**.

Lianlian Wu, Yuqi Wen, Dongjin Leng, Qinglong Zhang, Chong Dai, Zhongming Wang, Ziqi Liu, Bowei Yan, Yixin Zhang, Jing Wang, Song He and Xiaochen Bo. *Machine learning methods, databases and tools for drug combination prediction.* **Briefings in Bioinformatics, 23(1), 2022, 1-21.**

Shuyu Zheng, Jehad Aldahdooh, Tolou Shadbahr, Yinyin Wang, Dalal Aldahdooh, Jie Bao, Wenyu Wang, Jing Tang. *DrugComb update: a more comprehensive drug sensitivity data repository and analysis portal.* **Nucleic Acids Research, Volume 49, Issue W1, 2 July 2021, Pages W174-W184.**