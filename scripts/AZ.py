from acda.method_functions import *
from acda.plot_functions import *
from acda.general_functions import *
from preparation_functions import *
from sklearn.linear_model import LinearRegression
import pandas as pd

df = pd.read_pickle('df_AZ_lung.pklz')

v = MonteCarloCrossValidation(df, n=1, clf_for_CDA=LinearRegression())
print(v[0])


