import numpy as np
import pandas as pd
from statsmodels.formula.api import ols*

def regression_multivariable(s,x,y,z):
  df = pd.DataFrame({'s':S, 'x':X, 'y':Y, 'z':Z})
  model = ols("s ~z + x + y", df).fit() #trouve le modèle

  print(model.summary)
