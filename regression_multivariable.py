import numpy as np
import pandas as pd
from statsmodels.formula.api import ols*

def regression_multivariable(s,x,y,z):
  df = pd.DataFrame({'x':X, 'y':Y, 'z':Z})
  model = ols("z ~ x + y", df).fit() #trouve le modèle

  print(model.summary)
