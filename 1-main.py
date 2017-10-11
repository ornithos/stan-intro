import pystan
import sys

import numpy as np
from sklearn.datasets import load_breast_cancer

data_wc    = load_breast_cancer()
with open("stan-logreg.txt") as f:
    stan_code = f.read()

stan_data_dict = {'N': data_wc.target.size, 'p': data_wc.data.shape[1],
                  'X': data_wc.data, 'y': data_wc.target}

sm = pystan.StanModel(model_code=stan_code)
smoptim    = sm.optimizing(data=stan_data_dict, iter=50000);


if False:
    stan_fit   = pystan.stan(model_code=stan_code, data=stan_data_dict, iter=1000, chains=2)

    # ================================================

    import stan
    from pystan import StanModel
    import sys

    import numpy as np
    from sklearn.datasets import load_breast_cancer

    data_wc    = load_breast_cancer()
    with open("stan-logreg.txt") as f:
        stan_code = f.read()

    stan_data_dict = {'N': data_wc.target.size, 'p': data_wc.data.shape[1],
                      'X': data_wc.data, 'y': data_wc.target}

    bayes_logreg = StanModel(model_code=stan_code)
    stan_fit2    = bayes_logreg.sampling(data=stan_data_dict, iter=1000, chains=2)