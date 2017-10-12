from pystan import StanModel

import numpy as np
from sklearn.datasets import load_breast_cancer

data_wc    = load_breast_cancer()
with open("stan-logreg3.txt") as f:
    stan_code = f.read()

stan_data_dict = {'N': data_wc.target.size, 'p': data_wc.data.shape[1],
                  'X': data_wc.data, 'y': data_wc.target}

initBeta = [{'beta': np.arange(-3,3,0.2)*0.05}]

tmpHeader = "my-dotprod.hpp"
with open("tmp.hpp", "w") as f:
    with open(tmpHeader, "r") as f_cur:
        tmpHeaderTxt = f_cur.read()
    f.write(tmpHeaderTxt)

initBeta = [{'beta': np.arange(-3,3,0.2)*0.05}]

if True:
    bayes_logreg3_native = StanModel(model_code=stan_code)
    # stan_fit_native    = bayes_logreg2_native.sampling(data=stan_data_dict, iter=2000, chains=2);
    stan_optim3          = bayes_logreg3_native.optimizing(data=stan_data_dict, iter=50000, init=initBeta);

if True:
    #print("optimum value of original (sq) NLL is {}".format(stan_optim2))
    print("=================================")
    print("Now onto the customised version...")
    with open("stan_logreg_3my.txt") as f:
        stan_code_my  = f.read()
    bayes_logreg3_my  = StanModel(model_code=stan_code_my, allow_undefined=True,
                                 includes=['tmp.hpp'], include_dirs=["."])
    stan_optim3_my    = bayes_logreg3_my.optimizing(data=stan_data_dict, iter=50000, init=initBeta);