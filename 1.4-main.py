from pystan import StanModel

import numpy as np
from sklearn.datasets import load_breast_cancer

data_wc    = load_breast_cancer()


stan_data_dict = {'N': data_wc.target.size, 'p': 1,
                  'X': data_wc.data[:, [0]], 'y': data_wc.target}


tmpHeader = "my-dotsad4.hpp"
with open("tmp.hpp", "w") as f:
    with open(tmpHeader, "r") as f_cur:
        tmpHeaderTxt = f_cur.read()
    f.write(tmpHeaderTxt)


if True:
    #print("optimum value of original (sq) NLL is {}".format(stan_optim2))
    print("=================================")
    print("Now onto the customised version...")
    with open("stan_logreg_4my.txt") as f:
        stan_code_my  = f.read()
    if False:
        bayes_logreg3_my = StanModel(model_code=stan_code_my, verbose=True)
    if True:
        bayes_logreg3_my  = StanModel(model_code=stan_code_my, allow_undefined=True,
                                     includes=['tmp.hpp'], include_dirs=["."], verbose=True)
    # stan_optim3_my    = bayes_logreg3_my.optimizing(data=stan_data_dict, iter=50000, init=initBeta);