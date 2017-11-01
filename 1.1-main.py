
from pystan import StanModel

import numpy as np
from sklearn.datasets import load_breast_cancer

data_wc    = load_breast_cancer()
with open("stan-logreg_2.txt") as f:
    stan_code = f.read()

stan_data_dict = {'N': data_wc.target.size, 'p': 1, #data_wc.data.shape[1],
                  'X': data_wc.data[:,[0]], 'y': data_wc.target}

tmpHeader = "my-multiply.txt"
with open("tmp.hpp", "w") as f:
    with open(tmpHeader, "r") as f_cur:
        tmpHeaderTxt = f_cur.read()
    f.write(tmpHeaderTxt)


if False:
    bayes_logreg2_native = StanModel(model_code=stan_code)
    # stan_fit_native    = bayes_logreg2_native.sampling(data=stan_data_dict, iter=2000, chains=2);
    stan_optim2    = bayes_logreg2_native.optimizing(data=stan_data_dict, iter=50000);

#print("optimum value of original (sq) NLL is {}".format(stan_optim2))
print("=================================")
print("Now onto the customised version...")
with open("stan_logreg_2my.txt") as f:
    stan_code_my = f.read()
bayes_logreg2_my = StanModel(model_code=stan_code_my, allow_undefined=True,
                             includes=['tmp.hpp'], include_dirs=["."])
stan_optim2_my    = bayes_logreg2_my.optimizing(data=stan_data_dict, iter=50000);