from pystan import StanModel
import scipy.io as sio

# Load in example data and parameters from MATLAB implementation (pynamicalSystem)
ldsp   = sio.loadmat('exampleLDSforPython.mat')
# stan_data_dict['x0cov_mult'] = 1e-4;

with open("stan-lds.txt") as f:
    stan_code = f.read()

stan_data_dict = {'N': ldsp['N'][0][0], 'p': ldsp['p'][0][0], 'n': ldsp['n'][0][0], 'Y': ldsp['Y'],
                  'A': ldsp['A'], 'Q': ldsp['Q'], 'C': ldsp['C'], 'R': ldsp['R'], 'x0cov_mult': 1e-4}

lds_model  = StanModel(model_code=stan_code)
#lds_sample = lds_model.sampling(data=stan_data_dict, iter=1000, chains=2);