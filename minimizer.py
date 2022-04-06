from cobaya.model import get_model
from cobaya.run import run
import yaml
import os
import scipy.stats as st
import sys

# Read in the yaml file
config_fn = sys.argv[1]
with open(config_fn, "r") as fin:
    info = yaml.load(fin, Loader=yaml.FullLoader)

# Get the mean proposed in the yaml file for each parameter
p0 = {}
p_all = {}
for p in info['params']:
     if isinstance(info['params'][p], dict):
         if 'ref' in info['params'][p]:
             p0[p] = info['params'][p]['ref']['loc']
             p_all[p] = p0[p]
     else:
         p_all[p] = info['params'][p]
os.system('mkdir -p ' + info['output'])

print("params_dict = ", p0)

# Compute the likelihood
model = get_model(info)
l = model.likelihood['yxgxk_like.YxGxKLike']
loglikes, derived = model.loglikes(p0)
print("chi2 = ", -2 * loglikes[0])

# Minimize
updated_info, sampler = run(info)
bf = sampler.products()['minimum'].bestfit()
pf = {k: bf[k] for k in p0.keys()}
print("Final params: ")
print(pf)
loglikes, derived = model.loglikes(pf)
chi2 = -2*loglikes[0]
ndata = l.ndata
print("chi2 = ", chi2)
print("p = ", 1-st.chi2.cdf(chi2, ndata))
