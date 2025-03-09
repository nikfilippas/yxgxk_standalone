"""Interpolate the sigma8-redshift plane of Omega_grav using Battaglia."""
import pyccl as ccl
import numpy as np
from scipy.interpolate import RBFInterpolator
from scipy.stats import qmc
from pressure import BattagliaCalculator

def cosmo_s8(s8):
    return ccl.Cosmology(
        Omega_c=0.26066676, Omega_b=0.048974682, h=0.6766, n_s=0.9665,
        sigma8=s8, transfer_function="boltzmann_camb")

F_Omgr = BattagliaCalculator(base_model="gyksrA")._get_Omgr
lo_bounds, hi_bounds = [0.6, 0.05], [1.1, 0.40]

# Set up samples and sample the function.
sampler = qmc.LatinHypercube(d=2)
sample = sampler.random(n=1024)
sample = qmc.scale(sample, lo_bounds, hi_bounds)
d = np.array([F_Omgr(cosmo_s8(s8), z) for s8, z in sample])

# Build and save interpolator.
import pickle
rbf = RBFInterpolator(sample, d, neighbors=20)
with open("data/Omgr_Battaglia", "wb") as f:
    pickle.dump(rbf, f)
