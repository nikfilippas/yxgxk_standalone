"""
Verify that the NFW extrapolation of Duffy08 to 500c
is consistent with the cM relation in Ishiyama21.
"""
import numpy as np
# import matplotlib.pyplot as plt
import pyccl as ccl
from fit_duffy import ConcentrationDuffyFit500c, best_fit, used_fit

mass_def = ccl.halos.MassDef500c()
concentrations = [
    # ConcentrationDuffyFit500c(**dict(zip(["A", "B", "C"], best_fit()))),
    ConcentrationDuffyFit500c(**dict(zip(["A", "B", "C"], used_fit()))),
    ccl.halos.Concentration.from_name("Ishiyama21")(mass_def=mass_def)]
s8_arr = np.linspace(0.7, 1.1, 8)
ells = np.geomspace(6, 1024, 32)
hmc = ccl.halos.HMCalculator(mass_function="Tinker08",
                             halo_bias="Tinker10",
                             mass_def=mass_def)

for s8 in s8_arr:
    cosmo = ccl.Cosmology(Omega_c=0.25, Omega_b=0.05, h=0.67,
                          n_s=0.96, sigma8=s8)
    tr = ccl.CMBLensingTracer(cosmo)

    Cls = np.empty((2, ells.size))
    for i, c in enumerate(concentrations):
        prof_m = ccl.halos.HaloProfileNFW(c_m_relation=c)
        pka = ccl.halos.halomod_Pk2D(cosmo, hmc, prof_m, normprof=True)
        Cls[i] = cosmo.angular_cl(tr, tr, ell=ells, p_of_k_a=pka)

    dev = np.fabs(1 - np.divide(*Cls))
    # plt.figure()
    # plt.loglog(ells, dev)
    print(dev.max())
    # assert dev.max() < 5e-3
    # assert dev.max() < 5e-4
