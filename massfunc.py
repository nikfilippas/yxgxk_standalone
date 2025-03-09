"""
This is an explainer of the discrepancy we see in the posterior using the Despali16
mass function. D16 predicts an overabundance of high-mass haloes. Therefore, a smaller
value of sigma8 is needed to reconcile our observations.
Note, according to Press-Schecter, n(M) ~ exp(1/σ(M)), so it is *very* sensitive.
We observe a similar effect when B16 starts to deviate from our fiducial (T08).
After WISC-3, the posterior contours diverge.
"""
import pyccl as ccl
import numpy as np
import matplotlib.pyplot as plt

m_arr = np.logspace(10, 15, 128)
a_arr = np.linspace(1/1.4, 1, 16)

cosmo = ccl.CosmologyVanillaLCDM()


c_m = ccl.halos.ConcentrationIshiyama21()
mass_def = ccl.halos.MassDef500c(c_m=c_m)
massfunc_T08 = ccl.halos.MassFuncTinker08(cosmo, mass_def=mass_def)
massfunc_D16 = ccl.halos.MassFuncDespali16(cosmo, mass_def=mass_def)
massfunc_B16 = ccl.halos.MassFuncBocquet16(cosmo, mass_def=mass_def)

T08 = np.array([massfunc_T08.get_mass_function(cosmo, M=m_arr, a=a) for a in a_arr])
D16 = np.array([massfunc_D16.get_mass_function(cosmo, M=m_arr, a=a) for a in a_arr])
B16 = np.array([massfunc_B16.get_mass_function(cosmo, M=m_arr, a=a) for a in a_arr])

fig, ax = plt.subplots(1, 2, sharey=True, figsize=(12, 5))
ax[0].set_title("1-T08/D16")
ax[1].set_title("1-T08/B16")
ax[0].set_ylabel("scale factor")
ax[0].set_xlabel("log10(M)")
ax[1].set_xlabel("log10(M)")


kwargs = {"vmin": 0, "vmax": 1,
          "extent": [np.log10(m_arr[0]), np.log10(m_arr[-1]), a_arr[0], a_arr[-1]],
          "aspect": "auto"}

ax[0].imshow(np.abs(1-T08/D16), **kwargs)
ax[1].imshow(np.abs(1-T08/B16), **kwargs)
ax[1].colorbar()
