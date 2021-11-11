import numpy as np
import matplotlib.pyplot as plt
import getdist.plots as gplot
import getdist.mcsamples as gmc
import pyccl as ccl

burn = 0.3
model = "yxgxk_b08"
parname = "sigma8"
output = "figs"

z_arr = np.arange(0.07, 0.35, 0.05)
cosmo = ccl.Cosmology(Omega_c=0.25, Omega_b=0.05,
                      h=0.67, n_s=0.96, sigma8=0.81)
growth = np.array([cosmo.growth_factor(1/(1+z)) for z in z_arr])

fig, ax = plt.subplots()
ax.set_xlabel("z", fontsize=16)
ax.set_ylabel(f"{parname}", fontsize=16)

BF_arr = np.zeros((6, 3))
for ibin in range(6):
    fname = f"chains/{model}/{model}_{ibin}/cobaya"
    s = gmc.loadMCSamples(fname, settings={"ignore_rows": burn})
    dens = s.get1DDensity(f"{parname}")
    vmin, vmax = dens.getLimits(0.68)[:2]
    vbf = dens.getLimits(0.001)[0]
    summary = vbf, vbf-vmin, vmax-vbf
    BF_arr[ibin] = summary

BF_arr *= growth[:, None]  # NOTE: factor in growth
ax.errorbar(z_arr, BF_arr[:, 0], BF_arr[:, 1:].T, fmt="ko")
fig.savefig(f"{output}/{model}_{parname}.pdf", bbox_inches="tight")
