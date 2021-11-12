import numpy as np
import pyccl as ccl
from tqdm import tqdm
from scipy.interpolate import interp2d
import getdist.mcsamples as gmc
import matplotlib.pyplot as plt


class ChainCalculator(object):
    def __init__(self, *, new_interps=False):
        self.cosmo = ccl.CosmologyVanillaLCDM()
        self.cM = ccl.halos.ConcentrationDuffy08M500c()
        self.prof_g = ccl.halos.HaloProfileHOD(c_m_relation=self.cM)
        self.prof_k = ccl.halos.HaloProfileNFW(c_m_relation=self.cM)
        self.prof_y = ccl.halos.HaloProfilePressureGNFW()
        self.hmc = ccl.halos.HMCalculator(mass_function="Tinker08",
                                          halo_bias="Tinker10",
                                          mass_def="500c")

        # redshift distributions
        self.names = ["2mpz"] + ["wisc%d" % i for i in range(1, 6)]
        self.DIR = dict.fromkeys(self.names)
        self.get_z_arr()

        # interpolators
        self.interp_param_names = ["bPe", "Omth"]
        self.interpolators = self.get_interpolators(new_interps)

    def get_z_arr(self):
        fnames = [f"data/dndz/{name}_DIR.txt" for name in self.names]

        z_arr = []
        for name, fname in zip(self.names, fnames):
            z, nz = np.loadtxt(fname).T
            z_arr.append(np.average(z, weights=nz))
            self.DIR[name] = [np.column_stack((z, nz))]
        self.z_arr = np.array(z_arr)

    def update_parameters(self, *, lMmin_0=None, lM1_0=None,
                          mass_bias=None, sigma8=None):
        if sigma8 is not None:
            self.cosmo = ccl.Cosmology(
                Omega_c=0.26066676, Omega_b=0.048974682,
                h=0.6766, n_s=0.9665, sigma8=sigma8)
        if lMmin_0 is not None:
            self.prof_g.update_parameters(lMmin_0=lMmin_0, lM0_0=lMmin_0)
        if lM1_0 is not None:
            self.prof_g.update_parameters(lM1_0=lM1_0)
        if mass_bias is not None:
            self.prof_y.update_parameters(mass_bias=mass_bias)

    def get_interpolators(self, new_interps=False):
        import os
        if not os.path.isfile("interpolators.npy") or new_interps:
            I = dict.fromkeys(self.interp_param_names)
            for par in self.interp_param_names:
                I[par] = self.interpolate_param(par)
            np.save("interpolators.npy", I)
        else:
            I = np.load("interpolators.npy", allow_pickle=True).item()
        return I

    def calculate_bPe(self, z):
        bpe = ccl.halos.halomod_bias_1pt(
            self.cosmo, self.hmc,
            k=1e-3, a=1/(1+z),
            prof=self.prof_y, normprof=False)
        return 1e3 * bpe

    def calculate_Pe(self, z):
        pe = ccl.halos.halomod_mean_profile_1pt(
            self.cosmo ,self.hmc,
            k=1e-3, a=1/(1+z),
            prof=self.prof_y, normprof=False)
        return pe

    def calculate_Omth(self, z):
        pe = self.calculate_Pe(z)
        Y = 0.24
        prefac = (8-5*Y)/(4-2*Y)
        rho_th = pe*prefac/(1+z)**3
        # rho_critical in eV/cm^3
        rho_crit = 10537.0711*self.cosmo['h']**2
        return rho_th/rho_crit

    def interpolate_param(self, parname,
                          s8_min=0.2, s8_max=1.5, N_s8=16,
                          bH_min=0.005, bH_max=1.15, N_bH=16):
        # get derived parameter function
        func = getattr(self, f"calculate_{parname}")

        # define interpolation boundaries
        s8_arr = np.linspace(s8_min, s8_max, N_s8)
        bH_arr = np.linspace(bH_min, bH_max, N_bH)

        # loop through all redshifts within the boundaries
        F_interp = dict.fromkeys(self.names)
        for name, z in tqdm(zip(self.names, self.z_arr)):
            Arr = np.zeros((N_s8, N_bH))
            for i, s8 in enumerate(s8_arr):
                self.update_parameters(sigma8=s8)
                for j, bH in enumerate(bH_arr):
                    self.update_parameters(mass_bias=bH)
                    Arr[i, j] = func(z)

            F_interp[name] = interp2d(s8_arr, bH_arr, Arr, kind="cubic")

        return F_interp

    def from_chains(self, model, parname, latex=None):
        BF_arr = np.zeros((6, 3))
        for ibin, z in enumerate(tqdm(self.z_arr)):
            try:
                fname = f"chains/{model}/{model}_{ibin}/cobaya"
                s = gmc.loadMCSamples(fname, settings={'ignore_rows': 0.3})
            except ValueError:
                fname = f"chains/{model}/{model}_{ibin}_kmax1/cobaya"
                s = gmc.loadMCSamples(fname, settings={'ignore_rows': 0.3})
            p = s.getParams()

            if parname in self.interpolators:
                rel = self.interpolators[parname][self.names[ibin]]

                if hasattr(p, "sigma8"):
                    # free sigma8
                    bPe_chain = np.array([rel(s8, by) for s8, by
                        in zip(p.sigma8, p.ygk_mass_bias)]).squeeze()
                else:
                    # fixed sigma8
                    bPe_chain = np.array([rel(0.81, by)
                        for by in p.ygk_mass_bias]).squeeze()

                s.addDerived(bPe_chain, name=parname, label=latex)

            dens = s.get1DDensity(parname)
            vmin, vmax = dens.getLimits(0.68)[:2]
            vbf = dens.getLimits(0.001)[0]
            summary = vbf, vbf-vmin, vmax-vbf
            BF_arr[ibin] = summary

        # normalize by the growth factor for sigma8
        if parname == "sigma8":
            cosmo = ccl.CosmologyVanillaLCDM()
            BF_arr *= cosmo.growth_factor(1/(1+self.z_arr))[:, None]

        return BF_arr


def plot_tomo(models, parname, labels):
    if parname == "bPe":
        latex = "bP_e"
    elif parname == "Omth":
        latex = "Omega_th"
    else:
        latex = None

    fig, ax = plt.subplots()
    ax.set_xlabel("z", fontsize=16)
    ax.set_ylabel(latex, fontsize=16)
    fig.tight_layout()
    colors = ["k", "grey", "r", "brown", "orange"]

    for i, model in enumerate(models):
        BF = c.from_chains(model=model, parname=parname, latex=latex)
        ax.errorbar(c.z_arr+0.005*i, BF[:, 0], BF[:, 1:].T,
                    fmt="o", color=colors[i], label=labels[i])

    ax.legend(loc="upper right", fontsize=12)

c = ChainCalculator()

models = ["yxgxksig", "yxgxk_b08", "gxk", "gxk_kmax05"]
labels = models
plot_tomo(models, "sigma8", labels)

models = ["yxgxksig", "yxgxk"]
labels = models
plot_tomo(models, "ygk_mass_bias", labels)
plot_tomo(models, "Omth", labels)
