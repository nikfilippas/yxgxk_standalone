import os
import numpy as np
import pyccl as ccl
from tqdm import tqdm
from scipy.interpolate import RectBivariateSpline
import getdist.mcsamples as gmc
import getdist.plots as gplot
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True


class ChainCalculator(object):
    def __init__(self, *, new_interps=False):
        self.cosmo = ccl.Cosmology(
            Omega_c=0.26066676, Omega_b=0.048974682,
            h=0.6766, n_s=0.9665, sigma8=0.8102,
            mass_function="tinker")
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

        # MC samples
        self.mcsamples = np.empty(len(self.names), dtype=object)

        # interpolators
        self.interp_param_names = ["bPe", "Omth"]
        self.interpolators = self.get_interpolators(new_interps)

        # latex
        self.latex_names = {
            "bPe": "\\langle bP_e \\rangle\\ [\\mathrm{meV\\,cm^{-3}}]",
            "Omth": "\\Omega_{th}",
            "ygk_mass_bias": "1-b_{\\mathrm{H}}",
            "sigma8": "\\sigma_8"
            }
        self.latex_labels = {
            "yxgxksig": r"fiducial $3 \times 2\mathrm{pt}$ \,g,y,k",
            "yxgxk": r"$3 \times 2\mathrm{pt}$, \,g,y,k; fixed $\sigma_8=0.8122$",
            "yxgxk_b08": r"$3 \times 2\mathrm{pt}$ \,g,y,k; fixed $b_{\mathrm{H}}=0.80$",
            "yxgxk_b_uniform": r"$3 \times 2\mathrm{pt}$ \,g,y,k; $b_{\mathrm{H}} \sim U(0.60,0.90)$",
            "yxgxk_b_gauss": r"$3 \times 2\mathrm{pt}$ \,g,y,k; $b_{\mathrm{H}} \sim N(0.73,0.10)$",
            "gxk": r"$2 \times 2\mathrm{pt}$ \,g,k",
            "gxk_kmax05": r"$2 \times 2\mathrm{pt}$ \,g,k; $k_{\mathrm{max}}=0.5\,\mathrm{Mpc}^{-1}$",
            "yxg": r"Koukoufilippas et al., 2020",
            "dam_yxg": r"damonge Koukoufilippas et al., 2020"
            }

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

            F_interp[name] = RectBivariateSpline(s8_arr, bH_arr, Arr)

        return F_interp

    def get_summary(self, model, parname, latex=None):
        BF_arr = np.zeros((6, 3))
        for ibin, z in enumerate(tqdm(self.z_arr)):
            fname = f"chains/{model}/{model}_{ibin}/cobaya"

            s = gmc.loadMCSamples(fname, settings={'ignore_rows': 0.3})
            self.mcsamples[ibin] = s
            p = s.getParams()

            if parname in self.interpolators:
                rel = self.interpolators[parname][self.names[ibin]]

                if hasattr(p, "sigma8"):
                    # free sigma8
                    deriv_chain = np.array([rel(s8, by) for s8, by
                        in zip(p.sigma8, p.ygk_mass_bias)]).squeeze()
                else:
                    # fixed sigma8
                    deriv_chain = np.array([rel(0.81, by)
                        for by in p.ygk_mass_bias]).squeeze()

                s.addDerived(deriv_chain, name=parname, label=latex)

            dens = s.get1DDensity(parname)
            vmin, vmax = dens.getLimits(0.68)[:2]
            vbf = dens.getLimits(0.001)[0]
            summary = vbf, vbf-vmin, vmax-vbf
            BF_arr[ibin] = summary

        # normalize by the growth factor for sigma8
        if parname == "sigma8":
            BF_arr *= self.cosmo.growth_factor(1/(1+self.z_arr))[:, None]

        return BF_arr

    def theory_sigma8(self, ax):
        z_arr = np.linspace(0.01, 0.40, 64)
        s8 = self.cosmo.sigma8() * self.cosmo.growth_factor(1/(1+z_arr))
        ax.plot(z_arr, s8, "k--", lw=2, label="Planck 2018")

    def theory_bPe(self, ax):
        from battaglia import BattagliaSchockHeating
        z_arr = np.linspace(0.01, 0.40, 64)
        BSH = BattagliaSchockHeating(self.cosmo, z_arr).get_theory

        et2, et3, et5, etinf = [BSH(n_r) for n_r in [2, 3, 5, 20]]

        ax.plot(z_arr,et2,'-',label='$r_{\\rm max}=2\\,r_{200c}$',c='k')
        ax.plot(z_arr,et3,'--',label='$r_{\\rm max}=3\\,r_{200c}$',c='k')
        ax.plot(z_arr,et5,'-.',label='$r_{\\rm max}=5\\,r_{200c}$',c='k')
        ax.plot(z_arr,etinf,':',label='$r_{\\rm max}=\\infty$',c='k')

    def theory_Omth(self, ax):
        from battaglia import BattagliaSchockHeating
        z_arr = np.linspace(0.01, 0.40, 64)
        BSH = BattagliaSchockHeating(self.cosmo, z_arr).get_Om_th
        func = lambda nr: np.array([BSH(z, nr, 200) for z in z_arr])

        ot2, ot3, ot5, otinf = [func(n_r) for n_r in [2, 3, 5, 20]]

        ax.plot(z_arr,ot2,'-',label='$r_{\\rm max}=2\\,r_{200c}$',c='k')
        ax.plot(z_arr,ot3,'--',label='$r_{\\rm max}=3\\,r_{200c}$',c='k')
        ax.plot(z_arr,ot5,'-.',label='$r_{\\rm max}=5\\,r_{200c}$',c='k')
        ax.plot(z_arr,otinf,':',label='$r_{\\rm max}=\\infty$',c='k')

    def plot_tomo(self, models, parname, keep_on=False):
        latex = self.latex_names[parname]

        fig, ax = plt.subplots(figsize=(9,7))
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.tick_params(axis='both', which='minor', labelsize=12)
        ax.set_xlabel("$z$", fontsize=16)
        ax.set_ylabel(f"${latex}$", fontsize=16)
        ax.set_xlim(0.05, 0.40)
        fig.tight_layout()
        colors = ["k", "grey", "r", "brown", "orange",
                  "navy", "forestgreen", "crimson"]

        # First, plot theoretical models
        if parname in ["sigma8", "bPe"]:
            theory = getattr(self, f"theory_{parname}")
            theory(ax=ax)

        for i, model in enumerate(models):
            label = self.latex_labels[model]
            BF = self.get_summary(model=model, parname=parname, latex=latex)
            ax.errorbar(self.z_arr+0.005*i, BF[:, 0], BF[:, 1:].T,
                        fmt="o", color=colors[i], label=label)

        ax.legend(loc="upper right", fontsize=12, ncol=2, frameon=False)

        fname_out = f"figs/tomo_{parname}.pdf"
        if not os.path.isfile(fname_out):
            plt.savefig(fname_out, bbox_inches="tight")
        if not keep_on:
            plt.close()

    def plot_triangles(self, models, parnames="all", keep_on=False):
        for ibin, z in enumerate(self.z_arr):
            fnames = [f"chains/{model}/{model}_{ibin}/cobaya"
                      for model in models]
            s = [gmc.loadMCSamples(fname, settings={'ignore_rows': 0.3})
                 for fname in fnames]

            # list all free parameters
            if parnames == "all":
                exclude = "chi2"
                p = s[0].getParams()
                params = [par for par in p.__dict__ if exclude not in par]

            gdplot = gplot.get_subplot_plotter()
            gdplot.triangle_plot(s, params, filled=True,
                                 legend_labels=models)

            if len(models) == 1:
                fname_out = f"figs/triang_{models[0]}_{ibin}.pdf"
            else:
                fname_out = f"figs/triang_{ibin}.pdf"
            if not os.path.isfile(fname_out):
                plt.savefig(fname_out, bbox_inches="tight")
            if not keep_on:
                plt.close()
