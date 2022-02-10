import os
import sys
import warnings
from tqdm import tqdm

import sacc
import yaml
import numpy as np
import pyccl as ccl
from cobaya.model import get_model
from scipy.interpolate import RectBivariateSpline

import getdist.mcsamples as gmc
import getdist.plots as gplot
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

import _names


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
        self.hmc = None

        # redshift distributions
        self.names = ["2mpz"] + ["wisc%d" % i for i in range(1, 6)]
        self.DIR = dict.fromkeys(self.names)
        self.pz_old = dict.fromkeys(self.names)
        self.get_z_arr()

        # interpolators
        self.interp_param_names = ["bPe", "Omth"]
        self.interpolators = self.get_interpolators(new_interps)

        # latex
        self.latex_bins = _names.latex_bins
        self.latex_names = _names.latex_names
        self.latex_labels = _names.latex_labels

    def get_z_arr(self):
        fnames = [f"data/dndz/{name}_DIR.txt" for name in self.names]

        z_arr = []
        for name, fname in zip(self.names, fnames):
            z, nz = np.loadtxt(fname).T
            z_arr.append(np.average(z, weights=nz))
            self.DIR[name] = np.column_stack((z, nz))
        self.z_arr = np.array(z_arr)

        # get old dndz
        names = ["2MPZ_bin1"] + ["WISC_bin%d" % n for n in range(1, 6)]
        fnames = [f"data/dndz/{name}.txt" for name in names]
        for name, fname in zip(self.names, fnames):
            self.pz_old[name] = np.loadtxt(fname)

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
        if not os.path.isfile("data/interpolators.npy") or new_interps:
            I = dict.fromkeys(self.interp_param_names)
            for par in tqdm(self.interp_param_names):
                I[par] = self.interpolate_param(par)
            np.save("data/interpolators.npy", I)
        else:
            I = np.load("data/interpolators.npy", allow_pickle=True).item()
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
        mass_functions = ["Tinker08", "Despali16", "Tinker10", "Bocquet16"]
        s8_arr = np.linspace(s8_min, s8_max, N_s8)
        bH_arr = np.linspace(bH_min, bH_max, N_bH)

        # loop through all redshifts within the boundaries
        F_interp = dict.fromkeys(mass_functions)
        for mf in mass_functions:
            self.hmc = ccl.halos.HMCalculator(mass_function=mf,
                                              halo_bias="Tinker10",
                                              mass_def="500c")
            F_temp = dict.fromkeys(self.names)
            for name, z in tqdm(zip(self.names, self.z_arr)):
                Arr = np.zeros((N_s8, N_bH))
                for i, s8 in enumerate(s8_arr):
                    self.update_parameters(sigma8=s8)
                    for j, bH in enumerate(bH_arr):
                        self.update_parameters(mass_bias=bH)
                        Arr[i, j] = func(z)

                F_temp[name] = RectBivariateSpline(s8_arr, bH_arr, Arr)
            F_interp[mf] = F_temp

        return F_interp

    def get_summary(self, model, parname):
        path = model
        model = model.split("/")[-1]
        latex = self.latex_names[parname]

        BF_arr = np.zeros((6, 3))
        for ibin, z in enumerate(tqdm(self.z_arr)):
            # load set-up
            fname_setup = f"chains/{path}/{model}_{ibin}/cobaya.input.yaml"
            with open(fname_setup, "r") as f:
                config = yaml.safe_load(f)
            try:
                mf = config["likelihood"]["ygk_like.ygkLike"]["mf_name"]
            except KeyError:
                # old chains
                mf = config["likelihood"]["yxgxk_like.YxGxKLike"]["mf_name"]

            # load chain
            fname_chains = f"chains/{path}/{model}_{ibin}/cobaya"
            s = gmc.loadMCSamples(fname_chains, settings={'ignore_rows': 0.3})
            p = s.getParams()

            if parname in self.interpolators:
                rel = self.interpolators[parname][mf][self.names[ibin]]

                if hasattr(p, "sigma8"):
                    # free sigma8
                    deriv_chain = np.array([rel(s8, by) for s8, by
                        in zip(p.sigma8, p.ygk_mass_bias)]).squeeze()
                else:
                    # fixed sigma8
                    deriv_chain = np.array([rel(self.cosmo["sigma8"], by)
                        for by in p.ygk_mass_bias]).squeeze()

                s.addDerived(deriv_chain, name=parname, label=latex)

            dens = s.get1DDensity(parname)
            vmin, vmax = dens.getLimits(0.68)[:2]
            vbf = np.average(dens.getLimits(0.01)[:2])
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

    def plot_tomo(self, models, parname, keep_on=False, overwrite=False):
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
            print(model)
            BF = self.get_summary(model=model, parname=parname)
            model = model.split("/")[-1]
            label = self.latex_labels[model]
            ax.errorbar(self.z_arr+0.005*i, BF[:, 0], BF[:, 1:].T,
                        fmt="o", color=colors[i], label=label)

        ax.legend(loc="upper right", fontsize=12, ncol=2, frameon=False)

        hash_ = hash("".join(models)) + sys.maxsize + 1
        fname_out = f"figs/tomo_{parname}_{hash_}.pdf"
        if overwrite or not os.path.isfile(fname_out):
            fig.savefig(fname_out, bbox_inches="tight")
        if not keep_on:
            plt.close()

    def plot_triangles(self, models, parnames="all",
                       keep_on=False, overwrite=False):
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
                hash_ = hash(models[0]) + sys.maxsize + 1
                fname_out = f"figs/triang_{ibin}_{hash_}.pdf"
            else:
                hash_ = hash("".join(models) + sys.maxsize + 1)
                fname_out = f"figs/triang_{ibin}_{hash_}.pdf"
            if overwrite or not os.path.isfile(fname_out):
                plt.savefig(fname_out, bbox_inches="tight")
            if not keep_on:
                plt.close()

    def _get_MxN_axes(self, nrows, ncols):
        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

        figsize = (3*ncols, 4*nrows)
        fig = plt.figure(figsize=figsize)
        gs_main = GridSpec(nrows, ncols, figure=fig)

        axes = np.empty((ncols, nrows, 2), dtype=object)
        for col in range(ncols):
            for row in range(nrows):
                gs = GridSpecFromSubplotSpec(
                    2, 1, subplot_spec=gs_main[row, col],
                    height_ratios=[3, 1], hspace=0)
                for s in range(2):
                    sharey = None if col == 0 else axes[0, row, s]
                    sharex = None if row == 0 else axes[col, 0, 0]

                    ax = fig.add_subplot(gs[s], sharex=sharex, sharey=sharey)
                    axes[col, row, s] = ax

                    if (col, s) == (0, 0):
                        ax.set_ylabel(r"$C_{\ell}$", fontsize=16)
                    elif (col, s) == (0, 1):
                        ax.set_ylabel(r"$\Delta_{\ell}$", fontsize=16)
                    else:
                        ax.yaxis.set_visible(False)
                    if (row == nrows-1) and (s == 1):
                        ax.set_xlabel(r"$\ell$", fontsize=16)
                    else:
                        ax.xaxis.set_visible(False)
                    if s == 0:
                        ax.loglog()
                    elif s == 1:
                        ax.semilogx()
                        ax.axhline(0, c="k", ls=":", lw=3)

        axcol = axes[-1, :, 0]  # last column
        corr_names = ["g", "y", "\kappa"]  # correlation names
        for ax, name in zip(axcol, corr_names):
            ax.text(1.1, 0.33, r"$\mathbf{g \times %s}$" % name,
                    fontsize=16, rotation=-90,
                    horizontalalignment="center",
                    verticalalignment="center",
                    transform=ax.transAxes)

        fig.tight_layout(w_pad=0, h_pad=0)
        return fig, axes

    def plot_best_fit(self, model, keep_on=False, overwrite=False, out=False):
        fname_data = "cls_cov.fits"
        s_data = sacc.Sacc.load_fits(fname_data)
        fig, axes = self._get_MxN_axes(nrows=3, ncols=6)

        for ibin, _ in enumerate(tqdm(self.z_arr)):
            fname = f"chains/{model}/{model}_{ibin}/cobaya.input.yaml"
            with open(fname, "r") as stream:
                info = yaml.safe_load(stream)
            _ = [info.pop(key) for key in ["output", "sampler"]]
            mod = get_model(info)

            # get best fit
            fname = f"chains/{model}/{model}_{ibin}/cobaya"
            s = gmc.loadMCSamples(fname, settings={'ignore_rows': 0.3})
            p = s.getParams()
            p_bf = dict.fromkeys(
                par for par in p.__dict__.keys() if "chi2" not in par)

            argmin = p.chi2.argmin()
            for par in p_bf:
                p_bf[par] = getattr(p, par)[argmin]

            # get theory sacc object
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                loglikes, derived = mod.loglikes(p_bf)
            try:
                l = mod.likelihood['ygk_like.ygkLike']
            except KeyError:
                # old chains
                l = mod.likelihood['yxgxk_like.YxGxKLike']
            params = l.current_state['params'].copy()
            params.update(p_bf)
            s_pred = l.get_sacc_file(**params)

            # get and plot arrays
            chi2 = p.chi2.min()
            dof = 0
            xcorrs = s_pred.get_tracer_combinations()
            for nc, (t1, t2) in enumerate(xcorrs):
                l_t, cl_t = s_pred.get_ell_cl("cl_00", t1, t2)
                l_t_1h, cl_t_1h = s_pred.get_ell_cl("cl_1h", t1, t2)
                l_t_2h, cl_t_2h = s_pred.get_ell_cl("cl_2h", t1, t2)
                l_d, cl_d, cov = s_data.get_ell_cl(None, t1, t2,
                                                   return_cov=True)
                err = np.sqrt(np.diag(cov))

                # scale cuts
                lmin, lmax = l_t.min(), l_t.max()
                idx = np.where((lmax >= l_d) & (l_d >= lmin))[0]
                l_d, cl_d, err = l_d[idx], cl_d[idx], err[idx]

                # residuals
                res = (cl_d - cl_t) / err

                ax_cl, ax_dl = axes[ibin, nc]
                ax_cl.errorbar(l_d, cl_d, err, fmt="ro", ms=3)
                ax_cl.plot(l_t_1h, cl_t_1h, "darkred", alpha=0.3,
                           label=r"$\mathrm{1}$-$\mathrm{halo}$")
                ax_cl.plot(l_t_2h, cl_t_2h, "navy", alpha=0.3,
                           label=r"$\mathrm{2}$-$\mathrm{halo}$")
                ax_cl.plot(l_t, cl_t, "k-",
                           label=r"$\mathrm{1h+2h}$")
                ax_dl.errorbar(l_d, res, np.ones_like(res), fmt="ro", ms=3)
                handles, labels = ax_cl.get_legend_handles_labels()

                dof += len(cl_d)

            # display stats for this bin
            ax = axes[ibin, 0, 0]
            this_bin = self.latex_bins[ibin]
            this_stats = "$\\chi^2/N_{\\rm{d}}=%.2lf/%d$" % (chi2, dof)
            text = "\n".join([this_bin, this_stats])
            ax.text(0.02, 0.04, text, transform=ax.transAxes)

        fig.legend(handles, labels, fontsize="xx-large", ncol=3,
                   loc="center", bbox_to_anchor=(0.5, 1.0), frameon=False)

        fname_out = f"figs/bestfit_{model}.pdf"
        if overwrite or not os.path.isfile(fname_out):
            fig.savefig(fname_out, bbox_inches="tight")

        if not keep_on:
            plt.close()
        if out:
            return fig, axes

    def plot_nz(self, normed=False, compare=True, keep_on=False, overwrite=False):
        if not normed:
            import healpy as hp
            sizes = [476422, 3458260, 3851322, 4000017, 3412366, 1296810]
            fsky = np.mean(hp.read_map("data/maps/mask_v3.fits.gz",
                                       dtype=float))
            area = 4*np.pi*fsky*(180/np.pi)**2

        # set-up figure
        from matplotlib import cm
        cols = ['r']
        cm_wisc = cm.get_cmap('Blues')
        for i in np.arange(5) :
            cols.append(cm_wisc(0.2+((i+1.)/5.)*0.8))

        fig, ax = plt.subplots(tight_layout=True)
        ax.set_xlabel('$z$', fontsize=14)
        ax.set_ylabel('$dN/dz\\,d\\Omega\\,\\,[10^2\\,{\\rm deg}^{-2}]$',
                      fontsize=14)
        ax.tick_params(labelsize="large")

        # plot
        for i, name in enumerate(self.names):
            z, nz = self.DIR[name].T.copy()
            if not normed:
                nz *= sizes[i] / area
            ax.plot(z, nz/100, lw=2, c=cols[i], label=self.latex_bins[i])

            if compare:
                z, nz = self.pz_old[name].T.copy()
                if not normed:
                    nz *= sizes[i] / area
                ax.plot(z, nz/100, lw=1.5, c=cols[i], ls="--")

        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, )
        ax.legend(loc='upper right', ncol=1, frameon=False, fontsize=14)

        fname_out = f"figs/nzs_{compare}.pdf"
        if overwrite or not os.path.isfile(fname_out):
            fig.savefig(fname_out, bbox_inches="tight")

        if not keep_on:
            plt.close()
