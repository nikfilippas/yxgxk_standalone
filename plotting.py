import os
import copy
import numpy as np
import healpy as hp
from tqdm import tqdm
import getdist.mcsamples as gmc

import matplotlib.pyplot as plt
import getdist.plots as gplot
from matplotlib import cm
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.legend_handler import HandlerBase

import _names
from chanal import ChainCalculator
from pressure import ArnaudCalculator, BattagliaCalculator

from pyccl import Hashing
hashF = Hashing._hash_consistent

plt.rcParams["text.usetex"] = True


class AnyObjectHandler(HandlerBase):
    # Adapted from
    # https://matplotlib.org/users/legend_guide.html#legend-handlers
    def create_artists(self, legend, orig_handle,
                       x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0,y0+width], [0.7*height,0.7*height],
                        lw=3, linestyle="-", color="cadetblue")
        l2 = plt.Line2D([x0,y0+width], [0.3*height,0.3*height],
                        lw=3, linestyle="-", color="grey")
        return [l1, l2]


class Plotter(ChainCalculator):
    """Handles plotting."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._zarr = np.asarray(list(self.zmid.values()))
        self._zplot = np.linspace(0.01, 0.40, 64)
        self._setup_plot_names()

    def _setup_plot_names(self):
        self._corr_names = _names.gxg, _names.gxy, _names.gxk
        self._latex_bins = _names.latex_bins
        self._latex_names = _names.latex_names
        self._latex_labels = _names.latex_labels
        self._latex_short = _names.latex_labels_short
        self._colors = _names.colors
        self._markers = _names.markers

    def _overlay_sigma8(self, ax):
        s8 = self.cosmo.sigma8() * self.cosmo.growth_factor(1/(1+self._zplot))

        # get the Planck uncertainty band for sigma8
        from utils import COSMO_P18, CosmologyPlanck18
        par_min, par_max = COSMO_P18.copy(), COSMO_P18.copy()
        par_min["sigma8"] -= 0.0060
        par_max["sigma8"] += 0.0060
        cosmo_min = CosmologyPlanck18(**par_min)
        cosmo_max = CosmologyPlanck18(**par_max)
        s8_min = cosmo_min.sigma8() * cosmo_min.growth_factor(1/(1+self._zplot))
        s8_max = cosmo_max.sigma8() * cosmo_max.growth_factor(1/(1+self._zplot))

        ax.fill_between(self._zplot, s8_min, s8_max, color="k", alpha=0.3)
        ax.plot(self._zplot, s8, "k--", lw=2, label="Planck 2018")

    def _overlay_bPe(self, ax):
        bPe = BattagliaCalculator().get_bPe
        func = lambda n_r: bPe(self._zplot, n_r)

        et2, et3, et5, etinf = [func(n_r) for n_r in [2, 3, 5, 100]]

        ax.plot(self._zplot, et2, '-',label='$r_{\\rm max}=2\\,r_{200c}$', c='k')
        ax.plot(self._zplot, et3, '--',label='$r_{\\rm max}=3\\,r_{200c}$', c='k')
        ax.plot(self._zplot, et5, '-.',label='$r_{\\rm max}=5\\,r_{200c}$', c='k')
        ax.plot(self._zplot, etinf, ':',label='$r_{\\rm max}=\\infty$', c='k')

        # finally, we plot existing measurements of bPe
        from data.bpe_data_from_bibliography import bpe_data
        all_data = np.hstack(bpe_data)
        err = [all_data[1] - all_data[2], all_data[3] - all_data[1]]
        # plot error bars without markers
        ax.errorbar(all_data[0], all_data[1], yerr=err,
                    fmt="none", color="royalblue", alpha=0.2,
                    label="bibliography (see caption)")
        # plot name initials from bibliography
        biblio = ["Pandey19", "Chiang20", "Waerbeke13", "Yan20",
                  "Koukoufilippas20", "Vikram17"]
        for dat, bib in zip(bpe_data, biblio):
            ax.plot(dat[0], dat[1], "o", marker="$%s$" % bib[0],
                    color="royalblue", alpha=0.3, markersize=7)
        ax.set_ylim(0.07, 0.28)

    def _overlay_Omth(self, ax):
        Omth = ArnaudCalculator().get_Omth(self._zplot, n_r=100)
        # Omth = BattagliaCalculator().get_Omth(self._zplot, n_r=100)  # TODO: use Battaglia
        ax.plot(self._zplot, Omth, '-', label='GNFW profile', c='k')

        # Now we overlay the measurement from Fukugita & Peebles
        # (digitized from Chiang et al. 2020)
        # Multiply it by 1e8, and transpose it to the smallest redshift we use.
        ax.errorbar(
            self._zarr[0], [2.000], yerr=[[0.414], [0.490]],
            fmt="s", color="b", markersize=8,
            label=r"FP04, -$(2/3) f_b \Omega_{\mathrm{grav}, 0}$")


    def _overlay_ygk_mass_bias(self, ax):
        ax.axhline(0.72, ls=":", color="cadetblue")
        ax.axhspan(0.72-0.10, 0.72+0.10, color="cadetblue", alpha=0.3)
        ax.axhline(0.58, ls=":", color="grey")
        ax.axhspan(0.58-0.04, 0.58+0.06, color="grey", alpha=0.3)
        props = {"boxstyle": "round", "facecolor": "w", "alpha": 0.5}
        ax.text(0.396, 0.73, "CMB $\\kappa$ + N.C.",
                fontsize=10, fontweight="bold",
                horizontalalignment="right", verticalalignment="bottom",
                bbox=props, transform=ax.transData)
        ax.text(0.396, 0.59, "CMB + N.C.",
                fontsize=10, fontweight="bold",
                horizontalalignment="right", verticalalignment="bottom",
                bbox=props, transform=ax.transData)

    def tomographic(self, models, par, violins=[0],
                    keep_on=False, overwrite=False):
        fig, ax = plt.subplots(figsize=(9,7))
        ax.tick_params(labelsize=16)
        ax.set_xlabel("$z$", fontsize=20)
        ax.set_ylabel(f"${self._latex_names[par]}$", fontsize=20)
        ax.set_xlim(0.05, 0.40)
        fig.tight_layout()

        # First, overlay other plots
        func = f"_overlay_{par}"
        if hasattr(self, func):
            getattr(self, func)(ax)

        for i, model in enumerate(models):
            print(model)
            BF = self.get_summary(model, par).copy()
            if par == "Omth": BF *= 1e8
            label = self._latex_labels[model]
            ax.errorbar(self._zarr+0.005*i, BF[:, 0], BF[:, 1:].T,
                        fmt=self._markers[model], color=self._colors[model],
                        label=label)

        if violins is not None and par not in ["bPe", "Omth"]:
            vmods = [models[i] for i in violins]
            for ii, model in enumerate(vmods):
                chains = np.array(
                    [ibin[par] for ibin in self.get_chains(model).values()],
                    dtype=object)
                # normalize for y-scale if needed
                if par == "sigma8":
                    chains *= np.array([self.cosmo.growth_factor(1/(1+z))
                                        for z in self._zarr])
                elif par == "Omth": chains *= 1e8

                # plot violins
                vplot = ax.violinplot(
                    chains, self._zarr+0.005*violins[ii],
                    widths=0.01, showextrema=False)
                for pc in vplot["bodies"]:
                    pc.set_color(self._colors[model])
                    pc.set_alpha(0.2)

        # sophisticated legend
        handles, labels = ax.get_legend_handles_labels()
        handler_map = None
        if par == "ygk_mass_bias":
            handles = [object] + handles
            labels = ["Planck 2015"] + labels
            handler_map = {object: AnyObjectHandler()}
        ncol = 1 if len(handles) < 5 else 2
        kw = {"handles": handles, "labels": labels, "handler_map": handler_map,
              "loc": "best", "fontsize": 14, "ncol": ncol, "frameon": False}
        ax.legend(**kw)

        hash_ = hashF("".join(models))
        fname_out = f"figs/tomo_{par}_{hash_}.pdf"
        if overwrite or not os.path.isfile(fname_out):
            fig.savefig(fname_out, bbox_inches="tight")
        plt.show(block=False)
        if not keep_on:
            plt.close()

    def posterior(self, models, params="all", bins="all",
                  keep_on=False, overwrite=False):
        params_in = params
        tracers = self.tracers.copy()
        tracers = tracers if bins == "all" else [tracers[i] for i in bins]
        ibins = range(len(tracers)) if bins == "all" else bins
        for ibin in tqdm(ibins):
            fnames = [f"{self._chains_dir}/{model}/{model}_{ibin}/cobaya"
                      for model in models]
            s = [gmc.loadMCSamples(fname, settings={'ignore_rows': 0.3})
                 for fname in fnames]

            # list all free parameters
            if params_in == "all":
                exclude = "chi2"
                p = s[0].getParams()
                params = [par for par in p.__dict__ if exclude not in par]

            lbl = self._latex_labels if len(params) > 2 else self._latex_short
            gdplot = gplot.get_subplot_plotter()
            gdplot.triangle_plot(
                s, params, filled=True,
                legend_labels=[lbl[model] for model in models],
                contour_colors=[self._colors[model] for model in models],
                frameon=False)

            if len(models) == 1:
                hash_ = hashF(str(models[0]) + str(params))
                fname_out = f"figs/triang_{ibin}_{hash_}.pdf"
            else:
                hash_ = hashF("".join(models) + str(params))
                fname_out = f"figs/triang_{ibin}_{hash_}.pdf"
            if overwrite or not os.path.isfile(fname_out):
                plt.savefig(fname_out, bbox_inches="tight")
            plt.show(block=False)
            if not keep_on:
                plt.close()

    def _setup_MxN_axes(self, nrows, ncols):
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
        for ax, corr_name in zip(axcol, self._corr_names):
            ax.text(1.1, 0.33, rf"{corr_name}",
                    fontsize=16, rotation=-90,
                    horizontalalignment="center",
                    verticalalignment="center",
                    transform=ax.transAxes)

        fig.tight_layout(w_pad=0, h_pad=0)
        return fig, axes

    def best_fit(self, model, keep_on=False, overwrite=False, out=False):
        fig, axes = self._setup_MxN_axes(nrows=3, ncols=len(self.tracers))

        for ibin, tracer in enumerate(tqdm(self.tracers)):
            s_pred = self.get_sacc_file(model, tracer)

            chi2 = self.get_best_fit(model)[tracer]["chi2"]
            dof = 0
            xcorrs = s_pred.get_tracer_combinations()
            for nc, (t1, t2) in enumerate(xcorrs):
                l_t, cl_t = s_pred.get_ell_cl("cl_00", t1, t2)
                l_t_1h, cl_t_1h = s_pred.get_ell_cl("cl_1h", t1, t2)
                l_t_2h, cl_t_2h = s_pred.get_ell_cl("cl_2h", t1, t2)
                l_d, cl_d, cov = self._saccfile.get_ell_cl(None, t1, t2,
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
            this_bin = self._latex_bins[ibin]
            this_stats = "$\\chi^2/N_{\\rm{d}}=%.2lf/%d$" % (chi2, dof)
            text = "\n".join([this_bin, this_stats])
            ax.text(0.02, 0.04, text, transform=ax.transAxes)

        fig.legend(handles, labels, fontsize="xx-large", ncol=3,
                   loc="center", bbox_to_anchor=(0.5, 1.0), frameon=False)

        fname_out = f"figs/bestfit_{model}.pdf"
        if overwrite or not os.path.isfile(fname_out):
            fig.savefig(fname_out, bbox_inches="tight")

        plt.show(block=False)
        if not keep_on:
            plt.close()
        if out:
            return fig, axes

    def nz(self, normed=False, compare=True, keep_on=False, overwrite=False):
        if not normed:
            import healpy as hp
            sizes = [476422, 3458260, 3851322, 4000017, 3412366, 1296810]
            fsky = np.mean(hp.read_map("data/maps/mask_v3.fits.gz",
                                       dtype=float))
            area = 4*np.pi*fsky*(180/np.pi)**2

        # set-up figure
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
        for i, tracer in enumerate(self.tracers):
            z, nz = self.dndz[tracer].T.copy()
            if not normed:
                nz *= sizes[i] / area
            ax.plot(z, nz/100, lw=2, c=cols[i], label=self._latex_bins[i])

            if compare:
                z, nz = self.dndz_old[tracer].T.copy()
                if not normed:
                    nz *= sizes[i] / area
                ax.plot(z, nz/100, lw=1.5, c=cols[i], ls="--")

        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, )
        ax.legend(loc='upper right', ncol=1, frameon=False, fontsize=14)

        fname_out = f"figs/nzs_{compare}.pdf"
        if overwrite or not os.path.isfile(fname_out):
            fig.savefig(fname_out, bbox_inches="tight")

        plt.show(block=False)
        if not keep_on:
            plt.close()

    def _mpl_corr_block(self, tracer, overwrite=False, keep_on=True):
        corr = self.corrmats[tracer]

        fig, ax = plt.subplots()
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

        kw = lambda vert: {"fontsize": 24,
                           "horizontalalignment": "center",
                           "verticalalignment": "center",
                           "rotation": "vertical" if vert else "horizontal",
                           "transform": ax.transAxes}

        o = 0.05
        ax.text(0.17, 1+o, r"$g \times g$", **kw(False))
        ax.text(0.50, 1+o, r"$g \times y$", **kw(False))
        ax.text(0.83, 1+o, r"$g \times \kappa$", **kw(False))
        ax.text(-o, 0.83, r"$g \times g$", **kw(True))
        ax.text(-o, 0.50, r"$g \times y$", **kw(True))
        ax.text(-o, 0.17, r"$g \times \kappa$", **kw(True))

        ax.imshow(corr,
                  cmap=cm.gray, vmin=0, vmax=1,
                  interpolation="nearest", aspect="equal")

        fname_out = f"figs/corr_{tracer}.pdf"
        if overwrite or not os.path.isfile(fname_out):
            fig.savefig(fname_out, bbox_inches="tight")
        plt.show(block=False)
        if not keep_on:
            plt.close()

    def corr_matrices(self, overwrite=False, keep_on=True):
        for tracer in self.corrmats.keys():
            self._mpl_corr_block(tracer, overwrite=overwrite, keep_on=keep_on)

    def close_plots(self):
        plt.close("all")

    def table(self, model, params):
        samples = ["2MPZ"] + ["\\wisc-%d" % i for i in range(1, 6)]
        _ = [self.get_summary(model, par) for par in params]
        bf = copy.deepcopy(self.summary[model])

        s = ""
        for i, (sample, tracer) in enumerate(zip(samples, self.tracers)):
            s += sample + "  & "
            s += "%.2f" % self.zmid[tracer] + "  & "
            for par in params:
                vals = bf[par][i]
                if par == "Omth":
                    vals *= 1e8
                s += "\\vals(%.2f, %.2f, %.2f)" % (vals[0], vals[1], vals[2])
                s += "  & "
            s += "\\\\\n"
        return s

    @classmethod
    def _plot_map(cls, fname, nside=256, return_data=False):
        fname_use = os.path.join("data/maps", fname)
        data = hp.ud_grade(hp.read_map(fname_use, dtype=float), nside_out=nside)

        fig = plt.figure(num="mollview", figsize=(12, 7))
        kw = {"map"   : data,
              "fig"   : "mollview",
              "title" : None,
              "xsize" : 1200,
              "cmap"  : cm.gray,
              "cbar"  : False}

        hp.mollview(**kw)
        hp.graticule(c="maroon", alpha=0.4)
        if not return_data:
            return fig
        return fig, data

    @classmethod
    def plot_masks(cls, keep_on=False, overwrite=False):
        fnames = {"mask_g": 'mask_v3.fits.gz',
                  "mask_y": 'mask_planck60.fits.gz',
                  "mask_k": 'COM_Lensing_4096_R3.00_mask.fits.gz'}

        nside = 256
        combined = np.ones(hp.nside2npix(nside))  # total mask footprint
        individual = np.ones(3)  # individual mask footprint
        for i, (mask_name, fname) in enumerate(fnames.items()):
            fig, new_data = cls._plot_map(fname, nside=nside, return_data=True)
            combined *= new_data
            individual[i] = new_data.mean()

            fname_out = f"figs/{mask_name}.pdf"
            if overwrite or not os.path.isfile(fname_out):
                fig.savefig(fname_out)
            if not keep_on:
                plt.close()
        return individual.tolist() + [combined.mean()]
