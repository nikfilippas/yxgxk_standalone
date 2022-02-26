import numpy as np
import sacc
import pyccl as ccl
from itertools import product


class Container:
    """
    """

    def __init__(self,
                 fname_sacc="data/saccfiles/cls_cov.fits",
                 secondaries=["YMILCA", "KAPPA"]):
        self._secondaries = secondaries
        self._saccfile = sacc.Sacc.load_fits(fname_sacc)
        self._setup()

    def _setup(self):
        self.zbin_names = [f"LOWZ__{i}" for i in range(6)]
        self._names = ["2mpz"] + [f"wisc{i}" for i in range(1, 6)]
        self._get_dndz()  # defines self.dndz
        self._get_zmid()  # defines self.zmid
        self._get_lmax()  # defines self.lmax

    def _get_dndz(self):
        dic = dict.fromkeys(self.zbin_names)
        for name, store in zip(self._names, dic.keys()):
            fname = f"data/dndz/{name}_DIR.txt"
            dic[store] = np.loadtxt(fname)
        self.dndz = dic

    def _get_zmid(self):
        dic = self.dndz.copy()
        for name, value in dic.items():
            z, nz = value.T
            dic[name] = np.average(z, weights=nz)
        self.zmid = dic

    def _get_lmax(self):
        kmax = dict(zip(self.zbin_names, [0.5, 1., 1., 1., 1., 1.]))
        self._cosmo = ccl.Cosmology(
            Omega_c=0.25, Omega_b=0.045, h=0.67, n_s=0.96, sigma8=0.81)
        dic = self.zmid.copy()
        for name, zmid in self.zmid.items():
            chi = self._cosmo.comoving_radial_distance(1/(1+zmid))
            dic[name] = np.floor(kmax[name] * chi - 0.5)
        self.lmax = dic

    def _get_cov(self, tracer1, tracer2, tracer3, scalecut=False):
        ell, _, ind_AB = self._saccfile.get_ell_cl(
            None, tracer1, tracer2, return_ind=True)
        _, _, ind_CD = self._saccfile.get_ell_cl(
            None, tracer1, tracer3, return_ind=True)
        cov_ABCD = self._saccfile.covariance.covmat[ind_AB][:, ind_CD]
        # impose scalecut
        lmax = self.lmax[tracer1]
        idx = -1 if not scalecut else np.argmin(ell <= lmax)
        return cov_ABCD[:idx, :idx]

    def _build_cov_block(self, tracer):
        names = [tracer] + self._secondaries
        combs = list(product(*[names, names]))

        # collect all required blocks
        cov = [self._get_cov(tracer, tracer1, tracer2, scalecut=True)
               for tracer1, tracer2 in combs]
        # and now assemble them in a square
        n = len(names)
        cov = np.block([cov[n*i: n*(i+1)] for i in range(n)])
        return cov

    def _build_corr_block(self, tracer):
        cov = self._build_cov_block(tracer)
        diag = np.diag(cov)
        corr = cov/np.sqrt(diag[:, None] * diag[None, :])
        return corr

    def _mpl_corr_block(self, tracer, save=False, close=True):
        import matplotlib.pyplot as plt
        from matplotlib import cm
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

        corr = self._build_corr_block(tracer)
        ax.imshow(corr, cmap=cm.gray, interpolation="nearest", aspect="equal")

        if save:
            fname = f"figs/corr_{tracer}.pdf"
            fig.savefig(fname, bbox_inches="tight")
        if close:
            plt.close()

    def build_corr_matrices(self):
        for name in self.zbin_names:
            self._mpl_corr_block(name, save=True)

c = Container()
c.build_corr_matrices()
