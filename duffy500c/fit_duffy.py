import numpy as np
import pyccl as ccl
from pyccl.halos import Concentration, MassDef
from scipy.optimize import minimize
import matplotlib.pyplot as plt


class ConcentrationDuffyFit500c(Concentration):
    """
    Fit parameters A, B, C for Duffy08 with a mass definition of 500c.
    """

    def __init__(self, mass_def=None, **pars):
        self.A = pars.get("A")
        self.B = pars.get("B")
        self.C = pars.get("C")
        super().__init__(mass_def=mass_def)

    def _default_mass_def(self):
        self.mass_def = MassDef(500, 'critical')

    def _check_mass_def(self, mass_def):
        if (mass_def.Delta != 500) or (mass_def.rho_type != 'critical'):
            return True
        return False

    def _concentration(self, cosmo, M, a):
        M_pivot_inv = cosmo.cosmo.params.h * 5E-13
        return self.A * (M * M_pivot_inv)**self.B * a**(-self.C)


class Likelihood(object):
    """
    """

    def __init__(self):
        self.M_arr = np.logspace(7, 16, 64)
        self.a_arr = np.linspace(0.4, 1., 16)

        self.cosmo = ccl.CosmologyVanillaLCDM()
        self.mass_def = ccl.halos.MassDef.from_name("500c")()
        mfc = ccl.halos.MassFunc.from_name("Tinker08")
        mf = mfc(mass_def=self.mass_def)
        self.mass_function = np.array([
            mf.get_mass_function(self.cosmo, self.M_arr, a)
            for a in self.a_arr])

        self.c0 = ccl.halos.ConcentrationIshiyama21(mass_def=self.mass_def)
        self.pars = dict(A=3.67, B=-0.0903, C=-0.51)
        self.c1 = ConcentrationDuffyFit500c(mass_def=self.mass_def,
                                            **self.pars)

    def _init_duffy(self, **pars):
        if pars == self.pars:
            return
        self.pars = pars
        self.c1 = ConcentrationDuffyFit500c(mass_def=self.mass_def,
                                            **pars)

    def _get_statistic(self, con):
        stat = np.array([con.get_concentration(self.cosmo, self.M_arr, a)
                         for a in self.a_arr]) * self.mass_function
        return stat

    def chi2(self, **pars):
        self._init_duffy(**pars)
        v0 = self._get_statistic(self.c0)
        v1 = self._get_statistic(self.c1)
        chi2 = np.sum((1 - np.log(v1) / np.log(v0))**2)
        print(chi2)
        return chi2

    def minfunc(self, args):
        dic = dict(zip(["A", "B", "C"], args))
        return self.chi2(**dic)


def plot(lik):
    stat = lik._get_statistic(lik.c0) / lik._get_statistic(lik.c1)
    plt.imshow(np.log10(np.fabs(1 - stat)),
               extent=[np.log10(lik.M_arr.min()), np.log10(lik.M_arr.max()),
                       lik.a_arr.min(), lik.a_arr.max()],
               aspect="auto")
    plt.colorbar()
    plt.tight_layout()


def run_minimizer(lik):
    func = lik.minfunc
    p0 = list(lik.pars.values())
    res = minimize(func, p0)
    return res


def run_pipe():
    lik = Likelihood()
    res = run_minimizer(lik)
    plot(lik)
    print(res.x.tolist())


def best_fit():
    return [5.7220309114546595, -0.08053973326724174, -0.8526595438331771]
