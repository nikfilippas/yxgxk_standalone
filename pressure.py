import numpy as np
from scipy.integrate import quad
import pyccl as ccl
import pyccl.halos as hal

from utils import Container


class BattagliaCalculator(Container):
    """Theory wrapper for the Battaglia shock heating model."""

    def __init__(self):
        super().__init__()
        self.cosmo.compute_growth()

    def get_battaglia(self, m, z, Delta) :
        """Set all parameters needed to compute the Battaglia et al. profile.
        """
        cosmo = self.cosmo
        fb = cosmo["Omega_b"] / cosmo["Omega_m"]
        ez2 = (cosmo.h_over_h0(1/(1+z)))**2
        h = cosmo["h"]
        mr = m*1E-14
        p0 = 18.1 * mr**0.154 * (1+z)**(-0.758)
        rDelta = self.R_Delta(m, 1./(1+z), Delta=Delta) * (1+z)
        dic={'ups0': 0.518*p0*2.52200528E-19*Delta*m*h**2*ez2*fb*(1+z)/rDelta,
             'rDelta': rDelta,
             'xc': 0.497*(mr**(-0.00865))*((1+z)**0.731),
             'beta': 4.35*(mr**0.0393)*((1+z)**0.415)}
        return dic

    def ups_battaglia(self, x, bp) :
        """Battaglia pressure profile in units of pressure x = r/rDelta ."""
        xr = x/bp['xc']
        return bp['ups0'] * (xr**(-0.3)) * ((1+xr)**(-bp['beta']))

    def _integrate_profile(self, bp, n_r) :
        """Volume integral of the Battaglia pressure profile."""
        integrand = lambda x: x**2*self.ups_battaglia(x, bp)
        return 4*np.pi*(bp['rDelta'])**3 * quad(integrand, 0, n_r)[0]

    def _Eth_integrand(self, z, n_r, Delta):
        et = np.array([
            self._integrate_profile(
                self.get_battaglia(m, z, Delta), n_r)
            for m in self.hmc._mass])
        return et

    def _get_hmc_eff(self, Delta):
        """Effective CCL Halo Model calculator for given oDelta."""
        mass_def = hal.MassDef.from_name(f"{Delta}c")()
        mf_class = hal.MassFunc.from_name(self.hmc.mass_function.name)
        hb_class = hal.HaloBias.from_name(self.hmc.halo_bias.name)
        mass_function = mf_class(mass_def=mass_def)
        halo_bias = hb_class(mass_def=mass_def)
        return hal.HMCalculator(mass_function=mass_function,
                                halo_bias=halo_bias, mass_def=mass_def)

    def _get_bPe(self, z, n_r, Delta):
        hmc = self._get_hmc_eff(Delta)
        hmc._get_ingredients(1/(1+z), self.cosmo, True)
        et = self._Eth_integrand(z, n_r, Delta)
        return hmc._integrate_over_mbf(et)

    def get_bPe(self, z_arr, n_r, Delta):
        """Vectorized version of `_get_bPe`.
        Units: meV/cm^3.
        """
        return np.array([self._get_bPe(z, n_r, Delta) for z in z_arr])*1e6

    def _get_Pe(self, z, n_r, Delta):
        hmc = self._get_hmc_eff(Delta)
        hmc._get_ingredients(1/(1+z), self.cosmo, False)
        et = self._Eth_integrand(z, n_r, Delta)
        return hmc._integrate_over_mf(et)

    def get_Pe(self, z_arr, n_r, Delta):
        """Vectorized version of `_get_Eth`.
        Units: meV/cm^3.
        """
        return np.array([self._get_Pe(z, n_r, Delta) for z in z_arr])*1e6

    def get_Omth(self, z_arr, n_r, Delta):
        """
        Units: 1e8*Omth.
        Note: Not used in this pipeline.
        """
        E_th = self.get_Pe(z_arr, n_r, Delta)    # meV/cm^3
        rho_crit = 1.054e7 * self.cosmo["h"]**2  # meV/cm^3
        return 1e8*E_th/rho_crit

    def R_Delta(self, M, a, Delta=500):
        """Calculate the reference radius of a halo."""
        M, a = np.atleast_1d(M), np.atleast_1d(a)
        c1 = (self.cosmo["h"] * self.cosmo.h_over_h0(a))**2
        prefac = 1.16217766e12 * Delta * c1
        R = (M[..., None]/prefac)**(1/3)
        return R.squeeze()


class ArnaudCalculator(Container):
    """Theory wrapper for the Arnaud shock-heating model."""

    def __init__(self):
        super().__init__()
        self.cosmo.compute_growth()
        self.prof_y.update_parameters(mass_bias=0.73)  # sensible value

    def _get_bPe(self, z):
        bpe = hal.halomod_bias_1pt(self.cosmo, self.hmc,
                                   k=1e-3, a=1/(1+z),
                                   prof=self.prof_y, normprof=False)
        return bpe

    def get_bPe(self, z_arr):
        """Vectorized version of `_get_bPe`.
        Units: meV/cm^3.
        Note: Not used in this pipeline.
        """
        return np.array([self._get_bPe(z) for z in z_arr])*1e3

    def _get_Pe(self, z):
        pe = hal.halomod_mean_profile_1pt(self.cosmo ,self.hmc,
                                          k=1e-3, a=1/(1+z),
                                          prof=self.prof_y, normprof=False)
        return pe

    def get_Pe(self, z_arr):
        return np.array([self._get_Pe(z) for z in z_arr])

    def get_Omth(self, z_arr):
        """Units: 1e8*Omth."""
        pe = self.get_Pe(z_arr)
        Y = 0.24
        prefac = (8-5*Y)/(4-2*Y)
        rho_th = pe*prefac/(1+z_arr)**3
        # rho_critical in eV/cm^3
        rho_crit = 10537.0711*self.cosmo["h"]**2
        return 1e8*rho_th/rho_crit
