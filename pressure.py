import numpy as np
from scipy.integrate import simps
import pyccl.halos as hal
import pickle

from utils import Container


class BattagliaCalculator(Container):
    """Wrapper for the Battaglia profile post-predictions."""
    with open("data/Omgr_Battaglia", "rb") as f:
        rbfl = pickle.load(f)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.cosmo.compute_growth()
        # effective halo model stuff using 200c mass definition
        self._prof_y_b = hal.HaloProfilePressureBattaglia()
        cm = hal.ConcentrationDuffy08()
        self._prof_k_b = hal.HaloProfileNFW(c_m_relation=cm)
        self._hmc_b = hal.HMCalculator(
            mass_function=self.hmc.mass_function.name,
            halo_bias=self.hmc.halo_bias.name,
            mass_def="200c")

    def _get_bPe(self, z, n_r):
        self._hmc_b._get_ingredients(1/(1+z), self.cosmo, True)
        cprof = self._prof_y_b.profile_cumul_nr(
            n_r, self.cosmo, self._hmc_b._mass, 1/(1+z), self._hmc_b.mass_def)
        return self._hmc_b._integrate_over_mbf(cprof)

    def get_bPe(self, z_arr, n_r=100):
        """Vectorized version of `_get_bPe`."""
        return np.array([self._get_bPe(z, n_r) for z in z_arr])*1e3

    def _get_Omth(self, z, n_r):
        self._hmc_b._get_ingredients(1/(1+z), self.cosmo, False)
        cprof = self._prof_y_b.profile_cumul_nr(
            n_r, self.cosmo, self._hmc_b._mass, 1/(1+z), self._hmc_b.mass_def)
        pe = self._hmc_b._integrate_over_mf(cprof)
        Y = 0.24
        prefac = (8-5*Y) / (4-2*Y)
        rho_th = pe * prefac * (1/(1+z))**3
        # rho_critical in eV/cm^3
        rho_crit = 10537.0711 * self.cosmo['h']**2
        return rho_th / rho_crit

    def get_Omth(self, z_arr, n_r):
        """Vectorized version of `_get_Omth`."""
        return np.array([self._get_Omth(z, n_r) for z in z_arr])

    def _get_Omgr(self, cosmo, z):
        cosmo.compute_growth()
        ks = np.geomspace(1E-4, 100, 256)
        pk = hal.halomod_power_spectrum(
            cosmo, self._hmc_b, ks, 1/(1+z), self._prof_k_b,
            normprof=True, get_1h=True, get_2h=False)
        pkint = simps(ks*pk, x=np.log(ks), axis=-1)
        Om = cosmo["Omega_c"] + cosmo["Omega_b"]
        Obh2 = cosmo["Omega_b"] * cosmo["h"]**2
        return 3 * Om *Obh2 * 1.11265006E-7 * pkint * (1+z) / (16*np.pi**2)

    def get_Omgr(self, sigma8, z):
        """Use an RBF emulator with a clipped kernel to calculate Omega_grav."""
        points = np.atleast_2d(np.c_[sigma8, z])
        return self.rbfl(points)


class ArnaudCalculator(Container):
    """Theory wrapper for the Arnaud shock-heating model (not used anymore)."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.cosmo.compute_growth()
        self.update_parameters(mass_bias=0.73)  # sensible value
        self.Delta = 500

    def _get_hmc_eff(self, Delta):
        """Effective CCL Halo Model Calculator for given oDelta."""
        mass_def = hal.MassDef.from_name(f"{Delta}c")()
        mf_class = hal.MassFunc.from_name(self.hmc.mass_function.name)
        hb_class = hal.HaloBias.from_name(self.hmc.halo_bias.name)
        mass_function = mf_class(mass_def=mass_def)
        halo_bias = hb_class(mass_def=mass_def)
        return hal.HMCalculator(mass_function=mass_function,
                                halo_bias=halo_bias, mass_def=mass_def)

    def _get_bPe(self, z):
        hmc = self._get_hmc_eff(self.Delta)
        bpe = hal.halomod_bias_1pt(self.cosmo, hmc,
                                   k=1e-3, a=1/(1+z),
                                   prof=self.prof_y, normprof=False)
        return bpe

    def get_bPe(self, z_arr, n_r=100):
        """Vectorized version of `_get_bPe`.
        Units: meV/cm^3.
        Note: Not used in this pipeline.
        """
        self.update_parameters(x_out=n_r)
        return np.array([self._get_bPe(z) for z in z_arr])*1e3

    def _get_Pe(self, z):
        pe = hal.halomod_mean_profile_1pt(self.cosmo ,self.hmc,
                                          k=1e-3, a=1/(1+z),
                                          prof=self.prof_y, normprof=False)
        return pe

    def get_Pe(self, z_arr, n_r=100):
        self.update_parameters(x_out=n_r)
        return np.array([self._get_Pe(z) for z in z_arr])

    def get_Omth(self, z_arr, n_r=100):
        pe = self.get_Pe(z_arr, n_r=n_r)
        Y = 0.24
        prefac = (8-5*Y)/(4-2*Y)
        rho_th = pe*prefac/(1+z_arr)**3
        # rho_critical in eV/cm^3
        rho_crit = 10537.0711*self.cosmo["h"]**2
        return rho_th/rho_crit
