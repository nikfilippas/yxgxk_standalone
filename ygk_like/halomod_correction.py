import warnings
import pyccl as ccl
import numpy as np
from scipy.interpolate import interp2d


class HalomodCorrection(object):
    """Provides methods to estimate the correction to the halo
    model in the 1h - 2h transition regime.

    Args:
        mass_function, halo_bias, mass_def, concentration: halo model params
        k_range (list): range of k to use (in Mpc^-1).
        nlk (int): number of samples in log(k) to use.
        z_range (list): range of redshifts to use.
        nz (int): number of samples in redshift to use.
    """

    def __init__(self,
                 mass_function=None, halo_bias=None,
                 mass_def=None, concentration=None,
                 k_range=[1E-1, 5], nlk=20,
                 z_range=[0., 1.], nz=16):

        cosmo = ccl.CosmologyVanillaLCDM(transfer_function="bacco")
        lkarr = np.linspace(np.log10(k_range[0]),
                            np.log10(k_range[1]),
                            nlk)
        karr = 10.**lkarr
        zarr = np.linspace(z_range[0], z_range[1], nz)

        # halo model power spectrum
        if None in [mass_function, halo_bias, mass_def, concentration]:
            # CCL default
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")  # CCL deprecation warning
                pk_hm = np.array([ccl.halomodel_matter_power(cosmo, karr, a)
                                  for a in 1. / (1 + zarr)])
        else:
            hmc = ccl.halos.HMCalculator(mass_function=mass_function,
                                         halo_bias=halo_bias,
                                         mass_def=mass_def)
            hmd = ccl.halos.MassDef.from_name(mass_def)()
            cMc = ccl.halos.Concentration.from_name(concentration)
            cM = cMc(mass_def=hmd)
            prof = ccl.halos.HaloProfileNFW(c_m_relation=cM)
            pk_hm = np.array([ccl.halos.halomod_power_spectrum(
                cosmo, hmc, karr, a, prof, normprof=True)
                for a in 1/(1+zarr)])

        pk_hf = np.array([ccl.nonlin_matter_power(cosmo, karr, a)
                          for a in 1. / (1 + zarr)])
        ratio = pk_hf / pk_hm

        self.rk_func = interp2d(lkarr, 1/(1+zarr), ratio,
                                bounds_error=False, fill_value=1)

    def rk_interp(self, k, a):
        """
        Returns the halo model correction for an array of k
        values at a given redshift.

        Args:
            k (float or array): wavenumbers in units of Mpc^-1.
            a (float): value of the scale factor.
        """
        return self.rk_func(np.log10(k), a)-1
