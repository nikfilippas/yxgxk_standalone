import numpy as np
import pyccl as ccl
from tqdm import tqdm
from scipy.interpolate import interp2d
import getdist.mcsamples as gmc


class ArnaudCalculator(object):
    def __init__(self):
        self.cosmo = ccl.CosmologyVanillaLCDM()
        self.prof = ccl.halos.HaloProfilePressureGNFW()
        self.prof.update_parameters(mass_bias=0.8)
        self.hmc = ccl.halos.HMCalculator(mass_function="Tinker08",
                                          halo_bias="Tinker10",
                                          mass_def="500c")
        self.interpolators = {}
        self.z_arr = np.arange(0.07, 0.35, 0.05)

    def update_mass_bias(self, mass_bias):
        self.prof.update_parameters(mass_bias=mass_bias)

    def get_bPe(self, z):
        bpe = ccl.halos.halomod_bias_1pt(self.cosmo, self.hmc,
                                         1E-3, 1./(1+z),
                                         self.prof, normprof=False)
        return bpe*1E3

    def set_cosmo(self, sigma8):
        self.cosmo = ccl.Cosmology(Omega_c=0.26066676,
                                   Omega_b=0.048974682,
                                   h=0.6766,
                                   n_s=0.9665,
                                   sigma8=sigma8)

    def interpolate_bPe(self, z,
                        s8_min=0.2, s8_max=1.5, N_s8=16,
                        bH_min=0.005, bH_max=1.15, N_bH=16):
        s8_arr = np.linspace(s8_min, s8_max, N_s8)
        bH_arr = np.linspace(bH_min, bH_max, N_bH)
        Arr = np.zeros((N_s8, N_bH))
        for i, s8 in enumerate(s8_arr):
            self.set_cosmo(s8)
            for j, bH in enumerate(bH_arr):
                self.update_mass_bias(bH)
                Arr[i, j] = self.get_bPe(z)

        F_interp = interp2d(s8_arr, bH_arr, Arr, kind="cubic")
        self.interpolators[z] = F_interp

    def interpolate_all_redshifts(self):
        for z in tqdm(self.z_arr):
            self.interpolate_bPe(z)

    def from_chains(self, model):
        BF_arr = np.zeros((6, 3))
        for ibin, z in enumerate(tqdm(self.z_arr)):
            rel = self.interpolators[z]
            fname = f"chains/{model}/{model}_{ibin}/cobaya"
            s = gmc.loadMCSamples(fname, settings={'ignore_rows': 0.3})
            p = s.getParams()

            bPe_chain = np.array([rel(s8, by) for s8, by
                                  in zip(p.sigma8, p.ygk_mass_bias)])
            s.addDerived(bPe_chain, name="bPe", label="bP_e")

            dens = s.get1DDensity("bPe")
            vmin, vmax = dens.getLimits(0.68)[:2]
            vbf = dens.getLimits(0.001)[0]
            summary = vbf, vbf-vmin, vmax-vbf
            BF_arr[ibin] = summary
        return BF_arr


c = ArnaudCalculator()
c.interpolate_all_redshifts()
BF_yxgxksig = c.from_chains("yxgxksig")
