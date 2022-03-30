import os
import yaml
import numpy as np
import pyccl as ccl
import sacc
from itertools import product
from tqdm import tqdm
from scipy.interpolate import RectBivariateSpline


COSMO_P18 = {"Omega_c": 0.26066676,
             "Omega_b": 0.048974682,
             "h": 0.6766,
             "n_s": 0.9665,
             "sigma8": 0.8102}

class CosmologyPlanck18(ccl.Cosmology):
    """Planck 2018 cosmological parameters."""

    def __init__(self, **kwargs):
        super().__init__(**{**COSMO_P18, **kwargs})


class CosmoHaloModel:
    """Internally define all cosmology and halo model, given a setup file."""

    def __init__(self, base_model):
        fname = f"chains/{base_model}/{base_model}_0/cobaya.input.yaml"
        with open(fname, "r") as stream:
            info = yaml.safe_load(stream)
            like = next(iter(info["likelihood"].values()))
            theo = next(iter(info["theory"].values()))
            info = {**like, **theo}

        # setup halo model objects
        hal = ccl.halos
        cM_class = hal.Concentration.from_name(info["cm_name"])
        mf_class = hal.MassFunc.from_name(info["mf_name"])
        hb_class = hal.HaloBias.from_name(info["hb_name"])

        self.mass_def = hal.MassDef.from_name(info["mdef_name"])()
        self.c_m_relation = cM_class(mass_def=self.mass_def)
        self.prof_g = hal.HaloProfileHOD(c_m_relation=self.c_m_relation)
        self.prof_k = hal.HaloProfileNFW(c_m_relation=self.c_m_relation)
        self.prof_y = hal.HaloProfilePressureGNFW()
        self.mass_function = mf_class(mass_def=self.mass_def)
        self.halo_bias = hb_class(mass_def=self.mass_def)
        self.hmc = hal.HMCalculator(
            mass_function=self.mass_function,
            halo_bias=self.halo_bias,
            mass_def=self.mass_def)

        # set up cosmology
        tf = info["transfer_function"]
        self.cosmo = CosmologyPlanck18(transfer_function=tf)

    def update_parameters(self, *, sigma8=None, lMmin_0=None, lM1_0=None,
                          mass_bias=None, x_out=None, mass_function=None):
        if sigma8 is not None:
            self.cosmo = ccl.CosmologyPlanck18(sigma8=sigma8)
        if lMmin_0 is not None:
            self.prof_g.update_parameters(lMmin_0=lMmin_0, lM0_0=lMmin_0)
        if lM1_0 is not None:
            self.prof_g.update_parameters(lM1_0=lM1_0)
        if mass_bias is not None:
            self.prof_y.update_parameters(mass_bias=mass_bias)
        if x_out is not None:
            self.prof_y.update_parameters(x_out=x_out)
        if mass_function is not None:
            mf_class = ccl.halos.MassFunc.from_name(mass_function)
            self.mass_function = mf_class(c_m_relation=self.c_m_relation)
            self.hmc = ccl.halos.HMCalculator(
                mass_function=self.mass_function,
                halo_bias=self.halo_bias,
                mass_def=self.mass_def)


class Interpolator:
    """Enhanced dict type storing interpolators."""

    def __init__(self, dic):
        self.dic = dic
        self.parameters = set(self.dic.keys())
        self.mass_functions = set()
        self.tracers = set()

        for par in self.parameters:
            if not self.mass_functions:
                self.mass_functions = set(self.dic[par].keys())
            else:
                assert set(self.dic[par].keys()) == self.mass_functions
                for mf in self.mass_functions:
                    if not self.tracers:
                        self.tracers = set(self.dic[par][mf].keys())
                    else:
                        assert set(self.dic[par][mf].keys()) == self.tracers

    def validate(self, parameter=None, mass_function=None, tracer=None):
        """Validate that the input columns have been interpolated."""
        flag = True
        if parameter is not None:
            flag = flag and (set([parameter]) <= self.parameters)
        if mass_function is not None:
            flag = flag and (set([mass_function]) <= self.mass_functions)
        if tracer is not None:
            flag = flag and (set([tracer]) <= self.tracers)
        return flag


class Interpolate(CosmoHaloModel):
    """Handles interpolation of expensive parameters."""

    def __init__(self, base_model="gyksrA_T08"):
        super().__init__(base_model=base_model)
        self._save = "data/interpolators.npy"
        self._load()

    def _load(self):
        if os.path.isfile(self._save):
            dic = np.load(self._save, allow_pickle=True).item()
            self.interps = Interpolator(dic)
        else:
            pass  # TODO

    def _interp_bPe(self, z):
        bpe = ccl.halos.halomod_bias_1pt(
            self.cosmo, self.hmc,
            k=1e-3, a=1/(1+z),
            prof=self.prof_y, normprof=False)
        return 1e3 * bpe

    def _interp_Pe(self, z):
        pe = ccl.halos.halomod_mean_profile_1pt(
            self.cosmo ,self.hmc,
            k=1e-3, a=1/(1+z),
            prof=self.prof_y, normprof=False)
        return pe

    def _interp_Omth(self, z):
        pe = self.calculate_Pe(z)
        Y = 0.24
        prefac = (8-5*Y)/(4-2*Y)
        rho_th = pe*prefac/(1+z)**3
        # rho_critical in eV/cm^3
        rho_crit = 10537.0711*self.cosmo['h']**2
        return rho_th/rho_crit

    def _interpolate(self, par, mass_functions,
                    s8_min=0.2, s8_max=1.5, N_s8=16,
                    bH_min=0.005, bH_max=1.15, N_bH=16,):
        """Interpolate input parameter in the s8/bH grid over all input
        mass functions.
        """
        func = getattr(self, f"_interp_{par}")  # get derived par function
        s8_arr = np.linspace(s8_min, s8_max, N_s8)
        bH_arr = np.linspace(bH_min, bH_max, N_bH)

        # F_interp is the outer dictionary which will hold the interp grids
        # categorized by mass function
        F_interp = dict.fromkeys(mass_functions)
        for mf in mass_functions:
            self.update_parameters(mass_function=mf)

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

    def get_interpolator(self, parameter, mass_function, tracer):
        if self.interps.validate(parameter, mass_function, tracer):
            return self.interps.dic[parameter][mass_function][tracer]
        else:
            pass  # TODO


class Container(Interpolate):
    """Populate subclasses with useful parameters."""

    def __init__(self,
                 base_model="gyksrA_T08",
                 fname_sacc="data/saccfiles/cls_cov.fits",
                 secondaries=["YMILCA", "KAPPA"]):
        self._secondaries = secondaries
        self._saccfile = sacc.Sacc.load_fits(fname_sacc)

        self.cosmo = CosmologyPlanck18()
        self.tracers = [f"LOWZ__{i}" for i in range(6)]
        self._get_dndz()        # defines self.dndz
        self._get_zmid()        # defines self.zmid
        self._get_lmax()        # defines self.lmax
        self._build_corrmats()  # defines self.corrmats
        super().__init__(base_model=base_model)

    def _get_dndz(self):
        # new dndz
        self.dndz = dict.fromkeys(self.tracers)
        names = ["2mpz"] + [f"wisc{i}" for i in range(1, 6)]
        for name, store in zip(names, self.dndz.keys()):
            fname = f"data/dndz/{name}_DIR.txt"
            self.dndz[store] = np.loadtxt(fname)

        # old dndz
        self.dndz_old = dict.fromkeys(self.tracers)
        names = ["2MPZ_bin1"] + ["WISC_bin%d" % n for n in range(1, 6)]
        for name, store in zip(names, self.dndz_old.keys()):
            fname = f"data/dndz/{name}.txt"
            self.dndz_old[store] = np.loadtxt(fname)

    def _get_zmid(self):
        dic = self.dndz.copy()
        for name, value in dic.items():
            z, nz = value.T
            dic[name] = np.average(z, weights=nz)
        self.zmid = dic

    def _get_lmax(self):
        kmax = dict(zip(self.tracers, [0.5, 1., 1., 1., 1., 1.]))
        dic = self.zmid.copy()
        for name, zmid in self.zmid.items():
            chi = self.cosmo.comoving_radial_distance(1/(1+zmid))
            dic[name] = np.floor(kmax[name] * chi - 0.5)
        self.lmax = dic
        self.lmin = dict(zip(self.lmax.keys(), [0., 10., 10., 10., 10., 10.]))

    def _get_cov(self, tracer1, tracer2, tracer3, scalecut=False):
        ell, _, ind_AB = self._saccfile.get_ell_cl(
            None, tracer1, tracer2, return_ind=True)
        _, _, ind_CD = self._saccfile.get_ell_cl(
            None, tracer1, tracer3, return_ind=True)
        cov_ABCD = self._saccfile.covariance.covmat[ind_AB][:, ind_CD]
        # impose scalecut
        lmin, lmax = self.lmin[tracer1], self.lmax[tracer1]
        if scalecut:
            idx = np.where((lmin <= ell) & (ell <= lmax))[0]
            return cov_ABCD[idx[0]: idx[-1], idx[0]: idx[-1]]
        else:
            return cov_ABCD

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
        self.corrmats[tracer] = corr

    def _build_corrmats(self):
        self.corrmats = dict.fromkeys(self.tracers)
        for tracer in self.corrmats:
            self._build_corr_block(tracer)
