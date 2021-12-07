import warnings
import numpy as np
from scipy.interpolate import interp1d
import pyccl as ccl
import pyccl.nl_pt as pt
from cobaya.likelihood import Likelihood
from cobaya.log import LoggedError


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
        from scipy.interpolate import interp2d

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
            hmd = ccl.halos.mass_def_from_name(mass_def)()
            cMc = ccl.halos.concentration_from_name(concentration)
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


def beam_gaussian(ll, fwhm_amin):
    """
    Returns the SHT of a Gaussian beam.
    Args:
        l (float or array): multipoles.
        fwhm_amin (float): full-widht half-max in arcmins.
    Returns:
        float or array: beam sampled at `l`.
    """
    sigma_rad = np.radians(fwhm_amin / 2.355 / 60)
    return np.exp(-0.5 * ll * (ll + 1) * sigma_rad**2)


def beam_hpix(ll, ns):
    """
    Returns the SHT of the beam associated with a HEALPix
    pixel size.
    Args:
        l (float or array): multipoles.
        ns (int): HEALPix resolution parameter.
    Returns:
        float or array: beam sampled at `l`.
    """
    fwhm_hp_amin = 60 * 41.7 / ns
    return beam_gaussian(ll, fwhm_hp_amin)


class ConcentrationDuffy08M500c(ccl.halos.Concentration):
    """ Concentration-mass relation by Duffy et al. 2008
    (arXiv:0804.2486) extended to Delta = 500-critical.
    Args:
        mass_def (:class:`~pyccl.halos.massdef.MassDef`): a mass
            definition object that fixes
            the mass definition used by this c(M)
            parametrization.
    """
    name = 'Duffy08M500c'

    def __init__(self, mass_def=None):
        super(ConcentrationDuffy08M500c, self).__init__(mass_def=mass_def)

    def _default_mass_def(self):
        self.mass_def = ccl.halos.MassDef(500, 'critical')

    def _check_mass_def(self, mass_def):
        if (mass_def.Delta != 500) or (mass_def.rho_type != 'critical'):
            return True
        return False

    def _setup(self):
        self.A = 3.67
        self.B = -0.0903
        self.C = -0.51

    def _concentration(self, cosmo, M, a):
        M_pivot_inv = cosmo.cosmo.params.h * 5E-13
        return self.A * (M * M_pivot_inv)**self.B * a**(-self.C)


class YxGxKLike(Likelihood):
    # All parameters starting with this will be
    # identified as belonging to this stage.
    input_params_prefix: str = ""
    # Input sacc file
    input_file: str = ""
    # IA model
    ia_model: str = "IANone"
    # N(z) model name
    nz_model: str = "NzNone"
    # b(z) model name
    bz_model: str = "BzNone"
    # Shape systamatics
    shape_model: str = "ShapeNone"
    # Mass function name
    mf_name: str = "Tinker08"
    # Halo bias name
    hb_name: str = "Tinker10"
    # Mass definition
    mdef_name: str = "500c"
    # Concentration name
    cm_name: str = "Duffy08M500c"
    # zmax for HM calculator
    zmax_hmc: float = 0.6
    # #z for HM calculator
    nz_hmc: int = 16
    # #k for HM calculator
    nk_hmc: int = 256
    # k HM shot noise suppression
    k_1h_suppress: float = 0.01
    # List of bin names
    bins: list = []
    # List of default settings (currently only scale cuts)
    defaults: dict = {}
    # List of two-point functions that make up the data vector
    twopoints: list = []
    # Angular resolution
    nside: int = -1
    # M0-mode
    M0_track: bool = True
    # HM correction
    HM_correction: str = "HMCode"
    # allow satellites to form when no central is present
    ns_independent: bool = False

    def initialize(self):
        # Read SACC file
        self._read_data()
        # Ell sampling for interpolation
        self._get_ell_sampling()
        # Other global parameters
        self._init_globals()

    def _par(self, name, default, pars, altname=None, altdefault=None):
        """Helper method for quick parameter extraction."""
        prefix = self.input_params_prefix
        parname = "_".join(filter(None, [prefix, name]))
        par = pars.get(parname, default)
        if par is None:
            parname = "_".join(filter(None, [prefix, altname]))
            par = pars.get(parname, altdefault)
        return par

    def _init_globals(self):
        self.qabbr = {'galaxy_density': 'g',
                      'galaxy_shear': 'm',
                      'cmb_convergence': 'm',
                      'cmb_tSZ': 'y'}
        self.beam_pix = beam_hpix(self.l_sample, self.nside)

        if self.HM_correction == 'halofit':
            self.hmcorr = HalomodCorrection()
        else:
            # this includes all other options
            self.hmcorr = None

        # Initialize parameterless HM stuff
        if self.bz_model == 'HaloModel':
            self.massdef = ccl.halos.mass_def_from_name(self.mdef_name)()
            self.mfc = ccl.halos.mass_function_from_name(self.mf_name)
            self.hbc = ccl.halos.halo_bias_from_name(self.hb_name)
            cmc = ccl.halos.concentration_from_name(self.cm_name)
            self.cm = cmc(mass_def=self.massdef)
            hal = ccl.halos
            self.profs = {
                'galaxy_density': hal.HaloProfileHOD(
                    c_m_relation=self.cm,
                    ns_independent=self.ns_independent),
                'galaxy_shear': hal.HaloProfileNFW(c_m_relation=self.cm),
                'cmb_convergence': hal.HaloProfileNFW(c_m_relation=self.cm),
                'cmb_tSZ': hal.HaloProfilePressureGNFW()}
            self.p2pt_HOD = ccl.halos.Profile2ptHOD()

    def _read_data(self):
        """
        Reads sacc file
        Selects relevant data.
        Applies scale cuts
        Reads tracer metadata (N(z))
        Reads covariance
        """
        import sacc

        def get_suffix_for_tr(tr):
            q = tr.quantity
            if q in ['galaxy_density', 'cmb_convergence', 'cmb_tSZ']:
                return '0'
            elif q == 'galaxy_shear':
                return 'e'
            else:
                raise ValueError(f'dtype not found for quantity {q}')

        def get_lmax_from_kmax(cosmo, kmax, zmid):
            chi = ccl.comoving_radial_distance(cosmo, 1./(1+zmid))
            lmax = np.max([10., kmax * chi - 0.5])
            return lmax

        s = sacc.Sacc.load_fits(self.input_file)
        self.bin_properties = {}
        cosmo_lcdm = ccl.CosmologyVanillaLCDM(transfer_function="bacco")
        kmax_default = self.defaults.get('kmax', 0.1)
        for b in self.bins:
            if b['name'] not in s.tracers:
                raise LoggedError(self.log, "Unknown tracer %s" % b['name'])
            t = s.tracers[b['name']]
            if t.quantity in ['galaxy_density', 'galaxy_shear']:
                zmid = np.average(t.z, weights=t.nz)
                self.bin_properties[b['name']] = {'z_fid': t.z,
                                                  'nz_fid': t.nz,
                                                  'zmean_fid': zmid}
            else:
                self.bin_properties[b['name']] = {}

            # Ensure all tracers have ell_min
            if b['name'] not in self.defaults:
                self.defaults[b['name']] = {}
                self.defaults[b['name']]['lmin'] = self.defaults['lmin']

            # Give galaxy clustering an ell_max
            if t.quantity == 'galaxy_density':
                # Get lmax from kmax for galaxy clustering
                if 'kmax' in self.defaults[b['name']]:
                    kmax = self.defaults[b['name']]['kmax']
                else:
                    kmax = kmax_default
                lmax = get_lmax_from_kmax(cosmo_lcdm,
                                          kmax, zmid)
                self.defaults[b['name']]['lmax'] = lmax
            else:
                # Make sure everything else has an ell_max
                if 'lmax' not in self.defaults[b['name']]:
                    self.defaults[b['name']]['lmax'] = self.defaults['lmax']

        # First check which parts of the data vector to keep
        indices = []
        for cl in self.twopoints:
            tn1, tn2 = cl['bins']
            lmin = np.max([self.defaults[tn1].get('lmin', 2),
                           self.defaults[tn2].get('lmin', 2)])
            lmax = np.min([self.defaults[tn1].get('lmax', 1E30),
                           self.defaults[tn2].get('lmax', 1E30)])
            # Get the suffix for both tracers
            cl_name1 = get_suffix_for_tr(s.tracers[tn1])
            cl_name2 = get_suffix_for_tr(s.tracers[tn2])
            ind = s.indices('cl_%s%s' % (cl_name1, cl_name2), (tn1, tn2),
                            ell__gt=lmin, ell__lt=lmax)
            indices += list(ind)
        s.keep_indices(np.array(indices))

        # Now collect information about those
        # and put the C_ells in the right order
        indices = []
        self.cl_meta = []
        id_sofar = 0
        self.used_tracers = {}
        self.l_min_sample = 1E30
        self.l_max_sample = -1E30
        for cl in self.twopoints:
            # Get the suffix for both tracers
            tn1, tn2 = cl['bins']
            cl_name1 = get_suffix_for_tr(s.tracers[tn1])
            cl_name2 = get_suffix_for_tr(s.tracers[tn2])
            l, c_ell, cov, ind = s.get_ell_cl('cl_%s%s' % (cl_name1, cl_name2),
                                              tn1,
                                              tn2,
                                              return_cov=True,
                                              return_ind=True)
            if c_ell.size > 0:
                if tn1 not in self.used_tracers:
                    self.used_tracers[tn1] = s.tracers[tn1].quantity
                if tn2 not in self.used_tracers:
                    self.used_tracers[tn2] = s.tracers[tn2].quantity

            bpw = s.get_bandpower_windows(ind)
            if np.amin(bpw.values) < self.l_min_sample:
                self.l_min_sample = np.amin(bpw.values)
            if np.amax(bpw.values) > self.l_max_sample:
                self.l_max_sample = np.amax(bpw.values)

            self.cl_meta.append({'bin_1': tn1,
                                 'bin_2': tn2,
                                 'l_eff': l,
                                 'cl': c_ell,
                                 'cov': cov,
                                 'inds': (id_sofar +
                                          np.arange(c_ell.size,
                                                    dtype=int)),
                                 'l_bpw': bpw.values,
                                 'w_bpw': bpw.weight.T})
            indices += list(ind)
            id_sofar += c_ell.size
        indices = np.array(indices)
        # Reorder data vector and covariance
        self.data_vec = s.mean[indices]
        self.cov = s.covariance.covmat[indices][:, indices]
        # Check eigenvalues
        w, v = np.linalg.eigh(self.cov)
        if np.any(w <= 0):
            print(w)
            exit(1)
            #iw = 1./w
            #iw[w <= 0] = 0
            #self.inv_cov = np.dot(v, np.dot(np.diag(iw), v.T))
            self.inv_cov = np.linalg.inv(self.cov)
        else:
            self.inv_cov = np.linalg.inv(self.cov)
        self.ndata = len(self.data_vec)

    def _get_ell_sampling(self, nl_per_decade=30):
        # Selects ell sampling.
        # Ell max/min are set by the bandpower window ells.
        # It currently uses simple log-spacing.
        # nl_per_decade is currently fixed at 30
        if self.l_min_sample == 0:
            l_min_sample_here = 2
        else:
            l_min_sample_here = self.l_min_sample
        nl_sample = int(np.log10(self.l_max_sample / l_min_sample_here) *
                        nl_per_decade)
        l_sample = np.unique(
            np.geomspace(l_min_sample_here,
                         self.l_max_sample+1,
                         nl_sample).astype(int)).astype(float)

        if self.l_min_sample == 0:
            self.l_sample = np.concatenate((np.array([0.]), l_sample))
        else:
            self.l_sample = l_sample

    def _eval_interp_cl(self, cl_in, l_bpw, w_bpw):
        """ Interpolates C_ell, evaluates it at bandpower window
        ell values and convolves with window."""
        f = interp1d(self.l_sample, cl_in)
        cl_unbinned = f(l_bpw)
        cl_binned = np.dot(w_bpw, cl_unbinned)
        return cl_binned

    def _get_nz(self, cosmo, name, **pars):
        """ Get redshift distribution for a given tracer."""
        z = self.bin_properties[name]['z_fid']
        nz = self.bin_properties[name]['nz_fid']
        zm = self.bin_properties[name]['zmean_fid']
        dz = 0.
        wz = 1.
        if (self.nz_model == 'NzShift') or (self.nz_model == 'NzShiftWidth'):
            dz = self._par(f"{name}_dz", 0., pars)
        if (self.nz_model == 'NzShiftWidth') or (self.nz_model == 'NzWidth'):
            wz = self._par(f"{name}_wz", 1., pars)
        z = zm+dz+(z-zm)/wz
        msk = z >= 0
        z = z[msk]
        nz = nz[msk]
        return (z, nz)

    def _get_bz(self, cosmo, name, **pars):
        """ Get linear galaxy bias. Unless we're using a linear bias,
        this should be just 1."""
        z = self.bin_properties[name]['z_fid']
        zmean = self.bin_properties[name]['zmean_fid']
        bz = np.ones_like(z)
        if self.bz_model == 'Linear':
            b0 = self._par(f"{name}_b0", None, pars)
            bp = self._par(f"{name}_bp", 0., pars)
            bz = b0 + bp * (z - zmean)
        return (z, bz)

    def _get_ia_bias(self, cosmo, name, **pars):
        """ Intrinsic alignment amplitude.
        """
        if self.ia_model == 'IANone':
            return None
        else:
            z = self.bin_properties[name]['z_fid']
            if self.ia_model == 'IAPerBin':
                A = self._par(f"{name}_A_IA", 0., pars)
                A_IA = np.ones_like(z) * A
            elif self.ia_model == 'IADESY1':
                A0 = self._par("A_IA", 0., pars)
                eta = self._par("eta_IA", 0., pars)
                A_IA = A0 * ((1+z)/1.62)**eta
            else:
                raise LoggedError(self.log, "Unknown IA model %s" %
                                  self.ia_model)
            return (z, A_IA)

    def _get_tracers(self, cosmo, **pars):
        """ Transforms all used tracers into CCL tracers for the
        current set of parameters."""
        trs = {}
        for name, q in self.used_tracers.items():
            prefix = self.input_params_prefix + '_' + name
            if self.bz_model == 'HaloModel':
                prof = self.profs[q]
                normed = True

            if q == 'galaxy_density':
                nz = self._get_nz(cosmo, name, **pars)
                bz = self._get_bz(cosmo, name, **pars)
                t = ccl.NumberCountsTracer(cosmo, dndz=nz,
                                           bias=bz, has_rsd=False)
                if self.bz_model == 'EulerianPT':
                    z = self.bin_properties[name]['z_fid']
                    zmean = self.bin_properties[name]['zmean_fid']
                    b1 = pars[prefix + '_b1']
                    b1p = pars.get(prefix + '_b1p', 0.)
                    bz = b1 + b1p * (z - zmean)
                    b2 = pars.get(prefix + '_b2', 0.)
                    bs = pars.get(prefix + '_bs', 0.)
                    ptt = pt.PTNumberCountsTracer(b1=(z, bz), b2=b2, bs=bs)
                if self.bz_model == 'HaloModel':
                    hod_pars = {k: pars[prefix + '_' + k]
                                for k in ['lMmin_0', 'lM1_0']}
                    if self.M0_track:
                        hod_pars['lM0_0'] = hod_pars['lMmin_0']
                    else:
                        hod_pars['lM0_0'] = pars[prefix + '_lM0_0']
                    slM = pars.get(prefix + '_siglM_0', None)
                    if slM is None:
                        slM = self._par("siglM_0", None, pars)
                    hod_pars['siglM_0'] = slM
                    prof.update_parameters(**hod_pars)
            elif q == 'galaxy_shear':
                nz = self._get_nz(cosmo, name, **pars)
                ia = self._get_ia_bias(cosmo, name, **pars)
                t = ccl.WeakLensingTracer(cosmo, dndz=nz, ia_bias=ia)
                if self.bz_model == 'EulerianPT':
                    ptt = pt.PTMatterTracer()
            elif q == 'cmb_convergence':
                # B.H. TODO: pass z_source as parameter to the YAML file
                t = ccl.CMBLensingTracer(cosmo, z_source=1100)
                if self.bz_model == 'EulerianPT':
                    ptt = pt.PTMatterTracer()
            elif q == 'cmb_tSZ':
                t = ccl.tSZTracer(cosmo, z_max=3.)
                if self.bz_model == 'HaloModel':
                    o_m_b = self._par("mass_bias", 1., pars)
                    prof.update_parameters(mass_bias=o_m_b)
                    normed = False
                else:
                    raise NotImplementedError("Can't do tSZ without"
                                              " the halo model.")

            trs[name] = {}
            trs[name]['ccl_tracer'] = t
            if self.bz_model == 'EulerianPT':
                trs[name]['PT_tracer'] = ptt
            if self.bz_model == 'HaloModel':
                trs[name]['Profile'] = prof
                trs[name]['Normed'] = normed
        return trs

    def _get_pk_data(self, cosmo):
        """ Get all cosmology-dependent ingredients to create the
        different P(k)s needed for the C_ell calculation.
        For linear bias, this is just the matter power spectrum.
        """
        # Get P(k)s from CCL
        if self.bz_model == 'Linear':
            cosmo.compute_nonlin_power()
            pkmm = cosmo.get_nonlin_power(name='delta_matter:delta_matter')
            return {'pk_mm': pkmm}
        elif self.bz_model == 'EulerianPT':
            cosmo.compute_nonlin_power()
            pkmm = cosmo.get_nonlin_power(name='delta_matter:delta_matter')
            ptc = pt.PTCalculator(with_NC=True, with_IA=False,
                                  log10k_min=-4, log10k_max=2,
                                  nk_per_decade=20)
            pk_lin_z0 = ccl.linear_matter_power(cosmo, ptc.ks, 1.)
            ptc.update_pk(pk_lin_z0)
            return {'ptc': ptc, 'pk_mm': pkmm}
        elif self.bz_model == 'HaloModel':
            cosmo.compute_linear_power()
            cosmo.compute_nonlin_power()
            cosmo.compute_sigma()
            pkmm = cosmo.get_nonlin_power(name='delta_matter:delta_matter')
            mf = self.mfc(cosmo, mass_def=self.massdef)
            hb = self.hbc(cosmo, mass_def=self.massdef)
            hmc = ccl.halos.HMCalculator(mass_function=mf,
                                         halo_bias=hb,
                                         mass_def=self.massdef)
            return {'hmc': hmc, 'pk_mm': pkmm}
        else:
            raise LoggedError(self.log,
                              "Unknown bias model %s" % self.bz_model)

    def _get_pkxy(self, cosmo, clm, pkd, trs,
                  get_1h=True, get_2h=True,
                  **pars):
        """ Get the P(k) between two tracers. """
        q1 = self.used_tracers[clm['bin_1']]
        q2 = self.used_tracers[clm['bin_2']]

        if (self.bz_model == 'Linear') or (self.bz_model == 'BzNone'):
            if (q1 == 'galaxy_density') and (q2 == 'galaxy_density'):
                return pkd['pk_mm']  # galaxy-galaxy
            elif ((q1 != 'galaxy_density') and (q2 != 'galaxy_density')):
                return pkd['pk_mm']  # matter-matter
            else:
                return pkd['pk_mm']  # galaxy-matter
        elif (self.bz_model == 'EulerianPT'):
            if ((q1 != 'galaxy_density') and (q2 != 'galaxy_density')):
                return pkd['pk_mm']  # matter-matter
            else:
                ptt1 = trs[clm['bin_1']]['PT_tracer']
                ptt2 = trs[clm['bin_2']]['PT_tracer']
                pk_pt = pt.get_pt_pk2d(cosmo, ptt1, tracer2=ptt2,
                                       ptc=pkd['ptc'])
                return pk_pt
        elif self.bz_model == 'HaloModel':
            k_s = np.geomspace(1E-4, 1E2, self.nk_hmc)
            a_s = 1./(1+np.linspace(0., self.zmax_hmc, self.nz_hmc)[::-1])
            comb = self.qabbr[q1]+self.qabbr[q2]
            p1 = trs[clm['bin_1']]['Profile']
            p2 = trs[clm['bin_2']]['Profile']
            norm1 = trs[clm['bin_1']]['Normed']
            norm2 = trs[clm['bin_2']]['Normed']

            if q1 == q2 == 'galaxy_density':
                prof2pt = self.p2pt_HOD
            else:
                r = self._par(f"rho{comb}", 0., pars)
                prof2pt = ccl.halos.Profile2pt(r_corr=r)

            # halo model correction
            if not (get_1h and get_2h) or self.HM_correction != "HMCode":
                fsmooth = None
            elif self.HM_correction == "HMCode":
                alpha = self._par(f"alpha{comb}", None, pars, "alpha", 1.)
                fsmooth = lambda a: alpha

            # 1-halo suppression
            k_1h_suppress = self._par("k_1h_suppress", self.k_1h_suppress, pars)
            fsuppress = (lambda a: k_1h_suppress) if get_1h else None

            pkt = ccl.halos.halomod_power_spectrum(
                cosmo, pkd['hmc'], k_s, a_s,
                p1,
                prof_2pt=prof2pt, prof2=p2,
                normprof=norm1,
                normprof2=norm2,
                get_1h=get_1h, get_2h=get_2h,
                smooth_transition=fsmooth,
                suppress_1h=fsuppress)

            # halo model correction
            if get_1h and get_2h:
                if self.HM_correction == 'halofit':
                    A = self._par(f"Ahmc{comb}", None, pars, "Ahmc", 1.)
                    ratio = 1 + A*self.hmcorr.rk_interp(k_s, a_s)
                    pkt *= ratio

            pk = ccl.Pk2D(a_arr=a_s, lk_arr=np.log(k_s), pk_arr=np.log(pkt),
                          extrap_order_lok=1, extrap_order_hik=2,
                          cosmo=cosmo, is_logp=True)
            return pk
        else:
            raise LoggedError(self.log,
                              "Unknown bias model %s" % self.bz_model)

    def _get_pixel_window(self, clm):
        q1 = self.used_tracers[clm['bin_1']]
        q2 = self.used_tracers[clm['bin_2']]
        pixwin = np.ones(self.l_sample.size)
        # No window for CMB lensing or shear
        # (the latter only in the sampling limit)
        if q1 in ['galaxy_density', 'cmb_tSZ']:
            pixwin *= self.beam_pix
        if q2 in ['galaxy_density', 'cmb_tSZ']:
            pixwin *= self.beam_pix
        return pixwin

    def _get_cl_all(self, cosmo, pk, get_1h=True, get_2h=True, **pars):
        """ Compute all C_ells."""
        # Gather all tracers
        trs = self._get_tracers(cosmo, **pars)

        # Correlate all needed pairs of tracers
        cells = []
        for clm in self.cl_meta:
            pkxy = self._get_pkxy(cosmo, clm, pk, trs,
                                  get_1h=get_1h, get_2h=get_2h,
                                  **pars)
            cl = ccl.angular_cl(cosmo,
                                trs[clm['bin_1']]['ccl_tracer'],
                                trs[clm['bin_2']]['ccl_tracer'],
                                ell=self.l_sample, p_of_k_a=pkxy)
            # Pixel window function
            cl *= self._get_pixel_window(clm)
            clb = self._eval_interp_cl(cl, clm['l_bpw'], clm['w_bpw'])
            cells.append(clb)
        return cells

    def _apply_shape_systematics(self, cells, **pars):
        if self.shape_model == 'ShapeMultiplicative':
            # Multiplicative shear bias
            for i, clm in enumerate(self.cl_meta):
                q1 = self.used_tracers[clm['bin_1']]
                q2 = self.used_tracers[clm['bin_2']]
                if q1 == 'galaxy_shear':
                    m1 = self._par(f"{clm['bin_1']}_m", 0., pars)
                else:
                    m1 = 0.
                if q2 == 'galaxy_shear':
                    m2 = self._par(f"{clm['bin_2']}_m", 0., pars)
                else:
                    m2 = 0.
                prefac = (1+m1) * (1+m2)
                cells[i] *= prefac

    def get_cells_theory(self, get_1h=True, get_2h=True, **pars):
        # Get cosmological model
        res = self.provider.get_CCL()
        cosmo = res['cosmo']

        # First, gather all the necessary ingredients for the different P(k)
        pkd = res['pk_data']

        # Then pass them on to convert them into C_ells
        cells = self._get_cl_all(cosmo, pkd,
                               get_1h=get_1h, get_2h=get_2h,
                               **pars)

        # Multiplicative bias if needed
        self._apply_shape_systematics(cells, **pars)

        return cells

    def get_sacc_file(self, **pars):
        import sacc

        # Create empty file
        s = sacc.Sacc()

        # Add tracers
        for n, p in self.bin_properties.items():
            if n not in self.used_tracers:
                continue
            q = self.used_tracers[n]
            if q in ['galaxy_density', 'galaxy_shear']:
                s.add_tracer('NZ', n, quantity=q, spin=0,
                             z=p['z_fid'], nz=p['nz_fid'])
            else:
                s.add_tracer('Map', n, quantity=q, spin=0,
                             ell=np.arange(10), beam=np.ones(10))

        # Calculate power spectra
        cells = self.get_cells_theory(**pars)
        for clm, cl in zip(self.cl_meta, cells):
            s.add_ell_cl('cl_00', clm['bin_1'], clm['bin_2'], clm['l_eff'], cl)

        if self.bz_model == "HaloModel":
            # include the 1h/2h contributions
            cells_1h = self.get_cells_theory(get_2h=False, **pars)
            cells_2h = self.get_cells_theory(get_1h=False, **pars)
            for clm, cl1, cl2 in zip(self.cl_meta, cells_1h, cells_2h):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    s.add_ell_cl('cl_1h', clm['bin_1'], clm['bin_2'],
                                 clm["l_eff"], cl1)
                    s.add_ell_cl('cl_2h', clm['bin_1'], clm['bin_2'],
                                 clm["l_eff"], cl2)

        size = len(self.cov)
        cov = np.zeros((3*size, 3*size)) + 1e-16  # avoid zeros
        cov[:size, :size] = self.cov
        self.cov = cov
        s.add_covariance(self.cov)
        return s

    def _get_theory(self, **pars):
        """ Computes theory vector."""
        cells = self.get_cells_theory(**pars)

        # Flattening into a 1D array
        cl_out = np.zeros(self.ndata)
        for clm, cl in zip(self.cl_meta, cells):
            cl_out[clm['inds']] = cl

        return cl_out

    def get_requirements(self):
        # By selecting `self._get_pk_data` as a `method` of CCL here,
        # we make sure that this function is only run when the
        # cosmological parameters vary.
        return {'CCL': {'methods': {'pk_data': self._get_pk_data}}}

    def logp(self, **pars):
        """
        Simple Gaussian likelihood.
        """
        t = self._get_theory(**pars)
        r = t - self.data_vec
        chi2 = np.dot(r, np.dot(self.inv_cov, r))
        return -0.5*chi2
