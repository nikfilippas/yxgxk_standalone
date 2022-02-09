from pyccl.halos import Concentration, MassDef


class ConcentrationDuffy08M500c(Concentration):
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
        self.mass_def = MassDef(500, 'critical')

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
