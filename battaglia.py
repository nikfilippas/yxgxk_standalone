import numpy as np
from scipy.integrate import simps, quad
import pyccl as ccl


class BattagliaSchockHeating(object):
    """Theory wrapper for the Battaglia shock heating model.
    """

    def __init__(self, cosmo, z_arr):
        self.cosmo = cosmo
        self.z_arr = z_arr

    def get_theory(self, n_r):
        return np.array([self.get_bpe(z, n_r, 200) for z in self.z_arr])*1e6

    def get_battaglia(self, m,z,delta) :
        """Sets all parameters needed to compute the Battaglia et al. profile.
        """
        cosmo = self.cosmo
        fb=cosmo.cosmo.params.Omega_b/cosmo.cosmo.params.Omega_m
        ez2=(ccl.h_over_h0(cosmo,1/(1+z)))**2
        h=cosmo.cosmo.params.h
        mr=m*1E-14
        p0=18.1*mr**0.154*(1+z)**(-0.758)
        rDelta=self.R_Delta(m,1./(1+z),Delta=delta)*(1+z)
        dic={'ups0':0.518*p0*2.52200528E-19*delta*m*h**2*ez2*fb*(1+z)/rDelta,
             'rDelta':rDelta,
             'xc':0.497*(mr**(-0.00865))*((1+z)**0.731),
             'beta':4.35*(mr**0.0393)*((1+z)**0.415)}
        return dic

    def ups_battaglia(self, x, bp) :
        """Battaglia pressure profile in units of pressure x = r/rDelta ."""
        xr = x/bp['xc']
        return bp['ups0']*(xr**(-0.3))*((1+xr)**(-bp['beta']))

    def integrated_profile(self, bp, n_r) :
        """Volume integral of the Battaglia pressure profile."""
        integrand = lambda x: x**2*self.ups_battaglia(x,bp)
        return 4*np.pi*(bp['rDelta'])**3 * quad(integrand, 0, n_r)[0]

    def get_Eth_logmass(self, z, n_r, delta, nmass=256):
        lmarr = np.linspace(8.,16.,nmass)
        marr = 10.**lmarr
        et = np.array([
            self.integrated_profile(
                self.get_battaglia(m,z,delta),n_r)
            for m in marr])
        return et, lmarr

    def get_thermal_energy(self, z, n_r, delta):
        et = self.get_Eth_logmass(z, n_r, delta)
        return simps(et[0], x=et[1])

    def get_bpe(self, z, n_r, delta, nmass=256):
        cosmo = self.cosmo
        a = 1./(1+z)
        Dm = delta/ccl.omega_x(cosmo, a, "matter")  # CCL uses Delta_m

        et = self.get_Eth_logmass(z, n_r, delta, nmass=nmass)
        marr = 10**et[1]
        mfunc = ccl.massfunc(cosmo, marr, a, Dm)
        bh = ccl.halo_bias(cosmo, marr, a, Dm)
        return simps(et[0]*bh*mfunc,x=et[1])

    def get_Om_th(self, z, n_r, delta):
        E_th = self.get_thermal_energy(z, n_r, delta)
        rho_crit = 1.054e4 * self.cosmo["h"]**2
        return E_th/rho_crit

    def R_Delta(self, M, a, Delta=500):
        """Calculate the reference radius of a halo."""
        cosmo = self.cosmo
        # Input handling
        M, a = np.atleast_1d(M), np.atleast_1d(a)
        c1 = (cosmo["h"] * ccl.h_over_h0(cosmo, a))**2
        prefac = 1.16217766e12 * Delta * c1
        R = (M[..., None]/prefac)**(1/3)
        return R.squeeze()
