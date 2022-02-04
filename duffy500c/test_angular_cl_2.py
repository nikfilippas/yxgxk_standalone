from tqdm import tqdm
import numpy as np
import pyccl as ccl
from fit_duffy import ConcentrationDuffyFit500c, used_fit
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

def get_cosmo(sigma8):
    cosmo = ccl.Cosmology(
        Omega_c=0.25, Omega_b=0.05,
        h=0.67, n_s=0.96, sigma8=sigma8)
    cosmo.compute_sigma()
    return cosmo

def get_tracer(cosmo, z, nz, bz):
    return ccl.NumberCountsTracer(
        cosmo, dndz=(z,nz), bias=(z,bz), has_rsd=False)

def get_colors():
    from matplotlib import cm
    cols = ['r']
    cm_wisc = cm.get_cmap('Blues')
    for i in np.arange(5) :
        cols.append(cm_wisc(0.2+((i+1.)/5.)*0.8))
    return cols

def get_fname(index):
    dirdat = "/home/nick/Desktop/yxgxk_standalone/data/dndz/"
    zbins = ["2mpz"] + ["wisc%d" % i for i in range(1, 6)]
    fname = dirdat + zbins[index] + "_DIR.txt"
    return fname

def get_nz_bz(fname):
    z, nz = np.loadtxt(fname).T
    bz = np.ones_like(z)
    return z, nz, bz

def get_zmid(z, nz):
    return np.average(z, weights=nz)

def get_lmax_from_kmax(cosmo, kmax, zmid):
    chi = ccl.comoving_radial_distance(cosmo, 1./(1+zmid))
    lmax = kmax * chi - 0.5
    return lmax

def compute_cells(cosmo, hmc, tr, profs):
    Cls = np.empty((2, ells.size))
    for i, prof in enumerate(profs):
        pka = ccl.halos.halomod_Pk2D(cosmo, hmc, prof, normprof=True)
        Cls[i] = cosmo.angular_cl(tr, tr, ell=ells, p_of_k_a=pka)
    return Cls

def get_frac_diff(cells):
    return np.abs(1 - np.divide(*cells))


# define boundaries
ells = np.geomspace(6, 1400, 32)
s8_arr = np.linspace(0.7, 1.1, 2)
kmaxs = [0.5] + np.ones(5).tolist()

# cosmologies
cosmoP18 = get_cosmo(sigma8=0.81)
cosmos = [get_cosmo(s8) for s8 in s8_arr]

# halo model stuff
def get_hmc():
    mass_def = ccl.halos.MassDef500c()
    hmc = ccl.halos.HMCalculator(
        mass_function="Tinker08", halo_bias="Tinker10", mass_def=mass_def)
    return hmc

def get_profs():
    A, B, C = used_fit()
    pars = {"A": A, "B": B, "C": C}
    D500 = ConcentrationDuffyFit500c(**pars)
    I21 = ccl.halos.ConcentrationIshiyama21()
    profs = [ccl.halos.HaloProfileHOD(c_m_relation=c) for c in [D500, I21]]
    return profs

def get_devs():
    devs = np.empty((6, 2, ells.size))
    for i in tqdm(range(6)):
        fname = get_fname(i)
        z, nz, bz = get_nz_bz(fname)
        zmid = get_zmid(z, nz)
        lmax = get_lmax_from_kmax(cosmoP18, kmax=kmaxs[i], zmid=zmid)

        for j, s8 in enumerate(s8_arr):
            cosmo = cosmos[j]
            tr = get_tracer(cosmo, z, nz, bz)
            cells = compute_cells(cosmo, hmc, tr, profs)
            devs[i, j] = get_frac_diff(cells)
    return devs


hmc = get_hmc()
profs = get_profs()
devs = get_devs()

def plot(devs, save=True):
    fig, ax = plt.subplots(figsize=(9, 7))
    cols = get_colors()
    latex_bins = [r"$\mathrm{2MPZ}$"] + \
        [r"$\mathrm{WI \times SC}$ - $\mathrm{%d}$" % i
         for i in range(1, 6)]
    ax.axhline(1e-2, ls="--", lw=2, c="k")
    for i in range(6):
        fname = get_fname(i)
        z, nz, bz = get_nz_bz(fname)
        zmid = get_zmid(z, nz)
        lmax = get_lmax_from_kmax(cosmoP18, kmax=kmaxs[i], zmid=zmid)

        ax.fill_between(ells[ells < lmax],
                        devs[i, 0][ells < lmax],
                        devs[i, 1][ells < lmax],
                        color=cols[i],
                        alpha=0.25,
                        label=latex_bins[i])
        ax.loglog()

    ax.legend(loc="lower right", fontsize=17, ncol=2, frameon=False)
    ax.set_xlim(ells.min(), lmax)
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.set_xlabel(
        r"$\ell$",
        fontsize=32)
    ax.set_ylabel(
        r"$\left| 1 - \frac{C_{\ell}^{\rm D08,500c}}{C_{\ell}^{\rm I21}} \right|$",
        fontsize=32)
    fig.tight_layout()
    if save:
        fig.savefig("../figs/compare_angular_cl.pdf", bbox_inches="tight")

plot(devs)
