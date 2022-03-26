from plotting import Plotter
from _names import names

p = Plotter()
kwargs = {"keep_on": False, "overwrite": True}

# sigma8
m_s8 = [names[key] for key in [0, 2, 8, 3]]
p.tomographic(m_s8, "sigma8", **kwargs)

# gastrophysics
m_bh = [names[key] for key in [0, 8, 4]]
p.tomographic(m_bh, "ygk_mass_bias", **kwargs)
p.tomographic(m_bh, "bPe", **kwargs)
p.tomographic(m_bh, "Omth", **kwargs)
p.close_plots()

# gastrophysics mass functions
m_bh = [names[key] for key in [0, 10, 11, 7]]
p.tomographic(m_bh, "ygk_mass_bias", **kwargs)
p.tomographic(m_bh, "bPe", **kwargs)
p.tomographic(m_bh, "Omth", **kwargs)
p.close_plots()

# systematics
m_sys = [names[key] for key in [0, 1, 5, 6, 10]]
p.tomographic(m_sys, "sigma8", **kwargs)
p.tomographic(m_sys, "ygk_mass_bias", **kwargs)
p.tomographic(m_sys, "bPe", **kwargs)
p.close_plots()

# triangles fiducial
m_tri = [names[key] for key in [0, 8]]
p.posterior(m_tri, **kwargs)
p.posterior(m_tri, params=["ygk_mass_bias", "sigma8"], **kwargs)
p.close_plots()

# SZ-deprojected
m_tri = [names[key] for key in [0, 6]]
p.posterior(m_tri, **kwargs)
p.posterior(m_tri, params=["ygk_mass_bias", "sigma8"], **kwargs)
p.close_plots()

# triangles s8/bH mass functions
m_tri = [names[key] for key in [0, 10, 11, 7]]
p.posterior(m_tri, params=["ygk_mass_bias", "sigma8"], **kwargs)
p.close_plots()

# best fits
p.best_fit(names[0], **kwargs)
p.close_plots()

# dndz
p.nz(compare=True, **kwargs)
p.close_plots()

# correlation matrices
p.corr_matrices()

# results table
T = p.table(params=["z", "sigma8", "mass_bias", "bPe", "Omth", "chi2/Nd", "PTE"])
print(T)
