from plotting import Plotter
from _names import names

settings = {"base_model": "gyksrA",
            "chains_dir": "chains/chains_new",
            "fname_sacc": "data/saccfiles/cls_cov_new.fits"}
p = Plotter(**settings)
kwargs = {"keep_on": True, "overwrite": True}

# structure growth
m_s8 = [names[key] for key in [0, 1, 2, 3, 17]]
p.tomographic(m_s8, "sigma8", **kwargs)

# gas pressure
m_bh = [names[key] for key in [0, 2, 4, 17]]
p.tomographic(m_bh, "ygk_mass_bias", **kwargs)
p.tomographic(m_bh, "bPe", **kwargs)
p.tomographic(m_bh, "Omth", **kwargs)
p.close_plots()

# free Omega_m
m_om = [names[key] for key in [0, 2, 17]]
p.posterior(m_om, params=["ygk_mass_bias", "sigma8", "Omega_m"])
p.posterior(m_om, params=["all", "Omega_m"])

# # robustness checks
# m_sys = [names[key] for key in [0, 5, 6, 7, 8]]
# p.tomographic(m_sys, "sigma8", **kwargs)
# p.tomographic(m_sys, "ygk_mass_bias", **kwargs)
# p.tomographic(m_sys, "bPe", **kwargs)
# p.close_plots()

# # # mass functions
# mods = [names[key] for key in [0, 8, 9, 10]]
# p.tomographic(mods, "ygk_mass_bias", **kwargs)
# p.tomographic(mods, "bPe", **kwargs)
# p.tomographic(mods, "Omth", **kwargs)
# p.posterior(mods, params=["ygk_mass_bias", "sigma8"], bins=[3], **kwargs)
# p.close_plots()

# # triangles fiducial-gauss
# m_tri = [names[key] for key in [0, 2, 17]]
# p.posterior(m_tri, **kwargs)
# p.posterior(m_tri, params=["ygk_mass_bias", "sigma8"], **kwargs)
# p.close_plots()

# # triangles fiducial-SZdeproj
# m_tri = [names[key] for key in [0, 7]]
# p.posterior(m_tri, **kwargs)
# p.posterior(m_tri, params=["ygk_mass_bias", "sigma8"], **kwargs)
# p.close_plots()

# # best fits
# p.best_fit(names[0], **kwargs)
# p.close_plots()

# # dndz
# p.nz(compare=True, **kwargs)
# p.close_plots()

# # correlation matrices
# p.corr_matrices(**kwargs)

# # results table
# # tab_params = ["sigma8", "ygk_mass_bias", "bPe", "Omth"]
# T = p.table(model="gyksrA_T08", params=tab_params)
# print(T)

# # masks
# Plotter.plot_masks()
