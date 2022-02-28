from chanal import ChainCalculator
from chanal2 import Container
from _names import names_new as names

c = ChainCalculator()
keep_on = True
overwrite=True

# sigma8
m_s8 = [names[key] for key in [0, 2, 8, 3]]
c.plot_tomo(m_s8, "sigma8", keep_on=keep_on, overwrite=overwrite)

# gastrophysics
m_bh = [names[key] for key in [0, 4, 8]]
c.plot_tomo(m_bh, "ygk_mass_bias", keep_on=keep_on, overwrite=overwrite)
c.plot_tomo(m_bh, "bPe", keep_on=keep_on, overwrite=overwrite)
c.plot_tomo(m_bh, "Omth", keep_on=keep_on, overwrite=overwrite)
c.close_plots()

# gastrophysics mass functions
m_bh = [names[key] for key in [0, 10, 11, 7]]
c.plot_tomo(m_bh, "ygk_mass_bias", keep_on=keep_on, overwrite=overwrite)
c.plot_tomo(m_bh, "bPe", keep_on=keep_on, overwrite=overwrite)
c.plot_tomo(m_bh, "Omth", keep_on=keep_on, overwrite=overwrite)
c.close_plots()

# systematics
m_sys = [names[key] for key in [0, 1, 5, 6, 10]]
c.plot_tomo(m_sys, "sigma8", keep_on=keep_on, overwrite=overwrite)
c.plot_tomo(m_sys, "ygk_mass_bias", keep_on=keep_on, overwrite=overwrite)
c.plot_tomo(m_sys, "bPe", keep_on=keep_on, overwrite=overwrite)
c.close_plots()

# triangles
m_tri = [names[key] for key in [1, 7]]
c.plot_triangles(m_tri, keep_on=keep_on, overwrite=overwrite)
c.close_plots()

# best fits
c.plot_best_fit(names[0], keep_on=keep_on, overwrite=overwrite)
c.plot_best_fit(names[7], keep_on=keep_on, overwrite=overwrite)
c.plot_best_fit(names[3], keep_on=keep_on, overwrite=overwrite)
c.close_plots()

# dndz
c.plot_nz(compare=keep_on, keep_on=keep_on, overwrite=overwrite)
c.close_plots()

# correlation matrices
c = Container()
c.build_corr_matrices()
