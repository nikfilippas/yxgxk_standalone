from chanal import ChainCalculator
from _names import names_old as names

c = ChainCalculator()


# sigma8
m_s8 = [names[key] for key in [0, 2, 8, 3]]
c.plot_tomo(m_s8, "sigma8", keep_on=False)

# gastrophysics
m_bh = [names[key] for key in [0, 8, 7]]
c.plot_tomo(m_bh, "ygk_mass_bias", keep_on=False)
c.plot_tomo(m_bh, "bPe", keep_on=False)
c.plot_tomo(m_bh, "Omth", keep_on=False)

# systematics
m_sys = [names[key] for key in [0, 1, 5, 6, 10]]
c.plot_tomo(m_sys, "sigma8", keep_on=False)
c.plot_tomo(m_sys, "ygk_mass_bias", keep_on=False)
c.plot_tomo(m_sys, "bPe", keep_on=False)

# triangles
m_tri = [names[key] for key in [1, 7]]
c.plot_triangles(m_tri, keep_on=False)

# systematics
m_tri_sys = ["yxgxksig", "yxgxksig_mf_despali16"]
c.plot_triangles(m_tri_sys, keep_on=True)

# best fits
c.plot_best_fit(names[1], keep_on=True)
c.plot_best_fit(names[7], keep_on=True)
c.plot_best_fit(names[2], keep_on=True)
c.plot_best_fit(names[3], keep_on=True)

# dndz
c.plot_nz(compare=True, keep_on=True, overwrite=True)
