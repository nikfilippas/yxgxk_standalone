from chanal import ChainCalculator

c = ChainCalculator()

# sigma8
m_s8 = ["yxgxksig", "gxk", "yxgxk_b08"]
c.plot_tomo(m_s8, "sigma8", keep_on=True)

# gastrophysics
m_bh = ["yxgxksig", "yxgxk", "yxg", "yxgxk_b_gauss"]
c.plot_tomo(m_bh, "ygk_mass_bias", keep_on=True)
c.plot_tomo(m_bh, "bPe", keep_on=True)
c.plot_tomo(m_bh, "Omth", keep_on=True)

# systematics
m_sys = ["yxgxksig", "yxgxk_b_gauss",
         "yxgxksig_kmax05", "yxgxksig_hmc_hmcode",
         "yxgxksig_mf_despali16", "yxgxksig_ns_independent"]
c.plot_tomo(m_sys, "sigma8", keep_on=True)
c.plot_tomo(m_sys, "ygk_mass_bias", keep_on=True)
c.plot_tomo(m_sys, "bPe", keep_on=True)

# triangles
m_tri = ["yxgxksig", "yxgxk_b_uniform", "yxgxk_b_gauss"]
c.plot_triangles(m_tri, keep_on=True)

# systematics
m_tri_sys = ["yxgxksig", "yxgxksig_hmc_hmcode", "yxgxksig_mf_despali16"]
c.plot_triangles(m_tri_sys, keep_on=True)

# best fits
c.plot_best_fit("yxgxksig", keep_on=True)
c.plot_best_fit("yxgxksig_mf_despali16", keep_on=True)

# dndz
c.plot_nz(compare=True, keep_on=True, overwrite=True)
