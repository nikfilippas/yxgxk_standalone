from chanal import ChainCalculator

c = ChainCalculator()

# sigma8
m_s8 = ["yxgxksig", "yxgxk_b08",
        "yxgxk_b_uniform", "yxgxk_b_gauss",
        "gxk", "gxk_kmax05"]
c.plot_tomo(m_s8, "sigma8", keep_on=True)

# gastrophysics
m_bh = ["yxgxksig", "yxgxk", "yxgxk_b_uniform",
        "yxgxk_b_gauss", "yxg", "dam_yxg"]
c.plot_tomo(m_bh, "ygk_mass_bias", keep_on=True)
c.plot_tomo(m_bh, "bPe", keep_on=True)
c.plot_tomo(m_bh, "Omth", keep_on=True)

# systematics
m_sys = ["yxgxksig", "yxgxk_b_gauss", "gxk",
         "yxgxksig_hmc_hmcode", "yxgxksig_kmax05",
         "yxgxksig_mf_despali16", "yxgxgsig_ns_independent"]
c.plot_tomo(m_sys, "sigma8", keep_on=True)
c.plot_tomo(m_sys, "ygk_mass_bias", keep_on=True)

# triangles
m_tri = ["yxgxksig", "yxgxk_b_uniform", "yxgxk_b_gauss"]
c.plot_triangles(m_tri)

# best fits
c.plot_best_fit("yxgxksig", keep_on=True)
