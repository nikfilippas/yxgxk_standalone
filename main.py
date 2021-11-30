from chanal import ChainCalculator

c = ChainCalculator()

m_s8 = ["yxgxksig", "yxgxk_b08",
        "yxgxk_b_uniform", "yxgxk_b_gauss",
        "gxk", "gxk_kmax05", ]
c.plot_tomo(m_s8, "sigma8", keep_on=True)

m_bh = ["yxgxksig", "yxgxk", "yxgxk_b_uniform",
        "yxgxk_b_gauss", "yxg", "dam_yxg"]
c.plot_tomo(m_bh, "ygk_mass_bias", keep_on=True)
c.plot_tomo(m_bh, "bPe", keep_on=True)
c.plot_tomo(m_bh, "Omth", keep_on=True)

m_tri = ["yxgxksig", "yxgxk_b_uniform", "yxgxk_b_gauss"]
c.plot_triangles(m_tri)
