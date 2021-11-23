from chanal import ChainCalculator

c = ChainCalculator()

m_s8 = ["yxgxksig", "yxgxk_b08", "gxk",
        "gxk_kmax05", "yxgxk_b_uniform", "yxgxk_b_gauss"]
c.plot_tomo(m_s8, "sigma8", m_s8, keep_on=True)

m_bh = ["yxgxksig", "yxgxk", "yxgxk_b_uniform",
        "yxgxk_b_gauss", "yxg", "dam_yxg"]
c.plot_tomo(m_bh, "ygk_mass_bias", m_bh)
c.plot_tomo(m_bh, "bPe", m_bh)
c.plot_tomo(m_bh, "Omth", m_bh, keep_on=True)

m_tri = ["yxgxksig", "yxgxk_b_uniform", "yxgxk_b_gauss"]
c.plot_triangles(m_tri)
