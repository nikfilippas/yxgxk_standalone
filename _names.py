names_old = {
    0: "yxgxksig_hmc_Ahmc",       # Mmin, M1, 1-b, sigma8, rhogy, A
    1: "yxgxksig",                # Mmin, M1, 1-b, sigma8, rhogy
    2: "gxk",                     # Mmin, M1, sigma8
    3: "yxgxk_b08",               # Mmin, M1, sigma8, rhogy

    # robustness checks (include)
    4: "yxgxk",                   # Mmin, M1, 1-b, rhogy
    5: "yxgxksig_decouple_Ahmc",  # Mmin, M1, 1-b, sigma8, rhogy, Agg, Agy, Agm
    6: "yxgxk_szdeproj",          # SZdeproj : Mmin, M1, 1-b, sigma8, rhogy, A
    7: "yxgxk_mf_despali16",      # Despali16: Mmin, M1, 1-b, sigma8, rhogy
    8: "yxgxk_b_gauss",           # b_gauss : Mmin, M1, 1-b, sigma8, rhogy

    # robustness checks (do not include)
    9: "yxgxk_b_uniform",         # b_narrow: Mmin, M1, 1-b, sigma8, rhogy
    10: "yxgxk_ns_independent",   # ns_indep: Mmin, M1, 1-b, sigma8, rhogy

    # others
    11: "yxg",
    12: "dam_yxg",
    13: "gxk_kmax05",
    14: "yxgxksig_kmax05",
    15: "yxgxk_hmc_hmcode",
    16: "yxg_Ahmc",
}


names_new = {
    0: "gyksrA",        # Mmin, M1, 1-b, sigma8, rhogy, A
    1: "gyksr",         # Mmin, M1, 1-b, sigma8, rhogy
    2: "gksrA",         # Mmin, M1, sigma8, A
    3: "gyksrA_bf075",  # Mmin, M1, sigma8, rhogy, A

    # robustness checks (include)
    4: "gykrA",         # Mmin, M1, 1-b, rhogy, A
    5: "gyksrAAA",      # Mmin, M1, 1-b, sigma8, rhogy, Agg, Agy, Agm
    6: "gyksrA_SZ",     # SZdeproj : Mmin, M1, 1-b, sigma8, rhogy, A
    7: "gyksrA_D16",    # Despali16: Mmin, M1, 1-b, sigma8, rhogy, A
    8: "gyksrA_bG073",  # b_gauss  : Mmin, M1, 1-b, sigma8, rhogy, A

    # robustness checks (do not include)
    9: "gyksrA_bU075",  # b_narrow : Mmin, M1, 1-b, sigma8, rhogy, A
    10: "gyksrA_B16",   # Bocquet16: Mmin, M1, 1-b, sigma8, rhogy, A
    11: "gyksrA_Nsat",  # ns_indep : Mmin, M1, 1-b, sigma8, rhogy, A
}
