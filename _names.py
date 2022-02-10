names_old = {
    0: "yxgxksig_hmc_Ahmc",       # Mmin, M1, 1-b, sigma8, rhogy, A
    1: "yxgxksig",                # Mmin, M1, 1-b, sigma8, rhogy
    2: "gxk",                     # Mmin, M1, sigma8
    3: "yxgxk_b08",               # Mmin, M1, sigma8, rhogy

    # robustness checks (include)
    4: "yxgxk",                   # Mmin, M1, 1-b, rhogy
    5: "yxgxksig_decouple_Ahmc",  # Mmin, M1, 1-b, sigma8, rhogy, Agg, Agy, Agm
    6: "yxgxksig_szdeproj",       # SZdeproj : Mmin, M1, 1-b, sigma8, rhogy, A
    7: "yxgxksig_mf_despali16",   # Despali16: Mmin, M1, 1-b, sigma8, rhogy
    8: "yxgxk_b_gauss",           # b_gauss : Mmin, M1, 1-b, sigma8, rhogy

    # robustness checks (do not include)
    9: "yxgxk_b_uniform",         # b_narrow: Mmin, M1, 1-b, sigma8, rhogy
    10: "yxgxksig_ns_independent",# ns_indep: Mmin, M1, 1-b, sigma8, rhogy

    # others
    11: "yxg",
    12: "dam_yxg",
    13: "gxk_kmax05",
    14: "yxgxksig_kmax05",
    15: "yxgxk_hmc_hmcode",
    16: "yxg_Ahmc",
}
names_old = {k: f"chains_bak/{v}" for k, v in names_old.items()}


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


latex_bins = [r"$\mathrm{2MPZ}$"] + \
    [r"$\mathrm{WI \times SC}$ - $\mathrm{%d}$" % i for i in range(1, 6)]


latex_names = {
    "bPe": "\\langle bP_e \\rangle\\ [\\mathrm{meV\\,cm^{-3}}]",
    "Omth": "\\Omega_{\mathrm{th}}",
    "ygk_mass_bias": "1-b_{\\mathrm{H}}",
    "sigma8": "\\sigma_8"
    }


latex_labels = {
    "yxgxksig": r"fiducial $3 \times 2\mathrm{pt}$ \,g,y,k",
    "yxgxk": r"$3 \times 2\mathrm{pt}$, \,g,y,k; fixed $\sigma_8=0.8122$",
    "yxgxk_b08": r"$3 \times 2\mathrm{pt}$ \,g,y,k; fixed $1-b_{\mathrm{H}}=0.80$",
    "yxgxk_b_uniform": r"$3 \times 2\mathrm{pt}$ \,g,y,k; $1-b_{\mathrm{H}} \sim U(0.60,0.90)$",
    "yxgxk_b_gauss": r"$3 \times 2\mathrm{pt}$ \,g,y,k; $1-b_{\mathrm{H}} \sim N(0.73,0.10)$",
    "gxk": r"$2 \times 2\mathrm{pt}$ \,g,k",
    "gxk_kmax05": r"$2 \times 2\mathrm{pt}$ \,g,k; $k_{\mathrm{max}}=0.5\,\mathrm{Mpc}^{-1}$",
    "yxg": r"Koukoufilippas et al., 2020",
    "dam_yxg": r"damonge Koukoufilippas et al., 2020",
    "yxgxksig_hmc_Ahmc": r"$3 \times 2\mathrm{pt}$, \,g,y,k; HALOFIT 1h/2h transition",
    "yxgxksig_hmc_hmcode": r"$3 \times 2\mathrm{pt}$, \,g,y,k; HMCode 1h/2h transition",
    "yxgxksig_kmax05": r"$3 \times 2\mathrm{pt}$, \,g,y,k; $k_{\mathrm{max}}=0.5\,\mathrm{Mpc}^{-1}$",
    "yxgxksig_mf_despali16": r"$3 \times 2\mathrm{pt}$, \,g,y,k; Despali 2016 mass function",
    "yxgxksig_ns_independent": r"$3 \times 2\mathrm{pt}$, \,g,y,k; HOD $N_{\mathrm{sat}}$ independent of $N_{\mathrm{cen}}$",
    "yxgxksig_szdeproj": r"$3 \times 2\mathrm{pt}$ \,g,y,k,\,\textrm{SZ-deprojected}",
    "yxgxksig_decouple_Ahmc": r"$3 \times 2\mathrm{pt}$, \,g,y,k; decoupled HALOFIT 1h/2h transitions",
    }
