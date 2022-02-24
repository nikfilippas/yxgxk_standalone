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
    0: "gyksrA_T08",        # Mmin, M1, 1-b, sigma8, rhogy, A
    1: "gyksr_T08",         # Mmin, M1, 1-b, sigma8, rhogy
    2: "gksrA_T08",         # Mmin, M1, sigma8, A
    3: "gyksrA_bf075_T08",  # Mmin, M1, sigma8, rhogy, A

    # robustness checks (include)
    4: "gykrA_T08",         # Mmin, M1, 1-b, rhogy, A
    5: "gyksrAAA_T08",      # Mmin, M1, 1-b, sigma8, rhogy, Agg, Agy, Agm
    6: "gyksrA_SZ_T08",     # SZdeproj : Mmin, M1, 1-b, sigma8, rhogy, A
    7: "gyksrA_D16",        # Despali16: Mmin, M1, 1-b, sigma8, rhogy, A
    8: "gyksrA_bG073_T08",  # b_gauss  : Mmin, M1, 1-b, sigma8, rhogy, A

    # robustness checks (do not include)
    9: "gyksrA_bU075_T08",  # b_narrow : Mmin, M1, 1-b, sigma8, rhogy, A
    10: "gyksrA_B16",       # Bocquet16: Mmin, M1, 1-b, sigma8, rhogy, A
    11: "gyksrA",           # Tinker10 : Mmin, M1, 1-b, sigma8, rhogy, A
    12: "gyksrA_Nsat_T08",  # ns_indep : Mmin, M1, 1-b, sigma8, rhogy, A
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

gxg = r"$g \times g$"
gxy = r"$g \times y$"
gxk = r"$g \times \kappa$"

latex_labels_new = {
    "gyksrA_T08": rf"fiducial: {gxg}, {gxy}, {gxk}",
    "gksrA_T08" : rf"no tSZ: {gxg}, {gxk}",
    "gyksrA_bf075_T08": r"fixed $1-b_{\rm H} = 0.75$",
    "gyksrA_bG073_T08": r"$1-b_{\rm H} \sim G(0.73, 0.01)$",
    "gyksrA_bU075_T08": r"$1-b_{\rm H} \sum U(0.60, 0.90)$",
    "gykrA_T08" : r"fixed $\sigma_8$",
    "gyksr_T08" : r"fixed $A_{\rm HM}$",
    "gyksrA_B16": "Bocquet et al. 2016 mass function",
    "gyksrA"    : "Tinker et al. 2010 mass function",
    "gyksrA_D16": "Despali et al. 2016 mass function",
    "gyksrA_Nsat_T08": r"",
    "gyksrA_SZ_T08": "Planck 2018 SZ-deprojected map",
    "gyksrAAA_T08" : r"decoupled $A_{\rm HM}^{gg}$, $A_{\rm HM}^{gy}$, $A_{\rm HM}^{g\kappa}$",
    }


colors = {
    "gyksrA_T08": "k",
    "gksrA_T08" : "grey",
    "gyksrA_bf075_T08": "r",
    "gyksrA_bG073_T08": "brown",
    "gyksrA_bU075_T08": "orangered",
    "gykrA_T08" : "forestgreen",
    "gyksr_T08" : "crimson",
    "gyksrA_B16": "navy",
    "gyksrA"    : "royalblue",
    "gyksrA_D16": "deepskyblue",
    "gyksrA_Nsat_T08": "darkslategrey",
    "gyksrA_SZ_T08": "magenta",
    "gyksrAAA_T08" : "orange",
    }
