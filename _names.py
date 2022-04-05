names = {
    0: 'gyksrA_T08',        # Mmin, M1, 1-b, sigma8, rhogy, A
    1: 'gyksr_T08',         # Mmin, M1, 1-b, sigma8, rhogy
    2: 'gksrA_T08',         # Mmin, M1, sigma8, A
    3: 'gyksrA_bf075_T08',  # Mmin, M1, sigma8, rhogy, A

    # robustness checks (plotted)
    4: 'gykrA_T08',         # Mmin, M1, 1-b, rhogy, A
    5: 'gyksrAAA_T08',      # Mmin, M1, 1-b, sigma8, rhogy, Agg, Agy, Agm
    6: 'gyksrA_SZ_T08',     # SZdeproj : Mmin, M1, 1-b, sigma8, rhogy, A
    7: 'gyksrA_D16',        # Despali16: Mmin, M1, 1-b, sigma8, rhogy, A
    8: 'gyksrA_bG073_T08',  # b_gauss  : Mmin, M1, 1-b, sigma8, rhogy, A

    # robustness checks (not plotted)
    9: 'gyksrA_bU075_T08',  # b_narrow   : Mmin, M1, 1-b, sigma8, rhogy, A
    10: 'gyksrA_B16',       # Bocquet16  : Mmin, M1, 1-b, sigma8, rhogy, A
    11: 'gyksrA',           # Tinker10   : Mmin, M1, 1-b, sigma8, rhogy, A
    12: 'gyksrA_Nsat_T08',  # ns_indep   : Mmin, M1, 1-b, sigma8, rhogy, A
    13: 'gyksrA_T08_km4',   # kmax4=0.5  : Mmin, M1, 1-b, sigma8, rhogy, A
    14: 'gyksrA_T08_l40',   # lmin = 40  : Mmin, M1, 1-b, sigma8, rhogy, A
    15: 'gyksrA_T08_w4',    # w_4(z) ~ U(20%) : Mmin, M1, 1-b, sigma8, rhogy, A, w
    16: 'gyksrA_T08_km4_2', # kmax4=0.25 : Mmin, M1, 1-b, sigma8, rhogy, A
    17: 'gyrA_T08',         # gxg, gxy   : Mmin, M1, 1-b, rhogy, A
}

latex_bins = [r'$\mathrm{2MPZ}$'] + \
    [r'$\mathrm{WI \times SC}$ - $\mathrm{%d}$' % i for i in range(1, 6)]

latex_names = {
    'bPe': '\\langle bP_e \\rangle\\ [\\mathrm{meV\\,cm^{-3}}]',
    'Omth': '10^{-8}\,\\Omega_{\mathrm{th}}',
    'ygk_mass_bias': '1-b_{\\mathrm{H}}',
    'sigma8': '\\sigma_8'
}

gxg = r'$g \times g$'
gxy = r'$g \times y$'
gxk = r'$g \times \kappa$'

latex_labels = {
    'gyksrA_T08': rf'fiducial: {gxg}, {gxy}, {gxk}',
    'gksrA_T08' : rf'no tSZ: {gxg}, {gxk}',
    'gyksrA_bf075_T08': r'fixed $1-b_{\rm H} = 0.75$',
    'gyksrA_bG073_T08': r'$1-b_{\rm H} \sim G(0.73, 0.10)$',
    'gyksrA_bU075_T08': r'$1-b_{\rm H} \sum U(0.60, 0.90)$',
    'gykrA_T08' : r'fixed $\sigma_8$',
    'gyksr_T08' : r'fixed $A_{\rm HM}$',
    'gyksrA_B16': 'Bocquet et al. 2016 mass function',
    'gyksrA'    : 'Tinker et al. 2010 mass function',
    'gyksrA_D16': 'Despali et al. 2016 mass function',
    'gyksrA_Nsat_T08': r'$N_{rm sat}$ independent from $N_{\rm cen}$',
    'gyksrA_SZ_T08': 'Planck 2018 SZ-deprojected map',
    'gyksrAAA_T08' : r'decoupled $A_{\rm HM}^{gg}$, $A_{\rm HM}^{gy}$, $A_{\rm HM}^{g\kappa}$',
    'gyksrA_T08_km4': r'$k_{\rm max}^{\rm wisc\,4} = 0.5\,\rm Mpc^{-1}$',
    'gyksrA_T08_l40': r'$l_{\rm min} = 40$',
    'gyksrA_T08_w4' : r'$w^{\rm wisc\,4}(z) \sim U(0.80, 1.20)$',
    'gyksrA_T08_km4_2': r'$k_{\rm max}^{\rm wisc\,4} = 0.25\,\rm Mpc^{-1}$',
    'gyrA_T08': rf'no $\kappa$: {gxg}, {gxy}',
}

latex_labels_short = {
    'gyksrA_T08': 'fiducial',
    'gksrA_T08' : 'no tSZ',
    'gyksrA_bf075_T08': r'$1-b_{\rm H} = 0.75$',
    'gyksrA_bG073_T08': r'$1-b_{\rm H} \sim \rm Gauss$',
    'gyksrA_bU075_T08': r'$1-b_{\rm H} \sum \rm Uniform$',
    'gykrA_T08' : r'$\sigma_8 = 0.8102$',
    'gyksr_T08' : r'$A_{\rm HM} = 1$',
    'gyksrA_B16': 'Bocquet et al. 2016',
    'gyksrA'    : 'Tinker et al. 2010',
    'gyksrA_D16': 'Despali et al. 2016',
    'gyksrA_Nsat_T08': r'independent $N_{\rm sat}$',
    'gyksrA_SZ_T08': 'SZ-deprojected',
    'gyksrAAA_T08' : r'$A_{\rm HM}^{gg}$, $A_{\rm HM}^{gy}$, $A_{\rm HM}^{g\kappa}$',
    'gyksrA_T08_km4': r'$k_{\rm max}^{\rm wisc\,4} = 0.5\,\rm Mpc^{-1}$',
    'gyksrA_T08_km4_2': r'$k_{\rm max}^{\rm wisc\,4} = 0.25\,\rm Mpc^{-1}$',
    'gyksrA_T08_l40': r'$l_{\rm min} = 40$',
    'gyrA_T08': r'no $\kappa$',
}

colors = {
    'gyksrA_T08': 'k',
    'gksrA_T08' : 'grey',
    'gyksrA_bf075_T08': 'r',
    'gyksrA_bG073_T08': 'brown',
    'gyksrA_bU075_T08': 'orangered',
    'gykrA_T08' : 'forestgreen',
    'gyksr_T08' : 'crimson',
    'gyksrA_B16': 'navy',
    'gyksrA'    : 'royalblue',
    'gyksrA_D16': 'deepskyblue',
    'gyksrA_Nsat_T08': 'darkslategrey',
    'gyksrA_SZ_T08': 'magenta',
    'gyksrAAA_T08' : 'orange',
    'gyksrA_T08_km4': 'palegreen',
    'gyksrA_T08_km4_2': 'greenyellow',
    'gyksrA_T08_w4' : 'lime',
    'gyksrA_T08_l40': 'indigo',
    'gyrA_T08': 'darkslategrey',
}

markers = {
    'gyksrA_T08': 'o',        # fiducial
    'gksrA_T08' : 's',        # no tSZ
    'gyksrA_bf075_T08': 'D',  # fixed bH
    'gyksrA_bG073_T08': '^',  # Gauss
    'gykrA_T08' : 'p',        # fixed σ8
    'gyksr_T08' : 'd',        # fixed Ahmc
    'gyksrA_SZ_T08': 'P',     # SZ-deproj
    'gyksrAAA_T08' : 'x',     # Agg Agy Agk
    'gyksrA_B16': 'v',        # Bocquet16
    'gyksrA': 'v',            # Tinker10
    'gyksrA_D16': 'v',        # Despali16
    'gyksrA_T08_km4': '*',    # kmax4 = 0.50
    'gyksrA_T08_km4_2': 'H',  # kmax4 = 0.25
    'gyksrA_T08_w4' : '1',    # w_4(z) prior
    'gyksrA_T08_l40': '_',    # lmin = 40
    'gyksrA_bU075_T08': '^',  # Uniform
    'gyrA_T08': 's',          # no CMBκ
}
