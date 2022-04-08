
names_new = {
    0: 'gyksrA',            # fiducial : Mmin, M1, 1-b, sigma8, rhogy, A
    1: 'gksA',              # gxg, gxk : Mmin, M1, sigma8, A
    2: 'gyksrA_bG073',      # b_gauss  : Mmin, M1, 1-b, sigma8, rhogy, A
    3: 'gyksrA_bf075',      # b fixed  : Mmin, M1, sigma8, rhogy, A
    4: 'gykrA',             # s8 fixed : Mmin, M1, 1-b, rhogy, A

    # robustness checks (plotted)
    5: 'gyksr',             # A fixed  : Mmin, M1, 1-b, sigma8, rhogy
    6: 'gyksrAAA',          # A decoup : Mmin, M1, 1-b, sigma8, rhogy, Agg, Agy, Agm
    7: 'gyksrA_SZ',         # SZdeproj : Mmin, M1, 1-b, sigma8, rhogy, A

    # mass functions
    8: 'gyksrA_B16',        # Bocquet16 : Mmin, M1, 1-b, sigma8, rhogy, A
    9: 'gyksrA_T10',        # Tinker10  : Mmin, M1, 1-b, sigma8, rhogy, A
    10: 'gyksrA_D16',       # Despali16 : Mmin, M1, 1-b, sigma8, rhogy, A

    # robustness checks (not plotted)
    11: 'gyksrA_lmin40',    # lmin = 40  : Mmin, M1, 1-b, sigma8, rhogy, A
    12: 'gyksrA_kmax025',   # kmax = 0.25: Mmin, M1, 1-b, sigma8, rhogy, A
    13: 'gyksrA_bU075',     # b_narrow   : Mmin, M1, 1-b, sigma8, rhogy, A
    14: 'gyksrA_Nsat',      # ns_indep   : Mmin, M1, 1-b, sigma8, rhogy, A
    15: 'gyr',              # gxg, gxy   : Mmin, M1, 1-b, rhogy
    16: 'gyksrAw',          # w(z)       : Mmin, M1, 1-b, sigma8, rhogy, A, w
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
    'gyksrA'        : rf'fiducial: {gxg}, {gxy}, {gxk}',
    'gksA'          : rf'no tSZ: {gxg}, {gxk}',
    'gyksrA_bG073'  : r'$1-b_{\rm H} \sim G(0.73, 0.10)$',
    'gyksrA_bf075'  : r'fixed $1-b_{\rm H} = 0.75$',
    'gykrA'         : r'fixed $\sigma_8$',

    'gyksr'         : r'fixed $A_{\rm HM}$',
    'gyksrAAA'      : r'decoupled $A_{\rm HM}^{gg}$, $A_{\rm HM}^{gy}$, $A_{\rm HM}^{g\kappa}$',
    'gyksrA_SZ'     : 'Planck 2018 SZ-deprojected map',

    'gyksrA_B16'    : 'Bocquet et al. 2016 mass function',
    'gyksrA_T10'    : 'Tinker et al. 2010 mass function',
    'gyksrA_D16'    : 'Despali et al. 2016 mass function',

    'gyksrA_lmin40' : r'$l_{\rm min} = 40$',
    'gyksrA_kmax025': r'$k_{\rm max}^{\rm wisc\,4} = 0.25\,\rm Mpc^{-1}$',
    'gyksrA_bU075'  : r'$1-b_{\rm H} \sum U(0.60, 0.90)$',
    'gyksrA_Nsat'   : r'$N_{rm sat}$ independent from $N_{\rm cen}$',
    'gyr'           : r'no $\kappa$: %s, %s, fixed $A_{\rm hmc}$' % (gxg, gxy),
    'gyksrAw'       : r'$w^(z) \sim U(0.80, 1.20)$',
}

latex_labels_short = {
    'gyksrA'        : 'fiducial',
    'gksA'          : 'no tSZ',
    'gyksrA_bG073'  : r'$1-b_{\rm H} \sim \rm Gauss$',
    'gyksrA_bf075'  : r'$1-b_{\rm H} = 0.75$',
    'gykrA'         : r'$\sigma_8 = 0.8102$',

    'gyksr'         : r'$A_{\rm HM} = 1$',
    'gyksrAAA'      : r'$A_{\rm HM}^{gg}$, $A_{\rm HM}^{gy}$, $A_{\rm HM}^{g\kappa}$',
    'gyksrA_SZ'     : 'SZ-deprojected',

    'gyksrA_B16'    : 'Bocquet et al. 2016',
    'gyksrA_T10'    : 'Tinker et al. 2010',
    'gyksrA_D16'    : 'Despali et al. 2016',

    'gyksrA_lmin40' : r'$l_{\rm min} = 40$',
    'gyksrA_kmax025' : r'$k_{\rm max} = 0.25\,\rm Mpc^{-1}$',
    'gyksrA_bU075'  : r'$1-b_{\rm H} \sum \rm Uniform$',
    'gyksrA_Nsat'   : r'independent $N_{\rm sat}$',
    'gyr'           : r'g,y, fixed $A_{\rm hmc}$',
    'gyksrAw'       : r'$w^(z) \sim U(0.80, 1.20)$',
}

colors = {
    'gyksrA'        : 'k',
    'gksA'          : 'grey',
    'gyksrA_bG073'  : 'brown',
    'gyksrA_bf075'  : 'r',
    'gykrA'         : 'forestgreen',

    'gyksr'         : 'crimson',
    'gyksrAAA'      : 'orange',
    'gyksrA_SZ'     : 'magenta',

    'gyksrA_B16'    : 'navy',
    'gyksrA_T10'    : 'royalblue',
    'gyksrA_D16'    : 'deepskyblue',

    'gyksrA_lmin40' : 'indigo',
    'gyksrA_kmax025': 'greenyellow',
    'gyksrA_bU075'  : 'orangered',
    'gyksrA_Nsat'   : 'darkslategrey',
    'gyr'           : 'darkslategrey',
    'gyksrAw'       : 'lime',
}

markers = {
    'gyksrA'        : 'o',
    'gksA'          : 's',
    'gyksrA_bG073'  : '^',
    'gyksrA_bf075'  : 'D',
    'gykrA'         : 'p',

    'gyksr'         : 'd',
    'gyksrAAA'      : 'x',
    'gyksrA_SZ'     : 'P',

    'gyksrA_B16'    : 'v',
    'gyksrA_T10'    : 'v',
    'gyksrA_D16'    : 'v',

    'gyksrA_lmin40' : '_',
    'gyksrA_kmax025': 'H',
    'gyksrA_bU075'  : '^',
    'gyksrA_Nsat'   : '>',
    'gyr'           : 's',
    'gyksrAw'       : '1',
}
